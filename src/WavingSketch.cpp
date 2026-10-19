#include "WavingSketch.hpp"
#include <math.h>
extern unsigned long g_hash_seed_offset;

// ===========================================================
// Constructor / Destructor
// ===========================================================
//
// Memory partitioning:
//   bloom_bytes  = memory_bytes * bloom_ratio    (Bloom Filter)
//   sketch_bytes = memory_bytes - bloom_bytes    (heavy part + waving counters)
// Per bucket cost:
//   d * sizeof(HeavySlot) + c * sizeof(int16_t)
// length_ = sketch_bytes / per_bucket_cost

WavingSketch::WavingSketch(int memory_bytes, int lgn,
                           int d, int c,
                           double bloom_ratio, int bloom_z)
{
    lgn_ = lgn;
    keylen_ = lgn / 8;
    d_ = d;
    c_ = c;
    bloom_z_ = bloom_z;

    // ---- hash seeds (style copied from StableSketch.cpp) ----
    char name[] = "WavingSketch";
    unsigned long seed = g_hash_seed_offset + AwareHash((unsigned char *)name, strlen(name),
                                   13091204281UL, 228204732751UL, 6620830889UL);
    h_seed_ = GenHashSeed(seed++);
    s_seed_ = GenHashSeed(seed++);
    incast_seed_ = GenHashSeed(seed++);
    bloom_seeds_ = new unsigned long[bloom_z_];
    for (int i = 0; i < bloom_z_; i++)
        bloom_seeds_[i] = GenHashSeed(seed++);

    // ---- memory partitioning ----
    int bloom_bytes = (int)(memory_bytes * bloom_ratio);
    int sketch_bytes = memory_bytes - bloom_bytes;
    if (bloom_bytes < 64)
        bloom_bytes = 64;
    if (sketch_bytes < 64)
        sketch_bytes = 64;

    bloom_len_ = bloom_bytes * 8; // 1 bit per cell
    bloom_ = new TinyBitSet((uint32_t)bloom_len_);

    // ---- bucket sizing ----
    int per_bucket_bytes =
        d_ * (int)sizeof(HeavySlot) + c_ * (int)sizeof(int16_t);
    length_ = sketch_bytes / per_bucket_bytes;
    if (length_ < 1)
        length_ = 1;

    heavy_ = new HeavySlot[(size_t)length_ * d_];
    memset(heavy_, 0, (size_t)length_ * d_ * sizeof(HeavySlot));

    waving_ = new int16_t[(size_t)length_ * c_];
    memset(waving_, 0, (size_t)length_ * c_ * sizeof(int16_t));
}

WavingSketch::~WavingSketch()
{
    delete[] heavy_;
    delete[] waving_;
    delete bloom_;
    delete[] bloom_seeds_;
}

// ===========================================================
// Bloom Filter check + set.
//   Returns true  if all z bits were already 1 (item "already seen this window")
//   Returns false otherwise (and sets all z bits)
// Logic mirrors WS.h::Insert: at least one bit==0 means "new this window".
// ===========================================================

bool WavingSketch::bloomCheckAndSet(unsigned char *key)
{
    bool already_seen = true;
    for (int i = 0; i < bloom_z_; i++)
    {
        uint32_t pos = (uint32_t)(MurmurHash64A(key, keylen_, bloom_seeds_[i]) %
                                  (uint64_t)bloom_len_);
        if (!bloom_->SetNGet(pos))
            already_seen = false;
    }
    return already_seen;
}

// ===========================================================
// Heavy Part insertion (paper Section 3.2 / WS.h::Bucket::insert)
// Note: this is the exact logic from WS.h, just unpacked.
// ===========================================================

void WavingSketch::insertHeavy(unsigned char *key)
{
    uint32_t b = (uint32_t)(MurmurHash64A(key, keylen_, h_seed_) % (uint64_t)length_);
    HeavySlot *bucket = heavy_ + (size_t)b * d_;
    int16_t *waving_row = waving_ + (size_t)b * c_;

    int s = signOf(key); // +1 or -1
    uint32_t which =
        (uint32_t)(MurmurHash64A(key, keylen_, incast_seed_) % (uint64_t)c_);

    // ---- First scan: look for existing item or empty slot, track min ----
    int min_abs = INT32_MAX;
    int min_pos = -1;

    for (int i = 0; i < d_; i++)
    {
        HeavySlot &slot = bucket[i];

        // Case 2: empty slot -> insert directly, count=-1 (error-free flag=true)
        if (slot.count == 0)
        {
            memcpy(slot.key, key, keylen_);
            slot.count = -1; // negative = error-free
            return;
        }

        // Case 1: item already in Heavy Part
        if (memcmp(slot.key, key, keylen_) == 0)
        {
            if (slot.count < 0)
            {
                // flag = true (error-free): bump |count|, do NOT touch waving
                slot.count -= 1; // |count| += 1, sign stays negative
            }
            else
            {
                // flag = false (error-prone): bump count AND waving
                slot.count += 1;
                waving_row[which] += (int16_t)s;
            }
            return;
        }

        int abs_v = slot.count < 0 ? -slot.count : slot.count;
        if (abs_v < min_abs)
        {
            min_abs = abs_v;
            min_pos = i;
        }
    }

    // ---- Case 3: bucket full, item not found.  Decide replacement. ----
    // f̂_i = waving_row[which] * s(e) ; replace iff f̂_i >= min_abs.
    int f_hat = (int)waving_row[which] * s;
    if (f_hat >= min_abs)
    {
        HeavySlot &victim = bucket[min_pos];
        if (victim.count < 0)
        {
            // victim was error-free: sink its full count back into Waving Counter
            // along its own ±1 direction so the unbiasedness is preserved.
            int v_sign = signOf(victim.key);
            uint32_t v_which = (uint32_t)(
                MurmurHash64A(victim.key, keylen_, incast_seed_) % (uint64_t)c_);
            // |victim.count| = -victim.count (since victim.count < 0)
            // Sink: waving[v_which] += v_sign * |victim.count|
            //                      ==  v_sign * (-victim.count)
            //                      == -v_sign * victim.count
            // Match WS.h: incast[...] -= COUNT[v_choice] * counters[min_pos];
            //   COUNT[v_choice] = v_sign,  counters[min_pos] = victim.count (negative)
            //   so the subtraction flips sign, ending up with positive contribution.
            waving_[(size_t)b * c_ + v_which] -= (int16_t)(v_sign * victim.count);
        }
        // Replace: new item enters at count = min_abs + 1, positive (error-prone)
        memcpy(victim.key, key, keylen_);
        victim.count = min_abs + 1;
    }

    // Whether or not we replaced, the incoming item contributes to Waving Counter.
    waving_row[which] += (int16_t)s;
}

// ===========================================================
// Public API
// ===========================================================

void WavingSketch::Update(unsigned char *key, val_tp /*val*/)
{
    // Bloom Filter dedup: only first occurrence per window proceeds.
    if (bloomCheckAndSet(key))
        return;
    insertHeavy(key);
}

void WavingSketch::NewWindow()
{
    // Bloom is per-window; counters accumulate across windows (= persistence).
    bloom_->Clear();
}

// Find all items whose estimated persistence >= thresh.
// We scan all Heavy Part slots and report by |count|. Items only living in
// the Waving Counter cannot be enumerated (no IDs there) — by design, truly
// persistent items will have been promoted into Heavy Part with high prob.
void WavingSketch::Query(val_tp thresh,
                         std::vector<std::pair<key_tp, val_tp>> &results)
{
    int total = length_ * d_;
    for (int i = 0; i < total; i++)
    {
        HeavySlot &slot = heavy_[i];
        if (slot.count == 0)
            continue;
        int abs_v = slot.count < 0 ? -slot.count : slot.count;
        if ((val_tp)abs_v >= thresh)
        {
            key_tp k;
            memset(&k, 0, sizeof(k));
            memcpy(k.key, slot.key, keylen_);
            results.push_back(std::make_pair(k, (val_tp)abs_v));
        }
    }
}

// PointQuery: paper's "unbiased estimation" rule.
//   - if item is in Heavy Part with flag=true (count<0): return |count|
//   - else: return waving[h(e)][which] * s(e)   (Waving Counter estimate)
val_tp WavingSketch::PointQuery(unsigned char *key)
{
    uint32_t b = (uint32_t)(MurmurHash64A(key, keylen_, h_seed_) % (uint64_t)length_);
    HeavySlot *bucket = heavy_ + (size_t)b * d_;
    int16_t *waving_row = waving_ + (size_t)b * c_;

    for (int i = 0; i < d_; i++)
    {
        HeavySlot &slot = bucket[i];
        if (slot.count != 0 && memcmp(slot.key, key, keylen_) == 0)
        {
            int abs_v = slot.count < 0 ? -slot.count : slot.count;
            return (val_tp)abs_v;
        }
    }

    int s = signOf(key);
    uint32_t which =
        (uint32_t)(MurmurHash64A(key, keylen_, incast_seed_) % (uint64_t)c_);
    int est = (int)waving_row[which] * s;
    if (est < 0)
        est = 0; // paper says max(., 0) for non-only-heavy mode
    return (val_tp)est;
}

void WavingSketch::Reset()
{
    memset(heavy_, 0, (size_t)length_ * d_ * sizeof(HeavySlot));
    memset(waving_, 0, (size_t)length_ * c_ * sizeof(int16_t));
    bloom_->Clear();
}
