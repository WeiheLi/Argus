#include "OnOff.hpp"
#include <math.h>
extern unsigned long g_hash_seed_offset;

// ===========================================================
// Constructor / Destructor
// ===========================================================
//
// Memory accounting (per bucket):
//   - w KV pairs        : w * (keylen + 4) bytes
//   - 1 sketch counter  : 4 bytes
//   - w bucket flag bits: w bits
//   - 1 sketch flag bit : 1 bit

OnOff::OnOff(int memory_bytes, int lgn, int w_param)
{
    lgn_ = lgn;
    keylen_ = lgn / 8;
    w_ = w_param;

    // ---- hash seed (style copied from StableSketch.cpp) ----
    char name[] = "OnOffSketch";
    unsigned long seed = g_hash_seed_offset + AwareHash((unsigned char *)name, strlen(name),
                                   13091204281UL, 228204732751UL, 6620830889UL);
    hash_seed_ = GenHashSeed(seed++);

    // ---- compute length_ from total memory budget ----
    long long bits_per_bucket =
        8LL * ((long long)w_ * (keylen_ + 4) + 4) + (w_ + 1);
    long long total_bits = (long long)memory_bytes * 8;
    length_ = (int)(total_bits / bits_per_bucket);
    if (length_ < 1)
        length_ = 1;

    // ---- allocate ----
    kv_table_ = new KV[(size_t)length_ * w_];
    memset(kv_table_, 0, (size_t)length_ * w_ * sizeof(KV));

    sketch_ = new uint32_t[length_];
    memset(sketch_, 0, (size_t)length_ * sizeof(uint32_t));

    bucket_flag_ = new TinyBitSet((uint32_t)(length_ * w_));
    sketch_flag_ = new TinyBitSet((uint32_t)length_);
}

OnOff::~OnOff()
{
    delete[] kv_table_;
    delete[] sketch_;
    delete bucket_flag_;
    delete sketch_flag_;
}

// ===========================================================
// Insert — paper Section 3.3 original description
// (NOT the SIMD-friendly optimisation in Section 3.4 / Theorem 4.4).
//
//   (A) item already in bucket -> increment its KV counter (if unlocked)
//   (B) else, increment sketch[pos] (if unlocked)
//   (C) find min KV in bucket
//   (D) if sketch[pos] > B[pos]_min -> swap them, install item at min slot
//
// Note: replacement does NOT touch bucket_flag_ for the new occupant.
// ===========================================================

void OnOff::Update(unsigned char *key, val_tp /*val*/)
{
    uint32_t pos = (uint32_t)(MurmurHash64A(key, keylen_, hash_seed_) % (uint64_t)length_);
    KV *bucket = kv_table_ + (size_t)pos * w_;
    uint32_t bucket_bit_base = pos * w_;

    // (A) item already in bucket
    for (int i = 0; i < w_; i++)
    {
        if (bucket[i].counter > 0 && memcmp(bucket[i].key, key, keylen_) == 0)
        {
            if (!bucket_flag_->SetNGet(bucket_bit_base + i))
                bucket[i].counter += 1;
            return;
        }
    }

    // (B) increment sketch counter (at most once per window)
    if (!sketch_flag_->SetNGet(pos))
        sketch_[pos] += 1;

    // (C) find min KV in this bucket (linear scan; no Theorem 4.4 shortcut)
    int min_idx = 0;
    uint32_t min_val = bucket[0].counter;
    for (int i = 1; i < w_; i++)
    {
        if (bucket[i].counter < min_val)
        {
            min_val = bucket[i].counter;
            min_idx = i;
        }
    }

    // (D) if sketch[pos] > B[pos]_min, swap into the min slot
    if (sketch_[pos] > min_val)
    {
        memcpy(bucket[min_idx].key, key, keylen_);
        uint32_t tmp = sketch_[pos];
        sketch_[pos] = bucket[min_idx].counter;
        bucket[min_idx].counter = tmp;
        // intentionally do NOT lock bucket_flag_ for the swapped-in slot:
        // sticking to the paper's Section 3.3 textual description, which
        // omits this lock. The optimised version (Theorem 4.4) would lock it.
    }
}

// ===========================================================
// Query: enumerate all KV pairs whose counter >= thresh.
// ===========================================================

void OnOff::Query(val_tp thresh, std::vector<std::pair<key_tp, val_tp>> &results)
{
    int total = length_ * w_;
    for (int i = 0; i < total; i++)
    {
        KV &kv = kv_table_[i];
        if (kv.counter == 0)
            continue;
        if ((val_tp)kv.counter >= thresh)
        {
            key_tp k;
            memset(&k, 0, sizeof(k));
            memcpy(k.key, kv.key, keylen_);
            results.push_back(std::make_pair(k, (val_tp)kv.counter));
        }
    }
}

// ===========================================================
// PointQuery
//   - if item is in bucket  -> B[h1(e)][e]
//   - else                  -> sketch[h1(e)]
// ===========================================================

val_tp OnOff::PointQuery(unsigned char *key)
{
    uint32_t pos = (uint32_t)(MurmurHash64A(key, keylen_, hash_seed_) % (uint64_t)length_);
    KV *bucket = kv_table_ + (size_t)pos * w_;
    for (int i = 0; i < w_; i++)
    {
        if (bucket[i].counter > 0 && memcmp(bucket[i].key, key, keylen_) == 0)
            return (val_tp)bucket[i].counter;
    }
    return (val_tp)sketch_[pos];
}

// ===========================================================
// NewWindow: reset all flags ("set state to On"). Counters untouched.
// ===========================================================

void OnOff::NewWindow()
{
    bucket_flag_->Clear();
    sketch_flag_->Clear();
}

void OnOff::Reset()
{
    memset(kv_table_, 0, (size_t)length_ * w_ * sizeof(KV));
    memset(sketch_, 0, (size_t)length_ * sizeof(uint32_t));
    bucket_flag_->Clear();
    sketch_flag_->Clear();
}
