#include "Hypersistent.hpp"
#include <math.h>
extern unsigned long g_hash_seed_offset;

// ===========================================================
// Constructor / Destructor
// ===========================================================

Hypersistent::Hypersistent(int memory_bytes, int lgn,
                           double hot_ratio,
                           int l1_pct,
                           int delta1, int delta2,
                           int burst_kb, int burst_g,
                           int cold_d1, int cold_d2,
                           int hot_beta)
{
    lgn_ = lgn;
    keylen_ = lgn / 8;
    delta1_ = delta1;
    delta2_ = delta2;

    // ---- hash seeds (style copied from StableSketch.cpp) ----
    char name[] = "Hypersistent";
    unsigned long seed = g_hash_seed_offset + AwareHash((unsigned char *)name, strlen(name),
                                   13091204281UL, 228204732751UL, 6620830889UL);
    burst_seed_ = GenHashSeed(seed++);
    hot_seed_ = GenHashSeed(seed++);

    cold_d1_ = cold_d1;
    cold_d2_ = cold_d2;
    L1_seeds_ = new unsigned long[cold_d1_];
    L2_seeds_ = new unsigned long[cold_d2_];
    for (int i = 0; i < cold_d1_; i++)
        L1_seeds_[i] = GenHashSeed(seed++);
    for (int i = 0; i < cold_d2_; i++)
        L2_seeds_[i] = GenHashSeed(seed++);

    // ---- memory partitioning ----
    int hot_bytes = (int)(memory_bytes * hot_ratio);
    int burst_bytes = burst_kb * 1024;
    if (burst_bytes >= memory_bytes - hot_bytes)
        burst_bytes = (memory_bytes - hot_bytes) / 4; // safety: don't starve cold filter
    int cold_bytes = memory_bytes - hot_bytes - burst_bytes;
    if (cold_bytes < 256)
        cold_bytes = 256;

    int m1_bytes = cold_bytes * l1_pct / 100;
    int m2_bytes = cold_bytes - m1_bytes;

    // ---- Burst Filter ----
    burst_g_ = burst_g;
    burst_w_ = burst_bytes / (keylen_ * burst_g_);
    if (burst_w_ < 1)
        burst_w_ = 1;
    burst_keys_ = new unsigned char[(size_t)burst_w_ * burst_g_ * keylen_]();

    // ---- Cold Filter ----
    // each cell: 1 byte counter + 1 bit flag = 9 bits
    // cells per row = (m_bytes * 8 / 9) / d
    auto cold_row_cells = [](int bytes, int d) -> int
    {
        long long total_bits = (long long)bytes * 8;
        long long cells_total = total_bits / 9;
        long long per_row = cells_total / (long long)d;
        if (per_row < 1)
            per_row = 1;
        return (int)per_row;
    };
    m1_ = cold_row_cells(m1_bytes, cold_d1_);
    m2_ = cold_row_cells(m2_bytes, cold_d2_);

    L1_count_ = new uint8_t[(size_t)cold_d1_ * m1_]();
    L2_count_ = new uint8_t[(size_t)cold_d2_ * m2_]();
    L1_flag_ = new TinyBitSet((uint32_t)(cold_d1_ * m1_));
    L2_flag_ = new TinyBitSet((uint32_t)(cold_d2_ * m2_));

    // ---- Hot Part ----
    // entry ≈ keylen + 4 (per) bytes; flag stored separately in BitSet
    hot_beta_ = hot_beta;
    int entry_bytes = keylen_ + 4;
    int total_slots = hot_bytes / entry_bytes;
    hot_lambda_ = total_slots / (HOT_T * hot_beta_);
    if (hot_lambda_ < 1)
        hot_lambda_ = 1;

    int hot_total = HOT_T * hot_lambda_ * hot_beta_;
    hot_entries_ = new HotEntry[hot_total];
    memset(hot_entries_, 0, (size_t)hot_total * sizeof(HotEntry));
    hot_flag_ = new TinyBitSet((uint32_t)hot_total);
}

Hypersistent::~Hypersistent()
{
    delete[] burst_keys_;
    delete[] L1_count_;
    delete[] L2_count_;
    delete L1_flag_;
    delete L2_flag_;
    delete[] L1_seeds_;
    delete[] L2_seeds_;
    delete[] hot_entries_;
    delete hot_flag_;
}

// ===========================================================
// Stage 1: Burst Filter
//   - Per-window cache of unique IDs.
//   - Returns true  => absorbed  (already in cache, OR fit into empty cell)
//     Returns false => bucket full, item not present  -> caller pushes to Cold.
// ===========================================================

bool Hypersistent::insertBurst(unsigned char *key)
{
    uint32_t bucket_idx =
        (uint32_t)(MurmurHash64A(key, keylen_, burst_seed_) % (uint64_t)burst_w_);
    unsigned char *bucket =
        burst_keys_ + (size_t)bucket_idx * burst_g_ * keylen_;

    int empty_slot = -1;
    for (int i = 0; i < burst_g_; i++)
    {
        unsigned char *cell = bucket + (size_t)i * keylen_;
        if (memcmp(cell, key, keylen_) == 0)
        {
            // already there; if cell happens to be all-zero it'd match an
            // all-zero key — but real items shouldn't be all zero. This is the
            // same sentinel convention as the original Filter.h.
            if (!isZeroKey(cell, keylen_))
                return true;
        }
        if (empty_slot < 0 && isZeroKey(cell, keylen_))
            empty_slot = i;
    }
    if (empty_slot >= 0)
    {
        memcpy(bucket + (size_t)empty_slot * keylen_, key, keylen_);
        return true;
    }
    return false; // overflow -> caller falls through to Cold Filter
}

// ===========================================================
// Stage 2: Cold Filter, one layer (used for both L1 and L2)
//   Logic mirrors Filter.h::Insert in the paper's reference code.
//   Returns true  => absorbed  (already counted this window OR CU-updated)
//   Returns false => saturated, push to next layer.
// ===========================================================

bool Hypersistent::insertColdLayer(uint8_t *counts, TinyBitSet *flags,
                                   unsigned long *seeds, int d, int m,
                                   unsigned char *key, int threshold)
{
    // d is small (paper default 2); avoid heap alloc
    int positions[16];
    int min_count = INT32_MAX;
    bool any_unlocked = false;

    for (int i = 0; i < d; i++)
    {
        int pos = (int)(MurmurHash64A(key, keylen_, seeds[i]) % (uint64_t)m);
        positions[i] = pos;
        int row_pos = i * m + pos;
        int cur = counts[row_pos];
        if (cur < min_count)
            min_count = cur;
        if (!flags->Get((uint32_t)row_pos))
            any_unlocked = true;
    }

    if (!any_unlocked)
    {
        // all positions already locked this window -> already counted
        return true;
    }

    if (min_count >= threshold)
    {
        // saturated: lock all so we don't repeatedly push this window,
        // and signal caller to push to next layer
        for (int i = 0; i < d; i++)
            flags->Set((uint32_t)(i * m + positions[i]));
        return false;
    }

    // CU update: increment cells where counter == min AND flag is unlocked.
    // SetNGet returns OLD flag value: 0=unlocked -> increment by 1, 1=locked -> +0.
    for (int i = 0; i < d; i++)
    {
        int row_pos = i * m + positions[i];
        if ((int)counts[row_pos] == min_count)
        {
            counts[row_pos] +=
                (uint8_t)(!flags->SetNGet((uint32_t)row_pos));
        }
    }
    return true;
}

// ===========================================================
// Stage 3: Hot Part — 2-table d-left hashing, β entries per bucket
//   Logic mirrors HotStorage.h::Insert. Empty slots have per=0; the
//   probabilistic-replace formula  hash % (min_per + 1) == 0  naturally
//   makes empty slots (min_per=0) ALWAYS-replace.
// ===========================================================

void Hypersistent::insertHot(unsigned char *key)
{
    uint64_t h = MurmurHash64A(key, keylen_, hot_seed_);
    uint32_t hash0 = (uint32_t)(h % (uint64_t)hot_lambda_);
    uint32_t hash1 = (uint32_t)((h / (uint64_t)hot_lambda_) % (uint64_t)hot_lambda_);

    int min_per = INT32_MAX;
    int min_table = -1;
    int min_slot = -1;

    // ---- Search table 0 ----
    HotEntry *bucket0 = hot_entries_ + (size_t)(0 * hot_lambda_ + hash0) * hot_beta_;
    int base0 = (0 * hot_lambda_ + (int)hash0) * hot_beta_;
    for (int i = 0; i < hot_beta_; i++)
    {
        HotEntry &e = bucket0[i];
        if (e.per > 0 && memcmp(e.key, key, keylen_) == 0)
        {
            // found: increment iff not yet locked this window
            if (!hot_flag_->SetNGet((uint32_t)(base0 + i)))
                e.per += 1;
            return;
        }
        if ((int)e.per < min_per)
        {
            min_per = (int)e.per;
            min_table = 0;
            min_slot = i;
        }
    }

    // ---- Search table 1 ----
    HotEntry *bucket1 = hot_entries_ + (size_t)(1 * hot_lambda_ + hash1) * hot_beta_;
    int base1 = (1 * hot_lambda_ + (int)hash1) * hot_beta_;
    for (int i = 0; i < hot_beta_; i++)
    {
        HotEntry &e = bucket1[i];
        if (e.per > 0 && memcmp(e.key, key, keylen_) == 0)
        {
            if (!hot_flag_->SetNGet((uint32_t)(base1 + i)))
                e.per += 1;
            return;
        }
        if ((int)e.per < min_per)
        {
            min_per = (int)e.per;
            min_table = 1;
            min_slot = i;
        }
    }

    // ---- Not found: probabilistic replace at min slot ----
    if (min_table < 0)
        return; // shouldn't happen, but guard

    if ((h % ((uint64_t)min_per + 1ULL)) == 0)
    {
        uint32_t hash_chosen = (min_table == 0) ? hash0 : hash1;
        int base = (min_table * hot_lambda_ + (int)hash_chosen) * hot_beta_;
        HotEntry &e = hot_entries_[base + min_slot];
        memcpy(e.key, key, keylen_);
        // Increment iff slot's flag was unlocked this window.
        // NOTE: if the displaced slot was already locked (updated earlier this
        // window by some other item), the new occupant inherits per without +1.
        // This matches HotStorage.h's behaviour and keeps the at-most-once
        // invariant slot-local rather than item-local.
        if (!hot_flag_->SetNGet((uint32_t)(base + min_slot)))
            e.per += 1;
    }
}

// ===========================================================
// At end of window: drain Burst Filter into Cold/Hot, then clear flags.
// Paper III-D: "all items in the buckets are iterated over and inserted
//              into the Stage 2 Cold Filter. The Burst Filter is then cleared."
// We do the flush BEFORE clearing flags so the at-most-once invariant
// of the just-finished window still holds during the drain.
// ===========================================================

void Hypersistent::flushBurst()
{
    for (int b = 0; b < burst_w_; b++)
    {
        for (int c = 0; c < burst_g_; c++)
        {
            unsigned char *cell =
                burst_keys_ + (size_t)(b * burst_g_ + c) * keylen_;
            if (isZeroKey(cell, keylen_))
                continue;
            // Standard Cold->Hot cascade (skip Burst Filter on flush).
            if (!insertColdLayer(L1_count_, L1_flag_, L1_seeds_,
                                 cold_d1_, m1_, cell, delta1_))
            {
                if (!insertColdLayer(L2_count_, L2_flag_, L2_seeds_,
                                     cold_d2_, m2_, cell, delta2_))
                {
                    insertHot(cell);
                }
            }
            memset(cell, 0, keylen_);
        }
    }
}

// ===========================================================
// Public API
// ===========================================================

void Hypersistent::Update(unsigned char *key, val_tp /*val*/)
{
    // val is ignored: persistence is per-window-occurrence, not packet count.
    if (insertBurst(key))
        return;
    if (insertColdLayer(L1_count_, L1_flag_, L1_seeds_,
                        cold_d1_, m1_, key, delta1_))
        return;
    if (insertColdLayer(L2_count_, L2_flag_, L2_seeds_,
                        cold_d2_, m2_, key, delta2_))
        return;
    insertHot(key);
}

void Hypersistent::NewWindow()
{
    // 1. Clear all flags first so the upcoming flush can actually +1 the
    //    cold/hot cells (otherwise burst-overflow items from the just-
    //    finished window would have left those cells locked and the flush
    //    would be silently rejected).
    L1_flag_->Clear();
    L2_flag_->Clear();
    hot_flag_->Clear();

    // 2. Drain burst into cold/hot.  flushBurst() reuses the cold-filter
    //    insert path which sets flags on success; this prevents two burst
    //    IDs that hash-collide from double-incrementing the same cell.
    flushBurst();

    // 3. Clear flags again.  After flush, flags are partially saturated by
    //    flush itself.  If we left them set, the *new* window's first
    //    legitimate +1 on those cells would be wrongly rejected as
    //    "already counted this window".  Re-clearing makes the flag
    //    invariant hold: at the moment a new window opens, all cells are
    //    unlocked and ready to receive their first within-window +1.
    L1_flag_->Clear();
    L2_flag_->Clear();
    hot_flag_->Clear();
}

// Find all items whose estimated persistence >= thresh.
// Cold Filter cells don't store IDs, so by construction only Hot Part can
// supply candidates. This is fine as long as thresh > delta1 + delta2,
// which is the intended operating regime (paper sets delta1+delta2=115,
// while alpha*w is typically several hundred).
void Hypersistent::Query(val_tp thresh,
                         std::vector<std::pair<key_tp, val_tp>> &results)
{
    val_tp baseline = (val_tp)(delta1_ + delta2_);
    int total = HOT_T * hot_lambda_ * hot_beta_;
    for (int i = 0; i < total; i++)
    {
        HotEntry &e = hot_entries_[i];
        if (e.per == 0)
            continue;
        val_tp est = baseline + (val_tp)e.per;
        if (est >= thresh)
        {
            key_tp k;
            memset(&k, 0, sizeof(k));
            memcpy(k.key, e.key, keylen_);
            results.push_back(std::make_pair(k, est));
        }
    }
}

// Algorithm 5 (after-window mode — skip the Burst Filter check).
val_tp Hypersistent::PointQuery(unsigned char *key)
{
    // L1
    int v1 = INT32_MAX;
    for (int i = 0; i < cold_d1_; i++)
    {
        int pos = (int)(MurmurHash64A(key, keylen_, L1_seeds_[i]) % (uint64_t)m1_);
        int v = L1_count_[i * m1_ + pos];
        if (v < v1)
            v1 = v;
    }
    if (v1 < delta1_)
        return (val_tp)v1;

    // L2
    int v2 = INT32_MAX;
    for (int i = 0; i < cold_d2_; i++)
    {
        int pos = (int)(MurmurHash64A(key, keylen_, L2_seeds_[i]) % (uint64_t)m2_);
        int v = L2_count_[i * m2_ + pos];
        if (v < v2)
            v2 = v;
    }
    if (v2 < delta2_)
        return (val_tp)(delta1_ + v2);

    // Hot Part
    val_tp baseline = (val_tp)(delta1_ + delta2_);
    uint64_t h = MurmurHash64A(key, keylen_, hot_seed_);
    uint32_t hash0 = (uint32_t)(h % (uint64_t)hot_lambda_);
    uint32_t hash1 = (uint32_t)((h / (uint64_t)hot_lambda_) % (uint64_t)hot_lambda_);

    HotEntry *b0 = hot_entries_ + (size_t)(0 * hot_lambda_ + hash0) * hot_beta_;
    for (int i = 0; i < hot_beta_; i++)
    {
        if (b0[i].per > 0 && memcmp(b0[i].key, key, keylen_) == 0)
            return baseline + (val_tp)b0[i].per;
    }
    HotEntry *b1 = hot_entries_ + (size_t)(1 * hot_lambda_ + hash1) * hot_beta_;
    for (int i = 0; i < hot_beta_; i++)
    {
        if (b1[i].per > 0 && memcmp(b1[i].key, key, keylen_) == 0)
            return baseline + (val_tp)b1[i].per;
    }
    // Cold Filter saturated but not promoted to Hot Part — return baseline as
    // a conservative lower bound (per >= delta1+delta2).
    return baseline;
}

void Hypersistent::Reset()
{
    memset(burst_keys_, 0, (size_t)burst_w_ * burst_g_ * keylen_);
    memset(L1_count_, 0, (size_t)cold_d1_ * m1_);
    memset(L2_count_, 0, (size_t)cold_d2_ * m2_);
    L1_flag_->Clear();
    L2_flag_->Clear();
    int total = HOT_T * hot_lambda_ * hot_beta_;
    memset(hot_entries_, 0, (size_t)total * sizeof(HotEntry));
    hot_flag_->Clear();
}
