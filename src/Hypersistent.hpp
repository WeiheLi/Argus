#ifndef Hypersistent_H
#define Hypersistent_H

#include <vector>
#include <utility>
#include <cstring>
#include <cstdlib>
#include <cstdint>
#include <climits>
#include <iostream>
#include "datatypes.hpp"
extern "C"
{
#include "hash.h"
#include "util.h"
}

class Hypersistent
{
    // ---- Internal: tiny inline bitset ----
    struct TinyBitSet
    {
        uint64_t *bits;
        uint32_t n_words;
        uint32_t n_bits;

        explicit TinyBitSet(uint32_t n)
        {
            n_bits = n;
            n_words = (n + 63) / 64;
            bits = new uint64_t[n_words]();
        }
        ~TinyBitSet() { delete[] bits; }
        bool Get(uint32_t i) const
        {
            return ((bits[i >> 6] >> (i & 63)) & 1ULL) != 0;
        }
        void Set(uint32_t i)
        {
            bits[i >> 6] |= (1ULL << (i & 63));
        }
        // returns OLD value, sets to 1
        bool SetNGet(uint32_t i)
        {
            uint64_t mask = 1ULL << (i & 63);
            bool old = (bits[i >> 6] & mask) != 0;
            bits[i >> 6] |= mask;
            return old;
        }
        void Clear()
        {
            memset(bits, 0, (size_t)n_words * sizeof(uint64_t));
        }
    };

    // ---- Hot Part entry ----
    struct HotEntry
    {
        unsigned char key[LGN];
        uint32_t per; // persistence accumulated above (delta1 + delta2)
                      // per == 0  &&  key all-zero  =>  empty slot
    };

public:
    // memory_bytes: total memory budget in bytes (e.g. memory_kb * 1024)
    // lgn:          key width in bits (8 * keylen)
    // hot_ratio:    fraction of memory given to Hot Part        (paper finding-task default 0.4)
    // l1_pct:       L1 share within Cold Filter, in percent     (paper default 85, i.e. 17:3)
    // delta1, delta2:                                            paper defaults 15, 100
    // burst_kb:     Burst Filter size in KB                     (paper finding-task default 1)
    // burst_g:      cells per Burst Filter bucket               (paper default 16, SIMD-friendly)
    // cold_d1, cold_d2:  hash functions per Cold Filter layer   (paper default 2, 2)
    // hot_beta:     entries per Hot Part bucket                 (paper code SLOT=3)
    Hypersistent(int memory_bytes, int lgn,
                 double hot_ratio = 0.4,
                 int l1_pct = 85,
                 int delta1 = 15,
                 int delta2 = 100,
                 int burst_kb = 1,
                 int burst_g = 16,
                 int cold_d1 = 2,
                 int cold_d2 = 2,
                 int hot_beta = 3);

    ~Hypersistent();

    void Update(unsigned char *key, val_tp val); // val is ignored (persistence ≠ frequency)
    val_tp PointQuery(unsigned char *key);
    void Query(val_tp thresh, std::vector<std::pair<key_tp, val_tp>> &results);
    void NewWindow();
    void Reset();

private:
    // ---- config ----
    int lgn_, keylen_;
    int delta1_, delta2_;

    // ---- Stage 1: Burst Filter ----
    int burst_w_; // # buckets
    int burst_g_; // # cells per bucket
    unsigned char *burst_keys_; // burst_w_ * burst_g_ * keylen_ bytes (flat)
    unsigned long burst_seed_;

    // ---- Stage 2: Cold Filter (two layers) ----
    int cold_d1_, cold_d2_;
    int m1_, m2_;
    uint8_t *L1_count_;     // cold_d1_ * m1_  (1 byte per cell)
    uint8_t *L2_count_;     // cold_d2_ * m2_
    TinyBitSet *L1_flag_;   // 1 bit per cell
    TinyBitSet *L2_flag_;
    unsigned long *L1_seeds_;
    unsigned long *L2_seeds_;

    // ---- Stage 3: Hot Part (2-table d-left, β entries per bucket) ----
    static const int HOT_T = 2;
    int hot_lambda_; // # buckets per table
    int hot_beta_;   // # entries per bucket
    HotEntry *hot_entries_; // HOT_T * hot_lambda_ * hot_beta_
    TinyBitSet *hot_flag_;
    unsigned long hot_seed_;

    // ---- helpers ----
    bool insertBurst(unsigned char *key);
    bool insertColdLayer(uint8_t *counts, TinyBitSet *flags,
                         unsigned long *seeds, int d, int m,
                         unsigned char *key, int threshold);
    void insertHot(unsigned char *key);
    void flushBurst();

    static bool isZeroKey(const unsigned char *k, int n)
    {
        for (int i = 0; i < n; i++)
            if (k[i] != 0)
                return false;
        return true;
    }

    // disable copy
    Hypersistent(const Hypersistent &);
    Hypersistent &operator=(const Hypersistent &);
};

#endif // Hypersistent_H
