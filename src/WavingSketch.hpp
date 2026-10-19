#ifndef WavingSketch_H
#define WavingSketch_H

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

class WavingSketch
{
    // ---- Internal: tiny inline bitset (used by Bloom Filter) ----
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

    // ---- Heavy Part slot ----
    // count > 0  : "error-prone" (flag = false in paper) — counter from Waving promotion
    // count < 0  : "error-free"  (flag = true  in paper) — fresh insert / in-bucket update
    // count == 0 : empty slot
    // Real counter magnitude is |count|.
    struct HeavySlot
    {
        unsigned char key[LGN];
        int32_t count;
    };

public:
    // memory_bytes : total memory budget in bytes
    // lgn          : key width in bits
    // d            : Heavy Part slots per bucket (paper recommends 8~16)
    // c            : Waving Counters per bucket  (paper code uses several)
    // bloom_ratio  : fraction of memory given to Bloom Filter (paper: ~half)
    // bloom_z      : # hash functions in Bloom Filter (paper finding-persistent default 3)
    WavingSketch(int memory_bytes, int lgn,
                 int d = 16,
                 int c = 4,
                 double bloom_ratio = 0.5,
                 int bloom_z = 3);

    ~WavingSketch();

    void Update(unsigned char *key, val_tp val); // val ignored
    val_tp PointQuery(unsigned char *key);
    void Query(val_tp thresh, std::vector<std::pair<key_tp, val_tp>> &results);
    void NewWindow();
    void Reset();

private:
    int lgn_, keylen_;
    int d_, c_;       // d slots per Heavy Part, c Waving Counters per bucket
    int length_;      // # buckets

    HeavySlot *heavy_; // length_ * d_   (flat)
    int16_t *waving_;  // length_ * c_

    // Bloom Filter
    int bloom_z_;
    int bloom_len_;
    TinyBitSet *bloom_;
    unsigned long *bloom_seeds_;

    unsigned long h_seed_;       // bucket-locator hash
    unsigned long s_seed_;       // sign hash s(e) -> {+1,-1}
    unsigned long incast_seed_;  // chooses Waving Counter within a bucket

    // ---- helpers ----
    bool bloomCheckAndSet(unsigned char *key); // returns true if "already seen this window"
    void insertHeavy(unsigned char *key);
    inline int signOf(unsigned char *key) const
    {
        // s(e) -> {+1, -1} via low bit
        return ((MurmurHash64A(key, keylen_, s_seed_) & 1ULL) ? +1 : -1);
    }

    // disable copy
    WavingSketch(const WavingSketch &);
    WavingSketch &operator=(const WavingSketch &);
};

#endif // WavingSketch_H
