#ifndef OnOff_H
#define OnOff_H

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

class OnOff
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

    // ---- Per-bucket KV pair ----
    struct KV
    {
        unsigned char key[LGN];
        uint32_t counter;
    };

public:
    // memory_bytes : total memory budget in bytes
    // lgn          : key width in bits
    // w_param      : KV pairs per bucket. Paper Section 6.3.2 evaluates
    //                w in {2, 8, 32}; smaller w favours throughput at the
    //                cost of accuracy.
    OnOff(int memory_bytes, int lgn, int w_param = 2);

    ~OnOff();

    void Update(unsigned char *key, val_tp val); // val ignored
    val_tp PointQuery(unsigned char *key);
    void Query(val_tp thresh, std::vector<std::pair<key_tp, val_tp>> &results);
    void NewWindow();
    void Reset();

private:
    int lgn_, keylen_;
    int w_;
    int length_;

    KV *kv_table_;
    uint32_t *sketch_;
    TinyBitSet *bucket_flag_;
    TinyBitSet *sketch_flag_;

    unsigned long hash_seed_;

    // disable copy
    OnOff(const OnOff &);
    OnOff &operator=(const OnOff &);
};

#endif // OnOff_H
