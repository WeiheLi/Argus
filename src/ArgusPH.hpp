#ifndef ArgusPH_H
#define ArgusPH_H
#include <vector>
#include <utility>
#include <cstring>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "datatypes.hpp"
extern "C"
{
#include "hash.h"
#include "util.h"
}

/*
 * LCBSketch-PH: Persistent + Heavy Item Detection.
 *
 * Task: find items with V/T >= alpha AND C/V >= kappa.
 *
 * Confidence bounds via Hoeffding with failure probability delta:
 *   epsilon(dt) = sqrt(ln(2/delta) / (2 * dt))
 *   LCB = phat - epsilon,  UCB = phat + epsilon
 *
 * Struct-based bucket (16 bytes, no bit-packing, no hot_C):
 *   key            : 8 bytes
 *   V              : 2 bytes (uint16_t)
 *   C              : 2 bytes (uint16_t, capped at 65535)
 *   first_tracked  : 2 bytes (uint16_t)
 *   status         : 1 byte
 *   _pad           : 1 byte
 */
class ArgusPH
{
#pragma pack(push, 1)
    struct Bucket
    {
        unsigned char key[LGN]; // 8
        uint16_t V;             // 2
        uint16_t C;             // 2
        uint16_t first_tracked; // 2
        uint8_t  status;        // 1
    };                          // 15 bytes
#pragma pack(pop)

    struct lcb_type
    {
        val_tp sum;
        Bucket *counts;
        int depth, width, lgn;
        unsigned long *hash, *scale, *hardner;
    };

public:
    ArgusPH(int depth, int width, int lgn,
                double alpha = 0.3, double kappa = 10.0,
                double delta = 0.05, int n_windows = 2047);
    ~ArgusPH();

    void Update(unsigned char *key, val_tp value);
    void NewWindow();
    void QueryPH(std::vector<std::pair<key_tp, val_tp>> &results);

    struct PHResult { key_tp key; val_tp V; uint16_t C; };
    void QueryPHDetailed(std::vector<PHResult> &results);

    val_tp PointQuery(unsigned char *key);
    val_tp GetCount();
    void Reset();
    int GetTCurrent() const { return t_current; }
    static int BucketBytes() { return (int)sizeof(Bucket); }

private:
    inline void init_bucket(int idx, unsigned char *key, int keylen) {
        Bucket *b = &lcb_.counts[idx];
        memcpy(b->key, key, keylen);
        b->V = 1;
        b->C = 1;
        b->first_tracked = (uint16_t)t_current;
        b->status = 1;
    }

    lcb_type lcb_;
    int t_current;
    double alpha, kappa, delta;
    double *eps_table;   // eps_table[dt] = sqrt(ln(2/delta) / (2 * dt))
    int max_windows;
};

#endif
