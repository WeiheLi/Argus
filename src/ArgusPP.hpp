#ifndef ArgusPP_H
#define ArgusPP_H
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
 * LCBSketch-PP: Persistent + Priority Item Detection.
 *
 * Task: find items with V/T >= alpha AND priority >= pr_min.
 *
 * Confidence bounds via Hoeffding with failure probability delta:
 *   epsilon(dt) = sqrt(ln(2/delta) / (2 * dt))
 *   LCB = phat - epsilon,  UCB = phat + epsilon
 *
 * Struct-based bucket (14 bytes, no bit-packing):
 *   key            : 8 bytes
 *   V              : 2 bytes (uint16_t)
 *   first_tracked  : 2 bytes (uint16_t)
 *   status         : 1 byte
 *   priority       : 1 byte  (0-7)
 */
class ArgusPP
{
#pragma pack(push, 1)
    struct Bucket
    {
        unsigned char key[LGN]; // 8
        uint16_t V;             // 2
        uint16_t first_tracked; // 2
        uint8_t  status;        // 1
        uint8_t  priority;      // 1
    };                          // 14 bytes
#pragma pack(pop)

    struct lcb_type
    {
        val_tp sum;
        Bucket *counts;
        int depth, width, lgn;
        unsigned long *hash, *scale, *hardner;
    };

public:
    ArgusPP(int depth, int width, int lgn,
                double alpha = 0.3, unsigned char pr_min = 4,
                double delta = 0.05, int n_windows = 2047);
    ~ArgusPP();

    void Update(unsigned char *key, val_tp value, unsigned char priority);
    void NewWindow();
    void QueryPP(std::vector<std::pair<key_tp, val_tp>> &results);

    struct PPResult { key_tp key; val_tp V; unsigned char priority; };
    void QueryPPDetailed(std::vector<PPResult> &results);

    val_tp PointQuery(unsigned char *key);
    val_tp GetCount();
    void Reset();
    int GetTCurrent() const { return t_current; }
    static int BucketBytes() { return (int)sizeof(Bucket); }

private:
    inline void init_bucket(int idx, unsigned char *key, int keylen, unsigned char pr) {
        Bucket *b = &lcb_.counts[idx];
        memcpy(b->key, key, keylen);
        b->V = 1;
        b->first_tracked = (uint16_t)t_current;
        b->status = 1;
        b->priority = pr;
    }

    lcb_type lcb_;
    int t_current;
    double alpha, delta;
    unsigned char pr_min_;
    double *eps_table;   // eps_table[dt] = sqrt(ln(2/delta) / (2 * dt))
    int max_windows;
};

#endif
