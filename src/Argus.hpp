// Argus.hpp -- ARGUS, persistent-item detection on a full key (Algorithm 1).
//
// A bucket holds its occupant's key, the number of windows V in which that
// occupant has appeared, the window first_tracked in which it was installed,
// and a status bit set by an arrival in the current window. From V and the
// number of windows tracked the sketch forms a Hoeffding confidence interval
// for the occupant's persistence rate, and settles a collision by a three-way
// verdict over that interval: protect, evict or decay, or defer.
//
// ArgusPacked, included at the foot of this file, is the same algorithm on an
// 8-byte bucket; it is what runs when items are identified by a fingerprint.
#ifndef ARGUS_UNIFIED_H
#define ARGUS_UNIFIED_H
#include <vector>
#include <utility>
#include <cstring>
#include <cmath>
#include <cstdint>
#include "datatypes.hpp"
extern "C" {
#include "hash.h"
}
extern unsigned long g_hash_seed_offset;

// per-stage counters, reported in the output CSV
struct ArgusStats {
    uint64_t n_update = 0;
    uint64_t s1_empty = 0, s1_match = 0, s1_match_inc = 0;
    uint64_t s2_collision = 0;      // reached Stage 2
    uint64_t s3_drop = 0;           // victim active in the current window
    uint64_t protect = 0, defer = 0;
    uint64_t evict_ucb = 0;         // replacement via UCB < alpha at V=1
    uint64_t evict_override = 0;    // replacement via the V=1 override
    uint64_t decay = 0;             // V decremented (V>1)
};

class ArgusSketch
{
public:
#pragma pack(push, 1)
    struct Bucket {
        unsigned char key[LGN];
        uint16_t V;
        uint16_t first_tracked;
        uint8_t status;
    };
#pragma pack(pop)

    ArgusStats stats;

    ArgusSketch(int depth, int width, int lgn, double alpha, double delta, int n_windows)
        : depth_(depth), width_(width), alpha_(alpha)
    {
        keylen_ = lgn / 8;
        total_ = depth * width;
        counts_ = new Bucket[total_];
        memset(counts_, 0, sizeof(Bucket) * (size_t)total_);
        hardner_ = new unsigned long[depth];
        char name[] = "LCBSketch";
        unsigned long seed = g_hash_seed_offset +
            AwareHash((unsigned char *)name, strlen(name), 13091204281, 228204732751, 6620830889);
        for (int i = 0; i < depth; i++) GenHashSeed(seed++);            // hash[] (unused, keep order)
        for (int i = 0; i < depth; i++) GenHashSeed(seed++);            // scale[] (unused, keep order)
        for (int i = 0; i < depth; i++) hardner_[i] = GenHashSeed(seed++);
        t_current_ = 0;
        max_windows_ = n_windows + 2;
        // the interval half-width at T windows is c_delta / sqrt(T), tabulated once
        double c_delta = sqrt(log(2.0 / delta) / 2.0);
        eps_table_ = new double[max_windows_ + 2];
        eps_table_[0] = c_delta;
        for (int i = 1; i <= max_windows_ + 1; i++) eps_table_[i] = c_delta / sqrt((double)i);
    }
    ~ArgusSketch() { delete[] counts_; delete[] hardner_; delete[] eps_table_; }

    static int BucketBytes() { return (int)sizeof(Bucket); }
    int GetTCurrent() const { return t_current_; }

    // ---- bounds ---------------------------------------------------------
    // lo/hi are clipped to [0,1]; lo_u/hi_u are the raw values.
    inline void Bounds(int V, int dt, double &lo, double &hi, double &lo_u, double &hi_u) const
    {
        double phat = (double)V / (double)dt, e = eps_table_[dt];
        lo_u = phat - e; hi_u = phat + e;
        lo = lo_u < 0 ? 0 : lo_u; hi = hi_u > 1 ? 1 : hi_u;
    }

    // ---- Update (Algorithm 1) ------------------------------------------
    void Update(const unsigned char *key)
    {
        stats.n_update++;
        int idx[16];
        Bucket *b = nullptr;
        // ---- Stage 1: an empty bucket or the occupant itself ----
        for (int i = 0; i < depth_; i++) {
            unsigned long bucket = MurmurHash64A(key, keylen_, hardner_[i]) % width_;
            int index = i * width_ + (int)bucket;
            idx[i] = index;
            b = &counts_[index];
            if (b->V == 0) { init_bucket(index, key); stats.s1_empty++; return; }
            if (memcmp(key, b->key, keylen_) == 0) {
                stats.s1_match++;
                if (b->status == 0) { b->V += 1; b->status = 1; stats.s1_match_inc++; }
                return;
            }
        }
        stats.s2_collision++;

        // ---- Stage 2: the victim is the candidate with the lowest LCB ----
        int best = -1;
        double min_score = 1e18;
        for (int i = 0; i < depth_; i++) {
            b = &counts_[idx[i]];
            int dt = tracked(b);
            double lo, hi, lo_u, hi_u; Bounds(b->V, dt, lo, hi, lo_u, hi_u);
            if (lo_u < min_score) { min_score = lo_u; best = idx[i]; }
        }
        b = &counts_[best];

        // ---- Stage 3: within-window protection ----
        if (b->status == 1) { stats.s3_drop++; return; }

        // ---- Stage 4: the three-way verdict ----
        int v = b->V;
        int dt = tracked(b);
        double phat = (double)v / dt;
        double lo, hi, lo_u, hi_u; Bounds(v, dt, lo, hi, lo_u, hi_u);

        if (lo_u >= alpha_) { stats.protect++; return; }                 // protect
        if (hi < alpha_) {                                               // evict / decay via UCB
            if (v <= 1) { stats.evict_ucb++; init_bucket(best, key); }
            else { stats.decay++; b->V -= 1; }
            return;
        }
        if (v == 1 && phat < alpha_) {                                   // V=1 override
            stats.evict_override++; init_bucket(best, key); return;
        }
        stats.defer++;                                                   // defer
    }

    void NewWindow()
    {
        for (int i = 0; i < total_; i++) counts_[i].status = 0;
        t_current_++;
    }

    void QueryPersistent(int v_thresh, std::vector<std::pair<key_tp, val_tp>> &results) const
    {
        for (int i = 0; i < total_; i++) {
            const Bucket *b = &counts_[i];
            if (b->V > 0 && (int)b->V >= v_thresh) {
                key_tp k; memcpy(k.key, b->key, keylen_);
                results.push_back({k, (val_tp)b->V});
            }
        }
    }

    val_tp PointQuery(const unsigned char *key) const
    {
        for (int i = 0; i < depth_; i++) {
            unsigned long bucket = MurmurHash64A(key, keylen_, hardner_[i]) % width_;
            const Bucket *b = &counts_[i * width_ + (int)bucket];
            if (b->V > 0 && memcmp(b->key, key, keylen_) == 0) return b->V;
        }
        return 0;
    }

private:
    // the number of windows the occupant has been tracked for
    inline int tracked(const Bucket *b) const
    {
        int dt = t_current_ - (int)b->first_tracked + 1;
        if (dt < 1) dt = 1;
        if (dt > max_windows_) dt = max_windows_;
        return dt;
    }
    inline void init_bucket(int idx, const unsigned char *key)
    {
        Bucket *b = &counts_[idx];
        memcpy(b->key, key, keylen_);
        b->V = 1; b->first_tracked = (uint16_t)t_current_; b->status = 1;
    }

    int depth_, width_, keylen_, total_;
    double alpha_;
    Bucket *counts_;
    unsigned long *hardner_;
    int t_current_, max_windows_;
    double *eps_table_;
};

#include "ArgusPacked.inc"

#endif
