#ifndef TightSketchPH_H
#define TightSketchPH_H
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
 * Tight-Sketch (CIKM '23) for Persistent + Heavy Item Detection.
 *
 * Task (aligned with ArgusPH framework):
 *   V(e) / T   >= alpha   AND   C(e) / V(e)   >= kappa
 *   where V = persistence (# windows seen),
 *         C = frequency   (total occurrences).
 *
 * Per-bucket state (paper Section 3.4.4 — "significant item lookup"
 * data structure with double counters; we re-use it here because
 * the V+C double-counter design is exactly what's needed for the
 * V/T + C/V dual-ratio query):
 *   key  : item key
 *   fc   : frequency counter (= C)
 *   fa   : sustained arrival strength for frequency
 *   pc   : persistence counter (= V)
 *   pa   : sustained arrival strength for persistence,
 *          high bit doubles as the W (window-touched) flag.
 *
 * Update (paper Algorithm 1 + Section 3.4.4):
 *   Stage I:  probe rows; on match bump fc+fa unconditionally,
 *             bump pc+pa only if W=false (then set W=true).  On
 *             every non-matching row, decrement that row's fa & pa.
 *             Track victim = entry with smallest "significance"
 *             G = fc + pc within the probed rows.
 *   Stage II: probabilistic decay on the victim:
 *               G < M  ->  decay prob = 1/(G + 1)
 *               G >= M ->  decay prob = 1/(G * a_min + 1)
 *             A successful decay decrements both fc and pc by 1.
 *             When pc reaches 0 the victim is replaced (paper:
 *             "since persistence <= frequency, once pc decays to 0
 *             the new item can successfully replace").
 *
 * Bucket layout (struct-based; 16 bytes for LGN=8):
 *   key  : LGN bytes
 *   fc   : uint16_t
 *   fa   : uint16_t
 *   pc   : uint16_t
 *   pa   : uint16_t  (high bit = W flag, low 15 bits = strength)
 */
class TightSketchPH
{
    static const uint16_t W_MASK   = 0x8000;
    static const uint16_t PA_MASK  = 0x7FFF;
    static const int      CAP16    = 0x7FFE;
    static const int      M_DEFAULT = 10;

#pragma pack(push, 1)
    struct Bucket
    {
        unsigned char key[LGN]; // 8
        uint16_t fc;            // 2
        uint16_t fa;            // 2
        uint16_t pc;            // 2
        uint16_t pa;            // 2  (high bit = W)
    };                          // 16 bytes
#pragma pack(pop)

    struct ts_type
    {
        val_tp sum;
        Bucket *counts;
        int depth, width, lgn;
        unsigned long *hash, *hardner;
    };

public:
    // alpha   : persistence ratio threshold (V/T >= alpha)
    // kappa   : heavy ratio threshold       (C/V >= kappa)
    // M_thresh: paper M, decay-strategy switch
    TightSketchPH(int depth, int width, int lgn,
                  double alpha = 0.3, double kappa = 10.0,
                  int n_windows = 1600,
                  int M_thresh = M_DEFAULT);
    ~TightSketchPH();

    void Update(unsigned char *key, val_tp value);
    void NewWindow();

    // Query interface (compatible with LCB-PH / main.cpp framework):
    // returns items with V >= alpha*T AND C/V >= kappa.
    // val_tp in the pair carries V (persistence estimate).
    void QueryPH(std::vector<std::pair<key_tp, val_tp>> &results);

    struct PHResult { key_tp key; val_tp V; uint16_t C; };
    void QueryPHDetailed(std::vector<PHResult> &results);

    val_tp PointQuery(unsigned char *key);   // returns V (persistence)
    val_tp GetCount();
    void Reset();
    int GetTCurrent() const { return t_current; }
    static int BucketBytes() { return (int)sizeof(Bucket); }

private:
    static inline int  GetPA(uint16_t pa) { return (int)(pa & PA_MASK); }
    static inline bool GetW (uint16_t pa) { return (pa & W_MASK) != 0; }
    static inline uint16_t MakePA(int strength, bool w) {
        if (strength < 0) strength = 0;
        if (strength > (int)PA_MASK) strength = (int)PA_MASK;
        return (uint16_t)((w ? W_MASK : 0) | (strength & PA_MASK));
    }

    inline int significanceOf(const Bucket *b) const {
        return (int)b->fc + (int)b->pc;
    }

    inline void initBucket(int idx, unsigned char *key, int keylen) {
        Bucket *b = &ts_.counts[idx];
        memcpy(b->key, key, keylen);
        b->fc = 1;
        b->fa = 1;
        b->pc = 1;
        b->pa = MakePA(1, /*W=*/true);
    }

    ts_type ts_;
    int t_current;
    double alpha, kappa;
    int M_thresh;

    // disable copy
    TightSketchPH(const TightSketchPH &);
    TightSketchPH &operator=(const TightSketchPH &);
};

#endif
