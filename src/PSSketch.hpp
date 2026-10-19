#ifndef PSSketch_H
#define PSSketch_H
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
 * PSSketch: Persistent + Sparse flow detection (KDD '25).
 *
 * Two-layer architecture:
 *   - Competition Layer (CL): depth x width grid of compact entries.
 *     Each entry stores fingerprint, an 8-bit freq counter f, a 6-bit
 *     persistence counter p, and 2 status bits W (window-touched) and
 *     OF (overflowed-and-promoted).
 *   - Protection Layer (PL): flat array of <ID, f_of, p_of>; tracks
 *     overflow counts for promoted (persistent) flows.
 *
 * Reported density: d = (f_of * 2^L_f + f) / (p_of * 2^L_p + p)
 *
 * Bucket layout (struct-based, no bit-packing inside cell):
 *   key            : LGN bytes
 *   f              : 1 byte   (uint8_t)
 *   p_w_of         : 1 byte   (p:6 | W:1 | OF:1)
 */
class PSSketch
{
#pragma pack(push, 1)
    struct CLEntry
    {
        unsigned char key[LGN]; // fingerprint (full LGN bytes; FP collisions handled implicitly)
        uint8_t f;              // 8-bit freq, overflow at 256
        uint8_t p_w_of;         // bit layout: PPPPPP W OF
    };

    struct PLEntry
    {
        unsigned char key[LGN];
        uint16_t f_of; // f overflow count
        uint16_t p_of; // p overflow count
        uint8_t used;  // 0 = empty, 1 = occupied
    };
#pragma pack(pop)

    // ---- bit ops on p_w_of ----
    static inline int  GetP (uint8_t b) { return (b >> 2) & 0x3F; }
    static inline bool GetW (uint8_t b) { return (b >> 1) & 0x1; }
    static inline bool GetOF(uint8_t b) { return b & 0x1; }
    static inline uint8_t MakePWOF(int p, bool w, bool of) {
        return (uint8_t)(((p & 0x3F) << 2) | ((w & 1) << 1) | (of & 1));
    }
    static const int F_MAX = 256;          // 2^8
    static const int P_MAX = 64;           // 2^6
    static const int F_BITS = 8;
    static const int P_BITS = 6;

    struct cl_type
    {
        val_tp sum;
        CLEntry *cells;
        int depth, width, lgn;
        unsigned long *hash, *hardner;
    };

public:
    PSSketch(int depth, int width, int lgn,
             int pl_records = 4096,
             int p_threshold = 50,
             double d_threshold = 1.2);
    ~PSSketch();

    void Update(unsigned char *key, val_tp value);
    void NewWindow();

    // Query interface (compatible with main_hitter framework):
    // returns flows that satisfy persistence >= p_threshold AND density < d_threshold.
    // val_tp in the pair is the estimated persistence (V).
    void QueryPS(std::vector<std::pair<key_tp, val_tp>> &results);

    // Detailed query: also returns frequency and density estimates.
    struct PSResult { key_tp key; val_tp V; val_tp F; double D; };
    void QueryPSDetailed(std::vector<PSResult> &results);

    val_tp PointQuery(unsigned char *key);
    val_tp GetCount();
    void Reset();
    int GetTCurrent() const { return t_current; }
    static int BucketBytes() { return (int)sizeof(CLEntry); }

private:
    // ---- helpers ----
    void initCell(int idx, unsigned char *key, int keylen);
    int  findPLEntry(unsigned char *key, int keylen);   // returns idx or -1
    int  reservePLEntry();                              // returns idx or -1 (full)
    void clearPLEntry(int idx);
    bool reportToPL(unsigned char *key, int keylen, char mode);
    // mode = 'o' first promotion, 'f' f overflow, 'p' p overflow.
    // returns true on success, false if PL refuses (e.g. full / pruned).

    cl_type cl_;
    PLEntry *pl_;
    int pl_capacity_;

    int p_thr_;
    double d_thr_;
    int t_current;

    // disable copy
    PSSketch(const PSSketch &);
    PSSketch &operator=(const PSSketch &);
};

#endif
