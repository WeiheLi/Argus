// TinyLFU.hpp -- TinyLFU-style admission (Einziger, Friedman, Manes, ToS'17)
// adapted to persistent-item detection.
//
// Layout (memory split frac_table : 1-frac_table):
//   * a d x w bucket array identical to Argus's (key, V, status; 13 bytes),
//     victim = smallest V among the d candidates (LFU),
//   * a TinyLFU frequency filter: count-min sketch with `rows` rows of 8-bit
//     saturating counters that count *distinct windows* (a per-window presence
//     bit per counter deduplicates arrivals inside a window, as On-Off does),
//     halved every `reset_windows` windows (TinyLFU's aging/reset).
// Admission (TinyLFU rule): on a collision, the newcomer replaces the LFU
// victim iff its filter estimate is larger than the victim's count.
#ifndef TINYLFU_H
#define TINYLFU_H
#include <vector>
#include <utility>
#include <cstring>
#include <cstdint>
#include <cstdlib>
#include "datatypes.hpp"
extern "C" {
#include "hash.h"
}
extern unsigned long g_hash_seed_offset;

class TinyLFU
{
#pragma pack(push, 1)
    struct Bucket { unsigned char key[LGN]; uint16_t V; uint8_t status; };
#pragma pack(pop)
public:
    TinyLFU(int memory_bytes, int lgn, int depth = 4, double frac_table = 0.75,
            int rows = 4, int reset_windows = 128, int bucket_bytes = 0)
        : depth_(depth), rows_(rows), reset_windows_(reset_windows)
    {
        keylen_ = lgn / 8;
        int table_bytes = (int)(memory_bytes * frac_table);
        if (bucket_bytes <= 0) bucket_bytes = (int)sizeof(Bucket);
        width_ = table_bytes / (bucket_bytes * depth);
        if (width_ < 1) width_ = 1;
        total_ = depth * width_;
        counts_ = new Bucket[total_];
        memset(counts_, 0, sizeof(Bucket) * (size_t)total_);
        int filter_bytes = memory_bytes - table_bytes;
        // each cell = 8-bit counter + 1 presence bit = 9 bits
        m_ = (int)((filter_bytes * 8LL) / (9LL * rows));
        if (m_ < 16) m_ = 16;
        cm_ = new uint8_t[(size_t)rows * m_];
        memset(cm_, 0, (size_t)rows * m_);
        words_ = ((size_t)rows * m_ + 63) / 64;
        seen_ = new uint64_t[words_];
        memset(seen_, 0, words_ * 8);
        hardner_ = new unsigned long[depth + rows];
        char name[] = "TinyLFU";
        unsigned long seed = g_hash_seed_offset +
            AwareHash((unsigned char *)name, strlen(name), 13091204281, 228204732751, 6620830889);
        for (int i = 0; i < depth + rows; i++) hardner_[i] = GenHashSeed(seed++);
        windows_ = 0;
    }
    ~TinyLFU() { delete[] counts_; delete[] cm_; delete[] seen_; delete[] hardner_; }

    static int BucketBytes() { return (int)sizeof(Bucket); }

    inline int estimate_and_record(const unsigned char *key)
    {
        int est = 255;
        for (int r = 0; r < rows_; r++) {
            size_t pos = (size_t)r * m_ + (size_t)(MurmurHash64A(key, keylen_, hardner_[depth_ + r]) % (uint64_t)m_);
            uint64_t mask = 1ULL << (pos & 63);
            if (!(seen_[pos >> 6] & mask)) { seen_[pos >> 6] |= mask; if (cm_[pos] < 255) cm_[pos]++; }
            if (cm_[pos] < est) est = cm_[pos];
        }
        return est;
    }

    void Update(const unsigned char *key)
    {
        int est = estimate_and_record(key);   // filter sees every arrival (doorkeeper-free TinyLFU)
        int idx[16];
        Bucket *b;
        for (int i = 0; i < depth_; i++) {
            int index = i * width_ + (int)(MurmurHash64A(key, keylen_, hardner_[i]) % (uint64_t)width_);
            idx[i] = index; b = &counts_[index];
            if (b->V == 0) { memcpy(b->key, key, keylen_); b->V = 1; b->status = 1; return; }
            if (memcmp(key, b->key, keylen_) == 0) { if (b->status == 0) { b->V++; b->status = 1; } return; }
        }
        int best = idx[0]; int minV = counts_[idx[0]].V;
        for (int i = 1; i < depth_; i++) if (counts_[idx[i]].V < minV) { minV = counts_[idx[i]].V; best = idx[i]; }
        b = &counts_[best];
        if (est > (int)b->V) { memcpy(b->key, key, keylen_); b->V = 1; b->status = 1; }
    }

    void NewWindow()
    {
        for (int i = 0; i < total_; i++) counts_[i].status = 0;
        memset(seen_, 0, words_ * 8);
        windows_++;
        if (reset_windows_ > 0 && windows_ % reset_windows_ == 0)
            for (size_t i = 0; i < (size_t)rows_ * m_; i++) cm_[i] >>= 1;
    }

    void Query(int v_thresh, std::vector<std::pair<key_tp, val_tp>> &results) const
    {
        for (int i = 0; i < total_; i++) {
            const Bucket *b = &counts_[i];
            if (b->V > 0 && (int)b->V >= v_thresh) { key_tp k; memcpy(k.key, b->key, keylen_); results.push_back({k, (val_tp)b->V}); }
        }
    }

private:
    int depth_, rows_, reset_windows_, keylen_, width_, total_, m_;
    Bucket *counts_;
    uint8_t *cm_;
    uint64_t *seen_; size_t words_;
    unsigned long *hardner_;
    int windows_;
};
#endif
