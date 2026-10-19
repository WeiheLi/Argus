// PAComposite.hpp -- parallel-composition baseline for persistent-and-priority
// detection: a priority-aware layer in the spirit of PA-Sketch (Huang et al.)
// run side by side with an unmodified persistence sketch (Pandora or
// Stable-Sketch); the output is the intersection of the two.
//
// Priority layer: an array of cells (key, priority, count). An arrival with
// priority p probes p+1 cells (a higher priority hashes into more cells and is
// therefore harder to displace); it takes an empty cell, or displaces the
// probed cell with the lowest (priority, count) if that priority is strictly
// lower than p. The layer answers "is u present with priority >= pi_min".
#ifndef PACOMPOSITE_H
#define PACOMPOSITE_H
#include <vector>
#include <utility>
#include <cstring>
#include <cstdint>
#include "datatypes.hpp"
#include "Pandora.hpp"
#include "StableSketch.hpp"
extern "C" {
#include "hash.h"
}
extern unsigned long g_hash_seed_offset;

class PriorityLayer
{
#pragma pack(push, 1)
    struct Cell { unsigned char key[LGN]; uint8_t prio; uint16_t count; }; // 11 bytes
#pragma pack(pop)
public:
    PriorityLayer(int memory_bytes, int lgn, int max_prio = 8, int cell_bytes = 0)
    {
        keylen_ = lgn / 8;
        if (cell_bytes <= 0) cell_bytes = (int)sizeof(Cell);
        n_ = memory_bytes / cell_bytes; if (n_ < 8) n_ = 8;
        cells_ = new Cell[n_]; memset(cells_, 0, sizeof(Cell) * (size_t)n_);
        max_prio_ = max_prio;
        seeds_ = new unsigned long[max_prio + 1];
        char name[] = "PALayer";
        unsigned long seed = g_hash_seed_offset +
            AwareHash((unsigned char *)name, strlen(name), 13091204281, 228204732751, 6620830889);
        for (int i = 0; i <= max_prio; i++) seeds_[i] = GenHashSeed(seed++);
    }
    ~PriorityLayer() { delete[] cells_; delete[] seeds_; }
    void Update(const unsigned char *key, unsigned char prio)
    {
        int probes = (int)prio + 1; if (probes > max_prio_ + 1) probes = max_prio_ + 1;
        int weakest = -1; int wp = 1 << 30; int wc = 1 << 30;
        for (int i = 0; i < probes; i++) {
            int idx = (int)(MurmurHash64A(key, keylen_, seeds_[i]) % (uint64_t)n_);
            Cell *c = &cells_[idx];
            if (c->count == 0) { memcpy(c->key, key, keylen_); c->prio = prio; c->count = 1; return; }
            if (memcmp(c->key, key, keylen_) == 0) { if (c->count < 65535) c->count++; if (prio > c->prio) c->prio = prio; return; }
            if (c->prio < wp || (c->prio == wp && c->count < wc)) { wp = c->prio; wc = c->count; weakest = idx; }
        }
        if (weakest >= 0 && wp < (int)prio) { Cell *c = &cells_[weakest]; memcpy(c->key, key, keylen_); c->prio = prio; c->count = 1; }
    }
    int Priority(const unsigned char *key) const // -1 if absent
    {
        for (int i = 0; i <= max_prio_; i++) {
            int idx = (int)(MurmurHash64A(key, keylen_, seeds_[i]) % (uint64_t)n_);
            const Cell *c = &cells_[idx];
            if (c->count > 0 && memcmp(c->key, key, keylen_) == 0) return c->prio;
        }
        return -1;
    }
private:
    int keylen_, n_, max_prio_;
    Cell *cells_;
    unsigned long *seeds_;
};

template <class PersistSketch>
class PAComposite
{
public:
    // frac_pa of the budget goes to the priority layer, the rest to the persistence sketch
    PAComposite(int memory_bytes, int lgn, int depth, double frac_pa, int pers_bucket_bytes, int pa_cell_bytes = 0)
        : pa_((int)(memory_bytes * frac_pa), lgn, 8, pa_cell_bytes)
    {
        int pers_bytes = memory_bytes - (int)(memory_bytes * frac_pa);
        int width = pers_bytes / (pers_bucket_bytes * depth); if (width < 1) width = 1;
        pers_ = new PersistSketch(depth, width, lgn);
    }
    ~PAComposite() { delete pers_; }
    PriorityLayer pa_;
    PersistSketch *pers_;
};
#endif
