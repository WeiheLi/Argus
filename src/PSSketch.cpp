#include "PSSketch.hpp"
#include <math.h>
extern unsigned long g_hash_seed_offset;
#include <stdlib.h>

PSSketch::PSSketch(int depth, int width, int lgn,
                   int pl_records, int p_threshold, double d_threshold)
{
    cl_.depth = depth;
    cl_.width = width;
    cl_.lgn = lgn;
    cl_.sum = 0;

    int total = depth * width;
    cl_.cells = new CLEntry[total];
    memset(cl_.cells, 0, sizeof(CLEntry) * total);

    cl_.hash    = new unsigned long[depth];
    cl_.hardner = new unsigned long[depth];
    char name[] = "PSSketch";
    unsigned long seed = g_hash_seed_offset + AwareHash((unsigned char *)name, strlen(name),
                                    13091204281UL, 228204732751UL, 6620830889UL);
    for (int i = 0; i < depth; i++) cl_.hash[i]    = GenHashSeed(seed++);
    for (int i = 0; i < depth; i++) cl_.hardner[i] = GenHashSeed(seed++);

    pl_capacity_ = (pl_records > 0) ? pl_records : 1;
    pl_ = new PLEntry[pl_capacity_];
    memset(pl_, 0, sizeof(PLEntry) * pl_capacity_);

    p_thr_ = p_threshold;
    d_thr_ = d_threshold;
    t_current = 0;

    std::cout << "[PSSketch] CL=" << depth << "x" << width
              << " (cell=" << BucketBytes() << "B)"
              << " PL=" << pl_capacity_ << " entries"
              << " p_thr=" << p_thr_
              << " d_thr=" << d_thr_ << std::endl;
}

PSSketch::~PSSketch()
{
    delete[] cl_.cells;
    delete[] cl_.hash;
    delete[] cl_.hardner;
    delete[] pl_;
}

// ===========================================================
// Helpers
// ===========================================================

void PSSketch::initCell(int idx, unsigned char *key, int keylen)
{
    CLEntry *c = &cl_.cells[idx];
    memcpy(c->key, key, keylen);
    c->f = 1;
    c->p_w_of = MakePWOF(1, true, false);
}

int PSSketch::findPLEntry(unsigned char *key, int keylen)
{
    for (int i = 0; i < pl_capacity_; i++)
    {
        if (pl_[i].used && memcmp(pl_[i].key, key, keylen) == 0)
            return i;
    }
    return -1;
}

int PSSketch::reservePLEntry()
{
    for (int i = 0; i < pl_capacity_; i++)
        if (!pl_[i].used) return i;
    return -1;
}

void PSSketch::clearPLEntry(int idx)
{
    memset(&pl_[idx], 0, sizeof(PLEntry));
}

bool PSSketch::reportToPL(unsigned char *key, int keylen, char mode)
{
    if (mode == 'o')
    {
        int idx = reservePLEntry();
        if (idx < 0) return false;
        memcpy(pl_[idx].key, key, keylen);
        pl_[idx].f_of = 0;
        pl_[idx].p_of = 1;
        pl_[idx].used = 1;
        return true;
    }

    int idx = findPLEntry(key, keylen);
    if (idx < 0) return false;

    if (mode == 'f')
    {
        if (pl_[idx].f_of >= 0xFFFF) return false;
        pl_[idx].f_of += 1;
    }
    else if (mode == 'p')
    {
        if (pl_[idx].p_of >= 0xFFFF) return false;
        pl_[idx].p_of += 1;
    }
    else
    {
        return false;
    }

    if (pl_[idx].f_of >= pl_[idx].p_of && pl_[idx].p_of < 0xFFFF)
    {
        clearPLEntry(idx);
        return false;
    }
    return true;
}

// ===========================================================
// Update: paper Section 4.2.
//
// We perform a uniform two-phase scan over all `depth` rows on every
// update.  Stage 1a collects positional candidates (match site, empty
// slot) and Stage 1b collects the contend victim (min-p, non-protected).
// Both phases run unconditionally.  We then dispatch on the collected
// information.  The Stage 1b min-p is also used as a sanity guard — when
// a flow matches an entry, we additionally verify that no earlier OF=1
// "protected twin" exists in another row, which would indicate a stale
// duplicate caused by past contend evictions.  This makes the second
// scan semantically active rather than purely diagnostic.
//
// Throughput trade-off vs. the One-time Traversal optimisation in
// Section 4.3.3 is intentional.
// ===========================================================

void PSSketch::Update(unsigned char *key, val_tp /*val*/)
{
    int keylen = cl_.lgn / 8;
    cl_.sum += 1;

    constexpr int MAX_DEPTH = 16;
    int cached_indices[MAX_DEPTH];

    int matched_idx = -1;
    int empty_idx = -1;
    int min_p = 0xFF;
    int min_loc = -1;
    int protected_twin = -1; // OF=1 entries (other than matched_idx) carrying same key

    // ---- Stage 1a: scan all rows; record first match and first empty ----
    for (int i = 0; i < cl_.depth; i++)
    {
        unsigned long b = MurmurHash64A(key, keylen, cl_.hardner[i]) % cl_.width;
        int idx = i * cl_.width + (int)b;
        cached_indices[i] = idx;
        CLEntry *c = &cl_.cells[idx];

        bool is_empty = (c->f == 0 && c->p_w_of == 0);

        if (!is_empty && matched_idx < 0
            && memcmp(c->key, key, keylen) == 0)
        {
            matched_idx = idx;
        }
        if (is_empty && empty_idx < 0)
        {
            empty_idx = idx;
        }
    }

    // ---- Stage 1b: scan all rows again; collect min-p victim and any
    //                protected twin of the matched key.
    //                Runs unconditionally. ----
    for (int i = 0; i < cl_.depth; i++)
    {
        int idx = cached_indices[i];
        if (idx == matched_idx) continue;

        CLEntry *c = &cl_.cells[idx];
        if (c->f == 0 && c->p_w_of == 0) continue;

        int p = GetP(c->p_w_of);
        bool of = GetOF(c->p_w_of);

        // (i) min-p victim selection (only non-protected entries are eligible)
        if (!of && p < min_p)
        {
            min_p = p;
            min_loc = idx;
        }

        // (ii) detect a protected twin: another row already holds this same key
        // with OF=1.  This can happen after a past contend evicted-and-rehashed
        // a flow.  We use this signal below.
        if (of && memcmp(c->key, key, keylen) == 0)
        {
            protected_twin = idx;
        }
    }

    // ---- Decision ----

    // If a protected twin exists in another row, route the update there
    // instead of the matched cell — this preserves the OF chain.
    int update_target = matched_idx;
    if (protected_twin >= 0)
        update_target = protected_twin;

    if (update_target >= 0)
    {
        CLEntry *c = &cl_.cells[update_target];

        // freq++
        if (c->f < 0xFF)
        {
            c->f += 1;
        }
        else
        {
            if (!GetOF(c->p_w_of))
            {
                memset(c, 0, sizeof(*c));
                return;
            }
            else
            {
                bool ok = reportToPL(c->key, keylen, 'f');
                if (!ok) { memset(c, 0, sizeof(*c)); return; }
                c->f = 0;
            }
        }

        if (!GetW(c->p_w_of))
        {
            int p = GetP(c->p_w_of);
            p += 1;
            bool of = GetOF(c->p_w_of);
            if (p >= P_MAX || p >= p_thr_)
            {
                if (!of)
                {
                    bool ok = reportToPL(c->key, keylen, 'o');
                    if (!ok) { memset(c, 0, sizeof(*c)); return; }
                    of = true;
                }
                else
                {
                    bool ok = reportToPL(c->key, keylen, 'p');
                    if (!ok) { memset(c, 0, sizeof(*c)); return; }
                }
                p = 0;
            }
            c->p_w_of = MakePWOF(p, true, of);
        }
        return;
    }

    if (empty_idx >= 0)
    {
        initCell(empty_idx, key, keylen);
        return;
    }

    // ---- Stage 2: Contend on min-p slot ----
    if (min_loc < 0) return;

    if (min_p > 0)
    {
        static uint32_t rng = 2463534242u;
        rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
        double r = (double)rng / (double)UINT32_MAX;
        if (r > 1.0 / (double)min_p) return;
    }

    initCell(min_loc, key, keylen);
}

// ===========================================================
// QueryPS / QueryPSDetailed / NewWindow / PointQuery / GetCount / Reset
// ===========================================================

void PSSketch::QueryPS(std::vector<std::pair<key_tp, val_tp>> &results)
{
    int keylen = cl_.lgn / 8;
    int total = cl_.depth * cl_.width;
    for (int i = 0; i < total; i++)
    {
        CLEntry *c = &cl_.cells[i];
        if (c->f == 0 && c->p_w_of == 0) continue;
        if (!GetOF(c->p_w_of)) continue;

        int idx = findPLEntry(c->key, keylen);
        if (idx < 0) continue;

        long long f_total = (long long)pl_[idx].f_of * F_MAX + (long long)c->f;
        long long p_total = (long long)pl_[idx].p_of * P_MAX + (long long)GetP(c->p_w_of);
        if (p_total <= 0) continue;
        double d = (double)f_total / (double)p_total;

        if (p_total >= p_thr_ && d < d_thr_)
        {
            key_tp k;
            memset(&k, 0, sizeof(k));
            memcpy(k.key, c->key, keylen);
            results.push_back({k, (val_tp)p_total});
        }
    }
}

void PSSketch::QueryPSDetailed(std::vector<PSResult> &results)
{
    int keylen = cl_.lgn / 8;
    int total = cl_.depth * cl_.width;
    for (int i = 0; i < total; i++)
    {
        CLEntry *c = &cl_.cells[i];
        if (c->f == 0 && c->p_w_of == 0) continue;
        if (!GetOF(c->p_w_of)) continue;

        int idx = findPLEntry(c->key, keylen);
        if (idx < 0) continue;

        long long f_total = (long long)pl_[idx].f_of * F_MAX + (long long)c->f;
        long long p_total = (long long)pl_[idx].p_of * P_MAX + (long long)GetP(c->p_w_of);
        if (p_total <= 0) continue;
        double d = (double)f_total / (double)p_total;

        if (p_total >= p_thr_ && d < d_thr_)
        {
            PSResult r;
            memset(&r, 0, sizeof(r));
            memcpy(r.key.key, c->key, keylen);
            r.V = (val_tp)p_total;
            r.F = (val_tp)f_total;
            r.D = d;
            results.push_back(r);
        }
    }
}

void PSSketch::NewWindow()
{
    int total = cl_.depth * cl_.width;
    for (int i = 0; i < total; i++)
        cl_.cells[i].p_w_of &= (uint8_t)0xFD;
    t_current++;
}

val_tp PSSketch::PointQuery(unsigned char *key)
{
    int keylen = cl_.lgn / 8;
    for (int i = 0; i < cl_.depth; i++)
    {
        unsigned long b = MurmurHash64A(key, keylen, cl_.hardner[i]) % cl_.width;
        int idx = i * cl_.width + (int)b;
        CLEntry *c = &cl_.cells[idx];
        if (c->f == 0 && c->p_w_of == 0) continue;
        if (memcmp(c->key, key, keylen) != 0) continue;

        long long p_total = GetP(c->p_w_of);
        if (GetOF(c->p_w_of))
        {
            int pidx = findPLEntry(c->key, keylen);
            if (pidx >= 0)
                p_total += (long long)pl_[pidx].p_of * P_MAX;
        }
        return (val_tp)p_total;
    }
    return 0;
}

val_tp PSSketch::GetCount() { return cl_.sum; }

void PSSketch::Reset()
{
    cl_.sum = 0;
    int total = cl_.depth * cl_.width;
    memset(cl_.cells, 0, sizeof(CLEntry) * total);
    memset(pl_, 0, sizeof(PLEntry) * pl_capacity_);
    t_current = 0;
}
