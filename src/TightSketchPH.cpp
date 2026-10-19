#include "TightSketchPH.hpp"
#include <math.h>
extern unsigned long g_hash_seed_offset;
#include <stdlib.h>

TightSketchPH::TightSketchPH(int depth, int width, int lgn,
                             double alpha_val, double kappa_val,
                             int /*n_windows*/, int M_thresh_val)
    : alpha(alpha_val), kappa(kappa_val), M_thresh(M_thresh_val)
{
    ts_.depth = depth;
    ts_.width = width;
    ts_.lgn = lgn;
    ts_.sum = 0;

    int total_buckets = depth * width;
    ts_.counts = new Bucket[total_buckets];
    memset(ts_.counts, 0, sizeof(Bucket) * total_buckets);

    ts_.hash    = new unsigned long[depth];
    ts_.hardner = new unsigned long[depth];
    char name[] = "TightSketchPH";
    unsigned long seed = g_hash_seed_offset + AwareHash((unsigned char *)name, strlen(name),
                                    13091204281UL, 228204732751UL, 6620830889UL);
    for (int i = 0; i < depth; i++) ts_.hash[i]    = GenHashSeed(seed++);
    for (int i = 0; i < depth; i++) ts_.hardner[i] = GenHashSeed(seed++);

    t_current = 0;

    std::cout << "[TightSketchPH] " << depth << "x" << width
              << " (cell=" << BucketBytes() << "B)"
              << " alpha=" << alpha
              << " kappa=" << kappa
              << " M=" << M_thresh << std::endl;
}

TightSketchPH::~TightSketchPH()
{
    delete[] ts_.counts;
    delete[] ts_.hash;
    delete[] ts_.hardner;
}

// ===========================================================
// Update: paper Section 3.4.4.
//   Stage I:  probe each row. On match: bump fc+fa, bump pc+pa
//             only if W=false (then W=true). On non-match: pay
//             arrival-strength penalty (fa--, pa-- floor 0) and
//             track victim = entry with smallest significance G.
//   Stage II: probabilistic decay on victim. Decrementing both
//             fc & pc together; pc reaching 0 triggers replacement.
// ===========================================================

void TightSketchPH::Update(unsigned char *key, val_tp /*val*/)
{
    int keylen = ts_.lgn / 8;
    ts_.sum += 1;

    constexpr int MAX_DEPTH = 16;
    int cached_indices[MAX_DEPTH];

    int min_G_loc = -1;
    int min_G = INT32_MAX;

    // ---- Stage I ----
    for (int i = 0; i < ts_.depth; i++)
    {
        unsigned long h = MurmurHash64A(key, keylen, ts_.hardner[i]) % ts_.width;
        int idx = i * ts_.width + (int)h;
        cached_indices[i] = idx;
        Bucket *b = &ts_.counts[idx];

        bool is_empty = (b->fc == 0 && b->pc == 0 && b->fa == 0 && b->pa == 0);

        // (a) empty -> install
        if (is_empty)
        {
            initBucket(idx, key, keylen);
            return;
        }

        // (b) match -> update
        if (memcmp(b->key, key, keylen) == 0)
        {
            // freq always +1
            if (b->fc < CAP16) b->fc += 1;
            if (b->fa < CAP16) b->fa += 1;

            // persistence +1 only first hit per window
            if (!GetW(b->pa))
            {
                if (b->pc < CAP16) b->pc += 1;
                int pa_strength = GetPA(b->pa);
                if (pa_strength < (int)PA_MASK) pa_strength += 1;
                b->pa = MakePA(pa_strength, /*W=*/true);
            }
            return;
        }

        // (c) miss on this row: arrival-strength penalty, floor 0.
        // Paper Algorithm 1 line 11 — every conflicting row pays.
        if (b->fa > 0) b->fa -= 1;
        int pa_strength = GetPA(b->pa);
        if (pa_strength > 0) pa_strength -= 1;
        b->pa = MakePA(pa_strength, GetW(b->pa));

        // track victim by smallest G = fc + pc
        int G = significanceOf(b);
        if (G < min_G)
        {
            min_G = G;
            min_G_loc = idx;
        }
    }

    // ---- Stage II ----
    if (min_G_loc < 0) return;
    Bucket *v = &ts_.counts[min_G_loc];

    // Window-level protection: if victim was already updated this
    // window (W=true), don't aggressively decay it.
    if (GetW(v->pa)) return;

    // decay probability
    int a_min;
    {
        int pa_strength = GetPA(v->pa);
        a_min = (v->fa < pa_strength) ? (int)v->fa : pa_strength;
        if (a_min < 1) a_min = 1;
    }
    int G_v = significanceOf(v);
    if (G_v < 1) G_v = 1;

    double prob;
    if (G_v < M_thresh)
        prob = 1.0 / ((double)G_v + 1.0);
    else
        prob = 1.0 / ((double)G_v * (double)a_min + 1.0);

    static uint32_t rng = 0x9E3779B9u;
    rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
    double r = (double)rng / (double)UINT32_MAX;

    if (r < prob)
    {
        if (v->fc > 0) v->fc -= 1;
        if (v->pc > 0) v->pc -= 1;

        if (v->pc == 0)
        {
            initBucket(min_G_loc, key, keylen);
            return;
        }
    }
    // else: discard incoming.
}

// ===========================================================
// Query: items with V >= alpha*T AND C/V >= kappa.
// (Same predicate as ArgusPH::QueryPH — V <- pc, C <- fc.)
// ===========================================================

void TightSketchPH::QueryPH(std::vector<std::pair<key_tp, val_tp>> &results)
{
    int T = t_current + 1;
    int v_thresh = (int)(alpha * T);
    int total = ts_.width * ts_.depth;
    int keylen = ts_.lgn / 8;
    for (int i = 0; i < total; i++)
    {
        Bucket *b = &ts_.counts[i];
        if (b->pc == 0 || (int)b->pc < v_thresh) continue;
        if ((double)b->fc < kappa * (double)b->pc) continue;

        key_tp k;
        memset(&k, 0, sizeof(k));
        memcpy(k.key, b->key, keylen);
        results.push_back({k, (val_tp)b->pc});
    }
}

void TightSketchPH::QueryPHDetailed(std::vector<PHResult> &results)
{
    int T = t_current + 1;
    int v_thresh = (int)(alpha * T);
    int total = ts_.width * ts_.depth;
    int keylen = ts_.lgn / 8;
    for (int i = 0; i < total; i++)
    {
        Bucket *b = &ts_.counts[i];
        if (b->pc == 0 || (int)b->pc < v_thresh) continue;
        if ((double)b->fc < kappa * (double)b->pc) continue;

        PHResult r;
        memset(&r, 0, sizeof(r));
        memcpy(r.key.key, b->key, keylen);
        r.V = (val_tp)b->pc;
        r.C = b->fc;
        results.push_back(r);
    }
}

// ===========================================================
// NewWindow: clear W flag in every entry.
// ===========================================================

void TightSketchPH::NewWindow()
{
    int total = ts_.depth * ts_.width;
    for (int i = 0; i < total; i++)
        ts_.counts[i].pa &= PA_MASK;
    t_current++;
}

val_tp TightSketchPH::PointQuery(unsigned char *key)
{
    int keylen = ts_.lgn / 8;
    for (int i = 0; i < ts_.depth; i++)
    {
        unsigned long h = MurmurHash64A(key, keylen, ts_.hardner[i]) % ts_.width;
        int idx = i * ts_.width + (int)h;
        Bucket *b = &ts_.counts[idx];
        if (b->pc == 0) continue;
        if (memcmp(b->key, key, keylen) == 0)
            return (val_tp)b->pc;
    }
    return 0;
}

val_tp TightSketchPH::GetCount() { return ts_.sum; }

void TightSketchPH::Reset()
{
    ts_.sum = 0;
    int total = ts_.depth * ts_.width;
    memset(ts_.counts, 0, sizeof(Bucket) * total);
    t_current = 0;
}
