#include "ArgusPP.hpp"
#include <math.h>
extern unsigned long g_hash_seed_offset;

ArgusPP::ArgusPP(int depth, int width, int lgn,
                         double alpha_val, unsigned char pr_min,
                         double delta_val, int n_windows)
    : alpha(alpha_val), delta(delta_val), pr_min_(pr_min)
{
    lcb_.depth = depth;
    lcb_.width = width;
    lcb_.lgn = lgn;
    lcb_.sum = 0;

    int total_buckets = depth * width;
    lcb_.counts = new Bucket[total_buckets];
    memset(lcb_.counts, 0, sizeof(Bucket) * total_buckets);

    lcb_.hash = new unsigned long[depth];
    lcb_.scale = new unsigned long[depth];
    lcb_.hardner = new unsigned long[depth];
    char name[] = "LCBSketchPP";
    unsigned long seed = g_hash_seed_offset + AwareHash((unsigned char *)name, strlen(name),
                                    13091204281, 228204732751, 6620830889);
    for (int i = 0; i < depth; i++) lcb_.hash[i] = GenHashSeed(seed++);
    for (int i = 0; i < depth; i++) lcb_.scale[i] = GenHashSeed(seed++);
    for (int i = 0; i < depth; i++) lcb_.hardner[i] = GenHashSeed(seed++);

    t_current = 0;
    max_windows = n_windows + 2;

    // Hoeffding half-width: epsilon(dt) = sqrt(ln(2/delta) / (2 * dt))
    double c_delta = sqrt(log(2.0 / delta) / 2.0);
    eps_table = new double[max_windows + 2];
    eps_table[0] = c_delta;
    for (int i = 1; i <= max_windows + 1; i++)
        eps_table[i] = c_delta / sqrt((double)i);

    std::cout << "[ArgusPP] bucket=" << BucketBytes() << " bytes"
              << " alpha=" << alpha << " pr_min=" << (int)pr_min_
              << " delta=" << delta
              << " (c_delta=" << c_delta << ")" << std::endl;
}

ArgusPP::~ArgusPP()
{
    delete[] lcb_.counts;
    delete[] lcb_.hash;
    delete[] lcb_.scale;
    delete[] lcb_.hardner;
    delete[] eps_table;
}

void ArgusPP::Update(unsigned char *key, val_tp val, unsigned char priority)
{
    (void)val;
    int keylen = lcb_.lgn / 8;
    lcb_.sum += 1;

    Bucket *b = nullptr;
    int index;
    constexpr int MAX_DEPTH = 16;
    int cached_indices[MAX_DEPTH];

    // ============ Stage 1: sequential probe ============
    for (int i = 0; i < lcb_.depth; i++)
    {
        unsigned long bucket = MurmurHash64A(key, keylen, lcb_.hardner[i]) % lcb_.width;
        index = i * lcb_.width + bucket;
        cached_indices[i] = index;
        b = &lcb_.counts[index];

        if (b->V == 0)
        {
            init_bucket(index, key, keylen, priority);
            return;
        }
        if (memcmp(key, b->key, keylen) == 0)
        {
            if (priority > b->priority)
                b->priority = priority;
            if (b->status == 0)
            {
                b->V += 1;
                b->status = 1;
            }
            return;
        }
    }

    // ============ Stage 2: victim selection (PP score) ============
    // Score = (priority + 1) * (LCB + 1): multiplicative so priority and LCB
    // interact. Low if either factor is low; favors retention of occupants
    // that are both high-priority and likely persistent.
    int best_loc = -1;
    double min_score = 1e18;
    for (int i = 0; i < lcb_.depth; i++)
    {
        index = cached_indices[i];
        b = &lcb_.counts[index];
        int dt = t_current - (int)b->first_tracked + 1;
        if (dt < 1) dt = 1;
        if (dt > max_windows) dt = max_windows;
        double phat = (double)b->V / dt;
        double epsilon = eps_table[dt];
        double lcb = phat - epsilon;

        double score = (double)(b->priority + 1) * (lcb + 1.0);

        if (score < min_score) { min_score = score; best_loc = index; }
    }
    if (best_loc < 0) return;
    b = &lcb_.counts[best_loc];

    // ============ Stage 3: within-window protection ============
    if (b->status == 1) return;

    // ============ Stage 4: PP overlay on three-way decision ============
    int v = b->V;
    int dt = t_current - (int)b->first_tracked + 1;
    if (dt < 1) dt = 1;
    if (dt > max_windows) dt = max_windows;
    double phat = (double)v / dt;
    double epsilon = eps_table[dt];
    double lcb = phat - epsilon;
    double ucb = phat + epsilon;
    if (ucb > 1.0) ucb = 1.0;

    // Case 1 (PP-specific): task criterion violated (priority < pr_min) AND
    // LCB rules out persistence -> evict immediately.
    if ((int)b->priority < (int)pr_min_ && lcb < alpha)
    {
        init_bucket(best_loc, key, keylen, priority);
        return;
    }

    // Protect: LCB confirms p(v) >= alpha.
    if (lcb >= alpha) return;

    // Evict / decay: UCB rules out p(v) >= alpha, OR admission override at V=1.
    // Decay branch decrements only V (priority is an intrinsic attribute, not
    // statistical evidence, so it is not decayed).
    if (ucb < alpha || (v == 1 && phat < alpha))
    {
        if (v <= 1)
            init_bucket(best_loc, key, keylen, priority);   // replace
        else
            b->V -= 1;                                      // decay V only
        return;
    }

    // Defer: CI straddles alpha -- retain v, drop u.
}

void ArgusPP::QueryPP(std::vector<std::pair<key_tp, val_tp>> &results)
{
    int T = t_current + 1;
    int v_thresh = (int)(alpha * T);
    int total = lcb_.width * lcb_.depth;
    for (int i = 0; i < total; i++)
    {
        Bucket *b = &lcb_.counts[i];
        if (b->V == 0 || b->V < v_thresh) continue;
        if ((int)b->priority < (int)pr_min_) continue;
        key_tp reskey;
        memcpy(reskey.key, b->key, lcb_.lgn / 8);
        results.push_back({reskey, (val_tp)b->V});
    }
}

void ArgusPP::QueryPPDetailed(std::vector<PPResult> &results)
{
    int T = t_current + 1;
    int v_thresh = (int)(alpha * T);
    int total = lcb_.width * lcb_.depth;
    for (int i = 0; i < total; i++)
    {
        Bucket *b = &lcb_.counts[i];
        if (b->V == 0 || b->V < v_thresh) continue;
        if ((int)b->priority < (int)pr_min_) continue;
        PPResult r;
        memcpy(r.key.key, b->key, lcb_.lgn / 8);
        r.V = (val_tp)b->V;
        r.priority = b->priority;
        results.push_back(r);
    }
}

void ArgusPP::NewWindow()
{
    int total = lcb_.depth * lcb_.width;
    for (int i = 0; i < total; i++)
        lcb_.counts[i].status = 0;
    t_current++;
}

val_tp ArgusPP::PointQuery(unsigned char *key)
{
    int keylen = lcb_.lgn / 8;
    for (int i = 0; i < lcb_.depth; i++)
    {
        unsigned long bucket = MurmurHash64A(key, keylen, lcb_.hardner[i]) % lcb_.width;
        unsigned long index = i * lcb_.width + bucket;
        if (memcmp(lcb_.counts[index].key, key, keylen) == 0)
            return lcb_.counts[index].V;
    }
    return 0;
}

val_tp ArgusPP::GetCount() { return lcb_.sum; }

void ArgusPP::Reset()
{
    lcb_.sum = 0;
    int total = lcb_.depth * lcb_.width;
    memset(lcb_.counts, 0, sizeof(Bucket) * total);
    t_current = 0;
}
