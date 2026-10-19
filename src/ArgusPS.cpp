#include "ArgusPS.hpp"
#include <math.h>
extern unsigned long g_hash_seed_offset;

ArgusPS::ArgusPS(int depth, int width, int lgn,
                         double alpha_val, double theta_val,
                         double delta_val, int n_windows)
    : alpha(alpha_val), theta(theta_val), delta(delta_val)
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
    char name[] = "LCBSketchPS";
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

    std::cout << "[ArgusPS] bucket=" << BucketBytes() << " bytes"
              << " alpha=" << alpha << " theta=" << theta
              << " delta=" << delta
              << " (c_delta=" << c_delta << ")" << std::endl;
}

ArgusPS::~ArgusPS()
{
    delete[] lcb_.counts;
    delete[] lcb_.hash;
    delete[] lcb_.scale;
    delete[] lcb_.hardner;
    delete[] eps_table;
}

void ArgusPS::Update(unsigned char *key, val_tp val)
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
            init_bucket(index, key, keylen);
            return;
        }
        if (memcmp(key, b->key, keylen) == 0)
        {
            if (b->C < 0xFFFF) b->C += 1;  // cap at 65535
            if (b->status == 0)
            {
                b->V += 1;
                b->status = 1;
            }
            return;
        }
    }

    // ============ Stage 2: victim selection (PS score) ============
    // Score = LCB - penalty, where penalty = (C/V) / (theta + C/V) in [0,1).
    // Favors retention of sparse-looking occupants (small C/V).
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

        double ratio = (b->V > 0) ? (double)b->C / (double)b->V : 0;
        double penalty = ratio / (theta + ratio);
        double score = lcb - penalty;

        if (score < min_score) { min_score = score; best_loc = index; }
    }
    if (best_loc < 0) return;
    b = &lcb_.counts[best_loc];

    // ============ Stage 3: within-window protection ============
    if (b->status == 1) return;

    // ============ Stage 4: PS overlay on three-way decision ============
    int v = b->V;
    int c = b->C;
    int dt = t_current - (int)b->first_tracked + 1;
    if (dt < 1) dt = 1;
    if (dt > max_windows) dt = max_windows;
    double phat = (double)v / dt;
    double epsilon = eps_table[dt];
    double lcb = phat - epsilon;
    double ucb = phat + epsilon;
    if (ucb > 1.0) ucb = 1.0;

    // Case 1 (PS-specific): task criterion violated (C/V > theta) AND
    // LCB rules out persistence -> evict immediately.
    if ((double)c > theta * (double)v && lcb < alpha)
    {
        init_bucket(best_loc, key, keylen);
        return;
    }

    // Protect: LCB confirms p(v) >= alpha.
    if (lcb >= alpha) return;

    // Evict / decay: UCB rules out p(v) >= alpha, OR admission override at V=1.
    if (ucb < alpha || (v == 1 && phat < alpha))
    {
        if (v <= 1)
            init_bucket(best_loc, key, keylen);       // replace
        else {
            b->V -= 1;                                // decay V and C in step
            if (b->C > 0) b->C -= 1;
        }
        return;
    }

    // Defer: CI straddles alpha -- retain v, drop u.
}

void ArgusPS::QueryPS(std::vector<std::pair<key_tp, val_tp>> &results)
{
    int T = t_current + 1;
    int v_thresh = (int)(alpha * T);
    int total = lcb_.width * lcb_.depth;
    for (int i = 0; i < total; i++)
    {
        Bucket *b = &lcb_.counts[i];
        if (b->V == 0 || b->V < v_thresh) continue;
        if ((double)b->C > theta * (double)b->V) continue;
        key_tp reskey;
        memcpy(reskey.key, b->key, lcb_.lgn / 8);
        results.push_back({reskey, (val_tp)b->V});
    }
}

void ArgusPS::QueryPSDetailed(std::vector<PSResult> &results)
{
    int T = t_current + 1;
    int v_thresh = (int)(alpha * T);
    int total = lcb_.width * lcb_.depth;
    for (int i = 0; i < total; i++)
    {
        Bucket *b = &lcb_.counts[i];
        if (b->V == 0 || b->V < v_thresh) continue;
        if ((double)b->C > theta * (double)b->V) continue;
        PSResult r;
        memcpy(r.key.key, b->key, lcb_.lgn / 8);
        r.V = (val_tp)b->V;
        r.C = b->C;
        results.push_back(r);
    }
}

void ArgusPS::NewWindow()
{
    int total = lcb_.depth * lcb_.width;
    for (int i = 0; i < total; i++)
        lcb_.counts[i].status = 0;
    t_current++;
}

val_tp ArgusPS::PointQuery(unsigned char *key)
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

val_tp ArgusPS::GetCount() { return lcb_.sum; }

void ArgusPS::Reset()
{
    lcb_.sum = 0;
    int total = lcb_.depth * lcb_.width;
    memset(lcb_.counts, 0, sizeof(Bucket) * total);
    t_current = 0;
}
