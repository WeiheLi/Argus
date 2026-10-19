
// ---------------------------------------------------------------------------
// ARGUS on an 8-byte bucket:
//   32-bit key fingerprint | 16-bit V | 15-bit first_tracked | 1-bit status.
// The rule is the one above; only the layout differs, the window of the first
// arrival and the status bit sharing a single 16-bit field. Items are
// identified by their 32-bit fingerprint (queries return it in the first four
// key bytes) and the evaluation matches on the same fingerprint.
class ArgusPacked
{
public:
#pragma pack(push, 1)
    struct Bucket { uint32_t fp; uint16_t V; uint16_t t0s; };
#pragma pack(pop)
    ArgusStats stats;
    ArgusPacked(int depth, int width, int lgn, double alpha, double delta, int n_windows)
        : depth_(depth), width_(width), keylen_(lgn / 8), alpha_(alpha)
    {
        total_ = depth * width; counts_ = new Bucket[total_]; memset(counts_, 0, sizeof(Bucket) * (size_t)total_);
        hardner_ = new unsigned long[depth];
        char name[] = "LCBSketch";
        unsigned long seed = g_hash_seed_offset + AwareHash((unsigned char *)name, strlen(name), 13091204281, 228204732751, 6620830889);
        for (int i = 0; i < depth; i++) GenHashSeed(seed++);
        for (int i = 0; i < depth; i++) GenHashSeed(seed++);
        for (int i = 0; i < depth; i++) hardner_[i] = GenHashSeed(seed++);
        t_current_ = 0; max_windows_ = n_windows + 2;
        double c = sqrt(log(2.0 / delta) / 2.0);
        eps_table_ = new double[max_windows_ + 2]; eps_table_[0] = c;
        for (int i = 1; i <= max_windows_ + 1; i++) eps_table_[i] = c / sqrt((double)i);
    }
    ~ArgusPacked() { delete[] counts_; delete[] hardner_; delete[] eps_table_; }
    static int BucketBytes() { return (int)sizeof(Bucket); }
    static inline uint32_t Fingerprint(const unsigned char *key, int keylen) { return (uint32_t)(MurmurHash64A(key, keylen, 0x5bd1e995ULL) >> 16) | 1u; }
    static inline int t0(const Bucket *b) { return b->t0s & 0x7FFF; }
    static inline bool st(const Bucket *b) { return (b->t0s & 0x8000) != 0; }
    inline int tracked(const Bucket *b) const { int dt = t_current_ - t0(b) + 1; if (dt < 1) dt = 1; if (dt > max_windows_) dt = max_windows_; return dt; }
    inline void init_bucket(Bucket *b, uint32_t fp) { b->fp = fp; b->V = 1; b->t0s = (uint16_t)((t_current_ & 0x7FFF) | 0x8000); }
    void Update(const unsigned char *key)
    {
        stats.n_update++;
        uint32_t fp; if (keylen_ == 4) memcpy(&fp, key, 4); else fp = Fingerprint(key, keylen_);   // in fingerprint mode the key already is the identifier
        int idx[16]; Bucket *b = nullptr;
        for (int i = 0; i < depth_; i++) {
            int index = i * width_ + (int)(MurmurHash64A(key, keylen_, hardner_[i]) % width_);
            idx[i] = index; b = &counts_[index];
            if (b->V == 0) { init_bucket(b, fp); stats.s1_empty++; return; }
            if (b->fp == fp) { stats.s1_match++; if (!st(b)) { b->V += 1; b->t0s |= 0x8000; stats.s1_match_inc++; } return; }
        }
        stats.s2_collision++;
        int best = -1;
        double min_score = 1e18;
        for (int i = 0; i < depth_; i++) { b = &counts_[idx[i]]; int dt = tracked(b); double lo = (double)b->V / dt - eps_table_[dt]; if (lo < min_score) { min_score = lo; best = idx[i]; } }
        b = &counts_[best];
        if (st(b)) { stats.s3_drop++; return; }
        int v = b->V, dt = tracked(b); double phat = (double)v / dt, e = eps_table_[dt], lo = phat - e, hi = phat + e; if (hi > 1) hi = 1;
        if (lo >= alpha_) { stats.protect++; return; }
        if (hi < alpha_) { if (v <= 1) { stats.evict_ucb++; init_bucket(b, fp); } else { stats.decay++; b->V -= 1; } return; }
        if (v == 1 && phat < alpha_) { stats.evict_override++; init_bucket(b, fp); return; }
        stats.defer++;
    }
    void NewWindow() { for (int i = 0; i < total_; i++) counts_[i].t0s &= 0x7FFF; t_current_++; }
    void QueryPersistent(int v_thresh, std::vector<std::pair<key_tp, val_tp>> &results) const
    {
        for (int i = 0; i < total_; i++) { const Bucket *b = &counts_[i]; if (b->V > 0 && (int)b->V >= v_thresh) { key_tp k; memset(k.key, 0, LGN); memcpy(k.key, &b->fp, 4); results.push_back({k, (val_tp)b->V}); } }
    }
private:
    int depth_, width_, keylen_, total_; double alpha_; Bucket *counts_; unsigned long *hardner_; int t_current_, max_windows_; double *eps_table_;
};
