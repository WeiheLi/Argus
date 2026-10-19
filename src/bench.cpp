// bench.cpp -- unified benchmark for Argus and baselines on .dat streams.
//
//   bench --stream FILE --fmt 13|8|16 --name LABEL [--W 1600] [--alpha 0.3]
//         --task pi|ps|ph|pp [--theta 2.0] [--kappa 5] [--pimin 4]
//         --methods a,b,c --mems 10,20,30 --seeds 0,1,2 [--depth 4] [--delta 0.05]
//         --out results.csv [--fp] [--gtdir DIR] [--limit N] [--reps 1]
//
// Streams are loaded fully into memory as 8-byte keys. Windows are cut by
// arrival count into W equal parts. The ground truth (V, C per item) is cached
// next to the stream file.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <ctime>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <functional>

#include "datatypes.hpp"
#include "Argus.hpp"
#include "ArgusPS.hpp"
#include "ArgusPH.hpp"
#include "ArgusPP.hpp"
#include "PSketch.hpp"
#include "StableSketch.hpp"
#include "Pandora.hpp"
#include "WavingSketch.hpp"
#include "Hypersistent.hpp"
#include "OnOff.hpp"
#include "PSSketch.hpp"
#include "TightSketchPH.hpp"
#include "TinyLFU.hpp"
#include "PAComposite.hpp"

unsigned long g_hash_seed_offset = 0;
bool g_fp = false;
// the baselines' own parameters, at the settings the paper reports
int g_psketch_k = 18, g_pandora_def = 50;
static const int g_onoff_w = 4, g_waving_d = 16, g_waving_c = 4;
static const double g_hyper_hot = 0.4, g_tinylfu_frac = 0.75;

static inline double now_sec()
{
    struct timespec tv; clock_gettime(CLOCK_MONOTONIC, &tv);
    return tv.tv_sec + tv.tv_nsec * 1e-9;
}

// ---------------------------------------------------------------- options
struct Opt {
    std::string stream, name, task = "pi", out = "results.csv", gtdir;
    int fmt = 13, W = 1600, depth = 4, reps = 1;
    double alpha = 0.3, theta = 2.0, kappa = 5.0, delta = 0.05;
    int pimin = 4;
    std::vector<std::string> methods; std::vector<int> mems; std::vector<int> seeds;
    bool header = false; bool gtonly = false; uint64_t limit = 0;
};
static std::vector<std::string> split(const std::string &s, char c)
{
    std::vector<std::string> r; std::stringstream ss(s); std::string t;
    while (std::getline(ss, t, c)) if (!t.empty()) r.push_back(t);
    return r;
}
static Opt parse(int argc, char **argv)
{
    Opt o;
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        auto next = [&]() { if (i + 1 >= argc) { fprintf(stderr, "missing value for %s\n", a.c_str()); exit(1); } return std::string(argv[++i]); };
        if (a == "--stream") o.stream = next();
        else if (a == "--fmt") o.fmt = atoi(next().c_str());
        else if (a == "--name") o.name = next();
        else if (a == "--W") o.W = atoi(next().c_str());
        else if (a == "--alpha") o.alpha = atof(next().c_str());
        else if (a == "--task") o.task = next();
        else if (a == "--theta") o.theta = atof(next().c_str());
        else if (a == "--kappa") o.kappa = atof(next().c_str());
        else if (a == "--pimin") o.pimin = atoi(next().c_str());
        else if (a == "--methods") o.methods = split(next(), ',');
        else if (a == "--mems") { for (auto &x : split(next(), ',')) o.mems.push_back(atoi(x.c_str())); }
        else if (a == "--seeds") { for (auto &x : split(next(), ',')) o.seeds.push_back(atoi(x.c_str())); }
        else if (a == "--depth") o.depth = atoi(next().c_str());
        else if (a == "--delta") o.delta = atof(next().c_str());
        else if (a == "--out") o.out = next();
        else if (a == "--header") o.header = true;
        else if (a == "--gtonly") o.gtonly = true;
        else if (a == "--fp") g_fp = true;
        else if (a == "--limit") o.limit = strtoull(next().c_str(), nullptr, 10);
        else if (a == "--reps") o.reps = atoi(next().c_str());
        else if (a == "--gtdir") o.gtdir = next();
        else { fprintf(stderr, "unknown option %s\n", a.c_str()); exit(1); }
    }
    if (o.seeds.empty()) o.seeds.push_back(0);
    if (o.stream.empty() || (!o.gtonly && (o.methods.empty() || o.mems.empty()))) { fprintf(stderr, "usage: see header\n"); exit(1); }
    return o;
}

// ---------------------------------------------------------------- stream
static std::vector<uint64_t> load_stream(const std::string &path, int fmt, uint64_t limit)
{
    FILE *f = fopen(path.c_str(), "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", path.c_str()); exit(1); }
    fseek(f, 0, SEEK_END); uint64_t sz = ftell(f); fseek(f, 0, SEEK_SET);
    int rec = fmt, off = 0;
    if (fmt == 13) { rec = 13; off = 0; }
    else if (fmt == 8) { rec = 8; off = 0; }
    else if (fmt == 16) { rec = 16; off = 8; }
    else { fprintf(stderr, "bad fmt\n"); exit(1); }
    uint64_t n = sz / rec; if (limit && n > limit) n = limit;
    std::vector<uint64_t> keys(n);
    const size_t CH = 1 << 20; std::vector<unsigned char> buf(CH * rec);
    uint64_t done = 0;
    while (done < n) {
        size_t want = std::min<uint64_t>(CH, n - done);
        size_t got = fread(buf.data(), rec, want, f);
        if (got == 0) break;
        for (size_t i = 0; i < got; i++) memcpy(&keys[done + i], &buf[i * rec + off], 8);
        done += got;
    }
    fclose(f);
    keys.resize(done);
    return keys;
}

// priority class of an item: Zipf over classes 0..7 (higher classes rarer)
static unsigned char priority_of(uint64_t key)
{
    static double cdf[8]; static bool init = false;
    if (!init) { double H = 0; for (int k = 0; k < 8; k++) H += 1.0 / (k + 1); double acc = 0; for (int k = 0; k < 8; k++) { acc += 1.0 / (k + 1) / H; cdf[k] = acc; } init = true; }
    uint64_t h = MurmurHash64A(&key, 8, 0x9E3779B97F4A7C15ULL);
    double u = (double)(h % 1000003) / 1000003.0;
    for (int k = 0; k < 8; k++) if (u <= cdf[k]) return (unsigned char)k;
    return 7;
}

// ---------------------------------------------------------------- ground truth
struct GT { std::unordered_map<uint64_t, uint32_t> V, C; };

static void build_gt(const std::vector<uint64_t> &keys, int W, GT &gt)
{
    uint64_t N = keys.size();
    std::unordered_map<uint64_t, uint32_t> lastwin; lastwin.reserve(1 << 20);
    gt.V.reserve(1 << 20); gt.C.reserve(1 << 20);
    for (uint64_t i = 0; i < N; i++) {
        uint32_t w = (uint32_t)((i * (uint64_t)W) / N);
        uint64_t k = keys[i];
        gt.C[k]++;
        auto it = lastwin.find(k);
        if (it == lastwin.end()) { lastwin[k] = w; gt.V[k] = 1; }
        else if (it->second != w) { it->second = w; gt.V[k]++; }
    }
}
static bool load_gt(const std::string &path, GT &gt)
{
    FILE *f = fopen(path.c_str(), "rb"); if (!f) return false;
    uint64_t n; if (fread(&n, 8, 1, f) != 1) { fclose(f); return false; }
    gt.V.reserve(n); gt.C.reserve(n);
    for (uint64_t i = 0; i < n; i++) { uint64_t k; uint32_t v, c; if (fread(&k, 8, 1, f) != 1 || fread(&v, 4, 1, f) != 1 || fread(&c, 4, 1, f) != 1) { fclose(f); return false; } gt.V[k] = v; gt.C[k] = c; }
    fclose(f); return true;
}
static void save_gt(const std::string &path, const GT &gt)
{
    FILE *f = fopen(path.c_str(), "wb"); if (!f) return;
    uint64_t n = gt.V.size(); fwrite(&n, 8, 1, f);
    for (auto &kv : gt.V) { uint32_t c = gt.C.at(kv.first); fwrite(&kv.first, 8, 1, f); fwrite(&kv.second, 4, 1, f); fwrite(&c, 4, 1, f); }
    fclose(f);
}

// ---------------------------------------------------------------- method wrapper
struct Method {
    virtual ~Method() {}
    virtual void update(const unsigned char *key, unsigned char prio) = 0;
    virtual void window_end() {}                    // hooks before NewWindow (hotness etc.)
    virtual void new_window() = 0;
    virtual void query(int v_thresh, std::vector<std::pair<key_tp, val_tp>> &r) = 0;
    virtual ArgusStats *argus_stats() { return nullptr; }
};

#define KEYP(k) ((unsigned char *)(k))
#define LGNBITS (g_fp ? 32 : 8 * LGN)

struct MArgus : Method {
    ArgusSketch s;
    MArgus(int d, int w, double a, double dl, int W) : s(d, w, LGNBITS, a, dl, W) {}
    void update(const unsigned char *k, unsigned char) override { s.Update(k); }
    void new_window() override { s.NewWindow(); }
    void query(int v, std::vector<std::pair<key_tp, val_tp>> &r) override { s.QueryPersistent(v, r); }
    ArgusStats *argus_stats() override { return &s.stats; }
};
struct MArgusB : Method { ArgusPacked s;
    MArgusB(int d, int w, double a, double dl, int W) : s(d, w, LGNBITS, a, dl, W) {}
    void update(const unsigned char *k, unsigned char) override { s.Update(k); }
    void new_window() override { s.NewWindow(); }
    void query(int v, std::vector<std::pair<key_tp, val_tp>> &r) override { s.QueryPersistent(v, r); }
    ArgusStats *argus_stats() override { return &s.stats; } };
struct MPSketch : Method { PSketch s; MPSketch(int d, int w) : s(d, w, LGNBITS) {}
    void update(const unsigned char *k, unsigned char) override { s.Update(KEYP(k), 1); }
    void window_end() override { s.QueryHotness(); } void new_window() override { s.NewWindow(); }
    void query(int v, std::vector<std::pair<key_tp, val_tp>> &r) override { s.Query(v - 1, r); } };  // Query uses '>'
struct MPandora : Method { Pandora s; MPandora(int d, int w) : s(d, w, LGNBITS) {}
    void update(const unsigned char *k, unsigned char) override { s.Update(KEYP(k), 1); }
    void window_end() override { s.Inactivity(); } void new_window() override { s.NewWindow(); }
    void query(int v, std::vector<std::pair<key_tp, val_tp>> &r) override { s.Query(v - 1, r); } };
struct MWaving : Method { WavingSketch s; MWaving(int bytes) : s(bytes, LGNBITS, g_waving_d, g_waving_c) {}
    void update(const unsigned char *k, unsigned char) override { s.Update(KEYP(k), 1); }
    void new_window() override { s.NewWindow(); }
    void query(int v, std::vector<std::pair<key_tp, val_tp>> &r) override { s.Query(v, r); } };
struct MHyper : Method { Hypersistent s; MHyper(int bytes) : s(bytes, LGNBITS, g_hyper_hot, 85, 15, 100, 1, 16, 2, 2, 3) {}
    void update(const unsigned char *k, unsigned char) override { s.Update(KEYP(k), 1); }
    void new_window() override { s.NewWindow(); }
    void query(int v, std::vector<std::pair<key_tp, val_tp>> &r) override { s.Query(v, r); } };
struct MOnOff : Method { OnOff s; MOnOff(int bytes) : s(bytes, LGNBITS, g_onoff_w) {}
    void update(const unsigned char *k, unsigned char) override { s.Update(KEYP(k), 1); }
    void new_window() override { s.NewWindow(); }
    void query(int v, std::vector<std::pair<key_tp, val_tp>> &r) override { s.Query(v, r); } };
struct MTinyLFU : Method { TinyLFU s; MTinyLFU(int bytes, int d, int bb = 0) : s(bytes, g_fp ? 32 : 8 * LGN, d, g_tinylfu_frac, 4, 128, bb) {}
    void update(const unsigned char *k, unsigned char) override { s.Update(k); }
    void new_window() override { s.NewWindow(); }
    void query(int v, std::vector<std::pair<key_tp, val_tp>> &r) override { s.Query(v, r); } };
// compound tasks
struct MArgusPS : Method { ArgusPS s; MArgusPS(int d, int w, double a, double th, double dl, int W) : s(d, w, LGNBITS, a, th, dl, W) {}
    void update(const unsigned char *k, unsigned char) override { s.Update(KEYP(k), 1); }
    void new_window() override { s.NewWindow(); }
    void query(int, std::vector<std::pair<key_tp, val_tp>> &r) override { s.QueryPS(r); } };
struct MPSSketch : Method { PSSketch s; MPSSketch(int d, int w, int pl, int pth, double th) : s(d, w, LGNBITS, pl, pth, th) {}
    void update(const unsigned char *k, unsigned char) override { s.Update(KEYP(k), 1); }
    void new_window() override { s.NewWindow(); }
    void query(int, std::vector<std::pair<key_tp, val_tp>> &r) override { s.QueryPS(r); } };
struct MArgusPH : Method { ArgusPH s; MArgusPH(int d, int w, double a, double ka, double dl, int W) : s(d, w, LGNBITS, a, ka, dl, W) {}
    void update(const unsigned char *k, unsigned char) override { s.Update(KEYP(k), 1); }
    void new_window() override { s.NewWindow(); }
    void query(int, std::vector<std::pair<key_tp, val_tp>> &r) override { s.QueryPH(r); } };
struct MTightPH : Method { TightSketchPH s; MTightPH(int d, int w, double a, double ka, int W) : s(d, w, LGNBITS, a, ka, W, 10) {}
    void update(const unsigned char *k, unsigned char) override { s.Update(KEYP(k), 1); }
    void new_window() override { s.NewWindow(); }
    void query(int, std::vector<std::pair<key_tp, val_tp>> &r) override { s.QueryPH(r); } };
struct MArgusPP : Method { ArgusPP s; MArgusPP(int d, int w, double a, int pm, double dl, int W) : s(d, w, LGNBITS, a, (unsigned char)pm, dl, W) {}
    void update(const unsigned char *k, unsigned char p) override { s.Update(KEYP(k), 1, p); }
    void new_window() override { s.NewWindow(); }
    void query(int, std::vector<std::pair<key_tp, val_tp>> &r) override { s.QueryPP(r); } };
template <class P>
struct MPAComp : Method { PAComposite<P> c; int pimin; std::function<void(P &)> hook;
    MPAComp(int bytes, int d, int pm, int pers_bytes) : c(bytes, g_fp ? 32 : 8 * LGN, d, 0.35, pers_bytes, g_fp ? 7 : 0), pimin(pm) {}
    void update(const unsigned char *k, unsigned char p) override { c.pa_.Update(k, p); c.pers_->Update(KEYP(k), 1); }
    void window_end() override { if (hook) hook(*c.pers_); }
    void new_window() override { c.pers_->NewWindow(); }
    void query(int v, std::vector<std::pair<key_tp, val_tp>> &r) override {
        std::vector<std::pair<key_tp, val_tp>> tmp; c.pers_->Query(v - 1, tmp);
        for (auto &e : tmp) { int p = c.pa_.Priority(e.first.key); if (p >= pimin) r.push_back(e); }
    } };

static Method *make(const std::string &m, int mem_kb, const Opt &o, int v_thresh)
{
    int bytes = mem_kb * 1024, d = o.depth;
    int KB = g_fp ? 4 : 0;   // bytes saved per bucket by a 4-byte identifier
    auto w13 = [&](int bb) { int w = bytes / ((bb - KB) * d); return w < 1 ? 1 : w; };
    if (m == "argus") { if (g_fp) return new MArgusB(d, w13(8 + 4), o.alpha, o.delta, o.W); return new MArgus(d, w13(13), o.alpha, o.delta, o.W); }
    // each count-based baseline at the width its own fields need: an identifier,
    // two 16-bit counters and a status byte
    if (m == "psketch") return new MPSketch(d, w13(13));
    if (m == "pandora") return new MPandora(d, w13(13));
    if (m == "waving") return new MWaving(bytes);
    if (m == "hyper") return new MHyper(bytes);
    if (m == "onoff") return new MOnOff(bytes);
    if (m == "tinylfu") return new MTinyLFU(bytes, d, g_fp ? 7 : 0);
    if (m == "argus_ps") return new MArgusPS(d, w13(15), o.alpha, o.theta, o.delta, o.W);
    if (m == "pssketch") { int cl = bytes * 80 / 100, pl = bytes - cl; int w = cl / (PSSketch::BucketBytes() * d); if (w < 1) w = 1; int plr = pl / (LGN + 5); if (plr < 16) plr = 16; return new MPSSketch(d, w, plr, v_thresh, o.theta); }
    if (m == "argus_ph") return new MArgusPH(d, w13(15), o.alpha, o.kappa, o.delta, o.W);
    if (m == "tightph") return new MTightPH(d, w13(TightSketchPH::BucketBytes()), o.alpha, o.kappa, o.W);
    if (m == "argus_pp") return new MArgusPP(d, w13(14), o.alpha, o.pimin, o.delta, o.W);
    // the composite splits its budget, 35% to the priority layer and the rest to
    // the persistence sketch, which is charged the same bucket as on its own
    if (m == "pa_pandora") { auto *x = new MPAComp<Pandora>(bytes, d, o.pimin, 13 - KB); x->hook = [](Pandora &p) { p.Inactivity(); }; return x; }
    if (m == "pa_stable") return new MPAComp<StableSketch>(bytes, d, o.pimin, 13 - KB);
    fprintf(stderr, "unknown method %s\n", m.c_str()); exit(1);
}


int main(int argc, char **argv)
{
    Opt o = parse(argc, argv);
    std::vector<uint64_t> keys = load_stream(o.stream, o.fmt, o.limit);
    uint64_t N = keys.size();
    fprintf(stderr, "[%s] %llu arrivals\n", o.name.c_str(), (unsigned long long)N);

    // ground truth
    GT gt;
    std::string base = o.stream.substr(o.stream.find_last_of('/') + 1);
    std::string gtpath = (o.gtdir.empty() ? o.stream : (o.gtdir + "/" + base)) + ".W" + std::to_string(o.W) + ".L" + std::to_string(N) + ".gt";
    if (!load_gt(gtpath, gt)) { build_gt(keys, o.W, gt); save_gt(gtpath, gt); }
    if (o.gtonly) {
        printf("%s W=%d arrivals=%llu distinct=%zu", o.name.c_str(), o.W, (unsigned long long)N, gt.V.size());
        for (double al : {0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5}) { int th = (int)ceil(al * o.W - 1e-9); size_t c = 0; for (auto &kv : gt.V) if ((int)kv.second >= th) c++; printf("  a=%.2f:%zu", al, c); }
        printf("\n"); return 0;
    }
    int v_thresh = (int)ceil(o.alpha * o.W - 1e-9);
    std::unordered_set<uint64_t> truth; std::unordered_map<uint64_t, uint32_t> truthV;
    for (auto &kv : gt.V) {
        if ((int)kv.second < v_thresh) continue;
        uint32_t C = gt.C[kv.first]; double ratio = (double)C / kv.second;
        bool ok = true;
        if (o.task == "ps") ok = ratio <= o.theta;
        else if (o.task == "ph") ok = ratio >= o.kappa;
        else if (o.task == "pp") ok = priority_of(kv.first) >= o.pimin;
        if (ok) { truth.insert(kv.first); truthV[kv.first] = kv.second; }
    }
    fprintf(stderr, "[%s] distinct=%zu truth(%s)=%zu thresh=%d\n", o.name.c_str(), gt.V.size(), o.task.c_str(), truth.size(), v_thresh);

    std::vector<uint64_t> fpkeys;
    if (g_fp) { fpkeys.resize(N); for (uint64_t i = 0; i < N; i++) fpkeys[i] = ArgusPacked::Fingerprint((const unsigned char *)&keys[i], 8); }
    const std::vector<uint64_t> &ukeys = g_fp ? fpkeys : keys;
    // precompute per-arrival priorities for pp task
    std::vector<unsigned char> prio;
    if (o.task == "pp") { prio.resize(N); for (uint64_t i = 0; i < N; i++) prio[i] = priority_of(keys[i]); }


    std::ofstream out(o.out, std::ios::app);
    if (o.header) out << "stream,task,method,mem_kb,seed,depth,delta,W,alpha,precision,recall,f1,aae,are,n_truth,n_reported,tp,throughput_mops,"
                        "n_update,s1_empty,s1_match,s2_collision,s3_drop,protect,defer,evict_ucb,evict_override,decay\n";

    for (int seed : o.seeds) for (int mem : o.mems) for (auto &mname : o.methods) {
        g_hash_seed_offset = (unsigned long)seed * 7919ULL;
        srand(1000 + seed);
        Method *M = make(mname, mem, o, v_thresh);
        ArgusStats *st = M->argus_stats();

        double best_time = 1e18;
        for (int rep = 0; rep < o.reps; rep++) {
            if (rep > 0) { delete M; M = make(mname, mem, o, v_thresh); st = M->argus_stats(); }
            uint32_t cur = 0;
            double t0 = now_sec();
            for (uint64_t i = 0; i < N; i++) {
                uint32_t w = (uint32_t)((i * (uint64_t)o.W) / N);
                while (cur < w) {
                    M->window_end(); M->new_window(); cur++;
                }
                M->update((const unsigned char *)&ukeys[i], prio.empty() ? 0 : prio[i]);
            }
            double t1 = now_sec();
            if (t1 - t0 < best_time) best_time = t1 - t0;
        }
        double thr = N / best_time / 1e6;

        std::vector<std::pair<key_tp, val_tp>> res; M->query(v_thresh, res);
        std::unordered_map<uint64_t, val_tp> pred; pred.reserve(res.size() * 2);
        for (auto &r : res) { uint64_t k; memcpy(&k, r.first.key, 8); if (g_fp) k &= 0xFFFFFFFFULL; auto it = pred.find(k); if (it == pred.end() || r.second > it->second) pred[k] = r.second; }
        const std::unordered_map<uint64_t, uint32_t> *truthRef = &truthV;
        std::unordered_map<uint64_t, uint32_t> truthFP;
        if (g_fp) { for (auto &kv : truthV) { uint64_t f = ArgusPacked::Fingerprint((const unsigned char *)&kv.first, 8); auto it = truthFP.find(f); if (it == truthFP.end() || kv.second > it->second) truthFP[f] = kv.second; } truthRef = &truthFP; }
        uint64_t tp = 0; double aae = 0, are = 0;
        for (auto &kv : pred) { auto it = truthRef->find(kv.first); if (it != truthRef->end()) { tp++; double e = fabs((double)kv.second - it->second); aae += e; are += e / it->second; } }
        double prec = pred.empty() ? 0 : (double)tp / pred.size();
        double rec = truthRef->empty() ? 0 : (double)tp / truthRef->size();
        double f1 = (prec + rec > 0) ? 2 * prec * rec / (prec + rec) : 0;
        aae = tp ? aae / tp : 0; are = tp ? are / tp : 0;

        out << o.name << "," << o.task << "," << mname << (g_fp ? "@fp" : "") << "," << mem << "," << seed << "," << o.depth << "," << o.delta << "," << o.W << "," << o.alpha << ","
            << prec << "," << rec << "," << f1 << "," << aae << "," << are << "," << truth.size() << "," << pred.size() << "," << tp << "," << thr;
        if (st) out << "," << st->n_update << "," << st->s1_empty << "," << st->s1_match << "," << st->s2_collision << "," << st->s3_drop << "," << st->protect << "," << st->defer << "," << st->evict_ucb << "," << st->evict_override << "," << st->decay;
        else out << ",,,,,,,,,,";
        out << "\n"; out.flush();
        fprintf(stderr, "  %-18s mem=%4d seed=%d  P=%.4f R=%.4f F1=%.4f AAE=%.2f thr=%.2f Mops\n", mname.c_str(), mem, seed, prec, rec, f1, aae, thr);

        delete M;
    }

    return 0;
}
