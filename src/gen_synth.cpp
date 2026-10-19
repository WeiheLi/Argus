// gen_synth.cpp -- synthetic stream generators (8-byte key records).
//   gen_synth zipf   OUT N U s seed                 : N arrivals, universe U, Zipf exponent s
//   gen_synth markov OUT n_flows W avg_p rho seed pk : 2-state Markov per-window presence with
//                                                      stationary p_u (Zipf across flows), lag-1 corr rho
//   gen_synth lrd    OUT n_flows W avg_p a seed pk   : Pareto(a) on/off durations (long-range dependent
//                                                      for 1<a<2, Hurst H=(3-a)/2)
//   gen_synth onoff  OUT n_flows W avg_d L seed pk   : deterministic on/off period L (duty cycle Zipf)
// pk = packets emitted per (window, flow) firing. Arrivals inside a window are shuffled.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <string>

static void write_keys(const char *out, std::vector<uint64_t> &keys)
{
    FILE *f = fopen(out, "wb"); fwrite(keys.data(), 8, keys.size(), f); fclose(f);
    fprintf(stderr, "wrote %zu arrivals to %s\n", keys.size(), out);
}
static uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }

// per-flow marginal rates: Zipf-shaped across flows, scaled to mean avg_p
static std::vector<double> flow_rates(int n, double avg_p, std::mt19937_64 &rng)
{
    std::vector<double> p(n);
    double s = 1.0; double sum = 0;
    for (int i = 0; i < n; i++) { p[i] = 1.0 / pow(i + 1, s); sum += p[i]; }
    double scale = avg_p * n / sum;
    for (int i = 0; i < n; i++) { p[i] = std::min(1.0, p[i] * scale); }
    std::shuffle(p.begin(), p.end(), rng);
    return p;
}

int main(int argc, char **argv)
{
    if (argc < 3) { fprintf(stderr, "see header\n"); return 1; }
    std::string mode = argv[1]; const char *out = argv[2];
    if (mode == "zipf") {
        uint64_t N = strtoull(argv[3], 0, 10); uint64_t U = strtoull(argv[4], 0, 10); double s = atof(argv[5]); uint64_t seed = strtoull(argv[6], 0, 10);
        std::mt19937_64 rng(seed);
        std::vector<double> cdf(U); double acc = 0; for (uint64_t i = 0; i < U; i++) { acc += 1.0 / pow((double)(i + 1), s); cdf[i] = acc; }
        for (auto &c : cdf) c /= acc;
        std::uniform_real_distribution<double> uni(0, 1);
        std::vector<uint64_t> keys(N);
        for (uint64_t i = 0; i < N; i++) { double u = uni(rng); uint64_t r = std::lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin(); keys[i] = mix(r + 1) | 1; }
        write_keys(out, keys); return 0;
    }
    int n = atoi(argv[3]); int W = atoi(argv[4]); double avg = atof(argv[5]); double par = atof(argv[6]); uint64_t seed = strtoull(argv[7], 0, 10); int pk = argc > 8 ? atoi(argv[8]) : 2;
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uni(0, 1);
    std::vector<double> p = flow_rates(n, avg, rng);
    { std::string rp = std::string(out) + ".rates"; FILE *fr = fopen(rp.c_str(), "w");
      for (int u = 0; u < n; u++) fprintf(fr, "%016llx %.8f\n", (unsigned long long)(mix((uint64_t)u + 1) | 1), p[u]);
      fclose(fr); }
    std::vector<std::vector<uint32_t>> perwin(W);
    if (mode == "iid") {
        // independent per-window presence with probability p_u (one packet per firing), and every window
        // padded with one-off filler keys to exactly M = pk packets, so that count-based windows coincide
        // with the generator's windows and p_u is the per-window presence probability of the working model
        uint64_t M = (uint64_t)pk; uint64_t filler = (uint64_t)n + 1;
        std::vector<uint64_t> keys;
        for (int w = 0; w < W; w++) {
            std::vector<uint64_t> win; win.reserve(M);
            for (int u = 0; u < n; u++) if (uni(rng) < p[u]) win.push_back(mix((uint64_t)u + 1) | 1);
            if (win.size() > M) { fprintf(stderr, "window %d has %zu > M=%llu firings\n", w, win.size(), (unsigned long long)M); return 1; }
            while (win.size() < M) win.push_back(mix(filler++) | 1);
            std::shuffle(win.begin(), win.end(), rng);
            keys.insert(keys.end(), win.begin(), win.end());
        }
        write_keys(out, keys); return 0;
    }
    if (mode == "markov") {
        double rho = par;
        for (int u = 0; u < n; u++) {
            double pu = p[u]; double p01 = pu * (1 - rho), p11 = pu + (1 - pu) * rho;
            int x = uni(rng) < pu ? 1 : 0;
            for (int w = 0; w < W; w++) { if (x) perwin[w].push_back(u); double t = x ? p11 : p01; x = uni(rng) < t ? 1 : 0; }
        }
    } else if (mode == "lrd") {
        // alternating ON/OFF periods with Pareto(a, xmin) lengths; mean ON/OFF chosen to give marginal p_u
        double a = par;
        for (int u = 0; u < n; u++) {
            double pu = std::max(1e-4, std::min(0.999, p[u]));
            // mean of Pareto with xmin: a*xmin/(a-1). choose xmin_on, xmin_off so that mean_on/(mean_on+mean_off)=pu, with mean_on+mean_off = 40 windows
            double cycle = 40.0; double mean_on = pu * cycle, mean_off = cycle - mean_on;
            double xon = mean_on * (a - 1) / a, xoff = mean_off * (a - 1) / a;
            auto pareto = [&](double xm) { double u2 = uni(rng); double x = xm / pow(1 - u2, 1.0 / a); return x; };
            int x = uni(rng) < pu ? 1 : 0; double w = -uni(rng) * (x ? pareto(xon) : pareto(xoff)); // random phase
            for (int t = 0; t < W; ) {
                double len = x ? pareto(xon) : pareto(xoff); int L = (int)ceil(len + w); w = 0; if (L < 1) L = 1;
                for (int k = 0; k < L && t < W; k++, t++) if (x) perwin[t].push_back(u);
                x ^= 1;
            }
        }
    } else if (mode == "onoff") {
        int L = (int)par;
        for (int u = 0; u < n; u++) {
            double d = std::min(1.0, p[u]); int Lon = std::max(1, (int)round(d * L)); int phase = (int)(uni(rng) * L);
            for (int w = 0; w < W; w++) if (((w + phase) % L) < Lon) perwin[w].push_back(u);
        }
    } else { fprintf(stderr, "bad mode\n"); return 1; }
    std::vector<uint64_t> keys;
    for (int w = 0; w < W; w++) {
        std::vector<uint64_t> win; win.reserve(perwin[w].size() * pk);
        for (uint32_t u : perwin[w]) for (int k = 0; k < pk; k++) win.push_back(mix((uint64_t)u + 1) | 1);
        std::shuffle(win.begin(), win.end(), rng);
        keys.insert(keys.end(), win.begin(), win.end());
    }
    write_keys(out, keys);
    return 0;
}
