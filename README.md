# ARGUS

This repo is the implementation of **ARGUS**, a sketch framework for finding
persistent items in data streams. Persistence sketches decide, at a hash
collision, which of the colliding items keeps the bucket, and every such rule
we are aware of decides on an accumulated count. A count says nothing about how
long the item behind it has been watched: a count of five may have been earned
over five windows or over five hundred.

ARGUS evicts on the **persistence rate** and on the **observation length behind
it** instead. Each bucket holds a confidence interval for its occupant's rate,
whose width shrinks with the number of windows the occupant has been tracked,
and a collision is settled by a three-way verdict over that interval: *protect*
when the interval clears the persistence threshold, *evict* when it falls below
it, and *defer* when it straddles the threshold and the evidence is not yet
conclusive. The same interval carries over to three compound tasks, where an
item must be persistent *and* satisfy a second predicate:

| variant | task | predicate |
|---|---|---|
| `argus_ps` | persistent and sparse | at most `theta` arrivals per active window |
| `argus_ph` | persistent and heavy | at least `kappa` arrivals per active window |
| `argus_pp` | persistent and priority | a priority tag of at least `pimin` |

`technical-report.pdf` accompanies the paper: it carries the proofs of the
guarantees, the additional statements the paper states in prose, and the
figures and settings the paper has no room to print.

## Build

Requires a C++17 compiler. On Linux:

```
cd src
make
```

This produces two binaries: `bench`, the benchmark harness, and `gen_synth`,
the synthetic stream generator.

## Running

```
./bench --stream FILE --fmt 13|8|16 --name LABEL \
        --task pi|ps|ph|pp [--theta 2.0] [--kappa 5] [--pimin 4] \
        --methods argus,psketch,pandora,... --mems 10,20,30 --seeds 0,1,2 \
        [--W 1600] [--alpha 0.3] [--depth 4] [--delta 0.05] [--fp] \
        --gtdir DIR --out results.csv
```

A stream is a binary file of fixed-width records; `--fmt` gives the record
width in bytes and the first 8 bytes of each record are taken as the key.
`--W` cuts the stream by arrival count into that many equal windows, and
`--alpha` is the fraction of windows an item must appear in to count as
persistent. The ground truth is computed once per `(stream, W)` and cached in
`--gtdir`. `--mems` is the total memory budget in KB, charged to every method
alike. `--fp` identifies items by a 32-bit fingerprint of the key rather than
by the key itself, for every method and for the ground truth, so that no method
gains from the choice; this is the configuration the paper reports.

One line per `(method, budget, seed)` is appended to `--out`, with precision,
recall, F1, the average absolute error of the reported persistence, and
throughput. An ARGUS line also carries the per-stage counters: how many
arrivals took an empty bucket or matched an occupant, how many reached a
collision, and how each collision was settled, which is where the shares of
protect, defer and evict verdicts come from.

Example:

```
./bench --stream caida2018.dat --fmt 13 --name caida2018 --task pi \
        --methods argus,psketch,pandora,onoff,hyper,waving,tinylfu \
        --mems 10,20,30,40,50,60 --seeds 0,1,2,3,4 --fp \
        --gtdir ./data --out main.csv
```

Synthetic streams:

```
./gen_synth zipf   OUT N U s seed                  # Zipf popularity
./gen_synth markov OUT n_flows W avg_p rho seed pk # Markov-modulated presence
./gen_synth lrd    OUT n_flows W avg_p a seed pk   # Pareto on/off, Hurst (3-a)/2
./gen_synth onoff  OUT n_flows W avg_d L seed pk   # deterministic on/off cycle
```

## Methods

`--methods` takes a comma-separated list of:

| name | sketch |
|---|---|
| `argus` | ARGUS |
| `argus_ps`, `argus_ph`, `argus_pp` | the three compound-task variants |
| `psketch` | P-Sketch |
| `pandora` | Pandora |
| `onoff` | On-Off Sketch |
| `hyper` | Hypersistent |
| `waving` | Waving-Sketch, in its persistence mode |
| `tinylfu` | a TinyLFU admission policy over the same bucket array |
| `pssketch` | PSSketch, the persistent-and-sparse baseline |
| `tightph` | Tight-Sketch, the persistent-and-heavy baseline |
| `pa_pandora`, `pa_stable` | priority-adaptive composites, the persistent-and-priority baselines |

## Layout

```
README.md
technical-report.pdf
src/
  Argus.hpp, ArgusPacked.inc          ARGUS: the rule and its packed layout
  ArgusPS/PH/PP.{hpp,cpp}             the three compound-task variants
  PSketch, Pandora, OnOff,            the baselines
  Hypersistent, WavingSketch,
  TinyLFU, PSSketch, TightSketchPH,
  StableSketch, PAComposite
  bench.cpp                           benchmark harness
  gen_synth.cpp                       synthetic stream generator
  datatypes.hpp, hash.{c,h}, util.h   shared code
  Makefile
```
