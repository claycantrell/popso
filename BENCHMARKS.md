# Benchmarks

How the solver compares with the open-source solvers [HiGHS](https://highs.dev) and
[SCIP](https://www.scipopt.org). Every answer is checked independently: the returned
solution must be feasible, and its objective must match the known optimum.

_Last updated 2026-11-05 for popso v0.6.0._

## Setup

- Machine: Apple M4 Pro (12 cores, 24 GB), macOS (the 240-instance run: see its section).
- Versions: HiGHS 1.15.1 (highspy) and SCIP 10.0 (PySCIPOpt), with default settings.
- Time limit: 120 s per instance. Memory cap: 2–2.8 GB per solver process.
- Scoring:
  - **Solved** means the solver proved optimality, infeasibility or unboundedness within the
    time limit, and the answer passed the independent check.
  - **Shifted geometric mean** of the run times uses a shift of 1 s; an unsolved instance
    counts as 120 s. Lower is better.
- Gap tolerance: popso runs at a relative gap of 1e-6; the reference solvers run at 1e-9.
- Caveats:
  - The popso MIPLIB numbers come from a single seed.
  - The v0.6.0 numbers were measured with the release candidate. The release adds fixes for time-limit
    overruns on very large models and for two problems found by the 240-instance run (see the
    [CHANGELOG](CHANGELOG.md)); the instances whose results they change were rerun.
  - The reference times were measured separately, partly while other jobs were running.
  - Read small differences as noise.

## MIPLIB 2017 subset

63 instances from [MIPLIB 2017](https://miplib.zib.de), listed in
[`benchmarks/sets/miplib_wide_core.txt`](benchmarks/sets/miplib_wide_core.txt). popso and
HiGHS use 4 threads; SCIP's branch-and-bound is single-threaded.

**How the 63 instances were chosen.** They are a filtered subset of the official MIPLIB 2017
*benchmark* set: instances classified "easy" (known optimum) with at most about 200k nonzeros,
keeping those that at least one of HiGHS, SCIP or popso solves within 120 s, so that every
instance gives a meaningful time. This makes it a practical development set rather than an
unbiased sample, and popso has been tuned on it. The [full 240-instance benchmark
set](#full-miplib-2017-benchmark-set-240-instances) is below.

![MIPLIB: instances solved vs. time](docs/miplib_core.svg)

| solver | solved | shifted geomean (s) |
|---|---|---|
| popso | 59/63 | 9.0 |
| HiGHS | 53/63 | 18.9 |
| SCIP | 50/63 | 21.3 |

popso is fastest on 39 of the 63 instances. The table is seed 0. Over seeds 0, 1 and 2 popso
solves 59 / 61 / 59 instances with a shifted geometric mean of 9.0 / 9.2 / 9.0 s (v0.5.0: 58 / 59
/ 61 and 10.7 / 10.8 / 9.9 s). Single runs of the instances close to the time limit vary a lot
between seeds.

<details><summary>Per instance (seconds; TL = time limit; bold = fastest)</summary>

| instance | popso | HiGHS | SCIP |
|---|---|---|---|
| 30n20b8 | TL | TL | **116.6** |
| CMS750_4 | TL | TL | **77.0** |
| air05 | **15.7** | 16.8 | 21.8 |
| app1-1 | 12.9 | 14.7 | **3.8** |
| beasleyC3 | 21.7 | 12.3 | **9.1** |
| binkar10_1 | **2.2** | 23.0 | 16.6 |
| cbs-cta | **8.8** | 31.7 | 13.5 |
| cod105 | **1.0** | 62.9 | 32.7 |
| dano3_3 | **20.8** | 38.4 | 58.8 |
| dano3_5 | **42.9** | TL | TL |
| decomp2 | 1.1 | 2.1 | **1.0** |
| drayage-100-23 | 1.6 | **0.3** | 3.9 |
| drayage-25-23 | 5.4 | **0.9** | 8.7 |
| eil33-2 | **39.8** | 54.2 | 47.8 |
| enlight_hard | 1.7 | 6.4 | **0.0** |
| exp-1-500-5-5 | **0.3** | 2.5 | 1.3 |
| fiball | **2.2** | 57.1 | TL |
| h80x6320d | 27.9 | **8.5** | 71.9 |
| irp | **3.4** | 13.8 | 10.8 |
| mas76 | **12.4** | 59.2 | 34.9 |
| mc11 | 13.5 | 22.9 | **9.6** |
| mik-250-20-75-4 | **3.0** | 8.6 | 17.1 |
| mzzv11 | **10.9** | 44.4 | TL |
| mzzv42z | **3.7** | 26.1 | 69.3 |
| n2seq36q | **31.6** | 50.7 | 78.2 |
| n5-3 | 15.7 | **11.0** | 29.4 |
| neos-1122047 | **1.1** | 2.9 | 1.8 |
| neos-1171448 | **4.7** | 79.4 | 6.0 |
| neos-1445765 | **18.3** | 43.3 | 20.0 |
| neos-1582420 | **6.5** | 14.5 | 16.9 |
| neos-2978193-inde | **11.1** | 21.5 | TL |
| neos-3004026-krka | **0.5** | TL | 35.1 |
| neos-3083819-nubu | 4.7 | 5.3 | **3.9** |
| neos-3381206-awhea | 2.8 | **0.4** | 1.4 |
| neos-3627168-kasai | **10.2** | 48.7 | 74.4 |
| neos-4738912-atrato | **69.5** | TL | TL |
| neos-827175 | 6.3 | **2.0** | 2.8 |
| neos-911970 | **1.6** | 4.8 | 6.1 |
| neos-950242 | **1.1** | 19.0 | TL |
| neos-960392 | 13.3 | **4.8** | TL |
| neos17 | **1.4** | 4.7 | 3.5 |
| neos5 | **19.5** | 45.4 | TL |
| neos859080 | **0.3** | 0.7 | 0.3 |
| ns1208400 | **8.3** | 37.7 | TL |
| ns1830653 | **32.4** | TL | 72.0 |
| nu25-pr12 | 2.7 | **1.9** | 2.1 |
| p200x1188c | 3.3 | **0.3** | 2.4 |
| pg | 4.1 | **3.0** | 11.9 |
| pg5_34 | TL | TL | TL |
| piperout-08 | **2.3** | 10.4 | 7.7 |
| piperout-27 | **3.4** | 15.5 | 6.4 |
| pk1 | **5.1** | TL | 57.8 |
| qap10 | 44.0 | **35.5** | 68.3 |
| rmatr100-p10 | **12.7** | 41.8 | 87.8 |
| rocI-4-11 | TL | TL | **46.3** |
| rococoC10-001000 | **65.6** | TL | TL |
| roll3000 | 28.1 | **19.5** | 28.6 |
| seymour1 | **16.1** | 79.0 | 33.7 |
| sp150x300d | 0.1 | **0.0** | 0.2 |
| swath1 | **7.3** | 8.8 | 8.6 |
| swath3 | 89.4 | **60.4** | TL |
| timtab1 | **5.0** | 35.1 | 27.9 |
| unitcal_7 | **61.6** | 102.6 | TL |

</details>

## Hard MIPLIB instances (10-minute limit)

The 63 instances above were selected from a wider 106-instance set (MIPLIB 2017 benchmark set,
"easy", at most about 200k nonzeros). The remaining 43 are the ones that none of popso, HiGHS or
SCIP solved within 120 s. They are listed in
[`benchmarks/sets/miplib_hard.txt`](benchmarks/sets/miplib_hard.txt) and were run with a
600 s limit, otherwise with the same settings as above. All have known optimal values, so every
answer is checked.

![MIPLIB hard instances: instances solved vs. time](docs/miplib_hard.svg)

| solver | solved | shifted geomean (s) |
|---|---|---|
| popso | 18/43 | 281.5 |
| HiGHS | 12/43 | 487.5 |
| SCIP | 6/43 | 564.7 |

Most of these instances stay unsolved by all three solvers at 10 minutes; the plot shows only the
range where anything is solved. HiGHS hit its memory limit on supportcase18; popso stays within
2 GB on all 43. Earlier versions solved 9 (v0.3.0, shifted geometric mean 468.9 s) and 15
(v0.5.0, 324.9 s) of these instances. v0.6.0 solves app1-2, graphdraw-domain, icir97_tension,
neos-3216931-puriri and supportcase18 for the first time, and glass4 takes 22.8 s instead of
534 s. neos-1171737 and neos-1456979, solved by v0.5.0 in 184 s and 301 s, are unsolved at
600 s. Unlike the 63-instance set, this set was not used for tuning until v0.5.0; since then
several of its instances have been development targets.
HiGHS solves four that popso does not (blp-ic98, csched008, neos-1456979, sct2).

<details><summary>Per instance (seconds; TL = time limit; bold = fastest)</summary>

| instance | popso | HiGHS | SCIP |
|---|---|---|---|
| 50v-10 | TL | TL | TL |
| app1-2 | **327.6** | TL | TL |
| blp-ic98 | TL | **571.3** | TL |
| bppc4-08 | TL | TL | TL |
| chromaticindex512-7 | TL | TL | TL |
| cost266-UUE | TL | TL | TL |
| csched007 | **308.6** | 566.6 | TL |
| csched008 | TL | **152.9** | 585.9 |
| fhnw-binpack4-48 | **0.3** | TL | TL |
| gen-ip054 | **138.9** | 395.7 | TL |
| germanrr | TL | TL | TL |
| glass4 | **22.8** | TL | TL |
| gmu-35-40 | TL | TL | TL |
| gmu-35-50 | TL | TL | TL |
| graphdraw-domain | **130.7** | 490.7 | 409.7 |
| ic97_potential | TL | TL | TL |
| icir97_tension | **383.9** | TL | TL |
| leo1 | TL | TL | TL |
| lotsize | TL | TL | TL |
| mas74 | **73.9** | TL | TL |
| mcsched | **83.8** | 248.3 | 149.9 |
| milo-v12-6-r2-40-1 | TL | TL | TL |
| neos-1171737 | TL | TL | TL |
| neos-1456979 | TL | **325.2** | 425.5 |
| neos-3216931-puriri | **112.2** | 191.3 | TL |
| neos-4338804-snowy | TL | TL | TL |
| neos-4387871-tavua | TL | TL | TL |
| neos-4954672-berkel | TL | TL | TL |
| neos-5107597-kakapo | TL | TL | TL |
| neos-933966 | **30.1** | TL | 577.9 |
| nexp-150-20-8-5 | **17.3** | 127.4 | TL |
| ran14x18-disj-8 | **57.2** | 344.1 | TL |
| reblock115 | TL | TL | TL |
| rmatr200-p5 | TL | TL | TL |
| rococoB10-011000 | TL | TL | TL |
| sct2 | TL | **270.6** | TL |
| supportcase18 | **231.7** | TL | TL |
| supportcase26 | **471.7** | TL | TL |
| supportcase40 | **146.0** | TL | TL |
| tr12-30 | 229.9 | **162.7** | 386.4 |
| trento1 | TL | TL | TL |
| uct-subprob | TL | TL | TL |
| wachplan | **418.5** | TL | TL |

</details>

## Full MIPLIB 2017 benchmark set (240 instances)

All 240 instances of the official MIPLIB 2017 benchmark set, without selection
([`benchmarks/sets/miplib_benchmark.txt`](benchmarks/sets/miplib_benchmark.txt)). The reference
values are MIPLIB's own (232 optimal, 7 infeasible, 1 best known), so every answer is checked.
This run used a second machine:

- Machine: Intel Core i7-7700K (4 cores, 8 hardware threads, 64 GB), Linux. Other jobs ran on it
  at times.
- Versions and scoring as above; 120 s per instance, 4 threads (SCIP's branch-and-bound is
  single-threaded), memory cap 8 GB per solver process.
- The HiGHS and SCIP results are in
  [`benchmarks/reference/miplib_benchmark_x86.json`](benchmarks/reference/miplib_benchmark_x86.json).
- popso ran the release candidate. The run found a wrong "infeasible" answer on
  neos-5104907-jarama (a probing bug) and a bug in the solution checker on irish-electricity, both
  fixed in the release; the candidate also crashed once on neos-2657525-crna and ran past the time
  limit on bab2. Those four instances and four more whose result the release could change
  (neos-4387871-tavua, neos-3402454-bohle, square41, square47) were rerun with the release build;
  all eight stay unsolved at 120 s.

![MIPLIB 2017 benchmark set: instances solved vs. time](docs/miplib_benchmark.svg)

| solver | solved | shifted geomean (s) |
|---|---|---|
| popso | 75/240 | 62.5 |
| HiGHS | 63/240 | 76.2 |
| SCIP | 53/240 | 77.5 |

No solver solves 158 of the 240 instances within 2 minutes. popso is fastest on 50 instances
(HiGHS 17, SCIP 15) and solves 10 that neither reference solver does (comp07-2idx, dano3_5,
fhnw-binpack4-48, markshare_4_0, neos-3402294-bobin, neos-933966, netdiversion, nexp-150-20-8-5,
pk1, rococoC10-001000). HiGHS or SCIP solve 7 that popso does not: eil33-2, fhnw-binpack4-4,
neos-4413714-turia, neos-787933, ns1116954, rocI-4-11 and triptim1.

<details><summary>Per instance (seconds; TL = time limit; bold = fastest)</summary>

| instance | popso | HiGHS | SCIP |
|---|---|---|---|
| 30n20b8 | TL | TL | TL |
| 50v-10 | TL | TL | TL |
| CMS750_4 | TL | TL | TL |
| academictimetablesmall | TL | TL | TL |
| air05 | 33.9 | 29.8 | **28.8** |
| app1-1 | 34.0 | 19.1 | **6.7** |
| app1-2 | TL | TL | TL |
| assign1-5-8 | TL | TL | TL |
| atlanta-ip | TL | TL | TL |
| b1c1s1 | TL | TL | TL |
| bab2 | TL | TL | TL |
| bab6 | TL | TL | TL |
| beasleyC3 | 26.9 | 22.3 | **16.2** |
| binkar10_1 | **5.0** | 26.1 | 31.3 |
| blp-ar98 | TL | TL | TL |
| blp-ic98 | TL | TL | TL |
| bnatt400 | TL | TL | TL |
| bnatt500 | TL | TL | TL |
| bppc4-08 | TL | TL | TL |
| brazil3 | TL | TL | TL |
| buildingenergy | TL | TL | TL |
| cbs-cta | **16.2** | 46.6 | 21.3 |
| chromaticindex1024-7 | TL | TL | TL |
| chromaticindex512-7 | TL | TL | TL |
| cmflsp50-24-8-8 | TL | TL | TL |
| co-100 | TL | TL | TL |
| cod105 | **1.7** | 93.4 | 47.4 |
| comp07-2idx | **27.8** | TL | TL |
| comp21-2idx | TL | TL | TL |
| cost266-UUE | TL | TL | TL |
| cryptanalysiskb128n5obj14 | TL | TL | TL |
| cryptanalysiskb128n5obj16 | TL | TL | TL |
| csched007 | TL | TL | TL |
| csched008 | TL | TL | TL |
| cvs16r128-89 | TL | TL | TL |
| dano3_3 | **35.8** | 77.3 | 84.9 |
| dano3_5 | **72.7** | TL | TL |
| decomp2 | **2.1** | 19.0 | 2.1 |
| drayage-100-23 | 6.3 | **0.5** | 7.7 |
| drayage-25-23 | 9.2 | **1.6** | 15.5 |
| dws008-01 | TL | TL | TL |
| eil33-2 | TL | **74.9** | 78.1 |
| eilA101-2 | TL | TL | TL |
| enlight_hard | 4.6 | 11.5 | **0.0** |
| ex10 | **8.5** | 27.9 | TL |
| ex9 | **3.1** | 3.8 | 26.7 |
| exp-1-500-5-5 | **0.5** | 3.3 | 2.4 |
| fast0507 | TL | TL | TL |
| fastxgemm-n2r6s0t2 | TL | TL | TL |
| fhnw-binpack4-4 | TL | TL | **17.4** |
| fhnw-binpack4-48 | **0.5** | TL | TL |
| fiball | **4.3** | 93.2 | TL |
| gen-ip002 | TL | TL | TL |
| gen-ip054 | TL | TL | TL |
| germanrr | TL | TL | TL |
| gfd-schedulen180f7d50m30k18 | TL | TL | TL |
| glass-sc | TL | TL | TL |
| glass4 | TL | TL | TL |
| gmu-35-40 | TL | TL | TL |
| gmu-35-50 | TL | TL | TL |
| graph20-20-1rand | 72.0 | TL | **5.5** |
| graphdraw-domain | TL | TL | TL |
| h80x6320d | 51.5 | **16.6** | TL |
| highschool1-aigio | TL | TL | TL |
| hypothyroid-k1 | 80.9 | 13.9 | **9.3** |
| ic97_potential | TL | TL | TL |
| icir97_tension | TL | TL | TL |
| irish-electricity | TL | TL | TL |
| irp | **5.7** | 22.4 | 20.9 |
| istanbul-no-cutoff | **92.7** | TL | 93.3 |
| k1mushroom | TL | TL | TL |
| lectsched-5-obj | TL | TL | TL |
| leo1 | TL | TL | TL |
| leo2 | TL | TL | TL |
| lotsize | TL | TL | TL |
| mad | TL | TL | TL |
| map10 | TL | TL | TL |
| map16715-04 | TL | TL | TL |
| markshare2 | TL | TL | TL |
| markshare_4_0 | **18.8** | TL | TL |
| mas74 | TL | TL | TL |
| mas76 | **31.6** | 101.8 | 73.7 |
| mc11 | 33.2 | 48.2 | **17.6** |
| mcsched | TL | TL | TL |
| mik-250-20-75-4 | **4.9** | 14.5 | 36.7 |
| milo-v12-6-r2-40-1 | TL | TL | TL |
| momentum1 | TL | TL | TL |
| mushroom-best | TL | TL | TL |
| mzzv11 | **22.2** | 107.0 | TL |
| mzzv42z | **7.8** | 41.9 | TL |
| n2seq36q | **48.5** | 95.2 | TL |
| n3div36 | TL | TL | TL |
| n5-3 | **25.3** | 40.6 | 48.3 |
| neos-1122047 | **3.0** | 4.6 | 4.2 |
| neos-1171448 | 12.3 | TL | **10.7** |
| neos-1171737 | TL | TL | TL |
| neos-1354092 | TL | TL | TL |
| neos-1445765 | **29.1** | 80.5 | 39.1 |
| neos-1456979 | TL | TL | TL |
| neos-1582420 | **24.4** | 44.1 | 28.6 |
| neos-2075418-temuka | TL | TL | TL |
| neos-2657525-crna | TL | TL | TL |
| neos-2746589-doon | TL | TL | TL |
| neos-2978193-inde | **19.0** | 22.2 | TL |
| neos-2987310-joes | 7.5 | **6.0** | 17.5 |
| neos-3004026-krka | **0.8** | TL | 73.0 |
| neos-3024952-loue | TL | TL | TL |
| neos-3046615-murg | TL | TL | TL |
| neos-3083819-nubu | **5.8** | 10.6 | 7.9 |
| neos-3216931-puriri | TL | TL | TL |
| neos-3381206-awhea | 8.7 | 6.8 | **2.4** |
| neos-3402294-bobin | **44.6** | TL | TL |
| neos-3402454-bohle | TL | TL | TL |
| neos-3555904-turama | TL | TL | TL |
| neos-3627168-kasai | **22.4** | 76.5 | TL |
| neos-3656078-kumeu | TL | TL | TL |
| neos-3754480-nidda | TL | TL | TL |
| neos-3988577-wolgan | TL | TL | TL |
| neos-4300652-rahue | TL | TL | TL |
| neos-4338804-snowy | TL | TL | TL |
| neos-4387871-tavua | TL | TL | TL |
| neos-4413714-turia | TL | **27.0** | TL |
| neos-4532248-waihi | TL | TL | TL |
| neos-4647030-tutaki | TL | TL | TL |
| neos-4722843-widden | TL | TL | TL |
| neos-4738912-atrato | TL | TL | TL |
| neos-4763324-toguru | TL | TL | TL |
| neos-4954672-berkel | TL | TL | TL |
| neos-5049753-cuanza | TL | TL | TL |
| neos-5052403-cygnet | TL | TL | TL |
| neos-5093327-huahum | TL | TL | TL |
| neos-5104907-jarama | TL | TL | TL |
| neos-5107597-kakapo | TL | TL | TL |
| neos-5114902-kasavu | TL | TL | TL |
| neos-5188808-nattai | TL | TL | TL |
| neos-5195221-niemur | TL | TL | TL |
| neos-631710 | TL | TL | TL |
| neos-662469 | TL | TL | TL |
| neos-787933 | TL | TL | **1.2** |
| neos-827175 | 12.5 | **3.7** | 4.5 |
| neos-848589 | TL | TL | TL |
| neos-860300 | 67.1 | **15.9** | 23.8 |
| neos-873061 | TL | TL | TL |
| neos-911970 | **5.8** | 24.0 | 11.2 |
| neos-933966 | **22.3** | TL | TL |
| neos-950242 | **2.4** | 22.0 | TL |
| neos-957323 | **56.1** | 90.0 | TL |
| neos-960392 | **12.8** | 13.3 | TL |
| neos17 | **2.7** | 6.9 | 6.6 |
| neos5 | **74.5** | 105.3 | TL |
| neos8 | 5.9 | 4.0 | **1.6** |
| neos859080 | 1.2 | 1.5 | **0.6** |
| net12 | TL | TL | TL |
| netdiversion | **82.4** | TL | TL |
| nexp-150-20-8-5 | **20.9** | TL | TL |
| ns1116954 | TL | **2.5** | TL |
| ns1208400 | **23.6** | 90.7 | TL |
| ns1644855 | TL | TL | TL |
| ns1760995 | TL | TL | TL |
| ns1830653 | **64.5** | TL | 116.8 |
| ns1952667 | TL | TL | TL |
| nu25-pr12 | 8.3 | 5.7 | **4.0** |
| nursesched-medium-hint03 | TL | TL | TL |
| nursesched-sprint02 | 21.8 | **21.2** | 56.4 |
| nw04 | **19.9** | 54.8 | 46.6 |
| opm2-z10-s4 | TL | TL | TL |
| p200x1188c | 1.7 | **0.7** | 4.4 |
| peg-solitaire-a3 | TL | TL | TL |
| pg | 8.1 | **5.3** | 19.6 |
| pg5_34 | TL | TL | TL |
| physiciansched3-3 | TL | TL | TL |
| physiciansched6-2 | 12.0 | **9.8** | TL |
| piperout-08 | **4.5** | 23.9 | 18.0 |
| piperout-27 | **14.0** | 39.5 | 14.0 |
| pk1 | **12.8** | TL | TL |
| proteindesign121hz512p9 | TL | TL | TL |
| proteindesign122trx11p8 | TL | TL | TL |
| qap10 | **52.6** | 111.7 | 103.4 |
| radiationm18-12-05 | TL | TL | TL |
| radiationm40-10-02 | TL | TL | TL |
| rail01 | TL | TL | TL |
| rail02 | TL | TL | TL |
| rail507 | TL | TL | TL |
| ran14x18-disj-8 | TL | TL | TL |
| rd-rplusc-21 | TL | TL | TL |
| reblock115 | TL | TL | TL |
| rmatr100-p10 | **20.5** | 89.7 | TL |
| rmatr200-p5 | TL | TL | TL |
| rocI-4-11 | TL | TL | **101.4** |
| rocII-5-11 | TL | TL | TL |
| rococoB10-011000 | TL | TL | TL |
| rococoC10-001000 | **97.2** | TL | TL |
| roi2alpha3n4 | TL | TL | TL |
| roi5alpha10n8 | TL | TL | TL |
| roll3000 | 63.0 | **41.9** | 48.9 |
| s100 | TL | TL | TL |
| s250r10 | TL | TL | TL |
| satellites2-40 | TL | TL | TL |
| satellites2-60-fs | TL | TL | TL |
| savsched1 | TL | TL | TL |
| sct2 | TL | TL | TL |
| seymour | TL | TL | TL |
| seymour1 | **27.5** | TL | 50.4 |
| sing326 | TL | TL | TL |
| sing44 | TL | TL | TL |
| snp-02-004-104 | TL | TL | TL |
| sorrell3 | TL | TL | TL |
| sp150x300d | 0.1 | **0.1** | 0.3 |
| sp97ar | TL | TL | TL |
| sp98ar | TL | TL | TL |
| splice1k1 | TL | TL | TL |
| square41 | TL | TL | TL |
| square47 | TL | TL | TL |
| supportcase10 | TL | TL | TL |
| supportcase12 | TL | TL | TL |
| supportcase18 | TL | TL | TL |
| supportcase19 | TL | TL | TL |
| supportcase22 | TL | TL | TL |
| supportcase26 | TL | TL | TL |
| supportcase33 | TL | TL | TL |
| supportcase40 | TL | TL | TL |
| supportcase42 | TL | TL | TL |
| supportcase6 | TL | TL | TL |
| supportcase7 | 96.5 | **95.4** | TL |
| swath1 | **12.5** | 16.3 | 17.3 |
| swath3 | TL | TL | TL |
| tbfp-network | TL | TL | TL |
| thor50dday | TL | TL | TL |
| timtab1 | **16.5** | 53.4 | 55.2 |
| tr12-30 | TL | TL | TL |
| traininstance2 | TL | TL | TL |
| traininstance6 | TL | TL | TL |
| trento1 | TL | TL | TL |
| triptim1 | TL | **100.4** | TL |
| uccase12 | TL | TL | TL |
| uccase9 | TL | TL | TL |
| uct-subprob | TL | TL | TL |
| unitcal_7 | TL | TL | TL |
| var-smallemery-m6j6 | TL | TL | TL |
| wachplan | TL | TL | TL |

</details>

## Generated instances

148 models from random families of increasing size
([`benchmarks/instances/`](benchmarks/instances), produced by
[`benchmarks/generators/`](benchmarks/generators)):
- scheduling and lot-sizing (job shop, lot sizing, time-indexed);
- routing and networks (TSP, CVRP, fixed-charge network design, transport);
- combinatorial models (knapsack, multi-dimensional knapsack, set cover and packing, bin
  packing, generalized assignment, subset sum, pigeonhole);
- facility location (CFLP);
- random LPs and MILPs;
- models with logical constraints: minimum lot sizes with semi-continuous variables
  (`semilot`), networks with piecewise-linear costs, plant closing and single sourcing via
  SOS1/SOS2, indicator and semi-continuous constraints (`pwlnet`), and single-machine
  scheduling with indicator disjunctions (`indsched`);
- convex quadratic programs: portfolio optimization (`portfolio`), constrained least squares
  (`lsq`), and the MIQPs cardinality-constrained portfolio (`miqpport`) and random
  convex MIQPs (`miqprand`);
- hand-written LP/MIP edge cases (infeasible, unbounded, degenerate, empty rows, free
  variables, objective sense).

All solvers use 1 thread. HiGHS cannot read the 8 models with SOS or indicator constraints and has
no MIQP solver; its counts are over the other 130 models.

![Generated instances: instances solved vs. time](docs/generated.svg)

| solver | solved | shifted geomean (s) |
|---|---|---|
| popso | 146/148 | 0.8 |
| HiGHS | 125/130 | 1.3 |
| SCIP | 134/148 | 2.0 |

popso leaves gap_hard_1 and binpack_hard_1 unsolved. binpack_hard_1 is a lottery between seeds:
the bound of 18 bins is proved at once, and finding an 18-bin packing within the limit depends on
the seed (seeds 0, 1, 2 at 240 s: v0.5.0 TL / 6.8 s / TL, v0.4.2 23.7 s / TL / TL).

<details><summary>Per instance (seconds; TL = time limit; bold = fastest)</summary>

| instance | popso | HiGHS | SCIP |
|---|---|---|---|
| binpack_hard_1 | TL | 26.8 | **3.1** |
| binpack_large_1 | **0.1** | 0.1 | 0.3 |
| binpack_medium_1 | 0.0 | **0.0** | 0.0 |
| binpack_small_1 | 0.0 | 0.0 | **0.0** |
| binpack_tiny_1 | 0.0 | 0.0 | **0.0** |
| cflp_hard_1 | **40.6** | 47.7 | TL |
| cflp_large_1 | 92.4 | **25.7** | TL |
| cflp_medium_1 | **4.4** | 6.0 | 8.3 |
| cflp_small_1 | **0.1** | 0.8 | 0.6 |
| cflp_tiny_1 | 0.0 | 0.0 | **0.0** |
| cvrp_large_1 | **5.5** | 9.0 | 10.5 |
| cvrp_medium_1 | 1.4 | **1.2** | 2.1 |
| cvrp_small_1 | 0.1 | 0.1 | **0.1** |
| cvrp_tiny_1 | 0.0 | 0.0 | **0.0** |
| gap_hard_1 | TL | TL | TL |
| gap_large_1 | **10.4** | 12.6 | 82.4 |
| gap_medium_1 | **1.6** | 12.1 | 5.1 |
| gap_small_1 | **0.1** | 0.6 | 0.4 |
| gap_tiny_1 | 0.0 | 0.0 | **0.0** |
| indsched_large_1 | 0.3 | - | **0.0** |
| indsched_medium_1 | 0.0 | - | **0.0** |
| indsched_small_1 | 0.0 | - | **0.0** |
| indsched_tiny_1 | 0.0 | - | **0.0** |
| jobshop_hard_1 | **5.4** | 17.2 | 37.1 |
| jobshop_large_1 | **9.5** | 25.4 | 28.2 |
| jobshop_medium_1 | 1.1 | **0.3** | 8.0 |
| jobshop_small_1 | 0.1 | **0.0** | 0.5 |
| jobshop_tiny_1 | 0.0 | 0.0 | **0.0** |
| knapsack_large_1 | **0.0** | 0.1 | 0.3 |
| knapsack_medium_1 | 0.0 | 0.0 | **0.0** |
| knapsack_small_1 | 0.0 | 0.0 | **0.0** |
| knapsack_tiny_1 | 0.0 | 0.0 | **0.0** |
| lotsize_hard_1 | **49.9** | 92.7 | 82.1 |
| lotsize_large_1 | **9.8** | 11.9 | 14.9 |
| lotsize_medium_1 | **6.0** | 6.6 | 6.9 |
| lotsize_small_1 | **0.1** | 0.3 | 0.5 |
| lotsize_tiny_1 | 0.0 | 0.0 | **0.0** |
| lpsanity_beale_cycling_0 | 0.0 | **0.0** | 0.0 |
| lpsanity_degenerate_0 | 0.0 | **0.0** | 0.0 |
| lpsanity_dual_degenerate_0 | 0.0 | 0.0 | **0.0** |
| lpsanity_empty_col_0 | 0.0 | 0.0 | **0.0** |
| lpsanity_empty_row_0 | 0.0 | **0.0** | 0.0 |
| lpsanity_empty_row_infeasible_0 | 0.0 | **0.0** | 0.0 |
| lpsanity_extra_free_row_0 | 0.0 | 0.0 | **0.0** |
| lpsanity_feasible_0 | 0.0 | **0.0** | 0.0 |
| lpsanity_fixed_0 | 0.0 | **0.0** | 0.0 |
| lpsanity_fixed_format_0 | 0.0 | 0.0 | **0.0** |
| lpsanity_freevars_0 | 0.0 | 0.0 | **0.0** |
| lpsanity_illscaled_0 | 0.0 | **0.0** | 0.0 |
| lpsanity_inconsistent_eq_0 | 0.0 | **0.0** | 0.0 |
| lpsanity_infeasible_0 | 0.0 | **0.0** | 0.0 |
| lpsanity_infeasible_bounds_0 | 0.0 | **0.0** | 0.0 |
| lpsanity_klee_minty_0 | 0.0 | **0.0** | 0.0 |
| lpsanity_mip_equality_0 | 0.0 | 0.0 | **0.0** |
| lpsanity_mip_knapsack_0 | 0.0 | 0.0 | **0.0** |
| lpsanity_mip_marker_0 | 0.0 | 0.0 | **0.0** |
| lpsanity_mip_marker_defaults_0 | 0.0 | 0.0 | **0.0** |
| lpsanity_mip_mixed_bounds_0 | 0.0 | 0.0 | **0.0** |
| lpsanity_mip_negative_ub_0 | **0.0** | TL | TL |
| lpsanity_mip_parity_infeasible_0 | 0.0 | 0.0 | **0.0** |
| lpsanity_mip_unbounded_0 | 0.0 | **0.0** | 0.1 |
| lpsanity_negative_lb_0 | 0.0 | **0.0** | 0.0 |
| lpsanity_negative_ub_0 | **0.0** | TL | TL |
| lpsanity_objoffset_0 | 0.0 | 0.0 | **0.0** |
| lpsanity_objsense_inline_0 | 0.0 | **0.0** | 0.0 |
| lpsanity_objsense_maximize_0 | 0.0 | 0.0 | **0.0** |
| lpsanity_ranges_0 | 0.0 | **0.0** | 0.0 |
| lpsanity_redundant_eq_0 | 0.0 | 0.0 | **0.0** |
| lpsanity_unbounded_0 | 0.0 | 0.0 | **0.0** |
| lpsanity_unbounded_free_0 | 0.0 | 0.0 | **0.0** |
| lpsanity_upper_only_0 | 0.0 | 0.0 | **0.0** |
| lsq_large_1 | **1.4** | 2.0 | TL |
| lsq_medium_1 | **0.1** | 0.1 | TL |
| lsq_small_1 | 0.0 | **0.0** | 0.3 |
| lsq_tiny_1 | 0.0 | **0.0** | 0.0 |
| miqpport_hard_1 | **0.9** | - | TL |
| miqpport_large_1 | **0.2** | - | 13.4 |
| miqpport_medium_1 | **0.1** | - | 1.2 |
| miqpport_small_1 | **0.0** | - | 0.6 |
| miqpport_tiny_1 | **0.0** | - | 0.1 |
| miqprand_hard_1 | **2.5** | - | TL |
| miqprand_large_1 | **0.3** | - | TL |
| miqprand_medium_1 | **0.0** | - | 35.6 |
| miqprand_small_1 | **0.0** | - | TL |
| miqprand_tiny_1 | **0.0** | - | 0.1 |
| mkp_hard_1 | **5.6** | 43.9 | 22.6 |
| mkp_large_1 | **1.0** | 4.5 | 3.1 |
| mkp_medium_1 | **0.7** | 2.3 | 4.7 |
| mkp_small_1 | **0.0** | 0.1 | 0.1 |
| mkp_tiny_1 | 0.0 | 0.0 | **0.0** |
| netdesign_hard_1 | **19.0** | 30.9 | 36.4 |
| netdesign_large_1 | **0.8** | 5.4 | 4.1 |
| netdesign_medium_1 | **0.2** | 0.4 | 1.2 |
| netdesign_small_1 | **0.0** | 0.0 | 0.0 |
| netdesign_tiny_1 | 0.0 | 0.0 | **0.0** |
| pigeon_large_1 | 0.0 | **0.0** | 0.2 |
| pigeon_medium_1 | 0.0 | **0.0** | 0.1 |
| pigeon_small_1 | 0.0 | **0.0** | 0.0 |
| pigeon_tiny_1 | 0.0 | 0.0 | **0.0** |
| portfolio_large_1 | **0.0** | 0.0 | 1.2 |
| portfolio_medium_1 | 0.0 | **0.0** | 7.9 |
| portfolio_small_1 | 0.0 | **0.0** | 0.0 |
| portfolio_tiny_1 | 0.0 | **0.0** | 0.0 |
| pwlnet_large_1 | **1.4** | - | 5.1 |
| pwlnet_medium_1 | **0.1** | - | 0.2 |
| pwlnet_small_1 | 0.0 | - | **0.0** |
| pwlnet_tiny_1 | 0.0 | - | **0.0** |
| randlp_large_1 | **19.5** | 39.9 | 65.1 |
| randlp_medium_1 | 0.1 | **0.0** | 0.5 |
| randlp_small_1 | 0.0 | **0.0** | 0.0 |
| randmilp_hard_1 | **54.2** | TL | 71.9 |
| randmilp_large_1 | **12.2** | 25.3 | 21.2 |
| randmilp_medium_1 | **0.2** | 0.2 | 0.2 |
| randmilp_small_1 | 0.0 | 0.0 | **0.0** |
| randmilp_tiny_1 | 0.0 | 0.0 | **0.0** |
| semilot_large_1 | 5.4 | **1.3** | TL |
| semilot_medium_1 | 0.1 | 0.1 | **0.0** |
| semilot_small_1 | 0.0 | 0.0 | **0.0** |
| semilot_tiny_1 | 0.0 | 0.0 | **0.0** |
| setcover_hard_1 | **0.5** | 2.7 | 3.8 |
| setcover_large_1 | **0.7** | 2.4 | 4.2 |
| setcover_medium_1 | **0.1** | 1.1 | 0.7 |
| setcover_small_1 | 0.0 | 0.0 | **0.0** |
| setcover_tiny_1 | 0.0 | 0.0 | **0.0** |
| setpack_hard_1 | **56.6** | TL | TL |
| setpack_large_1 | **3.0** | 9.9 | 8.5 |
| setpack_medium_1 | **0.2** | 1.1 | 1.0 |
| setpack_small_1 | 0.0 | 0.0 | **0.0** |
| setpack_tiny_1 | 0.0 | 0.0 | **0.0** |
| subsetsum_medium_1 | **1.3** | 3.6 | 3.9 |
| subsetsum_small_1 | 0.0 | 0.0 | **0.0** |
| timeidx_hard_1 | 20.0 | **18.0** | 34.3 |
| timeidx_large_1 | **6.1** | 11.9 | 27.0 |
| timeidx_medium_1 | 0.5 | 1.1 | **0.3** |
| timeidx_small_1 | 0.1 | 0.0 | **0.0** |
| timeidx_tiny_1 | 0.0 | 0.0 | **0.0** |
| transinf_medium_1 | 0.0 | 0.0 | **0.0** |
| transinf_small_1 | 0.0 | 0.0 | **0.0** |
| transport_large_1 | **0.1** | 0.4 | 2.7 |
| transport_medium_1 | 0.0 | **0.0** | 0.0 |
| transport_small_1 | 0.0 | **0.0** | 0.0 |
| tsp_hard_1 | **39.9** | 106.1 | TL |
| tsp_large_1 | 98.7 | **68.9** | 115.1 |
| tsp_medium_1 | **3.0** | 6.3 | 3.7 |
| tsp_small_1 | **0.3** | 0.4 | 0.3 |
| tsp_tiny_1 | 0.0 | 0.0 | **0.0** |
| unbdmilp_medium_1 | 0.0 | **0.0** | 0.2 |
| unbdmilp_small_1 | 0.0 | **0.0** | 0.0 |

</details>

## popso over time

The same 63 MIPLIB instances and settings, for every released version (each measured once,
partly with other jobs running, so differences of a few percent are noise).

![popso versions: instances solved vs. time](docs/versions.svg)

| version | solved | shifted geomean (s) |
|---|---|---|
| v0.1.0 | 59/63 | 13.2 |
| v0.2.0 | 60/63 | 12.4 |
| v0.2.1 | 60/63 | 11.3 |
| v0.3.0 | 60/63 | 11.4 |
| v0.4.0 | 59/63 | 10.5 |
| v0.5.0 | 58/63 | 10.7 |
| v0.6.0 | 59/63 | 9.0 |

v0.3.0 adds new features (barrier LP, logical constraints) rather than MILP speed, so its
MIPLIB numbers match v0.2.1. v0.4.0 (new: QP / MIQP, better c-MIR aggregation, symmetry-breaking
reductions) is a few percent faster over three seeds; it solves one instance fewer in this run
because two instances close to the 120 s limit just miss it (both finish in 113 s and 157 s
with a longer limit). v0.5.0 (root strong branching, network cut-set cuts, a local-search
heuristic and earlier root heuristics, conflict analysis, presolve of long equations and implied
integers) is about 4% faster than v0.4.0 over three seeds (shifted geometric mean of all 189 runs
10.5 s against 10.9 s, 178 against 176 solved); on this set its single seed-0 run is no faster.
Most of its gains are on the hard instances. v0.6.0 (repeated root restarts, probing liftings and
sparsify in presolve, cheaper roots, budgets for the heuristics' LPs, a column-scaled cost
perturbation in the dual simplex, compact open-node storage) is about 14% faster than v0.5.0 over
three seeds (shifted geometric mean of all 189 runs 9.0 s against 10.5 s, 179 against 178 solved)
and solves 18 of the 43 hard instances against 15. The plot shows the last four versions.

## Commercial solvers

This page compares popso with open-source solvers only, so that anyone can reproduce every
number. Commercial MILP solvers are still clearly faster on MIPLIB; for independent
comparisons that include commercial solvers see Hans Mittelmann's
[benchmarks for optimization software](https://plato.asu.edu/bench.html).

## Reproduce

```sh
python3 -m venv .venv && .venv/bin/pip install -r benchmarks/requirements.txt
.venv/bin/python benchmarks/miplib/download.py --wide          # MIPLIB files (not in the repo)
.venv/bin/python benchmarks/run.py --solver build/popso --set miplib_wide_core -t 120 --threads 4 --mem-limit 2048
.venv/bin/python benchmarks/run.py --solver build/popso --set miplib_hard -t 600 --threads 4 --mem-limit 2048
.venv/bin/python benchmarks/miplib/download.py --benchmark     # all 240 (~640 MB)
.venv/bin/python benchmarks/run.py --solver build/popso --set miplib_benchmark -t 120 --threads 4 --mem-limit 8192 \
    --reference benchmarks/reference/miplib_benchmark_x86.json
.venv/bin/python benchmarks/run.py --solver build/popso -t 120 --threads 1 --mem-limit 2048
# tables and plots from a result file
.venv/bin/python benchmarks/tools/readme_table.py benchmarks/results/<run>.json --solvers popso,highs,scip
.venv/bin/python benchmarks/tools/plot_profile.py benchmarks/results/<run>.json --title "..." --out docs/<name>.svg
```

The reference results (HiGHS, SCIP) are produced by `benchmarks/reference.py`; see
[benchmarks/README.md](benchmarks/README.md).
