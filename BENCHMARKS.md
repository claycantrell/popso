# Benchmarks

How the solver compares with the open-source solvers [HiGHS](https://highs.dev) and
[SCIP](https://www.scipopt.org). Every answer is checked independently: the returned
solution must be feasible, and its objective must match the known optimum.

_Last updated 2026-11-07 for popso v0.7.0._

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
  - The v0.7.0 numbers were measured with the release code (only the version string differs).
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
| popso | 61/63 | 9.4 |
| HiGHS | 53/63 | 18.9 |
| SCIP | 50/63 | 21.3 |

popso is fastest on 41 of the 63 instances. The table is seed 0. Over seeds 0, 1 and 2 popso
solves 61 / 60 / 63 instances with a shifted geometric mean of 9.4 / 8.7 / 7.9 s (v0.6.0: 59 / 61
/ 59 and 9.0 / 9.2 / 9.0 s); over all 189 runs 184 solved and 8.7 s against 179 and 9.0 s. Single runs
of the instances close to the time limit vary a lot between seeds.

<details><summary>Per instance (seconds; TL = time limit; bold = fastest)</summary>

| instance | popso | HiGHS | SCIP |
|---|---|---|---|
| 30n20b8 | **59.5** | TL | 116.6 |
| CMS750_4 | 114.3 | TL | **77.0** |
| air05 | **15.8** | 16.8 | 21.8 |
| app1-1 | 24.1 | 14.7 | **3.8** |
| beasleyC3 | 18.2 | 12.3 | **9.1** |
| binkar10_1 | **2.6** | 23.0 | 16.6 |
| cbs-cta | **2.9** | 31.7 | 13.5 |
| cod105 | **1.0** | 62.9 | 32.7 |
| dano3_3 | **22.4** | 38.4 | 58.8 |
| dano3_5 | **47.5** | TL | TL |
| decomp2 | **0.3** | 2.1 | 1.0 |
| drayage-100-23 | 2.3 | **0.3** | 3.9 |
| drayage-25-23 | 19.5 | **0.9** | 8.7 |
| eil33-2 | 82.6 | 54.2 | **47.8** |
| enlight_hard | 1.4 | 6.4 | **0.0** |
| exp-1-500-5-5 | **0.3** | 2.5 | 1.3 |
| fiball | **2.2** | 57.1 | TL |
| h80x6320d | 34.4 | **8.5** | 71.9 |
| irp | **2.7** | 13.8 | 10.8 |
| mas76 | **12.6** | 59.2 | 34.9 |
| mc11 | 21.7 | 22.9 | **9.6** |
| mik-250-20-75-4 | **3.0** | 8.6 | 17.1 |
| mzzv11 | **13.0** | 44.4 | TL |
| mzzv42z | **3.8** | 26.1 | 69.3 |
| n2seq36q | **22.6** | 50.7 | 78.2 |
| n5-3 | 13.2 | **11.0** | 29.4 |
| neos-1122047 | **1.3** | 2.9 | 1.8 |
| neos-1171448 | **5.3** | 79.4 | 6.0 |
| neos-1445765 | **17.3** | 43.3 | 20.0 |
| neos-1582420 | **6.1** | 14.5 | 16.9 |
| neos-2978193-inde | **11.5** | 21.5 | TL |
| neos-3004026-krka | **0.4** | TL | 35.1 |
| neos-3083819-nubu | **2.6** | 5.3 | 3.9 |
| neos-3381206-awhea | 3.0 | **0.4** | 1.4 |
| neos-3627168-kasai | **21.3** | 48.7 | 74.4 |
| neos-4738912-atrato | **67.2** | TL | TL |
| neos-827175 | 4.5 | **2.0** | 2.8 |
| neos-911970 | **2.6** | 4.8 | 6.1 |
| neos-950242 | **1.3** | 19.0 | TL |
| neos-960392 | 7.9 | **4.8** | TL |
| neos17 | **1.0** | 4.7 | 3.5 |
| neos5 | **20.6** | 45.4 | TL |
| neos859080 | **0.3** | 0.7 | 0.3 |
| ns1208400 | **7.4** | 37.7 | TL |
| ns1830653 | **38.8** | TL | 72.0 |
| nu25-pr12 | 9.2 | **1.9** | 2.1 |
| p200x1188c | 4.5 | **0.3** | 2.4 |
| pg | 4.8 | **3.0** | 11.9 |
| pg5_34 | TL | TL | TL |
| piperout-08 | **2.1** | 10.4 | 7.7 |
| piperout-27 | **6.2** | 15.5 | 6.4 |
| pk1 | **8.5** | TL | 57.8 |
| qap10 | 37.9 | **35.5** | 68.3 |
| rmatr100-p10 | **12.2** | 41.8 | 87.8 |
| rocI-4-11 | TL | TL | **46.3** |
| rococoC10-001000 | **60.9** | TL | TL |
| roll3000 | 24.2 | **19.5** | 28.6 |
| seymour1 | **16.8** | 79.0 | 33.7 |
| sp150x300d | 0.1 | **0.0** | 0.2 |
| swath1 | **7.6** | 8.8 | 8.6 |
| swath3 | 88.7 | **60.4** | TL |
| timtab1 | **8.3** | 35.1 | 27.9 |
| unitcal_7 | **48.6** | 102.6 | TL |

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
| popso | 16/43 | 298.5 |
| HiGHS | 12/43 | 487.5 |
| SCIP | 6/43 | 564.7 |

Most of these instances stay unsolved by all three solvers at 10 minutes. HiGHS hit its memory limit on
supportcase18; popso stays within 2 GB on all 43. Earlier versions solved 9 (v0.3.0, shifted geometric
mean 468.9 s), 15 (v0.5.0, 324.9 s) and 18 (v0.6.0, 281.5 s) of these instances. v0.7.0 solves
neos-1171737 and neos-1456979 again (163 s and 170 s), but not app1-2, csched007, icir97_tension and
supportcase26, which v0.6.0 solved in 308-472 s in its run. Those four depend on the seed: with seed 1,
v0.6.3 and v0.7.0 each solve one of them (csched007, 260 s and 243 s) and none of the other three.
Unlike the 63-instance set, this set was not used for tuning until v0.5.0; since then several of its
instances have been development targets. popso solves 7 that neither reference solver does
(fhnw-binpack4-48, glass4, mas74, neos-1171737, supportcase18, supportcase40, wachplan); HiGHS solves
four that popso does not (blp-ic98, csched007, csched008, sct2), SCIP one (csched008).

<details><summary>Per instance (seconds; TL = time limit; bold = fastest)</summary>

| instance | popso | HiGHS | SCIP |
|---|---|---|---|
| 50v-10 | TL | TL | TL |
| app1-2 | TL | TL | TL |
| blp-ic98 | TL | **571.3** | TL |
| bppc4-08 | TL | TL | TL |
| chromaticindex512-7 | TL | TL | TL |
| cost266-UUE | TL | TL | TL |
| csched007 | TL | **566.6** | TL |
| csched008 | TL | **152.9** | 585.9 |
| fhnw-binpack4-48 | **0.1** | TL | TL |
| gen-ip054 | 434.7 | **395.7** | TL |
| germanrr | TL | TL | TL |
| glass4 | **48.9** | TL | TL |
| gmu-35-40 | TL | TL | TL |
| gmu-35-50 | TL | TL | TL |
| graphdraw-domain | **205.5** | 490.7 | 409.7 |
| ic97_potential | TL | TL | TL |
| icir97_tension | TL | TL | TL |
| leo1 | TL | TL | TL |
| lotsize | TL | TL | TL |
| mas74 | **105.0** | TL | TL |
| mcsched | **80.5** | 248.3 | 149.9 |
| milo-v12-6-r2-40-1 | TL | TL | TL |
| neos-1171737 | **163.5** | TL | TL |
| neos-1456979 | **169.6** | 325.2 | 425.5 |
| neos-3216931-puriri | **94.5** | 191.3 | TL |
| neos-4338804-snowy | TL | TL | TL |
| neos-4387871-tavua | TL | TL | TL |
| neos-4954672-berkel | TL | TL | TL |
| neos-5107597-kakapo | TL | TL | TL |
| neos-933966 | **40.6** | TL | 577.9 |
| nexp-150-20-8-5 | **13.1** | 127.4 | TL |
| ran14x18-disj-8 | **67.4** | 344.1 | TL |
| reblock115 | TL | TL | TL |
| rmatr200-p5 | TL | TL | TL |
| rococoB10-011000 | TL | TL | TL |
| sct2 | TL | **270.6** | TL |
| supportcase18 | **569.8** | TL | TL |
| supportcase26 | TL | TL | TL |
| supportcase40 | **147.3** | TL | TL |
| tr12-30 | 220.2 | **162.7** | 386.4 |
| trento1 | TL | TL | TL |
| uct-subprob | TL | TL | TL |
| wachplan | **334.2** | TL | TL |

</details>

## Full MIPLIB 2017 benchmark set (240 instances)

All 240 instances of the official MIPLIB 2017 benchmark set, without selection
([`benchmarks/sets/miplib_benchmark.txt`](benchmarks/sets/miplib_benchmark.txt)). The reference
values are MIPLIB's own (232 optimal, 7 infeasible, 1 best known), so every answer is checked.
This run used a second machine:

- Machine: AMD EPYC 4465P (12 cores, 24 hardware threads, 62 GB), Linux; six solver runs at a time.
- Versions and scoring as above; 120 s per instance, 4 threads (SCIP's branch-and-bound is
  single-threaded), memory cap 8 GB per solver process.
- The HiGHS and SCIP results were measured on the same machine with the same settings, in
  [`benchmarks/reference/miplib_benchmark_epyc.json`](benchmarks/reference/miplib_benchmark_epyc.json).
  v0.6.0's run used an Intel i7-7700K (75 solved, HiGHS 63, SCIP 53 there), so the two are not directly
  comparable.
- popso ran the v0.7.0 release code. neos-5104907-jarama went past the memory cap (counted unsolved).

![MIPLIB 2017 benchmark set: instances solved vs. time](docs/miplib_benchmark.svg)

| solver | solved | shifted geomean (s) |
|---|---|---|
| popso | 89/240 | 56.0 |
| HiGHS | 74/240 | 66.3 |
| SCIP | 65/240 | 69.1 |

No solver solves 144 of the 240 instances within 2 minutes. popso is fastest on 53 instances
(HiGHS 25, SCIP 18) and solves 11 that neither reference solver does (30n20b8, comp07-2idx, dano3_5,
fhnw-binpack4-48, mcsched, mushroom-best, neos-4738912-atrato, neos-933966, netdiversion,
nexp-150-20-8-5, ran14x18-disj-8). HiGHS or SCIP solve 7 that popso does not: csched008, fast0507,
neos-4722843-widden, rail507, rocI-4-11, supportcase33 and unitcal_7.

<details><summary>Per instance (seconds; TL = time limit; bold = fastest)</summary>

| instance | popso | HiGHS | SCIP |
|---|---|---|---|
| 30n20b8 | **93.5** | TL | TL |
| 50v-10 | TL | TL | TL |
| CMS750_4 | **40.3** | TL | 93.5 |
| academictimetablesmall | TL | TL | TL |
| air05 | 20.3 | **17.0** | 19.3 |
| app1-1 | 40.2 | 10.8 | **4.2** |
| app1-2 | TL | TL | TL |
| assign1-5-8 | TL | TL | TL |
| atlanta-ip | TL | TL | TL |
| b1c1s1 | TL | TL | TL |
| bab2 | TL | TL | TL |
| bab6 | TL | TL | TL |
| beasleyC3 | 24.5 | 13.1 | **10.4** |
| binkar10_1 | **3.7** | 15.9 | 18.9 |
| blp-ar98 | TL | TL | TL |
| blp-ic98 | TL | TL | TL |
| bnatt400 | TL | TL | TL |
| bnatt500 | TL | TL | TL |
| bppc4-08 | TL | TL | TL |
| brazil3 | TL | TL | TL |
| buildingenergy | TL | TL | TL |
| cbs-cta | **3.8** | 29.6 | 14.9 |
| chromaticindex1024-7 | TL | TL | TL |
| chromaticindex512-7 | TL | TL | TL |
| cmflsp50-24-8-8 | TL | TL | TL |
| co-100 | TL | TL | TL |
| cod105 | **1.5** | 58.3 | 33.3 |
| comp07-2idx | **9.7** | TL | TL |
| comp21-2idx | TL | TL | TL |
| cost266-UUE | TL | TL | TL |
| cryptanalysiskb128n5obj14 | TL | TL | TL |
| cryptanalysiskb128n5obj16 | TL | TL | TL |
| csched007 | TL | TL | TL |
| csched008 | TL | **97.8** | TL |
| cvs16r128-89 | TL | TL | TL |
| dano3_3 | **31.4** | 45.4 | 56.0 |
| dano3_5 | **68.6** | TL | TL |
| decomp2 | **0.5** | 11.4 | 0.9 |
| drayage-100-23 | 4.4 | **0.3** | 3.9 |
| drayage-25-23 | 11.3 | **0.9** | 8.9 |
| dws008-01 | TL | TL | TL |
| eil33-2 | 94.7 | **48.6** | 52.1 |
| eilA101-2 | TL | TL | TL |
| enlight_hard | 1.6 | 7.1 | **0.0** |
| ex10 | **9.0** | 22.4 | TL |
| ex9 | 3.0 | **2.6** | 19.0 |
| exp-1-500-5-5 | **0.5** | 2.1 | 1.5 |
| fast0507 | TL | TL | **116.0** |
| fastxgemm-n2r6s0t2 | TL | TL | TL |
| fhnw-binpack4-4 | **7.1** | TL | 10.2 |
| fhnw-binpack4-48 | **0.2** | TL | TL |
| fiball | **3.7** | 51.4 | TL |
| gen-ip002 | TL | TL | TL |
| gen-ip054 | TL | TL | TL |
| germanrr | TL | TL | TL |
| gfd-schedulen180f7d50m30k18 | TL | TL | TL |
| glass-sc | TL | TL | TL |
| glass4 | TL | TL | TL |
| gmu-35-40 | TL | TL | TL |
| gmu-35-50 | TL | TL | TL |
| graph20-20-1rand | 98.3 | TL | **3.3** |
| graphdraw-domain | TL | TL | TL |
| h80x6320d | 41.2 | **9.2** | 67.4 |
| highschool1-aigio | TL | TL | TL |
| hypothyroid-k1 | 11.6 | 6.4 | **3.7** |
| ic97_potential | TL | TL | TL |
| icir97_tension | TL | TL | TL |
| irish-electricity | TL | TL | TL |
| irp | **5.6** | 12.7 | 10.7 |
| istanbul-no-cutoff | 93.0 | TL | **56.7** |
| k1mushroom | TL | TL | TL |
| lectsched-5-obj | TL | TL | TL |
| leo1 | TL | TL | TL |
| leo2 | TL | TL | TL |
| lotsize | TL | TL | TL |
| mad | TL | TL | TL |
| map10 | TL | TL | TL |
| map16715-04 | TL | TL | TL |
| markshare2 | TL | TL | TL |
| markshare_4_0 | **17.4** | 77.6 | 71.2 |
| mas74 | TL | TL | TL |
| mas76 | **17.1** | 56.0 | 40.1 |
| mc11 | 38.5 | 28.0 | **10.4** |
| mcsched | **115.3** | TL | TL |
| mik-250-20-75-4 | **5.4** | 8.3 | 22.2 |
| milo-v12-6-r2-40-1 | TL | TL | TL |
| momentum1 | TL | TL | TL |
| mushroom-best | **118.8** | TL | TL |
| mzzv11 | **22.3** | 63.9 | TL |
| mzzv42z | **7.6** | 25.8 | 73.7 |
| n2seq36q | **32.3** | 47.6 | 76.1 |
| n3div36 | TL | TL | TL |
| n5-3 | **19.9** | 20.3 | 28.8 |
| neos-1122047 | 2.9 | **2.4** | 2.7 |
| neos-1171448 | 10.3 | 110.0 | **6.3** |
| neos-1171737 | TL | TL | TL |
| neos-1354092 | TL | TL | TL |
| neos-1445765 | 26.8 | 47.5 | **21.5** |
| neos-1456979 | TL | TL | TL |
| neos-1582420 | **11.6** | 23.2 | 16.5 |
| neos-2075418-temuka | TL | TL | TL |
| neos-2657525-crna | TL | TL | TL |
| neos-2746589-doon | TL | TL | TL |
| neos-2978193-inde | **9.7** | 11.4 | TL |
| neos-2987310-joes | 7.3 | **3.3** | 11.7 |
| neos-3004026-krka | **0.9** | TL | 42.3 |
| neos-3024952-loue | TL | TL | TL |
| neos-3046615-murg | TL | TL | TL |
| neos-3083819-nubu | **1.6** | 6.0 | 4.2 |
| neos-3216931-puriri | TL | TL | TL |
| neos-3381206-awhea | 3.9 | 3.6 | **1.3** |
| neos-3402294-bobin | **35.4** | 109.8 | TL |
| neos-3402454-bohle | TL | TL | TL |
| neos-3555904-turama | TL | TL | TL |
| neos-3627168-kasai | **21.2** | 46.4 | 89.2 |
| neos-3656078-kumeu | TL | TL | TL |
| neos-3754480-nidda | TL | TL | TL |
| neos-3988577-wolgan | TL | TL | TL |
| neos-4300652-rahue | TL | TL | TL |
| neos-4338804-snowy | TL | TL | TL |
| neos-4387871-tavua | TL | TL | TL |
| neos-4413714-turia | 89.4 | **13.0** | TL |
| neos-4532248-waihi | TL | TL | TL |
| neos-4647030-tutaki | TL | TL | TL |
| neos-4722843-widden | TL | **82.9** | TL |
| neos-4738912-atrato | **72.4** | TL | TL |
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
| neos-787933 | 3.0 | TL | **0.9** |
| neos-827175 | 8.1 | **2.1** | 2.8 |
| neos-848589 | TL | TL | TL |
| neos-860300 | 64.4 | **9.0** | 16.0 |
| neos-873061 | TL | TL | TL |
| neos-911970 | **5.2** | 15.0 | 7.6 |
| neos-933966 | **59.7** | TL | TL |
| neos-950242 | **3.1** | 12.5 | TL |
| neos-957323 | **28.0** | 43.4 | 97.6 |
| neos-960392 | 10.3 | **6.8** | TL |
| neos17 | **1.8** | 4.0 | 3.7 |
| neos5 | **53.7** | 61.2 | TL |
| neos8 | 13.1 | 2.1 | **0.8** |
| neos859080 | **0.3** | 0.8 | 0.4 |
| net12 | TL | TL | TL |
| netdiversion | **88.0** | TL | TL |
| nexp-150-20-8-5 | **18.1** | TL | TL |
| ns1116954 | 11.7 | **1.3** | TL |
| ns1208400 | **12.4** | 48.3 | TL |
| ns1644855 | TL | TL | TL |
| ns1760995 | TL | TL | TL |
| ns1830653 | **73.4** | TL | 79.9 |
| ns1952667 | TL | TL | TL |
| nu25-pr12 | 13.2 | 3.5 | **2.5** |
| nursesched-medium-hint03 | TL | TL | TL |
| nursesched-sprint02 | 15.2 | **13.9** | 36.2 |
| nw04 | **16.1** | 31.5 | 27.9 |
| opm2-z10-s4 | TL | TL | TL |
| p200x1188c | 7.3 | **0.4** | 2.7 |
| peg-solitaire-a3 | TL | TL | TL |
| pg | 7.0 | **3.2** | 12.0 |
| pg5_34 | TL | TL | TL |
| physiciansched3-3 | TL | TL | TL |
| physiciansched6-2 | 10.1 | **5.6** | 100.5 |
| piperout-08 | **4.3** | 16.7 | 8.9 |
| piperout-27 | 15.7 | 27.1 | **7.5** |
| pk1 | **12.2** | TL | 68.3 |
| proteindesign121hz512p9 | TL | TL | TL |
| proteindesign122trx11p8 | TL | TL | TL |
| qap10 | **49.1** | 64.0 | 65.3 |
| radiationm18-12-05 | TL | TL | TL |
| radiationm40-10-02 | TL | TL | TL |
| rail01 | TL | TL | TL |
| rail02 | TL | TL | TL |
| rail507 | TL | TL | **109.3** |
| ran14x18-disj-8 | **106.4** | TL | TL |
| rd-rplusc-21 | TL | TL | TL |
| reblock115 | TL | TL | TL |
| rmatr100-p10 | **17.4** | 46.3 | 87.1 |
| rmatr200-p5 | TL | TL | TL |
| rocI-4-11 | TL | TL | **58.8** |
| rocII-5-11 | TL | TL | TL |
| rococoB10-011000 | TL | TL | TL |
| rococoC10-001000 | **84.1** | 115.4 | TL |
| roi2alpha3n4 | TL | TL | TL |
| roi5alpha10n8 | TL | TL | TL |
| roll3000 | 51.3 | **25.8** | 30.8 |
| s100 | TL | TL | TL |
| s250r10 | TL | TL | TL |
| satellites2-40 | TL | TL | TL |
| satellites2-60-fs | TL | TL | TL |
| savsched1 | TL | TL | TL |
| sct2 | TL | TL | TL |
| seymour | TL | TL | TL |
| seymour1 | **21.8** | 75.7 | 29.6 |
| sing326 | TL | TL | TL |
| sing44 | TL | TL | TL |
| snp-02-004-104 | TL | TL | TL |
| sorrell3 | TL | TL | TL |
| sp150x300d | 0.2 | **0.0** | 0.2 |
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
| supportcase33 | TL | **91.8** | TL |
| supportcase40 | TL | TL | TL |
| supportcase42 | TL | TL | TL |
| supportcase6 | TL | TL | TL |
| supportcase7 | 65.0 | **52.1** | TL |
| swath1 | 10.8 | 9.9 | **9.5** |
| swath3 | **49.3** | 78.7 | TL |
| tbfp-network | **75.2** | 114.8 | TL |
| thor50dday | TL | TL | TL |
| timtab1 | **13.9** | 33.9 | 35.0 |
| tr12-30 | TL | TL | TL |
| traininstance2 | TL | TL | TL |
| traininstance6 | TL | TL | TL |
| trento1 | TL | TL | TL |
| triptim1 | 82.8 | **64.2** | TL |
| uccase12 | TL | TL | TL |
| uccase9 | TL | TL | TL |
| uct-subprob | TL | TL | TL |
| unitcal_7 | TL | **82.8** | TL |
| var-smallemery-m6j6 | TL | TL | TL |
| wachplan | TL | TL | TL |

</details>

## Generated instances

153 models from random families of increasing size
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
  (`lsq`), and the MIQPs cardinality-constrained portfolio (`miqpport`), random
  convex MIQPs (`miqprand`) and production targets in squares of integers (`sqtrack`);
- hand-written LP/MIP edge cases (infeasible, unbounded, degenerate, empty rows, free
  variables, objective sense).

All solvers use 1 thread. HiGHS cannot read the 8 models with SOS or indicator constraints and has
no MIQP solver; its counts are over the other 130 models. Three MIQPs have no proven reference
optimum; popso's optimal answers there are consistent with the known bounds and count as solved.

![Generated instances: instances solved vs. time](docs/generated.svg)

| solver | solved | shifted geomean (s) |
|---|---|---|
| popso | 152/153 | 0.8 |
| HiGHS | 125/130 | 1.3 |
| SCIP | 139/153 | 2.0 |

popso leaves gap_hard_1 unsolved (v0.6.0 also binpack_hard_1, now 10.5 s; that one depends on the seed:
the bound of 18 bins is proved at once, and finding an 18-bin packing within the limit is a matter of
luck). The `sqtrack` models take 0.02-0.14 s each (the squares as a MILP with tangent rows; the MIQP
search did not finish the large and hard ones in 60 s).

<details><summary>Per instance (seconds; TL = time limit; bold = fastest)</summary>

| instance | popso | HiGHS | SCIP |
|---|---|---|---|
| binpack_hard_1 | 10.5 | 26.8 | **3.1** |
| binpack_large_1 | **0.0** | 0.1 | 0.3 |
| binpack_medium_1 | 0.0 | **0.0** | 0.0 |
| binpack_small_1 | 0.0 | 0.0 | **0.0** |
| binpack_tiny_1 | 0.0 | 0.0 | **0.0** |
| cflp_hard_1 | 55.2 | **47.7** | TL |
| cflp_large_1 | 34.5 | **25.7** | TL |
| cflp_medium_1 | 7.0 | **6.0** | 8.3 |
| cflp_small_1 | **0.2** | 0.8 | 0.6 |
| cflp_tiny_1 | 0.0 | 0.0 | **0.0** |
| cvrp_large_1 | **4.8** | 9.0 | 10.5 |
| cvrp_medium_1 | 1.8 | **1.2** | 2.1 |
| cvrp_small_1 | 0.2 | 0.1 | **0.1** |
| cvrp_tiny_1 | 0.0 | 0.0 | **0.0** |
| gap_hard_1 | TL | TL | TL |
| gap_large_1 | **8.1** | 12.6 | 82.4 |
| gap_medium_1 | **1.6** | 12.1 | 5.1 |
| gap_small_1 | **0.2** | 0.6 | 0.4 |
| gap_tiny_1 | 0.0 | 0.0 | **0.0** |
| indsched_large_1 | 0.3 | - | **0.0** |
| indsched_medium_1 | 0.0 | - | **0.0** |
| indsched_small_1 | 0.0 | - | **0.0** |
| indsched_tiny_1 | 0.0 | - | **0.0** |
| jobshop_hard_1 | **5.3** | 17.2 | 37.1 |
| jobshop_large_1 | **14.8** | 25.4 | 28.2 |
| jobshop_medium_1 | 2.0 | **0.3** | 8.0 |
| jobshop_small_1 | 0.1 | **0.0** | 0.5 |
| jobshop_tiny_1 | 0.0 | 0.0 | **0.0** |
| knapsack_large_1 | **0.0** | 0.1 | 0.3 |
| knapsack_medium_1 | 0.0 | 0.0 | **0.0** |
| knapsack_small_1 | 0.0 | 0.0 | **0.0** |
| knapsack_tiny_1 | 0.0 | 0.0 | **0.0** |
| lotsize_hard_1 | **43.6** | 92.7 | 82.1 |
| lotsize_large_1 | **10.6** | 11.9 | 14.9 |
| lotsize_medium_1 | **5.5** | 6.6 | 6.9 |
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
| lsq_large_1 | **1.3** | 2.0 | TL |
| lsq_medium_1 | **0.1** | 0.1 | TL |
| lsq_small_1 | 0.0 | **0.0** | 0.3 |
| lsq_tiny_1 | 0.0 | **0.0** | 0.0 |
| miqpport_hard_1 | **0.8** | - | TL |
| miqpport_large_1 | **0.2** | - | 13.4 |
| miqpport_medium_1 | **0.1** | - | 1.2 |
| miqpport_small_1 | **0.0** | - | 0.6 |
| miqpport_tiny_1 | **0.0** | - | 0.1 |
| miqprand_hard_1 | **2.4** | - | TL |
| miqprand_large_1 | **0.3** | - | TL |
| miqprand_medium_1 | **0.0** | - | 35.6 |
| miqprand_small_1 | **0.0** | - | TL |
| miqprand_tiny_1 | **0.0** | - | 0.1 |
| mkp_hard_1 | **5.9** | 43.9 | 22.6 |
| mkp_large_1 | **1.0** | 4.5 | 3.1 |
| mkp_medium_1 | **0.8** | 2.3 | 4.7 |
| mkp_small_1 | **0.0** | 0.1 | 0.1 |
| mkp_tiny_1 | 0.0 | 0.0 | **0.0** |
| netdesign_hard_1 | **18.9** | 30.9 | 36.4 |
| netdesign_large_1 | **0.6** | 5.4 | 4.1 |
| netdesign_medium_1 | **0.2** | 0.4 | 1.2 |
| netdesign_small_1 | 0.0 | 0.0 | **0.0** |
| netdesign_tiny_1 | 0.0 | 0.0 | **0.0** |
| pigeon_large_1 | 0.0 | **0.0** | 0.2 |
| pigeon_medium_1 | 0.0 | **0.0** | 0.1 |
| pigeon_small_1 | 0.0 | **0.0** | 0.0 |
| pigeon_tiny_1 | 0.0 | 0.0 | **0.0** |
| portfolio_large_1 | 0.0 | **0.0** | 1.2 |
| portfolio_medium_1 | 0.0 | **0.0** | 7.9 |
| portfolio_small_1 | 0.0 | **0.0** | 0.0 |
| portfolio_tiny_1 | 0.0 | **0.0** | 0.0 |
| pwlnet_large_1 | **1.5** | - | 5.1 |
| pwlnet_medium_1 | **0.1** | - | 0.2 |
| pwlnet_small_1 | 0.1 | - | **0.0** |
| pwlnet_tiny_1 | 0.0 | - | **0.0** |
| randlp_large_1 | **19.6** | 39.9 | 65.1 |
| randlp_medium_1 | 0.1 | **0.0** | 0.5 |
| randlp_small_1 | 0.0 | **0.0** | 0.0 |
| randmilp_hard_1 | **62.5** | TL | 71.9 |
| randmilp_large_1 | **10.9** | 25.3 | 21.2 |
| randmilp_medium_1 | **0.1** | 0.2 | 0.2 |
| randmilp_small_1 | 0.0 | 0.0 | **0.0** |
| randmilp_tiny_1 | 0.0 | 0.0 | **0.0** |
| semilot_large_1 | 5.1 | **1.3** | TL |
| semilot_medium_1 | 0.1 | 0.1 | **0.0** |
| semilot_small_1 | 0.0 | 0.0 | **0.0** |
| semilot_tiny_1 | 0.0 | 0.0 | **0.0** |
| setcover_hard_1 | **0.5** | 2.7 | 3.8 |
| setcover_large_1 | **0.8** | 2.4 | 4.2 |
| setcover_medium_1 | **0.1** | 1.1 | 0.7 |
| setcover_small_1 | 0.0 | 0.0 | **0.0** |
| setcover_tiny_1 | 0.0 | 0.0 | **0.0** |
| setpack_hard_1 | **56.6** | TL | TL |
| setpack_large_1 | **3.0** | 9.9 | 8.5 |
| setpack_medium_1 | **0.2** | 1.1 | 1.0 |
| setpack_small_1 | 0.0 | 0.0 | **0.0** |
| setpack_tiny_1 | 0.0 | 0.0 | **0.0** |
| sqtrack_hard_1 | **0.1** | - | 1.8 |
| sqtrack_large_1 | **0.0** | - | 0.3 |
| sqtrack_medium_1 | **0.0** | - | 0.0 |
| sqtrack_small_1 | **0.0** | - | 1.7 |
| sqtrack_tiny_1 | 0.0 | - | **0.0** |
| subsetsum_medium_1 | **1.2** | 3.6 | 3.9 |
| subsetsum_small_1 | 0.0 | 0.0 | **0.0** |
| timeidx_hard_1 | 19.6 | **18.0** | 34.3 |
| timeidx_large_1 | **7.7** | 11.9 | 27.0 |
| timeidx_medium_1 | 0.6 | 1.1 | **0.3** |
| timeidx_small_1 | 0.1 | 0.0 | **0.0** |
| timeidx_tiny_1 | 0.0 | 0.0 | **0.0** |
| transinf_medium_1 | 0.0 | 0.0 | **0.0** |
| transinf_small_1 | 0.0 | 0.0 | **0.0** |
| transport_large_1 | **0.1** | 0.4 | 2.7 |
| transport_medium_1 | 0.0 | **0.0** | 0.0 |
| transport_small_1 | 0.0 | **0.0** | 0.0 |
| tsp_hard_1 | **61.7** | 106.1 | TL |
| tsp_large_1 | 98.4 | **68.9** | 115.1 |
| tsp_medium_1 | 4.3 | 6.3 | **3.7** |
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
| v0.7.0 | 61/63 | 9.4 |

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
and solves 18 of the 43 hard instances against 15. v0.7.0 (objective propagation with clique
implications, a root cut loop that backs off on wide models, faster clique and c-MIR separation, the trivial
objective bound, presolve substitution with implied bounds, a concurrent barrier on large roots) solves
two more in this seed-0 run at a slightly higher geometric mean; over three seeds it is about 4% faster
than v0.6.0 (all 189 runs: 8.7 s against 9.0 s, 184 against 179 solved). The plot shows the last four
versions.

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
