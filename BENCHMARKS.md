# Benchmarks

How the solver compares with the open-source solvers [HiGHS](https://highs.dev) and
[SCIP](https://www.scipopt.org). Every answer is checked independently: the returned
solution must be feasible, and its objective must match the known optimum.

_Last updated 2026-11-02 for popso v0.3.0._

## Setup

- Machine: Apple M4 Pro (12 cores, 24 GB), macOS.
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
unbiased sample, and popso has been tuned on it. A run on the full 240-instance MIPLIB 2017
benchmark set is planned.

![MIPLIB: instances solved vs. time](docs/miplib_core.svg)

| solver | solved | shifted geomean (s) |
|---|---|---|
| popso | 59/63 | 10.5 |
| HiGHS | 53/63 | 18.9 |
| SCIP | 50/63 | 21.3 |

popso is fastest on 38 of the 63 instances. Over three seeds the shifted geometric mean is
10.5 / 10.9 / 11.4 s (v0.3.0: 11.4 / 11.4 / 11.1 s); single runs of the instances close to the
time limit vary a lot between seeds.

<details><summary>Per instance (seconds; TL = time limit; bold = fastest)</summary>

| instance | popso | HiGHS | SCIP |
|---|---|---|---|
| 30n20b8 | TL | TL | **116.6** |
| CMS750_4 | TL | TL | **77.0** |
| air05 | **16.7** | 16.8 | 21.8 |
| app1-1 | 9.6 | 14.7 | **3.8** |
| beasleyC3 | 25.6 | 12.3 | **9.1** |
| binkar10_1 | **2.6** | 23.0 | 16.6 |
| cbs-cta | **6.4** | 31.7 | 13.5 |
| cod105 | **0.8** | 62.9 | 32.7 |
| dano3_3 | **34.2** | 38.4 | 58.8 |
| dano3_5 | **43.9** | TL | TL |
| decomp2 | 1.4 | 2.1 | **1.0** |
| drayage-100-23 | 4.2 | **0.3** | 3.9 |
| drayage-25-23 | 7.9 | **0.9** | 8.7 |
| eil33-2 | 70.6 | 54.2 | **47.8** |
| enlight_hard | 0.9 | 6.4 | **0.0** |
| exp-1-500-5-5 | 1.8 | 2.5 | **1.3** |
| fiball | **14.3** | 57.1 | TL |
| h80x6320d | 28.7 | **8.5** | 71.9 |
| irp | **5.7** | 13.8 | 10.8 |
| mas76 | **7.9** | 59.2 | 34.9 |
| mc11 | 27.5 | 22.9 | **9.6** |
| mik-250-20-75-4 | **3.9** | 8.6 | 17.1 |
| mzzv11 | **13.4** | 44.4 | TL |
| mzzv42z | **3.2** | 26.1 | 69.3 |
| n2seq36q | **22.3** | 50.7 | 78.2 |
| n5-3 | **10.1** | 11.0 | 29.4 |
| neos-1122047 | **1.1** | 2.9 | 1.8 |
| neos-1171448 | **4.4** | 79.4 | 6.0 |
| neos-1445765 | 21.0 | 43.3 | **20.0** |
| neos-1582420 | **5.2** | 14.5 | 16.9 |
| neos-2978193-inde | **16.8** | 21.5 | TL |
| neos-3004026-krka | **1.0** | TL | 35.1 |
| neos-3083819-nubu | 6.6 | 5.3 | **3.9** |
| neos-3381206-awhea | 1.8 | **0.4** | 1.4 |
| neos-3627168-kasai | **14.3** | 48.7 | 74.4 |
| neos-4738912-atrato | **71.5** | TL | TL |
| neos-827175 | 4.6 | **2.0** | 2.8 |
| neos-911970 | **1.8** | 4.8 | 6.1 |
| neos-950242 | **7.1** | 19.0 | TL |
| neos-960392 | 18.3 | **4.8** | TL |
| neos17 | **0.7** | 4.7 | 3.5 |
| neos5 | **27.9** | 45.4 | TL |
| neos859080 | **0.1** | 0.7 | 0.3 |
| ns1208400 | **25.5** | 37.7 | TL |
| ns1830653 | **71.0** | TL | 72.0 |
| nu25-pr12 | 4.1 | **1.9** | 2.1 |
| p200x1188c | 1.8 | **0.3** | 2.4 |
| pg | **2.7** | 3.0 | 11.9 |
| pg5_34 | **103.1** | TL | TL |
| piperout-08 | **2.6** | 10.4 | 7.7 |
| piperout-27 | 8.3 | 15.5 | **6.4** |
| pk1 | **4.5** | TL | 57.8 |
| qap10 | 49.2 | **35.5** | 68.3 |
| rmatr100-p10 | **10.3** | 41.8 | 87.8 |
| rocI-4-11 | TL | TL | **46.3** |
| rococoC10-001000 | **61.0** | TL | TL |
| roll3000 | 32.5 | **19.5** | 28.6 |
| seymour1 | **16.5** | 79.0 | 33.7 |
| sp150x300d | 0.1 | **0.0** | 0.2 |
| swath1 | **5.8** | 8.8 | 8.6 |
| swath3 | TL | **60.4** | TL |
| timtab1 | **6.3** | 35.1 | 27.9 |
| unitcal_7 | **73.4** | 102.6 | TL |

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
| popso | 9/43 | 468.9 |
| HiGHS | 12/43 | 487.5 |
| SCIP | 6/43 | 564.7 |

Most of these instances stay unsolved by all three solvers at 10 minutes; the plot shows only the
range where anything is solved. popso hit its 2 GB memory limit on graphdraw-domain, HiGHS on
supportcase18. popso has not been tuned on this set.

<details><summary>Per instance (seconds; TL = time limit; bold = fastest)</summary>

| instance | popso | HiGHS | SCIP |
|---|---|---|---|
| 50v-10 | TL | TL | TL |
| app1-2 | TL | TL | TL |
| blp-ic98 | TL | **571.3** | TL |
| bppc4-08 | TL | TL | TL |
| chromaticindex512-7 | TL | TL | TL |
| cost266-UUE | TL | TL | TL |
| csched007 | **326.2** | 566.6 | TL |
| csched008 | TL | **152.9** | 585.9 |
| fhnw-binpack4-48 | TL | TL | TL |
| gen-ip054 | 541.3 | **395.7** | TL |
| germanrr | TL | TL | TL |
| glass4 | **355.8** | TL | TL |
| gmu-35-40 | TL | TL | TL |
| gmu-35-50 | TL | TL | TL |
| graphdraw-domain | TL | 490.7 | **409.7** |
| ic97_potential | TL | TL | TL |
| icir97_tension | TL | TL | TL |
| leo1 | TL | TL | TL |
| lotsize | TL | TL | TL |
| mas74 | **101.1** | TL | TL |
| mcsched | **118.3** | 248.3 | 149.9 |
| milo-v12-6-r2-40-1 | TL | TL | TL |
| neos-1171737 | TL | TL | TL |
| neos-1456979 | TL | **325.2** | 425.5 |
| neos-3216931-puriri | TL | **191.3** | TL |
| neos-4338804-snowy | TL | TL | TL |
| neos-4387871-tavua | TL | TL | TL |
| neos-4954672-berkel | TL | TL | TL |
| neos-5107597-kakapo | TL | TL | TL |
| neos-933966 | **50.0** | TL | 577.9 |
| nexp-150-20-8-5 | TL | **127.4** | TL |
| ran14x18-disj-8 | **67.4** | 344.1 | TL |
| reblock115 | TL | TL | TL |
| rmatr200-p5 | TL | TL | TL |
| rococoB10-011000 | TL | TL | TL |
| sct2 | TL | **270.6** | TL |
| supportcase18 | TL | TL | TL |
| supportcase26 | TL | TL | TL |
| supportcase40 | **176.6** | TL | TL |
| tr12-30 | TL | **162.7** | 386.4 |
| trento1 | TL | TL | TL |
| uct-subprob | TL | TL | TL |
| wachplan | **542.4** | TL | TL |

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
| popso | 147/148 | 0.8 |
| HiGHS | 125/130 | 1.3 |
| SCIP | 134/148 | 2.0 |

<details><summary>Per instance (seconds; TL = time limit; bold = fastest)</summary>

| instance | popso | HiGHS | SCIP |
|---|---|---|---|
| binpack_hard_1 | 6.8 | 26.8 | **3.1** |
| binpack_large_1 | 0.3 | **0.1** | 0.3 |
| binpack_medium_1 | 0.0 | **0.0** | 0.0 |
| binpack_small_1 | 0.0 | 0.0 | **0.0** |
| binpack_tiny_1 | 0.0 | 0.0 | **0.0** |
| cflp_hard_1 | 64.9 | **47.7** | TL |
| cflp_large_1 | 63.8 | **25.7** | TL |
| cflp_medium_1 | 9.0 | **6.0** | 8.3 |
| cflp_small_1 | **0.2** | 0.8 | 0.6 |
| cflp_tiny_1 | 0.0 | 0.0 | **0.0** |
| cvrp_large_1 | **7.5** | 9.0 | 10.5 |
| cvrp_medium_1 | **1.1** | 1.2 | 2.1 |
| cvrp_small_1 | 0.1 | 0.1 | **0.1** |
| cvrp_tiny_1 | 0.0 | 0.0 | **0.0** |
| gap_hard_1 | TL | TL | TL |
| gap_large_1 | 17.5 | **12.6** | 82.4 |
| gap_medium_1 | **1.0** | 12.1 | 5.1 |
| gap_small_1 | **0.1** | 0.6 | 0.4 |
| gap_tiny_1 | 0.0 | 0.0 | **0.0** |
| indsched_large_1 | 0.3 | - | **0.0** |
| indsched_medium_1 | **0.0** | - | 0.0 |
| indsched_small_1 | **0.0** | - | 0.0 |
| indsched_tiny_1 | 0.0 | - | **0.0** |
| jobshop_hard_1 | 25.9 | **17.2** | 37.1 |
| jobshop_large_1 | **16.6** | 25.4 | 28.2 |
| jobshop_medium_1 | 1.2 | **0.3** | 8.0 |
| jobshop_small_1 | 0.1 | **0.0** | 0.5 |
| jobshop_tiny_1 | 0.0 | 0.0 | **0.0** |
| knapsack_large_1 | **0.0** | 0.1 | 0.3 |
| knapsack_medium_1 | 0.0 | 0.0 | **0.0** |
| knapsack_small_1 | 0.0 | 0.0 | **0.0** |
| knapsack_tiny_1 | 0.0 | 0.0 | **0.0** |
| lotsize_hard_1 | **58.9** | 92.7 | 82.1 |
| lotsize_large_1 | 14.4 | **11.9** | 14.9 |
| lotsize_medium_1 | 7.7 | **6.6** | 6.9 |
| lotsize_small_1 | **0.2** | 0.3 | 0.5 |
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
| lpsanity_ranges_0 | 0.1 | **0.0** | 0.0 |
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
| miqpport_medium_1 | **0.0** | - | 1.2 |
| miqpport_small_1 | **0.0** | - | 0.6 |
| miqpport_tiny_1 | **0.0** | - | 0.1 |
| miqprand_hard_1 | **2.4** | - | TL |
| miqprand_large_1 | **0.3** | - | TL |
| miqprand_medium_1 | **0.0** | - | 35.6 |
| miqprand_small_1 | **0.0** | - | TL |
| miqprand_tiny_1 | **0.0** | - | 0.1 |
| mkp_hard_1 | **5.3** | 43.9 | 22.6 |
| mkp_large_1 | **1.0** | 4.5 | 3.1 |
| mkp_medium_1 | **0.5** | 2.3 | 4.7 |
| mkp_small_1 | **0.0** | 0.1 | 0.1 |
| mkp_tiny_1 | 0.0 | 0.0 | **0.0** |
| netdesign_hard_1 | **16.0** | 30.9 | 36.4 |
| netdesign_large_1 | **1.9** | 5.4 | 4.1 |
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
| pwlnet_large_1 | **3.8** | - | 5.1 |
| pwlnet_medium_1 | 0.3 | - | **0.2** |
| pwlnet_small_1 | 0.0 | - | **0.0** |
| pwlnet_tiny_1 | 0.0 | - | **0.0** |
| randlp_large_1 | **17.6** | 39.9 | 65.1 |
| randlp_medium_1 | 0.1 | **0.0** | 0.5 |
| randlp_small_1 | 0.0 | **0.0** | 0.0 |
| randmilp_hard_1 | 84.8 | TL | **71.9** |
| randmilp_large_1 | **13.0** | 25.3 | 21.2 |
| randmilp_medium_1 | **0.1** | 0.2 | 0.2 |
| randmilp_small_1 | 0.0 | 0.0 | **0.0** |
| randmilp_tiny_1 | 0.0 | 0.0 | **0.0** |
| semilot_large_1 | 3.4 | **1.3** | TL |
| semilot_medium_1 | 0.1 | 0.1 | **0.0** |
| semilot_small_1 | 0.1 | 0.0 | **0.0** |
| semilot_tiny_1 | 0.0 | 0.0 | **0.0** |
| setcover_hard_1 | **0.5** | 2.7 | 3.8 |
| setcover_large_1 | **0.8** | 2.4 | 4.2 |
| setcover_medium_1 | **0.1** | 1.1 | 0.7 |
| setcover_small_1 | 0.0 | 0.0 | **0.0** |
| setcover_tiny_1 | 0.0 | 0.0 | **0.0** |
| setpack_hard_1 | **44.4** | TL | TL |
| setpack_large_1 | **2.6** | 9.9 | 8.5 |
| setpack_medium_1 | **0.2** | 1.1 | 1.0 |
| setpack_small_1 | 0.0 | 0.0 | **0.0** |
| setpack_tiny_1 | 0.0 | 0.0 | **0.0** |
| subsetsum_medium_1 | **0.9** | 3.6 | 3.9 |
| subsetsum_small_1 | **0.0** | 0.0 | 0.0 |
| timeidx_hard_1 | **17.1** | 18.0 | 34.3 |
| timeidx_large_1 | **5.2** | 11.9 | 27.0 |
| timeidx_medium_1 | 0.5 | 1.1 | **0.3** |
| timeidx_small_1 | 0.1 | 0.0 | **0.0** |
| timeidx_tiny_1 | 0.0 | 0.0 | **0.0** |
| transinf_medium_1 | 0.0 | 0.0 | **0.0** |
| transinf_small_1 | 0.0 | 0.0 | **0.0** |
| transport_large_1 | **0.1** | 0.4 | 2.7 |
| transport_medium_1 | 0.0 | **0.0** | 0.0 |
| transport_small_1 | 0.0 | **0.0** | 0.0 |
| tsp_hard_1 | **47.8** | 106.1 | TL |
| tsp_large_1 | **54.8** | 68.9 | 115.1 |
| tsp_medium_1 | **2.7** | 6.3 | 3.7 |
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

v0.3.0 adds new features (barrier LP, logical constraints) rather than MILP speed, so its
MIPLIB numbers match v0.2.1. v0.4.0 (new: QP / MIQP, better c-MIR aggregation, symmetry-breaking
reductions) is a few percent faster over three seeds; it solves one instance fewer in this run
because two instances close to the 120 s limit just miss it (both finish in 113 s and 157 s
with a longer limit). The plot leaves out v0.2.0 (12.4 s, 60/63) for readability.

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
.venv/bin/python benchmarks/run.py --solver build/popso -t 120 --threads 1 --mem-limit 2048
# tables and plots from a result file
.venv/bin/python benchmarks/tools/readme_table.py benchmarks/results/<run>.json --solvers popso,highs,scip
.venv/bin/python benchmarks/tools/plot_profile.py benchmarks/results/<run>.json --title "..." --out docs/<name>.svg
```

The reference results (HiGHS, SCIP) are produced by `benchmarks/reference.py`; see
[benchmarks/README.md](benchmarks/README.md).
