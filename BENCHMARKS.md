# Benchmarks

How the solver compares with the open-source solvers [HiGHS](https://highs.dev) and
[SCIP](https://www.scipopt.org). Every answer is checked independently: the returned
solution must be feasible, and its objective must match the known optimum.

_Last updated 2026-11-02 for popso v0.2.0._

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
| popso | 60/63 | 12.4 |
| HiGHS | 53/63 | 18.9 |
| SCIP | 50/63 | 21.3 |

popso is fastest on 28 of the 63 instances.

<details><summary>Per instance (seconds; TL = time limit; bold = fastest)</summary>

| instance | popso | HiGHS | SCIP |
|---|---|---|---|
| 30n20b8 | 119.4 | TL | **116.6** |
| CMS750_4 | **58.5** | TL | 77.0 |
| air05 | 20.0 | **16.8** | 21.8 |
| app1-1 | 8.6 | 14.7 | **3.8** |
| beasleyC3 | 22.7 | 12.3 | **9.1** |
| binkar10_1 | **2.6** | 23.0 | 16.6 |
| cbs-cta | **7.3** | 31.7 | 13.5 |
| cod105 | **12.7** | 62.9 | 32.7 |
| dano3_3 | 38.7 | **38.4** | 58.8 |
| dano3_5 | **61.1** | TL | TL |
| decomp2 | 1.4 | 2.1 | **1.0** |
| drayage-100-23 | 4.1 | **0.3** | 3.9 |
| drayage-25-23 | 9.4 | **0.9** | 8.7 |
| eil33-2 | 103.4 | 54.2 | **47.8** |
| enlight_hard | 1.0 | 6.4 | **0.0** |
| exp-1-500-5-5 | 2.1 | 2.5 | **1.3** |
| fiball | 95.5 | **57.1** | TL |
| h80x6320d | 27.2 | **8.5** | 71.9 |
| irp | **4.8** | 13.8 | 10.8 |
| mas76 | **7.8** | 59.2 | 34.9 |
| mc11 | 25.7 | 22.9 | **9.6** |
| mik-250-20-75-4 | **1.4** | 8.6 | 17.1 |
| mzzv11 | 58.1 | **44.4** | TL |
| mzzv42z | **11.2** | 26.1 | 69.3 |
| n2seq36q | **24.9** | 50.7 | 78.2 |
| n5-3 | 13.8 | **11.0** | 29.4 |
| neos-1122047 | **1.2** | 2.9 | 1.8 |
| neos-1171448 | **4.4** | 79.4 | 6.0 |
| neos-1445765 | 21.4 | 43.3 | **20.0** |
| neos-1582420 | **6.3** | 14.5 | 16.9 |
| neos-2978193-inde | 43.5 | **21.5** | TL |
| neos-3004026-krka | **1.2** | TL | 35.1 |
| neos-3083819-nubu | 5.4 | 5.3 | **3.9** |
| neos-3381206-awhea | 2.1 | **0.4** | 1.4 |
| neos-3627168-kasai | **13.7** | 48.7 | 74.4 |
| neos-4738912-atrato | **42.5** | TL | TL |
| neos-827175 | 4.9 | **2.0** | 2.8 |
| neos-911970 | **1.7** | 4.8 | 6.1 |
| neos-950242 | **6.0** | 19.0 | TL |
| neos-960392 | 14.2 | **4.8** | TL |
| neos17 | 4.7 | 4.7 | **3.5** |
| neos5 | **29.3** | 45.4 | TL |
| neos859080 | 0.4 | 0.7 | **0.3** |
| ns1208400 | **20.3** | 37.7 | TL |
| ns1830653 | 78.5 | TL | **72.0** |
| nu25-pr12 | 4.7 | **1.9** | 2.1 |
| p200x1188c | 1.8 | **0.3** | 2.4 |
| pg | **2.8** | 3.0 | 11.9 |
| pg5_34 | TL | TL | TL |
| piperout-08 | **2.6** | 10.4 | 7.7 |
| piperout-27 | 9.2 | 15.5 | **6.4** |
| pk1 | **5.9** | TL | 57.8 |
| qap10 | 48.3 | **35.5** | 68.3 |
| rmatr100-p10 | **10.7** | 41.8 | 87.8 |
| rocI-4-11 | TL | TL | **46.3** |
| rococoC10-001000 | **48.0** | TL | TL |
| roll3000 | 59.6 | **19.5** | 28.6 |
| seymour1 | **14.0** | 79.0 | 33.7 |
| sp150x300d | 0.2 | **0.0** | 0.2 |
| swath1 | 9.6 | 8.8 | **8.6** |
| swath3 | TL | **60.4** | TL |
| timtab1 | **13.4** | 35.1 | 27.9 |
| unitcal_7 | **50.8** | 102.6 | TL |

</details>

## Generated instances

118 models from random families of increasing size
([`benchmarks/instances/`](benchmarks/instances), produced by
[`benchmarks/generators/`](benchmarks/generators)):
- scheduling and lot-sizing (job shop, lot sizing, time-indexed);
- routing and networks (TSP, CVRP, fixed-charge network design, transport);
- combinatorial models (knapsack, multi-dimensional knapsack, set cover and packing, bin
  packing, generalized assignment, subset sum, pigeonhole);
- facility location (CFLP);
- random LPs and MILPs;
- hand-written LP/MIP edge cases (infeasible, unbounded, degenerate, empty rows, free
  variables, objective sense).

All solvers use 1 thread.

![Generated instances: instances solved vs. time](docs/generated.svg)

| solver | solved | shifted geomean (s) |
|---|---|---|
| popso | 116/118 | 1.0 |
| HiGHS | 113/118 | 1.5 |
| SCIP | 111/118 | 1.7 |

<details><summary>Per instance (seconds; TL = time limit; bold = fastest)</summary>

| instance | popso | HiGHS | SCIP |
|---|---|---|---|
| binpack_hard_1 | TL | 26.8 | **3.1** |
| binpack_large_1 | 0.4 | **0.1** | 0.3 |
| binpack_medium_1 | 0.0 | **0.0** | 0.0 |
| binpack_small_1 | 0.0 | 0.0 | **0.0** |
| binpack_tiny_1 | 0.0 | 0.0 | **0.0** |
| cflp_hard_1 | 64.2 | **47.7** | TL |
| cflp_large_1 | 63.8 | **25.7** | TL |
| cflp_medium_1 | 9.1 | **6.0** | 8.3 |
| cflp_small_1 | **0.2** | 0.8 | 0.6 |
| cflp_tiny_1 | 0.0 | 0.0 | **0.0** |
| cvrp_large_1 | **3.6** | 9.0 | 10.5 |
| cvrp_medium_1 | 1.2 | **1.2** | 2.1 |
| cvrp_small_1 | 0.1 | 0.1 | **0.1** |
| cvrp_tiny_1 | 0.0 | 0.0 | **0.0** |
| gap_hard_1 | TL | TL | TL |
| gap_large_1 | 23.8 | **12.6** | 82.4 |
| gap_medium_1 | **1.0** | 12.1 | 5.1 |
| gap_small_1 | **0.1** | 0.6 | 0.4 |
| gap_tiny_1 | 0.0 | 0.0 | **0.0** |
| jobshop_hard_1 | 21.2 | **17.2** | 37.1 |
| jobshop_large_1 | **14.8** | 25.4 | 28.2 |
| jobshop_medium_1 | 0.7 | **0.3** | 8.0 |
| jobshop_small_1 | 0.1 | **0.0** | 0.5 |
| jobshop_tiny_1 | **0.0** | 0.0 | 0.0 |
| knapsack_large_1 | **0.0** | 0.1 | 0.3 |
| knapsack_medium_1 | 0.0 | 0.0 | **0.0** |
| knapsack_small_1 | 0.0 | 0.0 | **0.0** |
| knapsack_tiny_1 | 0.0 | 0.0 | **0.0** |
| lotsize_hard_1 | **57.0** | 92.7 | 82.1 |
| lotsize_large_1 | 16.4 | **11.9** | 14.9 |
| lotsize_medium_1 | 7.0 | **6.6** | 6.9 |
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
| lpsanity_mip_unbounded_0 | **0.0** | 0.0 | 0.1 |
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
| mkp_hard_1 | **5.2** | 43.9 | 22.6 |
| mkp_large_1 | **1.0** | 4.5 | 3.1 |
| mkp_medium_1 | **0.5** | 2.3 | 4.7 |
| mkp_small_1 | **0.0** | 0.1 | 0.1 |
| mkp_tiny_1 | 0.0 | 0.0 | **0.0** |
| netdesign_hard_1 | **21.5** | 30.9 | 36.4 |
| netdesign_large_1 | **1.9** | 5.4 | 4.1 |
| netdesign_medium_1 | **0.1** | 0.4 | 1.2 |
| netdesign_small_1 | 0.1 | 0.0 | **0.0** |
| netdesign_tiny_1 | 0.0 | 0.0 | **0.0** |
| pigeon_large_1 | 0.0 | **0.0** | 0.2 |
| pigeon_medium_1 | **0.0** | 0.0 | 0.1 |
| pigeon_small_1 | 0.0 | **0.0** | 0.0 |
| pigeon_tiny_1 | 0.0 | 0.0 | **0.0** |
| randlp_large_1 | **17.7** | 39.9 | 65.1 |
| randlp_medium_1 | 0.0 | **0.0** | 0.5 |
| randlp_small_1 | 0.0 | **0.0** | 0.0 |
| randmilp_hard_1 | 89.6 | TL | **71.9** |
| randmilp_large_1 | **13.1** | 25.3 | 21.2 |
| randmilp_medium_1 | **0.2** | 0.2 | 0.2 |
| randmilp_small_1 | 0.0 | 0.0 | **0.0** |
| randmilp_tiny_1 | 0.0 | 0.0 | **0.0** |
| setcover_hard_1 | **0.5** | 2.7 | 3.8 |
| setcover_large_1 | **0.7** | 2.4 | 4.2 |
| setcover_medium_1 | **0.1** | 1.1 | 0.7 |
| setcover_small_1 | 0.0 | 0.0 | **0.0** |
| setcover_tiny_1 | 0.0 | 0.0 | **0.0** |
| setpack_hard_1 | **45.2** | TL | TL |
| setpack_large_1 | **2.4** | 9.9 | 8.5 |
| setpack_medium_1 | **0.2** | 1.1 | 1.0 |
| setpack_small_1 | 0.0 | 0.0 | **0.0** |
| setpack_tiny_1 | 0.0 | 0.0 | **0.0** |
| subsetsum_medium_1 | **1.0** | 3.6 | 3.9 |
| subsetsum_small_1 | **0.0** | 0.0 | 0.0 |
| timeidx_hard_1 | **17.7** | 18.0 | 34.3 |
| timeidx_large_1 | **5.3** | 11.9 | 27.0 |
| timeidx_medium_1 | 0.5 | 1.1 | **0.3** |
| timeidx_small_1 | 0.1 | 0.0 | **0.0** |
| timeidx_tiny_1 | 0.0 | 0.0 | **0.0** |
| transinf_medium_1 | 0.0 | 0.0 | **0.0** |
| transinf_small_1 | 0.0 | 0.0 | **0.0** |
| transport_large_1 | **0.1** | 0.4 | 2.7 |
| transport_medium_1 | 0.0 | **0.0** | 0.0 |
| transport_small_1 | 0.0 | **0.0** | 0.0 |
| tsp_hard_1 | **64.6** | 106.1 | TL |
| tsp_large_1 | **47.9** | 68.9 | 115.1 |
| tsp_medium_1 | **2.5** | 6.3 | 3.7 |
| tsp_small_1 | **0.2** | 0.4 | 0.3 |
| tsp_tiny_1 | 0.0 | 0.0 | **0.0** |
| unbdmilp_medium_1 | 0.0 | **0.0** | 0.2 |
| unbdmilp_small_1 | 0.0 | **0.0** | 0.0 |

</details>

## Reproduce

```sh
python3 -m venv .venv && .venv/bin/pip install -r benchmarks/requirements.txt
.venv/bin/python benchmarks/miplib/download.py --wide          # MIPLIB files (not in the repo)
.venv/bin/python benchmarks/run.py --solver build/popso --set miplib_wide_core -t 120 --threads 4 --mem-limit 2048
.venv/bin/python benchmarks/run.py --solver build/popso -t 120 --threads 1 --mem-limit 2048
# tables and plots from a result file
.venv/bin/python benchmarks/tools/readme_table.py benchmarks/results/<run>.json --solvers popso,highs,scip
.venv/bin/python benchmarks/tools/plot_profile.py benchmarks/results/<run>.json --title "..." --out docs/<name>.svg
```

The reference results (HiGHS, SCIP) are produced by `benchmarks/reference.py`; see
[benchmarks/README.md](benchmarks/README.md).
