# Changelog

## v0.6.2 — 2026-11

A correctness fix on top of v0.6.1. The benchmark numbers of v0.6.0 still apply.

**Fixes**
- Wrong optimal answer after two root restarts: when the first restarted search found nothing better
  than the incumbent it started from, the second restart was handed that start in the previous
  problem's columns; the second restarted root was pruned against the mismatched point and the old
  incumbent came back as optimal (MIPLIB drayage-25-23 on x86: 101344.27 reported optimal, optimum
  101282.65). The restart now passes the start in the restarted problem's columns, and a start whose
  size does not match the problem is ignored with a warning.

## v0.6.1 — 2026-11

Correctness and time-limit fixes on top of v0.6.0. The benchmark numbers of v0.6.0 still apply.

**Fixes**
- LP row scaling: entries below 1e-9 of a row's largest no longer set its geometric scale. A data-noise
  coefficient (MIPLIB var-smallemery-m6j6: 8.9e-16 next to 106.75) scaled the row by 2^22, so the dual
  simplex's 1e-7 tolerance on its logical accepted an unscaled reduced cost of 0.138: the root LP
  stopped at -139.80 instead of -156.06, and a run claimed OPTIMAL at -141.53 (optimum -149.38).
- Node tableau c-MIR cuts dropped tableau entries below 1e-11 of the row's largest without their
  bounds, which is not a relaxation where a column has a huge range (MIPLIB neos-4413714-turia: a
  budget row with right-hand side 9e8; cuts violated by the optimum, a reported bound 46.47 above the
  optimum 45.37). A tiny entry now leaves only with its bound.
- Conflict analysis could read conflicts already freed by the bounded per-worker conflict pool: adding
  a learned conflict could drop conflicts the node's bound trail still referred to (MIPLIB
  neos-2657525-crna: a crash in 1 of 24 runs at 4 threads; none in 48 runs now).
- Root cut loop: when the root cuts do not move the bound and are removed, the LP continues from the
  initial LP's optimal basis. Removing a cut row whose slack was nonbasic dropped the whole basis and
  the LP was re-solved from scratch (MIPLIB triptim1: that re-solve ran for the rest of a 300 s run).
- Time limits: root dives are capped by the time left when each dive starts, not when it is queued
  (`--deterministic`, 2 threads: ns1116954 ended at 235 s with a 120 s limit, now at 122 s), and
  feasibility jump checks the time limit before publishing each point (savsched1 ran 6 s past it).
- Deterministic runs give different node counts than v0.6.0 where the changed cuts or LP scaling apply.

## v0.6.0 — 2026-11

MILP: presolve sparsify, probing liftings over repeated passes with kept implications, repeated root
restarts, cheaper roots, budgets for the heuristics' LPs, a column-scaled cost perturbation in the dual
simplex, faster Forrest-Tomlin updates and compact open-node storage. Time limits are kept on very large
models.

**Benchmarks**
- Full MIPLIB 2017 benchmark set, 240 instances (Intel i7-7700K, Linux, 120 s, 4 threads), first run:
  75 solved, shifted geometric mean 62.5 s (HiGHS 63 / 76.2 s, SCIP 53 / 77.5 s); 10 solved by popso
  only, 7 by HiGHS or SCIP but not popso.
- 63-instance MIPLIB 2017 subset (120 s, 4 threads), seeds 0 / 1 / 2: 59 / 61 / 59 solved, shifted
  geometric mean 9.0 / 9.2 / 9.0 s (v0.5.0: 58 / 59 / 61, 10.7 / 10.8 / 9.9 s); about 14% faster over
  the three seeds (all 189 runs: 9.0 s against 10.5 s, 179 against 178 solved).
- 43 hard MIPLIB 2017 instances (600 s, 4 threads): 18 solved, shifted geometric mean 281.5 s (v0.5.0:
  15, 324.9 s; HiGHS 12, SCIP 6). Newly solved: app1-2, graphdraw-domain, icir97_tension,
  neos-3216931-puriri, supportcase18; glass4 534 -> 23 s. neos-1171737 and neos-1456979 are now
  unsolved at 600 s.
- Generated models (1 thread, 120 s): 146 of 148 solved, as v0.5.0 (binpack_hard_1 and gap_hard_1).
- Measured with the release candidate; the fixes since then are listed below (Time limits on large
  models, Fixes).

**Node-LP incumbents re-solved warm at a tight tolerance**
- When a node LP point is integral but, rounded, violates the model by more than 1e-9, the node LP is
  re-solved from its basis with the integer columns fixed and primal tolerance 1e-9. Its 1e-6 slack could
  grow past the model check after postsolve (dano3_5: 2.7e-6), and the cold fixed-integer re-solves that
  repaired such points could run out of time on large LPs (dano3_5: 30k iterations), so the solve was
  repeated without presolve and ended without a solution. `MILP_POLISH_LP=0` turns it off.
- miplib_wide_core (120 s, 4 threads, seeds 1-2): neutral, 1.026x (90% CI 0.994-1.060), solved 120 -> 121
  of 126.

**Presolve: sparsify**
- Equations are added to the rows of their shortest column, `r += lambda e`, when that cancels more
  entries of r than it adds (two-row nonzero cancellation, as in HiGHS's presolve). Equations of at
  most 50 entries, |lambda| <= 1e3, a row's coefficient range may not grow beyond 1e6; rows of integer
  columns with integral coefficients stay integral; set partitioning equations are not used and set
  packing / covering / partitioning rows are not changed (their cliques). It runs when the other
  reductions are exhausted, before the clique extension and probing (which find more on the sparser
  rows), up to 3 times; the following rounds remove the singleton and doubleton rows it creates. Rows
  only: no postsolve; LP and fixed-MIP duals are computed on the original rows as before.
- Root gap closed by the presolved LP (presolve only): app1-1 1.2% -> 9.4% (3603 x 1206 / 24.7k ->
  3332 x 1008 / 22.7k), app1-2 0.4% -> 16.5% (39555 x 13265 / 274k -> 38175 x 12342 / 236k; its
  presolve 2 -> 5 s, mostly probing), mzzv42z 39.2% -> 46.0%, mzzv11 34.4% -> 35.8%; fewer nonzeros on
  ran14x18-disj-8 (10.3k -> 8.0k), neos-4954672-berkel, trento1 (947 columns fewer), n5-3,
  neos-3083819-nubu, dano3_3 / dano3_5.
- The hard instances it was made for (300 s, 4 threads, one seed): neos-3216931-puriri 116 -> 70 s
  (202 -> 70 nodes), trento1 gap 0.075% -> 0.028%, neos-4954672-berkel 8.39% -> 8.22%,
  ran14x18-disj-8 56 -> 54 s; app1-2's dual bound reaches the optimum -41 (main: -43; root bound -239
  -> -66), its incumbent is worse (-22 vs -30).
- miplib_wide_core (120 s, 4 threads, seeds 1-2): neutral, 1.010x (90% CI 0.974-1.050), solved 121 ->
  118 of 126 (near-limit instances whose presolved model it does not change). On the 10 instances
  whose presolved model it changes 1.07x: rows combined with equations lose knapsack / MIR structure
  for the cut separators (sp150x300d 0.1 -> 1 s, root bound 69 -> 65; n5-3; neos859080), against
  neos-3083819-nubu (2-5x faster) and dano3_5. `MILP_PRESOLVE_OFF=sparsify` turns it off.
- Equations with a continuous column that has a variable bound (an inequality with continuous columns
  of one sign and integer columns of the other: `x <= u y`, `sum of flows <= 7 y1 + 31 y2`) are not
  used: the cut separators aggregate through those columns (flow balances with variable upper bounds
  on sp150x300d, neos-4954672-berkel and ran14x18-disj-8; arc totals cancelled out of n5-3's capacity
  rows). Their root bounds come back (sp150x300d 65 -> 69, the optimum; n5-3, berkel, ran14x18-disj-8),
  the presolve gains above stay. miplib_wide_core against sparsify without it: 0.965x (90% CI
  0.931-0.996), solved 118 -> 119 of 126; sp150x300d 1.0 -> 0.1 s. `MILP_PRESOLVE_OFF=sparsify_protect`
  turns the rule off.
- Fixed on the way (also on main): two wrong answers of the presolve on small models: GCD tightening of
  a row whose upper side lay just below a multiple of the gcd by round-off (declared infeasible), and a
  round-off coefficient (1.4e-17 on a free integer column) taken as a substitution pivot.

**Heuristic LP budgets**
- Tree dives get a share of the search's LP iterations that follows their success, as SCIP's diving
  heuristics: 5% x (1 + 10 (new incumbents + 1) / (dives + 1)) per thread, at most the previous 10%
  (10% for the first 10 dives, 5.5% after 100 without a new incumbent), and a dive starts only when the
  budget left covers a whole dive (the mean of the thread's finished dives; after a dive that ended at its
  limit, twice its iterations): a dive cut off by the budget finds nothing (pg5_34: 13 of 42 tree dives).
  `MILP_DIVE_ADAPT=0`: the flat 10%.
- Without an incumbent the guided root dive is skipped: it repeated the fractional dive exactly (same LP,
  same choices), in every sub-MIP root and in roots without a start point. `MILP_DIVE_DEDUP=0`: both.
- Sub-MIPs that look for an improvement (RINS, mutation, crossover, local branching, RENS with an
  incumbent) no longer run the feasibility pump: it ignores their cutoff, and its points were nearly
  always worse than the incumbent (6 sub-MIP incumbents from 94 pumps on 13 core instances) while it took
  a quarter to a half of those sub-MIPs' LP iterations (pg5_34, pg, nu25-pr12). `MILP_SUB_FP=1`: the pump.
- A `Heuristics:` summary line: tree dives (calls, points, new incumbents, limit hits, LP iterations) and
  tree sub-MIPs (calls, new incumbents, LP iterations, of them in the sub-MIPs' own dives and pumps).
- miplib_wide_core (120 s, 4 threads, seeds 1-2): solved 121 / 121 of 126, shifted geomean 0.987x (90% CI
  0.962-1.013), first incumbents and final gaps unchanged; more nodes in the same time where both sides
  stop at the limit (rocI-4-11 +18% / +79%, swath3 +46%).

**Probing: liftings, repeated passes, kept implications**
- Probing lifts one-sided rows with what fixing a binary x implies: when the bounds that x = v
  implies leave a row `a y <= b` a slack delta it can never use, x enters the row (x = 1:
  `a_x += delta`; x = 0: `a_x -= delta`, `b -= delta`). At most one lifting per row per pass, never
  equations or ranged rows, new entries capped at 10% of the nonzeros per pass. E.g. `y <= 100 x`
  becomes `y <= 5 x` when x = 1 forces y <= 5 through other rows.
- Probing runs up to 3 passes while a pass fixes, substitutes or lifts something, each with half the
  previous work limit and at most the work the previous one used.
- Implications between binaries found by probing (x = v forces y = w) are kept and propagated
  together with the rows by later probes and passes, so they also work backwards (y = 1 - w forces
  x = 1 - v) where the rows alone do not show it.
- Binaries are no longer substituted out of equations before the first probing pass: probing may
  fix them (neos-3216931-puriri: a binary y = sum of 96 binaries cannot be 0, which turns five such
  rows into partitioning rows).
- Root gap closed by the presolved LP (presolve only): neos-3216931-puriri 1.0% -> 12.3% (LP 58940
  -> 60354), mzzv11 19.2% -> 34.4%, mzzv42z 18.4% -> 39.2% (8769 x 9116 / 126k -> 6250 x 6444 / 74k);
  app1-1, app1-2 and milo-v12-6-r2-40-1 unchanged. At the root (2 threads) the bound after cuts on
  puriri goes 59231 -> 62214 (optimum 71320); mzzv42z is solved at the root (31 -> 19 s), mzzv11 in
  2 nodes. The probing work grows (mzzv42z 179M -> 444M units, presolve about 2x; puriri 21M -> 46M).
- miplib_wide_core (120 s, 4 threads, seeds 1-2, against the same main): solved 119 -> 118 of 126,
  shifted geomean 0.954x (90% CI 0.873-1.034); faster in every run on neos-950242 (8 -> 1-2 s),
  piperout-27 (6 -> 2 s), neos-3381206-awhea and drayage-100-23, slower on cbs-cta (2.4 -> 6.0 s: the
  heuristics find its objective-0 solution later).
- Knobs: `MILP_PRESOLVE_OFF=probing_lift`, `probing_passes` (one pass), `probing_implications`,
  `probing_keep` (all four: the previous presolve); `MILP_PROBING_PASSES=<n>`; opt-in
  `MILP_PRESOLVE_ON=probing_cliques` (the implications x = 1 => y = 0 also extend set-packing rows
  in the clique extension; no better than the liftings alone, and more nonzeros).

**Cheaper roots**
- A root whose restart is already due (enough integer columns fixed by reduced-cost fixing before the
  probes) skips root strong branching; the restarted root probes the smaller problem. The test is the
  restart rule itself, with the root's cost so far (a restart that would be refused keeps the probes), and
  it applies only to roots with some LP work (`MILP_RESTART_SKIP_SB_WORK`, 2e7 iterations x (rows +
  columns)): quick roots keep their probes, which can settle the root (exp-1-500-5-5). h80x6320d (2
  threads): 34.1 / 42.1 / 34.9 s -> 29.9 / 40.4 / 33.4 s; irp (4 threads) 6.5 / 4.6 -> 2.9 / 2.3 s.
  `MILP_RESTART_SKIP_SB=0`: probes in every root.
- With an incumbent known before the root heuristics (from local MIP or FJ), RENS now runs beside the
  dives with that incumbent as its cutoff, as without one, instead of after them with their best point:
  when the dives find the optimum, RENS can only fail and ran to its cap (h80x6320d: 10 s). h80x6320d with
  the above: 29.7 s shifted geomean over 3 seeds (34.3 s without). `MILP_ROOT_MERGE_INC=0`: the old order.
- Restarted roots from the second restart on run their dives, RENS and root LNS only if the search before
  improved the incumbent; an improvement below the gap tolerance (copies of the same point) no longer
  counts (h80x6320d `--deterministic`: 45.2 -> 32.7 s). `MILP_RESTART_IMPROVED_GAP=0`: the old test.
- The root-LP sub-MIP (runs without an incumbent; the root waits for it) gets a budget relative to the
  root: at most max(0.5 s, the root's time) (was max(2 s, ...)), deterministic max(5e7, 3x the root LP
  work) work units (was twice that with a 2e8 floor), and it stops after a first sub-MIP that is
  infeasible, where the LP accepted the fixings and lower fixing rates are guesses. enlight_hard: root 1.15
  -> 0.37 s (`--deterministic`, 2.9 -> 2.1 s). `MILP_LP_RENS_REL=0`: the old caps.
- A root strong-branching pass stops after 4 probes in a row that learn nothing (a child at the iteration
  limit, no fixing): dano3_3 26k probe iterations in 3 passes -> 3-4k (2 threads, 3 seeds: 42.5 -> 35.7 s
  shifted geomean); qap10's root 45 -> 40 s. `MILP_ROOT_SB_ABORT=k`: after k (0: off).
- The root cut loop's cap on flat rounds (the bound not moving on a dual degenerate LP) is now in work units
  instead of wall-clock time, the same in `--deterministic` mode (which had no cap): LP work plus an
  estimate of the separation, at most max(2e8, 3x the initial LP's work). The rounds no longer depend on
  the machine's load (neos-3381206-awhea needs 7 flat rounds, which a loaded machine could cut off), and
  neos-960392's 12 flat rounds in `--deterministic` mode become 1 (124 -> 89 s). `MILP_FLAT_WORK=0`: the
  time cap; `MILP_FLAT_WORK_MIN`: the floor.
- miplib_wide_core (120 s, 4 threads, seeds 1-2, against the code before these changes): solved 120 ->
  122 of 126, shifted geometric mean 10.48 -> 9.25 s (0.892x, 90% CI 0.850-0.935), on a busy machine: part
  of it is neos-3381206-awhea (22.8 / 8.7 -> 3.7 / 3.6 s), whose flat cut rounds the time cap had cut off.
  `MILP_FJ_FIRST=1` stays opt-in (an earlier combined A/B with it: 0.991x; no gain of its own here).

**Dual simplex: column-scaled cost perturbation**
- The dual simplex's cost perturbation, 5e-7 (1 + |c_j|)(1 + random) times a cost scale, now uses the
  column's largest entry (after row scaling) in place of the 1, which is what HiGHS's perturbation of
  its column-scaled LP amounts to. On columns with large entries the old perturbation sat at the level
  of the dual tolerance, where the ratio test's Harris groups tie: csched008's cold root LP crawled at
  its optimal value for ~4900 iterations (HiGHS 1458 in all). The DSE weights now also follow the
  primal simplex pivots that clean up after the perturbation is removed, instead of being reset to 1
  (every warm start from such a basis had unit weights).
- miplib_wide_core (120 s, 4 threads, 4 seeds): 241 -> 245 of 252 solved, shifted geomean 0.977x (90%
  CI 0.946-1.011); root LPs of long degenerate models cheaper (csched008 -24%, dano3 -16%); warm
  re-solves from the same vertex about as before (+2% iterations). `MILP_LP_PERTURB_COLSCALE=0`,
  `MILP_LP_PRIMAL_DSE=0` restore the old behaviour. Measured and left off: perturbing the logicals,
  an objective-stagnation trigger for the cold-solve boost, a signed Harris bound (knobs listed in
  PLAN.md).

**Repeated root restarts (pseudocost carry-over opt-in)**
- The root restarts (presolve again with the root's global fixings and cut rows) now chain: after
  the first one (at 10% of the integer columns newly fixed at the root, as before, at least 10
  columns, was 20) a restarted root restarts again at 10%, up to 4 times (`MILP_MAX_RESTARTS`, was
  1). From the second on a restart needs the presolve after the previous one to have removed at
  least 5% of the integer columns, and the restarted roots so far plus the last one's cost again
  must stay within the first root's cost: its LP work in `--deterministic` mode (so the decision is
  reproducible), its time otherwise (at least 5% of the time limit). Restarted roots from the
  second restart on are cheaper: at most 5 cut rounds (the earlier roots' cuts are carried as
  rows), and the root dives, RENS and root LNS run only if the search before the restart improved
  the incumbent.
- miplib_wide_core + mas74 (120 s, 4 threads, seeds 1-2, both sides on the same main): 122 / 122 of
  128 solved, shifted geomean 1.000x (90% CI 0.968-1.034), measured with the first restart at 5%;
  with it back at 10% the three quick instances that restarted only at 5% (binkar10_1, neos-911970,
  neos17: 0.4-1.4 s slower in both seeds) do not restart, about 0.986x. Chained restarts:
  drayage-25-23 8.8 / 7.6 -> 6.0 / 5.7 s, swath3 solved on one seed (62 s, both TL before),
  h80x6320d 31 / 38 -> 32 / 40 s, swath1 6.8 / 6.2 -> 8.7 / 6.9 s.
- Opt-in `MILP_RESTART_PC`: restarted searches start from the pseudocosts (and inference history)
  of the search before, mapped through the new presolve: `=1` as history with counts capped at the
  reliability threshold (`MILP_RESTART_PC_CAP`), `=2` as estimates only for columns without
  observations of their own. Off by default: on core instances that restart, both needed ~15-25%
  more LP iterations (deterministic runs; mas76 1.60M -> 2.51M / 2.70M, drayage-25-23 46k -> 63k /
  97k).
- Knobs: `MILP_RESTART_FRAC`, `MILP_RESTART_FRAC2`, `MILP_RESTART_MIN`, `MILP_RESTART_MINRED`,
  `MILP_RESTART_WORK`, `MILP_RESTART_ROUNDS`; the old rule is `MILP_MAX_RESTARTS=1
  MILP_RESTART_MIN=20` (identical search). The command line prints the number of restarts;
  `SolveStats` has `restarts` and `tree_restarts`.

**Interrupts and the time limit during the root cut probing**
- The probing that collects implied bounds and literal conflicts before the root cut loop now
  checks the interrupt flag (every 64 probes): an interrupt landing there took up to 1.2-2.2 s to
  stop the solve on cflp_hard_1 (2 threads), now ~0.02 s.
- It also stops at the time limit: its work budget alone could take 6.5 s (fiball, 1 thread), so
  short time limits ran over (fiball at a 4 s limit now ends at 4.04 s).

**Compact open-node storage, gradual memory drain**
- The queue of open nodes keeps a packed 88-byte record per node (was a 112-byte node) with 20-byte
  heap keys that grow by 1.25x (were 24 bytes, doubling), chunked slots that compact under memory
  pressure, warm-start bases stored as single allocations with a 16-byte header (were ~144 bytes for
  a small basis), and the branching bound of an open node inline instead of a bound-path segment of
  its own. gen-ip054 at ~490k open nodes: 493 -> ~265 bytes per open node (334 -> ~220 without
  stored bases); in drain mode at a 300 MB cap 375 -> 200 bytes. Node selection and the applied
  bounds are unchanged: --deterministic runs that stay clear of the memory limit give the same
  nodes, LP iterations and results as before.
- Under --mem-limit: the memory states are decided on the bytes the open nodes use (the hard stop on
  what is allocated); once the limit gets close, queued nodes that an improved incumbent has pruned
  are dropped; drain mode alternates deepest-first pops with best-bound pops (k : 1, k adapted so that
  the tree shrinks slowly; all deepest first close to the stop) and ends below 70% of the limit
  (was 60%), so the bound keeps moving while the tree drains. gen-ip054 at --mem-limit 300 (2
  threads, 300 s): the bound was frozen at 6794.70 from 122 s on; now 6801.41 at the end, still
  rising in drain mode.

**Faster Forrest-Tomlin updates on large sparse bases**
- When the entering column was solved densely, the simplex copied the whole m-row result into the
  update's spike and the update then tested every row. The dense solve now copies only blocks of 8
  rows that hold a nonzero and lists them, and the update reads only those (same entries in the same
  order, so the search is unchanged). Same deterministic work, 2 threads: rmatr100-p10 35.0 ->
  29.8 s, unitcal_7 (160 nodes) 99 -> 91 s, ns1208400 20.0 -> 19.5 s, h80x6320d unchanged. Profiling
  the strong-branching probes found no other per-probe overhead (state restore 0.1-1.2% of probe
  time; a probe iteration costs what an iteration of the root's late cut-round re-solves costs).

**Strong-branching control and score corrections**
- Dynamic strong-branching reliability, on by default (`MILP_SB_DYN=2`): reliability falls from 4
  towards 1 as strong-branching iterations exceed half of the node LP iterations. miplib_wide_core
  (120 s, 4 threads, seeds 1-2, x86-64 Linux machine): solved 109 -> 114 of 126, shifted geomean
  18.45 -> 17.96 s (0.975x); 30n20b8 now solved on both seeds, neos-3381206-awhea 41.7 -> 16.5 s;
  slower neos-3627168-kasai (1.39x). `MILP_SB_DYN=0` turns it off, `=3` uses only the probe
  iteration limit of twice the average node LP, `=1` both.
- Opt-in (off by default): `MILP_SD` (bits), score corrections after Shah and Dey for reliability and root strong
  branching: 1 = child gains capped at the node's gap to the incumbent, 2 = `q_min^0.3 q_max^0.7`
  instead of the product, 4 = an asymmetry factor from the direction of the last branching of the
  infeasible / integral leaves (`MILP_SD_ETA`, default 0.15).
- Opt-in: `MILP_PL` (1 exponential, 2 Pareto fit), probabilistic lookahead after Mexi et al.; once the best
  candidate has not changed for 0.6x the lookahead, probing stops when one more candidate is not
  expected to shrink the abstract tree implied by the gap and the gains seen at the node.
- Opt-in: `MILP_APS` (discount factor), ancestral pseudocosts, used at nodes where every candidate is
  reliable in both its own and its ancestral pseudocosts.
- Deterministic runs on miplib_wide_core (2 threads, node limit 10000, 90 s) against the v0.5.0 defaults:
  `MILP_SD=7` and `MILP_PL=1` change nodes by -1.4% / +1.4% and LP iterations by +6% / +0.6% on the
  instances both solve, with gaps on the unsolved ones better on about as many as worse;
  `MILP_SB_DYN=2` cuts strong branching from 35% to 28% of the LP iterations for 11% more nodes at
  the same total iterations and closes more gap on unsolved instances (7 better, 3 worse);
  `MILP_SB_DYN=1` loses on unitcal_7 and mc11 (longer probes).

**Feasibility jump before local MIP (opt-in, off by default)**
- `MILP_FJ_FIRST=1`: the helper thread runs its first feasibility-jump attempt before the local-MIP
  search instead of after it. Local MIP without a solution restarts until its whole budget is used
  (3e9 effort, 5-10 s, longer on a loaded machine), so on instances where FJ finds points at once and
  local MIP does not, the root ran without an incumbent: neos-3381206-awhea, FJ's first point at
  0.2-0.6 s in v0.4.0, none from local MIP in 30 s on a loaded machine (FJ never ran before the
  thread's deadline). With the knob (2 threads): 4 seeds, shifted geomean 3.86 -> 3.03 s (v0.4.0 3.12 s);
  on a loaded machine 30.4 -> 9.7 s (v0.4.0 15.8 s). Off: the same search as before.

**Time limits on large models**
- Phases that ran past the time limit on large MIPLIB 2017 models (bab2 was stopped by the benchmark runner
  after 144 s at a 120 s limit):
  - Packing-orbitope detection: one work budget for all symmetry components, checked per hashed column (was
    per component, between passes): bab2 several minutes -> 0.7 s.
  - Symmetry-breaking (SST) reductions: the work budget is also checked while the orbits are scored.
  - Sub-MIP presolve stops with the sub-MIP (its stop flag and time limit; the feasibility-jump repair on
    supportcase19, 1.4M columns, ended 4 s late).
  - Probing liftings that add a new entry append it without searching the column: ns1760995's probing 37.6
    -> 2.8 s (~200k new entries per pass), the same presolved model.
- Deterministic runs are unchanged.

**Tree-size-estimation restarts (opt-in, off by default)**
- `MILP_EST_RESTART=1`: the search counts its closed leaves per depth (tree weight: the sum of
  2^-depth, independent of the order in which leaves close, so `--deterministic` reads the same value)
  and estimates the tree size from it with HiGHS's windowed secant. When the estimate stays more than
  20 r times the nodes so far (r from the share of columns fixed) and the gap hardly shrinks, the
  search restarts once (presolve again with what the tree fixed), within the first half of the time
  limit. Verbosity 2 prints the tree weight. miplib_wide_core plus 5 larger instances (120 s, 4 threads,
  seeds 1-2): solved 125 -> 121 of 136, 1.021x (90% CI 0.990-1.055): its restarts are mostly on big
  trees (mas74) that a restart only sets back. Knobs `MILP_EST_FACTOR`, `MILP_EST_GATE`,
  `MILP_EST_MIN_NODES`, `MILP_EST_TRACE`.

**Cut management knobs (opt-in, off by default)**
- Diagnostics: `MILP_CUT_TRACE=1` logs the root cut rounds (candidates, rejections as too parallel or
  over the round limit, re-added and re-separated cuts, rows removed by age); verbosity 2 adds each
  search's node time split (LP solves, cut-round re-solves, separation per separator, pool scans) and
  pool churn; `MILP_POOL_TRACE=1` per progress line.
- Knobs: `MILP_CUT_AGE=k` (root cut rows stay until k + 1 re-solves in a row found them slack),
  `MILP_POOL_DEDUP=1` (pool cuts with the coefficients of a live cut map to it),
  `MILP_POOL_AGE_NEW=k` (pool rows that never bound leave after k slack node LPs), `MILP_POOL_PAR`,
  `MILP_POOL_SPARSE`, `MILP_NODE_SEP_WORK`. miplib_wide_core (120 s, 4 threads, seeds 1-2): dedup +
  aging of new rows 1.002x (90% CI 0.968-1.040), root row aging 1.016x (0.961-1.072); unset, the search
  is unchanged.

**Fixes**
- Probing could prove a feasible model infeasible: a continuous column whose implied bounds crossed
  within the tolerance was fixed at the midpoint, which could lie outside its current domain, and over a
  long row such crossed bounds added up (neos-5104907-jarama: INFEASIBLE after 3 s, optimum 935, in
  the second probing pass). The column is now fixed inside its domain.
- The benchmark verifier (`benchmarks/lib/mps.py`) took indented lines starting with `*` for comments
  and dropped irish-electricity's right-hand sides (its RHS vector is named `*RHS*`), so correct
  solutions were reported infeasible. Only a `*` in column 1 starts a comment now.

## v0.5.0 — 2026-11

MILP: root strong branching, network cut-set cuts, a local-search primal heuristic and earlier
root heuristics, conflict analysis by generalized resolution, presolve of long equations and
implied integers. Also LP and fixed-MIP dual values, interrupts, an incumbent callback and gzipped
input.

**Benchmarks**
- 43 hard MIPLIB 2017 instances (600 s, 4 threads): 15 solved, shifted geometric mean 324.9 s
  (v0.3.0: 9, 468.9 s; HiGHS 12, SCIP 6). Newly solved: nexp-150-20-8-5 (16 s), tr12-30,
  supportcase26, neos-1171737, neos-1456979, fhnw-binpack4-48; all nine solved before still are
  (csched007 326 -> 171 s, gen-ip054 541 -> 217 s; glass4 356 -> 534 s).
- 63-instance MIPLIB 2017 subset (120 s, 4 threads), seeds 0 / 1 / 2: 58 / 59 / 61 solved, shifted
  geometric mean 10.7 / 10.8 / 9.9 s (v0.4.0: 59 / 59 / 58, 10.5 / 10.9 / 11.4 s); about 4% faster
  over the three seeds.
- Generated models (1 thread, 120 s): 146 of 148 solved (v0.4.0: 147); the difference is
  binpack_hard_1, which depends on the seed in both versions.

**Root strong branching rounds**
- After the root cut loop and the root heuristics, the root strong-branches the fractional
  candidates in pseudocost order (dual simplex warm starts from the root basis, at most
  root LP iterations / 10 clamped to [100, 500] iterations per child, stop after 8 candidates
  without a better score, total budget 0.5x the root LP iterations + 5000). A candidate whose one
  child is infeasible or cut off is fixed to the other side globally; the root is then processed
  again (propagation of the fixings, LP, up to 5 cut rounds that keep the earlier cut rows,
  another pass) while passes fix and budget remains. The probes initialise the pseudocosts, the
  root branches on the best probed candidate with its children's exact bounds, and
  max over candidates of the smaller child bound is a global bound. Not in sub-MIP heuristics.
  Knobs: `MILP_ROOT_SB=0` off, `MILP_ROOT_SB_WORK`, `_OFS`, `_ITERS`, `_CANDS`, `_LOOKAHEAD`,
  `_PASSES`, `_CUTS`, `_SUBMIP=1`.
- Root gap (miplib_wide_core, mean share of the root gap closed): 66.8 -> 67.3% (SCIP 64.6%), app1-1
  78 -> 92.5% (SCIP 95%), seymour1 43 -> 48%; root time shifted geomean 3.86 -> 4.06 s.
- End to end (miplib_wide_core, 120 s, 4 threads, A and B back to back): on the 35 instances
  where root strong branching fixes columns or moves the root bound, 0.968x shifted geomean over seeds
  1, 2 (1.010 / 0.928; solved 34 -> 33 / 33 -> 34), swath3 78 -> 31 s, timtab1 14 -> 6 s, ns1830653
  94 -> 49 s; slower h80x6320d, mzzv11. Root time on the other 28: 1.045x. --deterministic
  cflp_medium_1: 7286 -> 3623 nodes.
- `CutLoopOptions::keep_rows`: a cut loop on an LP that already has cut rows only removes its own.

**Interfaces: dual values, interrupts, incumbent callback, gzipped input**
- LP duals: row duals, reduced costs and the optimal basis of LPs solved to optimality. Presolve
  only postsolves primal points, so the duals are recovered in the original space after the solve:
  a crossover (primal push only; `crossover()` now accepts a point without duals) from the
  postsolved optimal point gives a basis of the original LP, and the simplex finishes it (primal
  first). The duals of that optimal basis are complementary to the returned solution (any optimal
  primal-dual pair is), with or without presolve; the primal solution and all search statistics are
  unchanged. Cost: 0-75 cleanup iterations on the generated instances, < 0.01 s up to the medium
  ones, 0.9 s on randlp_large (24 s solve). New `solver/duals.hpp`.
- Fixed-MIP duals: the duals of the LP with the integer columns fixed at the incumbent, after the
  solve, on request (`SolverParams::duals` = 2; 1 = LPs only, 0 = none). Not for models with
  logical constraints (SOS, indicators, semi-continuous) or quadratic objectives.
- Exposed as before for QPs: CLI `--duals <file>` (also sets fixed-MIP duals for MILPs; a comment
  line marks them), C API `popso_get_duals` / `popso_get_reduced_costs` with the new parameter
  `duals` (default 1 in the C API and Python, 0 in the CLI without `--duals`), `popso_has_duals`
  returns 2 for fixed-MIP duals, new `popso_has_basis` / `popso_get_basis` (`POPSO_BASIS_BASIC`,
  `_LOWER`, `_UPPER`, `_ZERO`); Python `Result.duals`, `reduced_costs`, `fixed_mip_duals`,
  `var_basis`, `constr_basis` (`popso.BASIC`, `AT_LOWER`, `AT_UPPER`, `AT_ZERO`), `solve(duals=2)`.
- Interrupts: `SolverParams::interrupt` (an `std::atomic<bool>*`); once set, the solve stops with
  the existing status `Interrupted` (C `POPSO_STATUS_INTERRUPTED`) and the best solution so far. The
  flag is installed as a thread-local pointer (`core/interrupt.hpp`) that every solver thread
  inherits, and polled by the branch-and-bound loops (opportunistic and deterministic), the simplex
  iterations, the barrier iterations, its Cholesky factorization (per supernode and panel) and CG
  refinement, the cut loop, presolve rounds and probing, symmetry detection, feasibility jump,
  diving, the feasibility pump, sub-MIP retries and the MIQP search. Measured stop latency after
  the request: 3-40 ms (cflp_hard root and tree, randlp_large in the barrier, which took up to
  1.5 s before the factorization and CG checks). Separate models in different threads do not
  interrupt each other. C API `popso_interrupt(m)` (one lock-free atomic store: thread-safe and
  async-signal-safe; the request is cleared when `popso_solve` starts). CLI: SIGINT (`std::signal`)
  stops the solve and prints the normal summary; a second Ctrl-C ends the process. Python:
  `Model.interrupt()`; Ctrl-C during `solve()` returns the result with status `"INTERRUPTED"` (the
  native solve runs on a helper thread when called from the main thread).
- Incumbent callback: `SolverParams::on_incumbent(x, objective, bound)` / C
  `popso_set_incumbent_callback(m, fn, user_data)` with `int fn(void* user_data, double objective,
  double bound, const double* x, int num_vars)`, Python `Model.set_incumbent_callback(fn)` with
  `fn(objective, bound, x)`. The point is mapped back through presolve, symmetry reductions,
  restarts and logical reformulations to the model's columns; calls are serialized, strictly
  improving and, in deterministic mode, made at the deterministic merge points (identical sequences
  across runs); the final postsolved solution is reported too if it improves. A nonzero / true
  return (or a Python exception, re-raised by `solve()`) stops the solve like an interrupt; no calls
  follow. Not for quadratic objectives.
- Gzipped input: `.mps.gz` / `.lp.gz` (any file starting with the gzip magic bytes, whatever its
  name) are decompressed while reading by a from-scratch inflate (`io/gzip.hpp`, RFC 1951/1952:
  stored, fixed- and dynamic-Huffman blocks, all header flags, concatenated members, CRC-32 and
  length checks; streaming input and output with a 32 KB window; ~280 MB/s). MIPLIB `.mps.gz` files
  can be read directly (`benchmarks/miplib/download.py` still unpacks them).
- No change to search results: `--deterministic --threads 2` gives the same nodes / LP iterations
  (cflp_medium 7286 / 269906, gap_medium 8551 / 186860, lotsize_medium 7469 / 160109) at the same
  speed (time ratios 0.99-1.004).
- Tests: `lp_duals_kkt_random` (150 random LPs with every bound and row type and presolvable
  structure, minimization and maximization, presolve on and off: stationarity, sign conditions,
  complementarity <= 1e-7, strong duality, square basis), `lp_duals_generated_instances` (randlp,
  transport, 20 lpsanity LPs), `fixed_mip_duals_match_fixed_lp` (80 random MILPs: the fixed-MIP duals
  are optimal for the LP with the integers fixed, against its own solve), interrupt tests (preset
  flag, flag set from a thread during the search, LP / barrier / MIQP paths, scope inheritance),
  incumbent callback order / determinism / stop, gzip test vectors made with Python's zlib (levels
  0-9, Huffman-only strategy, all header flags, two members, 368 KB with window slides, byte-by-byte
  input), corrupt and truncated data, gzipped model files; C API and Python tests for all of it
  (including a simulated Ctrl-C).

**Primal heuristics: earlier incumbents**
- Local-MIP style local search (new `heuristics/lmip.hpp`, after Lin et al., "Local-MIP", 2024/25):
  tight moves, lift moves and breakthrough moves with weighted constraint scores, run first on the
  feasibility-jump thread (default effort 3e9, about 5-10 s; random restarts while no solution is
  known); in opportunistic mode the thread then restarts it from each new incumbent with the
  incumbent as cutoff, at most 5% of the elapsed time. Single-threaded root heuristic study: a first
  feasible point on 58 of 63 core instances (feasibility jump alone: 50) and 36 of 43 hard ones (28).
  Measured end to end (miplib_wide_core, 120 s, 4 threads, seeds 1-2, x86-64 Linux machine): solved 102 ->
  109 of 126, shifted geomean 20.2 -> 18.4 s. `MILP_LMIP=0` turns it off.
- Root-LP sub-MIP (`lp_rens` in `heuristics/submip`), run at the root when no incumbent exists:
  RENS with relaxed, checked fixing. Integer columns that are integral in the root LP are fixed in
  batches (largest |reduced cost| first, then the longest columns), each batch propagated, the LP
  re-solved and the restricted problem presolved (light probing); a batch the LP or presolve rejects
  is undone and split, a single rejected column stays free; with no integral column left the least
  fractional one is rounded. The sub-MIP on the remaining domain is retried with a bisection on the
  fixing rate (exhausted: fix less; limit hit: fix more), and when the two bounds meet the fixings
  between them are left free and the dive starts again from the target (0.7). Work budget in
  deterministic mode (2x the root sub-MIP iteration cap), time otherwise (max(2 s, root time),
  at most 10% of the remaining time). `MILP_NO_LP_RENS=1` turns it off.
- Root heuristics without an incumbent run the sub-MIP heuristics (RENS, feasibility pump,
  root-LP sub-MIP) and the four dives as one parallel batch (sub-MIPs first) instead of dives first;
  in opportunistic mode every point is published as soon as its task ends (blp-ic98: first
  incumbent 8 s instead of 16-21 s). `MILP_NO_ROOT_MERGE`, `MILP_NO_EARLY_PUBLISH`.
- Root LNS rounds get at least 3 s (deterministic: 4x the LP iterations) when the incumbent is more
  than 1% above the root bound (`MILP_ROOT_LNS_POOR`, seconds; 0 = off).
- Tree LNS: on slow nodes the sub-MIPs no longer wait 100 nodes between calls; 10 nodes and 50000
  search LP iterations since the last call suffice (`MILP_LNS_ITER_GAP`, 0 = off). The share of LP
  work they may use is unchanged.
- RINS caps its fixing rate: when the agreement between LP point and incumbent fixes > 60% and the
  sub-MIP is exhausted (often already in presolve: nexp-150-20-8-5 fixed 95%), the next RINS
  leaves a random subset of the agreeing columns free down to 10% below that rate (at least 50%);
  a limit hit raises the cap again (`MILP_NO_RINS_CAP`).
- Measured (60 instances of miplib_wide_core, 120 s, 4 threads, seeds 1-2, A and B back to
  back): solved 110 -> 113 of 120, shifted geomean 13.67 -> 12.17 s (ratio 0.898, 90% CI
  [0.832, 0.969]); blp-ic98 first incumbent 12.6-21 s -> 5.6-8 s. neos-3216931-puriri still has no
  incumbent in 120 s.

**Cutting planes and presolve**
- Network cut-set cuts (new `cuts/mcf.hpp`, separator "mcf", root cut loop): multi-commodity flow
  structure is detected (flow conservation rows, arc or undirected edge capacity rows with integer
  design columns, nodes shared by the commodities), and for single nodes and the node sets formed by
  contracting the non-binding arcs, the flow rows of the set plus the capacity rows of the crossing
  arcs give c-MIR / aggregated cover cuts (as SCIP's mcf separator). nexp-150-20-8-5: root bound
  92.4 -> 230.2 (optimum 231), solved in 22 s (4 threads; before: not solved in 600 s). `MILP_MCF=0`
  turns it off.
- Presolve no longer substitutes continuous columns with an objective coefficient out of an equation
  when the equation does not imply their bounds (the bounds became the row's sides and the columns'
  variable bounds were lost to the c-MIR path aggregation): tr12-30 root bound 127730 -> ~130000
  (optimum 130596), solved in 274 s (4 threads; before: not solved in 600 s).
  `MILP_PRESOLVE_ON=subst_bounded` restores the old behavior.

**Presolve: long equations, implied integers**
- Column singletons are substituted out of equations of any length (the 1000-entry limit now only
  applies to non-singleton pivots), implied free singletons also when their coefficient is small
  relative to the row (objective-defining rows z = c^T x: blp-ic98 z with 1 vs 1973), implied free
  ones first; a singleton that does not qualify no longer keeps the other columns from being tried.
  Presolved nonzeros: blp-ic98 191.7k -> 169.6k, leo1 114.2k -> 107.4k, swath1/3 31.7k -> 25.5k,
  mcsched 7.4k -> 4.5k, germanrr 168.6k -> 157.9k, sct2 20.1k -> 19.0k (HiGHS: 175k / 107k / 25.5k /
  4.6k / 158k / 18.8k).
- Dual bound strengthening keeps an integral continuous bound integral (csched008's makespan bound
  was 190.000000191: a relative 1e-9 margin on top of an exact 190).
- Implied integers (`presolve/implied_int.hpp`): continuous columns integral in every feasible point
  with integral integers (primal rule: an equation whose other columns are (implied) integer with
  integral ratios and side, iterated to a fixed point) or movable to an integral value without loss
  (dual rule: integral bounds and only such rows). Stored as `MipProblem::implied_int` (not
  `is_int`): never branched on, not required for feasibility; their bounds are rounded; used for
  objective integrality (replaces the one-level equation check of the search) and treated as integers
  by the c-MIR and Gomory separators (with integral bounds). csched008 212 of 212 continuous columns,
  csched007 292/292, uct-subprob 1445/1445, neos-3216931-puriri 159/159, trento1 668/683,
  neos-827175 10248/10248; csched007/008 and uct-subprob now have an integral objective.
  `MILP_PRESOLVE_OFF=implint` / `implint_dual` switch the detection / the dual rule off.
- Measured (A and B back to back, 4 threads, 2 seeds, 120 s): 46 of the 63 core instances
  presolve to byte-identical models (unchanged search, presolve time 9.68 vs 9.57 s); on the 17 others
  solved 15 / 13 (seed 1) and 13 / 15 (seed 2), geomean time ratio 0.957, WRONG 0 (irp 0.44, eil33-2
  0.58; swath3 / CMS750_4 lost a seed-1 solve). At 300 s: csched007 solved (288 s; was 2.9% gap),
  uct-subprob optimum 314 found and bound 295.4 -> 298, sct2 gap 0.025% -> 0.017%; csched008 bound still
  171; blp-ic98 35% fewer nodes (denser objective), bound 4454.7 -> 4450.5, incumbent 4508.4 -> 4504.9.
- Tests: long objective-row substitution, singleton blocking, integral dual bound, primal chain and
  dual (makespan) detection, a fractional start for an implied integer, 800 random implied-integer
  models (presolve on / off and an objective-scaled reference without objective integrality), LP
  duals with a substituted 1200-entry objective row.

**Conflict analysis**
- Cut-based conflict analysis by generalized resolution (new `search/genres.hpp`, after Mexi et al.,
  "Cut-based conflict analysis in mixed integer programming", IJOC 2025): an infeasible or cut-off
  node's dual proof is resolved against the reason constraints of the trail's bound changes, with
  c-MIR reduction of the intermediate constraints, and the learned linear constraints are used next
  to the existing conflict clauses. Validity tests by enumeration on random pure-integer and mixed
  models. Measured (miplib_wide_core, 120 s, 4 threads, 3 seeds, A and B back to back): seed 1
  solved 58 -> 60 of 63 at the same shifted geomean; seeds 2-3 solved 104 -> 105 of 126, shifted
  geomean 20.16 -> 19.77 s. `MILP_GENRES=0` turns it off (`2`: conflict clauses only where the
  resolution fails).

**Fixes**
- binkar10_1 (`--gap 1e-6`) sometimes ended OPTIMAL at a worse objective (6742.22-6743.43 for the
  optimum 6742.200024; reproducible with `--deterministic --threads 2 --seed 4`). Presolve
  substitutes its penalty columns (cost 1e6, column singletons of equations) out, so their bound
  becomes a row of the presolved problem; an LP point may violate that row within the feasibility
  tolerance, and postsolve, which clamps the substituted column into its bound, moved the violation
  into the objective (2.5e-8 x 1e6). The search's incumbent and proof were for the presolved value
  (6742.1993); the postsolved point was 0.026 worse. Now a postsolved point whose objective came out
  worse than the search's gets its continuous columns re-optimized (integers fixed, LP at primal
  tolerance 1e-9; 6742.200024 here), and the better point is kept (new `solver/postsolve.hpp`). The
  search is unchanged.
- A search stopped by a limit or an interrupt right after an in-tree restart reported its bound
  as -inf; it now keeps the bound from before the restart.
- Gzip input: the decompressed size is capped (default 16 GiB, `milp::set_gzip_max_output`), so a
  small high-ratio file can no longer exhaust memory; the error names the file and the limit.
- LP and fixed-MIP duals are dropped (with a log line) when the basis they come from does not
  reproduce the solution's objective (relative 1e-6; fixed-MIP: may be better than the incumbent,
  which is logged). Documented what fixed-MIP duals mean after an early stop (the fixed LP's
  optimum, not necessarily complementary to the incumbent's continuous part).
- The dual recovery (crossover and simplex cleanup) now honours interrupts and the time limit
  (at least 1 s) and gives up on the duals when stopped; the crossover used to ignore both.
- C++ API: a stop requested by the incumbent callback set the caller's `SolverParams::interrupt`
  flag, so reusing the params interrupted the next solve at once. The callback now uses an
  internal flag; the caller's flag is only read (`InterruptScope` takes two flags).
- Python: any exception (not only KeyboardInterrupt) escaping the wait in `solve()` now interrupts
  the native solve and waits for the helper thread, which used to keep running on a model that
  could be freed.
- Time limit inside the barrier: checked with the interrupt's fine-grained checks (Cholesky
  supernodes and panels, CG refinement) through a thread-local deadline (`DeadlineScope`,
  `interrupted_or_expired()`); `--time-limit 0.0001` on randlp_large stops after 0.5 s instead of
  7 s.

## v0.4.2 — 2026-11

**Fixes**
- Diving heuristics could loop forever (and grow memory until the limit) on an integer column
  whose LP value lies just outside an integer bound (x = ub + 1e-6, within the LP's tolerance):
  rounding it does not change its domain. Such columns are no longer chosen. MIPLIB
  graphdraw-domain: from the 2 GB memory limit to optimal in 293 s (4 threads, 357 MB).

## v0.4.1 — 2026-11

**Fixes**
- Windows (MSVC): C API calls that reject their arguments deep inside the library (for example
  `popso_add_q_entries` with an out-of-range variable index) aborted the process instead of
  returning `POPSO_ERROR_INVALID_ARGUMENT`. The shared library is now compiled with `/EHs`
  instead of MSVC's default `/EHsc`, which assumes that `extern "C"` functions never throw.

## v0.4.0 — 2026-11

New problem classes: convex quadratic programs (QP) and mixed-integer convex quadratic programs
(MIQP): a quadratic objective `c'x + 1/2 x'Qx` with linear constraints. Also better c-MIR cuts
and symmetry-breaking reductions for MILP.

**Benchmarks**
- 63-instance MIPLIB 2017 subset (120 s, 4 threads): 59 solved, shifted geometric mean 10.5 s
  (v0.3.0: 60, 11.4 s); over three seeds 10.5 / 10.9 / 11.4 s against 11.4 / 11.4 / 11.1 s.
  The instance count differs by two instances just above the time limit.
- Generated models, now 148 including 18 QPs and MIQPs (1 thread): 147 solved; SCIP solves 134,
  HiGHS 125 of the 130 it handles (no SOS, indicator or MIQP).
- Maros–Meszaros convex QP test set: all 97 convex problems of the 98 smallest, 32 of the 40
  larger ones (60 s).

**Quadratic objectives**
- Model: a symmetric Q stored as its lower triangle. Input: MPS `QUADOBJ` (one triangle),
  `QMATRIX` and CPLEX `QSECTION` (whole matrix); LP-format objective terms
  `[ 2 x ^ 2 + 3 x * y ] / 2`. The MPS writer writes `QUADOBJ`. Quadratic constraints (MPS
  `QCMATRIX` or a `QSECTION` of another row, brackets in LP constraints) are rejected with a
  "not supported yet" error.
- Convexity check before solving (sparse Cholesky of the unit-diagonal scaled Q plus a tiny
  shift; a negative pivot means a negative eigenvalue). A nonconvex objective (concave
  objective when maximizing) gets the new status `NONCONVEX` (C API `POPSO_STATUS_NONCONVEX`).
- Continuous QPs: the barrier method with Q in the Newton system (augmented system, or the normal
  equations when Q is diagonal), no crossover; infeasibility and unboundedness are confirmed
  by the simplex on auxiliary LPs. Once the dual regularization's error dominates the primal
  residual, the Newton systems are refined against the unregularized matrix (problems with
  large optimal duals such as YAO stalled otherwise).
- Dual values of continuous QPs: row duals and reduced costs (d = c + Qx - A'y) of the
  barrier's final point, in the model's objective sense. CLI `--duals <file>` (`row <name>
  <value>` / `col <name> <value>` lines), C API `popso_has_duals`, `popso_get_duals`,
  `popso_get_reduced_costs`, Python `Result.duals`, `Result.reduced_costs`, `Result.dual(c)`,
  `Result.reduced_cost(x)`. Not for MIQPs (LPs and MILPs: see "Interfaces" below).
- Presolve for quadratic models (QP and MIQP): fixed and empty columns (their quadratic terms
  folded into the costs and the objective constant), dual fixing of columns without quadratic
  terms, empty, free and singleton rows, redundant row sides and infeasibility from row
  activities, integer bounds from activities (MIQP). No column substitution, so the Hessian
  stays as it is; the duals are carried back (a singleton row whose bound is active takes the
  reduced cost of its column). The convexity test runs on the presolved Hessian (a cross term
  with a fixed column is linear).
- MIQPs: branch-and-bound (best bound with plunging, domain propagation of the rows),
  deterministic. The root relaxation is solved by the barrier; the nodes by a dual active-set
  method (Goldfarb–Idnani on the KKT system of the working set: dense Cholesky of the free
  columns' Hessian plus a Schur complement of the active rows, later working-set changes as
  borders of that factorization), warm-started from the parent's working set and, when it is
  still cached, its factorization; a singular Hessian gets a tiny proximal term. Node bounds
  are safe Lagrangian bounds (valid for any primal point and row multipliers by convexity),
  which also give reduced-cost tightening of integer bounds at the nodes and at the root.
  Reliability branching (strong branching by the active-set method on candidates with fewer
  than 4 pseudocost observations), plunging into the child with the smaller pseudocost
  estimate. Heuristics: rounding, a QP dive and feasibility jump on the objective linearized
  at the root solution. The barrier remains the fallback when the active-set method fails.
  Compared with barrier relaxations at every node: miqprand_large 4.9 s -> 0.3 s,
  miqprand_hard (new) 43 s -> 2.7 s, miqpport_hard (new) 13.4 s -> 0.8 s.
- C API: `popso_set_q`, `popso_add_q_entries`, `popso_num_q_nonzeros`, `popso_get_q`. Python:
  `popso.QuadExpr` (products and squares of variables and linear expressions) in
  `set_objective`, `Model.get_q()`.
- Benchmarks: new families `portfolio`, `lsq` (constrained least squares) and the MIQPs
  `miqpport` (cardinality-constrained portfolio) and `miqprand`; the solution checker
  recomputes quadratic objectives; `benchmarks/tools/qp_fuzz.py` compares random QPs and MIQPs
  with HiGHS and SCIP; `benchmarks/tools/maros_meszaros.py` fetches and converts the
  Maros–Meszaros QP test set. The 16 generated QP/MIQP instances are solved in a shifted
  geometric mean of 0.25 s (2 threads; references: HiGHS 1.5 s on 14, SCIP 5.6 s on 12).
- Maros–Meszaros (qpbenchmark's collection, 2 threads): the 98 problems up to 60 KB are all
  solved (97 convex, 1.8 s in total; VALUES is reported as nonconvex); of the 40 larger ones
  (60 s limit) 32 are solved (BOYD1/2, CONT-300, EXDATA among them); LISWET1, 7-10 and 12 end
  with a numerical error and CVXQP1_L / CVXQP3_L at the time limit.

**Symmetry**
- Symmetry-breaking reductions from a stabilizer chain (Schreier-Sims table / SST): after
  presolve, leaders whose orbit contains columns in conflict with them (a row forbids both at
  1) give x_leader >= x_follower rows and fix the conflicting followers to 0, then the search
  continues in (a subgroup of) the leader's stabilizer. Used when the fixings remove at least
  10% of the columns; the problem is presolved again and the remaining symmetry is handled as
  before (packing orbitopes, orbital branching). cod105: 1024 -> 341 columns, solved at the
  root (12 s -> 1-3 s).
- CLI `--symmetry on|off` (default on) switches all symmetry handling.

**Fixes**
- Barrier: inequality rows with fixed columns had their slack bounds shifted twice (the
  crossover repaired the point, so LP results were right, but the interior solution was not).

**Cutting planes**
- c-MIR path aggregation now eliminates a continuous column through the tight row with the largest
  LP dual (the rows that prove the LP bound) instead of the first tight row found (not for columns
  with a simple variable upper bound x <= u y, where the old order suits fixed-charge paths). On the
  63-instance MIPLIB subset the root node closes 65.2% of the root gap on average (63.8% before;
  SCIP 64.6% with node limit 1); neos17 goes from 18% to 94% of its root gap closed and solves in
  about 1 s instead of 4-6 s. End to end the subset is unchanged (60 s / 120 s A/Bs within noise).

## v0.3.1 — 2026-11

Benchmark update; the solver is the same as in v0.3.0.

- BENCHMARKS.md: new section on 43 hard MIPLIB 2017 instances (those of the 106-instance
  development set that no solver finishes within 120 s), run with a 10-minute limit. popso solves
  9, HiGHS 12, SCIP 6; set file `benchmarks/sets/miplib_hard.txt`.
- HiGHS and SCIP reference results for these instances at 600 s.
- `plot_profile.py`: `--tmin` and `--ymax` to zoom the time and count axes.

## v0.3.0 — 2026-11

New capabilities: an interior-point (barrier) LP solver and the logical modelling features real
models use — semi-continuous variables, SOS1/SOS2 and indicator constraints.

**Benchmarks**
- 63-instance MIPLIB 2017 subset (120 s, 4 threads): 60 solved, shifted geometric mean 11.4 s
  (unchanged from v0.2.1: this release adds features, not MILP speed).
- Generated models, now 130 including three families with logical constraints (1 thread):
  129 solved; SCIP solves 122 of 130, HiGHS 117 of the 122 it can read.
- LP relaxations of the MIPLIB instances: the new default `--lp-method auto` is about 27% faster
  than the dual simplex alone. BENCHMARKS.md has a new plot of popso's progress across versions.

**Barrier (interior-point) LP solver**
- New LP algorithm, written from scratch like the rest: Mehrotra's predictor-corrector
  primal-dual method with Gondzio's centrality correctors, on the presolved and scaled problem
  (geometric scaling; free variables, ranges and two-sided bounds handled directly).
- Linear algebra: an approximate minimum degree ordering (AMD) and a supernodal sparse Cholesky
  factorization, multithreaded for large dense blocks. Each solve picks the cheaper of the normal
  equations and the regularized augmented system (factored as a quasi-definite LDL^T, which suits
  LPs with far more rows than columns and is the form a future QP solver needs). Dense columns in
  the normal equations go through a Sherman–Morrison–Woodbury update; every solve is refined
  iteratively.
- Crossover to an optimal basis: dual and primal push phases, then a short simplex cleanup, so
  the answer is an exact vertex as with the simplex.
- `--lp-method auto|simplex|barrier` (C API parameter `lp_method`: 0, 1, 2; the Python bindings
  also take the names). The method is used for a pure LP and for the root LP of a MIP. `auto`
  (the default) runs the dual simplex and the barrier at the same time on large pure LPs and
  takes the first to finish (not in `--deterministic` mode, where it stays with the simplex);
  MIP root LPs stay with the simplex unless `barrier` is chosen.
- On LP relaxations of 69 MIPLIB instances (4 threads), `auto` takes 27% less time than the
  dual simplex alone in shifted geometric mean and 44% less in total; the barrier alone is about
  as fast as HiGHS' interior-point solver (shifted geometric mean 0.40 s vs 0.43 s, with
  crossover). Examples: chromaticindex512-7 20 s -> 1.6 s, dano3_3 6.2 s -> 2.7 s, qap10
  3.4 s -> 1.1 s, rmatr200-p5 5.9 s -> 2.9 s, brazil3 2.4 s -> 0.8 s.

**Model features**
- Semi-continuous and semi-integer variables (x = 0 or lb <= x <= ub), SOS1 and SOS2
  constraints, and indicator constraints (if a binary variable takes a given value, a linear
  constraint must hold).
- They are read from MPS files (`SC`/`SI` bounds, `SOS` section, `INDICATORS` section) and LP
  files (`Semi-Continuous` section, `SOS` section with `S1::`/`S2::` entries, and indicator
  constraints written `name: b = 1 -> x + y <= 5`). The MPS writer writes them, so models
  round-trip.
- C API: `POPSO_SEMICONTINUOUS` and `POPSO_SEMIINTEGER` variable types, `popso_add_sos`,
  `popso_add_indicator`, `popso_num_sos`, `popso_num_indicators`, and `popso_read_lp` for LP
  files. Python: `SEMICONTINUOUS`/`SEMIINTEGER` variable types, `Model.add_sos`,
  `Model.add_indicator(b, 1, x + y <= 5)` and `popso.read_lp`.
- How they are solved: the solver adds a binary variable for each semi-continuous variable, SOS
  member or SOS2 segment and links it with big-M rows; an indicator constraint becomes a big-M
  row. The big-M values come from the variable bounds, tightened by bound propagation over the
  model rows. Presolve, cuts, heuristics and branching then work as for any MILP, and the
  final solution is checked against the original constraints. When no finite bound exists (or
  it would exceed 1e6), the constraint is not relaxed by a big-M: the solver branches on its
  binary variable and enforces the constraint in that subtree. Presolve is skipped for such
  models.
- Testing: unit tests compare small models with optima found by enumeration, using both the
  big-M form and branching only, and `benchmarks/tools/logical_fuzz.py` compares random small
  models with SCIP. Three new generated benchmark families, `semilot` (minimum
  lot sizes), `pwlnet` (transportation with piecewise-linear concave costs, plant opening and
  single sourcing) and `indsched` (single-machine scheduling with indicator disjunctions and no
  big-M), have reference answers from SCIP (and from HiGHS for `semilot`). The benchmark
  checker verifies the new constraints.

## v0.2.1 — 2026-11

A patch release: faster root cuts on pure integer models and a robustness fix for deterministic
parallel mode under tight memory limits.

**Benchmarks**
- 63-instance MIPLIB 2017 subset (120 s, 4 threads): 60 solved, shifted geometric mean 11.3 s
  (v0.2.0: 12.4 s). See [BENCHMARKS.md](BENCHMARKS.md).

**Solver**
- Deterministic parallel mode (`--deterministic`) now drains the search tree reliably near the
  memory cap (`--mem-limit`). Before, a parallel deterministic run could stop with `MEMORY_LIMIT`
  where a 1-thread run under the same cap proved optimality. The work between two
  synchronization points is now sized so that the tree stays below the cap. Reruns stay
  identical, and runs that stay well below the cap are unchanged.

**Root cuts**
- c-MIR path aggregation continues through general integer columns on pure integer models (root cut loop,
  work-capped, not on dual degenerate root LPs): mzzv11 / mzzv42z are now solved at the root (4 threads:
  51 s -> 17 s and 15 s -> 3.3 s).
- Clique cuts are extended to maximal cliques by zero-valued literals.
- `benchmarks/tools/rootgap.py`: root gap closed per instance vs SCIP / HiGHS (node limit 1) and the known
  optimum; the solver prints a machine-readable `ROOTSTAT` line after each root.

## v0.2.0 — 2026-11

popso becomes usable as a library and on all three major platforms: a C API and Python bindings,
starting solutions, an LP-format reader, CI on Linux, macOS and Windows, a stronger presolve and a
more robust LP factorization.

**Benchmarks**
- 63-instance MIPLIB 2017 subset (120 s, 4 threads): 60 solved (v0.1.0: 59), shifted geometric
  mean 12.4 s (v0.1.0: 13.2 s); HiGHS solves 53, SCIP 50. 118 generated models (1 thread): 116
  solved. See [BENCHMARKS.md](BENCHMARKS.md).

**Solver**
- Stronger presolve:
  - clique extension: set-packing rows grow with conflicting binaries, dominated packing rows are
    dropped, and binaries that conflict with a whole equality row are fixed;
  - dual bound strengthening and GCD-based tightening of integer rows;
  - a set-cover rule for dominated columns;
  - faster probing that also tightens continuous bounds.

  On the MIPLIB subset the presolved models now close more of the root gap on average than
  SCIP's or HiGHS's presolve, and presolve takes about 30% less time than in v0.1.0.
- LP factorization: an update that would make the factorization inaccurate (a very large row
  multiplier in the Forrest–Tomlin update) now triggers a refactorization instead of carrying the
  error forward.

**Library API**
- C API in `include/popso.h`, built as the shared library `popso` (`libpopso.so`,
  `libpopso.dylib`, `popso.dll`). It can build models (variables, sparse rows, objective) or read
  MPS files, and it sets parameters, solves, and queries status, objective, bound, gap, node count,
  time and the solution. It can also write MPS and solution files. Calls return error codes, and
  `popso_last_error` gives a message. Separate models can be solved from different threads.
- Python bindings (`python/popso`), using `ctypes` with no compiled extension and no
  dependencies. They support linear expressions (`m.add_constr(x + 2 * y <= 5)`), a result object
  that can be indexed by variable or name, parameters, and starting solutions. The package can be
  installed with pip.

**Starting solutions**
- `--start <file>` reads a starting solution in the `name value` format that `--sol` writes. A
  complete feasible start becomes the first incumbent. A partial or infeasible start keeps its
  integer values and is completed by a short sub-MIP; if that fails it is ignored, and the log
  says so. Starts are also available in the C API and in Python.

**Input / output**
- CPLEX LP-format reader (`.lp` files, as written by CPLEX, HiGHS, SCIP and modelling
  tools): objective constant, ranged rows (`lo <= expr <= hi`), all bound forms
  (`free`, `-inf`/`+inf`, fixed, double-sided), General/Binary sections, comments and
  case-insensitive keywords. Quadratic terms, indicator constraints, semi-continuous
  variables and SOS are rejected with a clear "not supported yet" error and line number.
- The command line picks the reader by file extension (`.lp` vs `.mps`, case-insensitive);
  `--format mps|lp` overrides it.
- Testing: a pure-Python LP writer in the benchmark modeller and
  `benchmarks/tools/lp_crosscheck.py`, which checks that every benchmark instance read as
  LP gives exactly the same model as the MPS file, and solves LP files written by HiGHS
  and SCIP.

**Platforms**
- Continuous integration on GitHub Actions: every push and pull request is built and tested on
  Linux (gcc and clang), macOS (clang) and Windows (MSVC). Each build runs the unit tests, the
  C API and Python binding tests, and solves the tiny benchmark instances with independently
  verified answers. A separate Linux build runs the unit and C API tests under
  AddressSanitizer and UndefinedBehaviorSanitizer.
- Builds without warnings with gcc 15 and clang. The new CMake option `POPSO_WERROR` turns
  warnings into errors.
- Builds with gcc on macOS (the macOS system headers used for memory measurement did not
  compile with gcc).
- Windows/MSVC: `windows.h` no longer defines `min`/`max` macros, there are no warnings about
  standard C functions such as `getenv`, and the MSVC warning set matches the gcc/clang one.
- Benchmark runner: `benchmarks/run.py` needs only the Python standard library unless it runs
  the reference solvers. It finds `popso.exe` and `build/Release/` on Windows, works there
  without `ps`, and has a `--strict` flag that fails when an instance is not solved.

## v0.1.0 — 2026-11

First public release: a complete, multithreaded LP/MILP solver written from scratch in
C++20, using only the standard library.

**Solves**
- Linear programs (LP) and mixed-integer linear programs (MILP) with continuous, integer
  and binary variables, ranged rows and free variables.
- It proves optimality, infeasibility or unboundedness.
- On a 63-instance MIPLIB 2017 subset it solves 59 within 120 s on 4 threads, against 53
  for HiGHS and 50 for SCIP. On 118 generated models it solves 116 on 1 thread, against
  113 for HiGHS and 111 for SCIP. See [BENCHMARKS.md](BENCHMARKS.md).

**Input / output**
- MPS reader (free and fixed format, including RANGES, OBJSENSE and integer markers).
- Solution file output (`name value` per line), plus a summary of status, objective, bound,
  gap, nodes and time.
- Command-line interface with time, gap, node and memory limits, thread count, seed,
  deterministic mode and racing.

**LP engine**
- Bounded dual simplex.
- LU factorization with Forrest–Tomlin updates and hyper-sparse solves.
- Dual steepest-edge pricing, cost perturbation, and a primal simplex for cleanup.
- Warm starts from stored bases throughout the tree.

**Presolve**
- Removes singleton rows and columns, parallel rows and columns, and dominated columns.
- Tightens bounds and coefficients.
- Aggregates and substitutes variables.
- Probing with implications, and clique merging.
- Full postsolve; restarts after root fixings.

**Cutting planes**
- Gomory mixed-integer cuts and c-MIR cuts on aggregated rows (with variable-bound
  substitution).
- Lifted knapsack covers, also on aggregated rows.
- Clique, implied-bound and zero-half cuts.
- Tableau cuts at tree nodes.
- A dynamic cut pool with aging.

**Primal heuristics**
- Rounding, diving, local search.
- Feasibility pump and feasibility jump.
- Large-neighbourhood sub-MIPs (RINS, RENS, crossover) with a solution pool.

**Branch-and-bound**
- Branching: reliability branching on pseudocosts with strong branching.
- Pruning:
  - node-level domain propagation;
  - conflict analysis;
  - reduced-cost fixing.
- Symmetry: orbital fixing and packing orbitopes.
- Restarts at the root and in the tree.
- Node selection mixes best bound and best estimate.
- Recovery from numerical failures in node LPs.

**Parallelism and resources**
- Multithreaded tree search, opportunistic by default.
- A deterministic mode that gives the same result on every run with the same thread count.
- Search-tree memory stays within a configurable cap, using compact node storage; near the
  cap the solver finishes off open nodes rather than creating new ones.

**Known limitations**
- No quadratic or nonlinear models.
- No SOS or indicator constraints, and no semi-continuous variables.
- No LP-format reader, no library API, and no barrier/interior-point LP method.
