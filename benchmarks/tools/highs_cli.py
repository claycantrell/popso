#!/usr/bin/env python3
"""Stand-in solver implementing our CLI contract on top of HiGHS.

    highs_cli.py <model.mps> [--sol <file>] [--time-limit <sec>] [--threads <n>] [--mem-limit <MB>] [--gap <rel>]

Used to test benchmarks/run.py before (and alongside) the real C++ solver,
e.g.  run.py --solver "benchmarks/tools/highs_cli.py".

Set FAKE_BUG to inject a defect and check that run.py catches it loudly:
  FAKE_BUG=perturb      shift one variable value by +0.5 (infeasible/non-int)
  FAKE_BUG=objective    report objective + 1
  FAKE_BUG=infeasible   claim INFEASIBLE for every model
  FAKE_BUG=suboptimal   stop at 5% gap but claim OPTIMAL
FAKE_BUG_MATCH=<substring> restricts the bug to matching model paths.

HiGHS' MPS reader differs from the CPLEX convention on negative UP bounds
(see generators/lp_sanity.py HAND_DERIVED_ONLY), so this stand-in is expected to get
lpsanity_negative_ub_0 / lpsanity_mip_negative_ub_0 wrong.
"""

import argparse
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from lib.solfile import write_sol  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("model")
    ap.add_argument("--sol")
    ap.add_argument("--time-limit", type=float, default=math.inf)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--mem-limit", type=int, default=0, help="MB (accepted, not enforced by HiGHS)")
    ap.add_argument("--gap", type=float, default=1e-9, help="relative MIP gap")
    a = ap.parse_args()

    import highspy

    bug = os.environ.get("FAKE_BUG", "")
    if bug and os.environ.get("FAKE_BUG_MATCH", "") not in a.model:
        bug = ""
    h = highspy.Highs()
    h.setOptionValue("output_flag", True)
    if math.isfinite(a.time_limit):
        h.setOptionValue("time_limit", a.time_limit)
    h.setOptionValue("threads", a.threads)
    h.setOptionValue("mip_rel_gap", 0.05 if bug == "suboptimal" else a.gap)
    if h.readModel(a.model) == highspy.HighsStatus.kError:
        print("cannot read model", file=sys.stderr)
        return 3
    h.run()
    MS = highspy.HighsModelStatus
    ms = h.getModelStatus()
    if ms == MS.kUnboundedOrInfeasible:
        h.setOptionValue("presolve", "off")
        h.clearSolver()
        h.run()
        ms = h.getModelStatus()
    status = {MS.kOptimal: "OPTIMAL", MS.kInfeasible: "INFEASIBLE", MS.kUnbounded: "UNBOUNDED",
              MS.kUnboundedOrInfeasible: "INF_OR_UNBD", MS.kTimeLimit: "TIME_LIMIT"}.get(ms, "UNKNOWN")
    info = h.getInfo()
    name = Path(a.model).name.split(".")[0]
    values = None
    obj = None
    if info.primal_solution_status == 2 and status in ("OPTIMAL", "TIME_LIMIT"):
        lp = h.getLp()
        sol = h.getSolution()
        values = {h.getColName(j)[1]: v for j, v in enumerate(sol.col_value)}
        obj = info.objective_function_value
    if bug == "perturb" and values:
        k = next(iter(values))
        values[k] += 0.5
    if bug == "objective" and obj is not None:
        obj += 1
    if bug == "infeasible":
        status, values, obj = "INFEASIBLE", None, None
    if a.sol:
        write_sol(a.sol, name, status, obj, values)
    print(f"status {status} objective {obj}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
