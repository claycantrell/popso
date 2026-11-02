"""Thin wrappers around the reference solvers (HiGHS, SCIP).

Each `solve_<name>(path, time_limit, threads, want_solution)` returns a dict:
  status      one of lib.solfile statuses (optimal/infeasible/unbounded/
              inf_or_unbd/timelimit/unknown)
  objective   primal objective (None if no solution)
  bound       dual bound (None if unavailable)
  time        wall-clock seconds of the solve call (excl. reading)
  solution    {name: value} if want_solution and a solution exists
"""

from __future__ import annotations

import math
import time

from lib.solfile import INF_OR_UNBD, INFEASIBLE, OPTIMAL, TIMELIMIT, UNBOUNDED, UNKNOWN

# Reference solves aim at proven optimality: tight gaps.
REL_GAP = 1e-9
ABS_GAP = 1e-9


def _num(v):
    if v is None:
        return None
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) and abs(v) < 1e29 else None


def solve_highs(path, time_limit, threads=1, want_solution=False, mem_mb=None):
    import highspy

    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("time_limit", float(time_limit))
    h.setOptionValue("mip_rel_gap", REL_GAP)
    h.setOptionValue("mip_abs_gap", ABS_GAP)
    h.setOptionValue("threads", int(threads))
    st = h.readModel(str(path))
    if st == highspy.HighsStatus.kError:
        raise RuntimeError("HiGHS failed to read model")
    t0 = time.perf_counter()
    h.run()
    elapsed = time.perf_counter() - t0
    ms = h.getModelStatus()
    MS = highspy.HighsModelStatus
    if ms == MS.kUnboundedOrInfeasible:
        # disambiguate without presolve
        h.setOptionValue("presolve", "off")
        h.clearSolver()
        t1 = time.perf_counter()
        h.run()
        elapsed += time.perf_counter() - t1
        ms = h.getModelStatus()
    status = {
        MS.kOptimal: OPTIMAL,
        MS.kInfeasible: INFEASIBLE,
        MS.kUnbounded: UNBOUNDED,
        MS.kUnboundedOrInfeasible: INF_OR_UNBD,
        MS.kTimeLimit: TIMELIMIT,
    }.get(ms, UNKNOWN)
    info = h.getInfo()
    has_sol = info.primal_solution_status == 2  # kSolutionStatusFeasible
    lp = h.getLp()
    is_mip = any(t != highspy.HighsVarType.kContinuous for t in (lp.integrality_ or []))
    obj = _num(info.objective_function_value) if has_sol else None
    if status in (INFEASIBLE, UNBOUNDED, INF_OR_UNBD):
        obj = None
    if is_mip:
        bound = _num(info.mip_dual_bound)
    else:
        bound = obj if status == OPTIMAL else None
    res = {"status": status, "objective": obj, "bound": bound, "time": elapsed,
           "raw_status": h.modelStatusToString(ms)}
    if is_mip:
        res["nodes"] = int(info.mip_node_count)
    if want_solution and has_sol and obj is not None:
        sol = h.getSolution()
        names = [h.getColName(j)[1] for j in range(lp.num_col_)]
        res["solution"] = dict(zip(names, sol.col_value))
    return res


def solve_scip(path, time_limit, threads=1, want_solution=False, mem_mb=None):
    import pyscipopt

    m = pyscipopt.Model()
    m.hideOutput()
    m.readProblem(str(path))
    m.setParam("limits/time", float(time_limit))
    m.setParam("limits/gap", REL_GAP)
    m.setParam("limits/absgap", ABS_GAP)
    m.setParam("randomization/randomseedshift", 0)
    if mem_mb:
        m.setParam("limits/memory", float(mem_mb) * 0.8)
    t0 = time.perf_counter()
    m.optimize()
    elapsed = time.perf_counter() - t0
    raw = m.getStatus()
    status = {
        "optimal": OPTIMAL,
        "gaplimit": OPTIMAL,
        "infeasible": INFEASIBLE,
        "unbounded": UNBOUNDED,
        "inforunbd": INF_OR_UNBD,
        "timelimit": TIMELIMIT,
        "memlimit": "memlimit",
    }.get(raw, UNKNOWN)
    obj = None
    if m.getNSols() > 0 and status not in (INFEASIBLE, UNBOUNDED, INF_OR_UNBD):
        obj = _num(m.getObjVal())
    bound = None
    if status in (OPTIMAL, TIMELIMIT):
        bound = _num(m.getDualbound())
    res = {"status": status, "objective": obj, "bound": bound, "time": elapsed,
           "raw_status": raw, "nodes": int(m.getNNodes())}
    if want_solution and obj is not None:
        best = m.getBestSol()
        res["solution"] = {v.name: m.getSolVal(best, v) for v in m.getVars(transformed=False)}
    return res


SOLVERS = {"highs": solve_highs, "scip": solve_scip}


def available_solvers():
    out = []
    try:
        import highspy  # noqa: F401
        out.append("highs")
    except ImportError:
        pass
    try:
        import pyscipopt  # noqa: F401
        out.append("scip")
    except ImportError:
        pass
    return out
