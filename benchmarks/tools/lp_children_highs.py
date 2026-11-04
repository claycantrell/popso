#!/usr/bin/env python3
"""HiGHS counterpart of lp_children.cpp: root LP (integrality relaxed, dual simplex), the children listed in
the candidates file (each warm from the root basis) and the dive chain, with iteration counts and times.
Usage: lp_children_highs.py <model.mps> <candidates file>"""
import sys
import time

import highspy

h = highspy.Highs()
h.setOptionValue("output_flag", False)
h.setOptionValue("threads", 1)
h.readModel(sys.argv[1])
lp = h.getLp()
n = lp.num_col_
h.changeColsIntegrality(n, list(range(n)), [highspy.HighsVarType.kContinuous] * n)
h.setOptionValue("solver", "simplex")
h.setOptionValue("simplex_strategy", 1)  # dual
t = time.perf_counter()
h.run()
print(f"root: obj {h.getInfo().objective_function_value:.10g} its {h.getInfo().simplex_iteration_count} "
      f"time {time.perf_counter() - t:.2f}")
root = h.getBasis()
lb0, ub0 = list(lp.col_lower_), list(lp.col_upper_)
tot = cnt = 0
dive = []
for line in open(sys.argv[2]):
    f = line.split()
    if f[0] == "D":
        dive.append((int(f[1]), float(f[2]), float(f[3])))
        continue
    j, d, v = int(f[0]), int(f[1]), float(f[2])
    h.setBasis(root)
    if d == 0:
        h.changeColBounds(j, lb0[j], v)
    else:
        h.changeColBounds(j, v, ub0[j])
    t = time.perf_counter()
    h.run()
    it = h.getInfo().simplex_iteration_count
    print(f"child col {j} {'up' if d else 'down'}: obj {h.getInfo().objective_function_value:.10g} its {it} "
          f"time {time.perf_counter() - t:.3f}")
    tot += it
    cnt += 1
    h.changeColBounds(j, lb0[j], ub0[j])
if cnt:
    print(f"children: {cnt}, mean its {tot / cnt:.0f}")
h.setBasis(root)
h.run()
for k, (j, l, u) in enumerate(dive):
    h.changeColBounds(j, l, u)
    t = time.perf_counter()
    h.run()
    print(f"dive {k} col {j}: obj {h.getInfo().objective_function_value:.10g} its "
          f"{h.getInfo().simplex_iteration_count} time {time.perf_counter() - t:.3f}")
