"""Read/write "name value" .sol files (the CLI contract of our solver).

    # Solution for model <name>
    # Status = OPTIMAL
    # Objective value = <value>
    <varname> <value>
    ...
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

# Canonical status strings used throughout the benchmark tooling.
OPTIMAL = "optimal"
INFEASIBLE = "infeasible"
UNBOUNDED = "unbounded"
INF_OR_UNBD = "inf_or_unbd"
TIMELIMIT = "timelimit"
MEMLIMIT = "memlimit"
UNKNOWN = "unknown"

_STATUS_MAP = {
    "OPTIMAL": OPTIMAL,
    "INFEASIBLE": INFEASIBLE,
    "UNBOUNDED": UNBOUNDED,
    "INF_OR_UNBD": INF_OR_UNBD,
    "INFEASIBLE_OR_UNBOUNDED": INF_OR_UNBD,
    "TIME_LIMIT": TIMELIMIT,
    "TIMELIMIT": TIMELIMIT,
    "MEMORY_LIMIT": MEMLIMIT,
    "MEM_LIMIT": MEMLIMIT,
}


@dataclass
class Solution:
    status: str = UNKNOWN
    raw_status: str = ""
    objective: float | None = None
    model_name: str = ""
    values: dict = field(default_factory=dict)
    errors: list = field(default_factory=list)


def read_sol(path) -> Solution:
    s = Solution()
    with open(path, "rt") as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                body = line[1:].strip()
                low = body.lower()
                if low.startswith("solution for model"):
                    s.model_name = body[len("solution for model"):].strip()
                elif low.startswith("status"):
                    val = body.split("=", 1)[1].strip() if "=" in body else ""
                    s.raw_status = val
                    s.status = _STATUS_MAP.get(val.upper().replace(" ", "_"), UNKNOWN)
                elif low.startswith("objective value"):
                    try:
                        s.objective = float(body.split("=", 1)[1])
                    except (IndexError, ValueError):
                        s.errors.append(f"line {lineno}: bad objective line")
                continue
            toks = line.split()
            if len(toks) != 2:
                s.errors.append(f"line {lineno}: expected '<name> <value>'")
                continue
            try:
                v = float(toks[1])
            except ValueError:
                s.errors.append(f"line {lineno}: bad value {toks[1]!r}")
                continue
            if toks[0] in s.values:
                s.errors.append(f"line {lineno}: duplicate variable {toks[0]}")
            s.values[toks[0]] = v
    return s


def write_sol(path, model_name, status, objective=None, values=None):
    """Write a .sol file. `status` is a contract string such as OPTIMAL."""
    with open(path, "wt") as f:
        f.write(f"# Solution for model {model_name}\n")
        f.write(f"# Status = {status}\n")
        if objective is not None and math.isfinite(objective):
            f.write(f"# Objective value = {objective:.17g}\n")
        if values:
            for name, v in values.items():
                f.write(f"{name} {v:.17g}\n")
