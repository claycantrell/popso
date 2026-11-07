"""Primal gap and primal integral (Berthold, "Measuring the impact of primal heuristics", Oper. Res.
Letters 41 (2013)) of a run's incumbent trace, used by run.py --trace and ab.py --primal."""

from __future__ import annotations


def pgap(z, ref):
    """Primal gap of objective z to the reference z*: 0 at z*, 1 without a solution or when the signs
    differ, else |z - z*| / max(|z|, |z*|)."""
    if z is None or ref is None:
        return 1.0
    if abs(z - ref) <= 1e-9 * max(1.0, abs(ref)):
        return 0.0
    if z * ref < 0:
        return 1.0
    return abs(z - ref) / max(abs(z), abs(ref))


def trace_points(incumbents):
    """[[time, objective], ...] from a popso --profile "incumbents" list (an event's first and, for a run
    of improvements, its last point)."""
    pts = []
    for e in incumbents:
        pts.append([e["time"], e["objective"]])
        if e.get("count", 1) > 1 and e.get("last_time") is not None:
            pts.append([e["last_time"], e["last_objective"]])
    pts.sort(key=lambda p: p[0])
    return pts


def integral(points, ref, horizon):
    """Primal integral over [0, horizon] divided by horizon (the mean primal gap, in [0, 1]): gap 1 until
    the first incumbent, then the gap of the incumbent at each time. The gap never increases: incumbents
    only improve, and one beyond a reference that is not optimal keeps the gap before it."""
    if ref is None or horizon <= 0:
        return None
    total, t_prev, g_prev = 0.0, 0.0, 1.0
    for t, z in points:
        t = min(max(t, 0.0), horizon)
        total += g_prev * (t - t_prev)
        t_prev, g_prev = t, min(g_prev, pgap(z, ref))
    total += g_prev * (horizon - t_prev)
    return total / horizon
