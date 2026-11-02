#!/usr/bin/env python3
"""SVG plot "instances solved vs. time" (popso vs. the reference solvers).

One step line per solver: y = number of instances solved within t seconds,
x = t on a log scale.  Higher and further left is better; the right end of a
line is the solver's solved count.  Our times come from one run.py result
file, the reference times from reference/solutions.json (+ miplib_solu.json).
The SVG is self-contained (no dependencies) and follows the viewer's light /
dark preference.

    .venv/bin/python benchmarks/tools/plot_profile.py benchmarks/results/<run>.json \
        --title "MIPLIB 2017 subset (63 instances)" --out docs/miplib.svg
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH))

from run import load_references  # noqa: E402

SOLVED = ("optimal", "infeasible", "unbounded", "inf_or_unbd")
NAMES = {"ours": "popso", "highs": "HiGHS", "scip": "SCIP"}
# categorical slots 1-3 (validated light / dark), fixed per solver
COLORS = {"ours": ("#2a78d6", "#3987e5"), "highs": ("#eb6834", "#d95926"), "scip": ("#1baf7a", "#199e70")}


def solve_times(result_path, solvers, refs):
    d = json.loads(Path(result_path).read_text())
    recs = d["results"] if isinstance(d, dict) else d
    times = {s: [] for s in solvers}
    for r in recs:
        n = r["instance"]
        ref = refs.get(n, {})
        ok = r["class"] in ("ok", "noref")
        if r["class"].lower() == "wrong" and r.get("status") == "optimal" and ref.get("status") == "optimal":
            z, obj = ref.get("objective"), r.get("objective")
            ok = z is not None and obj is not None and abs(obj - z) <= 1e-6 * max(1.0, abs(z))
        times["ours"].append(r["time"] if ok else None)
        for s in solvers[1:]:
            t = ref.get("times", {}).get(s)
            times[s].append(t[0] if t and t[0] is not None and t[1] in SOLVED else None)
    return times, len(recs)


def svg(times, n, title, subtitle, tl):
    W, H = 720, 400
    left, right, top, bottom = 56, 118, 58, 46
    pw, ph = W - left - right, H - top - bottom
    tmin = 0.01
    lx0, lx1 = math.log10(tmin), math.log10(tl)

    def X(t):
        return left + (math.log10(max(t, tmin)) - lx0) / (lx1 - lx0) * pw

    def Y(k):
        return top + ph - k / n * ph

    o = []
    o.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
             f'role="img" aria-label="{title}: instances solved versus time">')
    css = ["svg{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif}",
           ".bg{fill:#fcfcfb}.t1{fill:#0b0b0b}.t2{fill:#52514e}.grid{stroke:#e4e3df}.axis{stroke:#a3a29c}"]
    css += [f".s-{s}{{stroke:{c[0]}}}.f-{s}{{fill:{c[0]}}}" for s, c in COLORS.items()]
    dark = [".bg{fill:#1a1a19}.t1{fill:#f0efec}.t2{fill:#c3c2b7}.grid{stroke:#333331}.axis{stroke:#6b6a64}"]
    dark += [f".s-{s}{{stroke:{c[1]}}}.f-{s}{{fill:{c[1]}}}" for s, c in COLORS.items()]
    o.append("<style>" + "".join(css) + "@media (prefers-color-scheme: dark){" + "".join(dark) + "}</style>")
    o.append(f'<rect class="bg" width="{W}" height="{H}" rx="8"/>')
    o.append(f'<text class="t1" x="{left}" y="24" font-size="15" font-weight="600">{title}</text>')
    o.append(f'<text class="t2" x="{left}" y="42" font-size="12">{subtitle}</text>')

    # grid + axes (recessive)
    step = 10 if n > 40 else 5
    for k in range(0, n + 1, step):
        y = Y(k)
        o.append(f'<line class="grid" x1="{left}" x2="{left + pw}" y1="{y:.1f}" y2="{y:.1f}" stroke-width="1"/>')
        o.append(f'<text class="t2" x="{left - 8}" y="{y + 4:.1f}" font-size="11" text-anchor="end">{k}</text>')
    for t, lab in ((0.01, "0.01"), (0.1, "0.1"), (1, "1"), (10, "10"), (tl, f"{tl:g}")):
        x = X(t)
        o.append(f'<line class="grid" x1="{x:.1f}" x2="{x:.1f}" y1="{top}" y2="{top + ph}" stroke-width="1"/>')
        o.append(f'<text class="t2" x="{x:.1f}" y="{top + ph + 18}" font-size="11" text-anchor="middle">{lab}</text>')
    o.append(f'<line class="axis" x1="{left}" x2="{left + pw}" y1="{top + ph}" y2="{top + ph}" stroke-width="1"/>')
    o.append(f'<text class="t2" x="{left + pw / 2}" y="{H - 8}" font-size="12" text-anchor="middle">'
             f'time limit per instance (seconds, log scale)</text>')
    o.append(f'<text class="t2" transform="translate(16 {top + ph / 2}) rotate(-90)" font-size="12" '
             f'text-anchor="middle">instances solved</text>')

    # step lines, direct labels at the right end (nudged apart to avoid collisions)
    ends = []
    for s, ts in times.items():
        done = sorted(t for t in ts if t is not None and t < tl)
        pts = [(X(tmin), Y(0))]
        k = 0
        for t in done:
            x = X(t)
            pts.append((x, Y(k)))
            k += 1
            pts.append((x, Y(k)))
        pts.append((X(tl), Y(k)))
        d = "M" + " L".join(f"{x:.1f} {y:.1f}" for x, y in pts)
        o.append(f'<path class="s-{s}" d="{d}" fill="none" stroke-width="2" stroke-linejoin="round">'
                 f'<title>{NAMES[s]}: {k}/{n} solved</title></path>')
        ends.append([Y(k), s, k])
    ends.sort()
    for i in range(1, len(ends)):
        if ends[i][0] - ends[i - 1][0] < 15:
            ends[i][0] = ends[i - 1][0] + 15
    for y, s, k in ends:
        x = X(tl) + 8
        o.append(f'<circle class="f-{s}" cx="{x + 4}" cy="{y:.1f}" r="4"/>')
        o.append(f'<text class="t1" x="{x + 13}" y="{y + 4:.1f}" font-size="12">{NAMES[s]} '
                 f'<tspan class="t2">{k}/{n}</tspan></text>')
    o.append("</svg>")
    return "\n".join(o) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("result", help="run.py result .json (our times)")
    ap.add_argument("--title", required=True)
    ap.add_argument("--subtitle", default="")
    ap.add_argument("--time-limit", type=float, default=120.0)
    ap.add_argument("--solvers", default="popso,highs,scip")
    ap.add_argument("--reference", default=str(BENCH / "reference" / "solutions.json"))
    ap.add_argument("--miplib-solu", default=str(BENCH / "reference" / "miplib_solu.json"))
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    solvers = ["ours" if s.strip() == "popso" else s.strip() for s in args.solvers.split(",") if s.strip()]
    refs = load_references(args.reference, args.miplib_solu)
    times, n = solve_times(args.result, solvers, refs)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(svg(times, n, args.title, args.subtitle, args.time_limit))
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
