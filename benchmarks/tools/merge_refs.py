#!/usr/bin/env python3
"""Merge reference results (benchmarks/reference.py output files, e.g. computed on another machine with --out)
into benchmarks/reference/solutions.json: each instance of the source replaces the destination's entry (its file
path rewritten to the repository-relative one of the destination when the instance exists there).

    .venv/bin/python benchmarks/tools/merge_refs.py benchmarks/results/remote-.../qc_scip_ref3.json [more.json ...]
        [--into benchmarks/reference/solutions.json] [--keep NAME,NAME]

--keep: instances whose destination entry stays (e.g. a reference the source solve could not reproduce).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parent.parent


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sources", nargs="+")
    ap.add_argument("--into", default=str(BENCH / "reference" / "solutions.json"))
    ap.add_argument("--keep", default="", help="comma list of instances not to replace")
    args = ap.parse_args(argv)
    keep = set(filter(None, args.keep.split(",")))
    dst = json.loads(Path(args.into).read_text())
    added = replaced = kept = 0
    for src in args.sources:
        s = json.loads(Path(src).read_text())
        for name, rec in s.get("instances", {}).items():
            if name in keep and name in dst["instances"]:
                kept += 1
                continue
            old = dst["instances"].get(name)
            if old is not None and "file" in old:
                rec["file"] = old["file"]
            replaced += old is not None
            added += old is None
            dst["instances"][name] = rec
    tmp = Path(args.into).with_suffix(".json.tmp")
    tmp.write_text(json.dumps(dst, indent=1, sort_keys=False) + "\n")
    os.replace(tmp, args.into)
    print(f"{args.into}: {added} added, {replaced} replaced, {kept} kept")
    return 0


if __name__ == "__main__":
    sys.exit(main())
