"""CLI tests of the user-model tooling: popso --write-anon / --anon-map / --anon-seed /
--anon-keep-order / --solve and --profile (ctest "cli_tooling").

usage: python3 test_cli_tooling.py <popso executable> <source dir>
"""
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest

POPSO = None
SOURCE = None


def run(*args):
    return subprocess.run([POPSO, *args], capture_output=True, text=True, timeout=120)


def read(path):
    with open(path) as f:
        return f.read()


def instance(name):
    return os.path.join(SOURCE, "benchmarks", "instances", name + ".mps")


def mps_names(path):
    """Row (other than N rows) and column names of an MPS file."""
    rows, cols, section = set(), set(), None
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("*"):
                continue
            if not line[0].isspace():
                section = line.split()[0]
                continue
            tok = line.split()
            if section == "ROWS" and tok[0].upper() != "N":
                rows.add(tok[1])
            elif section == "COLUMNS" and "'MARKER'" not in tok:
                cols.add(tok[0])
    return rows, cols


class CliTooling(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def path(self, name):
        return os.path.join(self.dir, name)

    def test_write_anon_exits_without_solving(self):
        src = instance("pwlnet_tiny_1")
        out, mp = self.path("a.mps"), self.path("a.map")
        r = run(src, "--write-anon", out, "--anon-map", mp)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertNotIn("Status:", r.stdout)
        text = read(out)
        tokens = set(text.split())
        rows, cols = mps_names(src)
        generic = re.compile(r"^[xcsi]\d+$")
        for name in rows | cols:
            if not generic.match(name):
                self.assertNotIn(name, tokens)
        self.assertTrue(text.startswith("NAME anon\n"))
        # The map: a comment line, then "<new> <original>" for every column, row, SOS set and indicator.
        lines = read(mp).splitlines()
        self.assertTrue(lines[0].startswith("#"))
        pairs = [line.split() for line in lines[1:]]
        self.assertTrue(all(len(p) == 2 for p in pairs))
        originals = {p[1] for p in pairs}
        self.assertTrue(rows <= originals)
        self.assertTrue(cols <= originals)
        new_cols = [p[0] for p in pairs if p[0].startswith("x")]
        self.assertEqual(len(new_cols), len(cols))

    def test_seed_and_keep_order(self):
        src = instance("setcover_tiny_1")
        a, b, c, k = (self.path(n) for n in ("a.mps", "b.mps", "c.mps", "k.mps"))
        self.assertEqual(run(src, "--write-anon", a, "--verbosity", "0").returncode, 0)
        self.assertEqual(run(src, "--write-anon", b, "--verbosity", "0").returncode, 0)
        self.assertEqual(run(src, "--write-anon", c, "--anon-seed", "17", "--verbosity", "0").returncode, 0)
        kmap = self.path("k.map")
        self.assertEqual(run(src, "--write-anon", k, "--anon-keep-order", "--anon-map", kmap).returncode, 0)
        self.assertEqual(read(a), read(b))
        self.assertNotEqual(read(a), read(c))
        # Keep order: x<j> is the j-th column of the original.
        _, cols = mps_names(src)
        with open(src) as f:
            order = []
            section = None
            for line in f:
                if line.strip() and not line[0].isspace():
                    section = line.split()[0]
                elif section == "COLUMNS" and "'MARKER'" not in line and line.split()[0] not in order:
                    order.append(line.split()[0])
        pairs = [line.split() for line in read(kmap).splitlines()[1:]]
        self.assertEqual([p[1] for p in pairs if p[0].startswith("x")], order)

    def test_option_rules(self):
        src = instance("knapsack_tiny_1")
        r = run(src, "--write-anon", self.path("a.mps"), "--profile", self.path("p.json"))
        self.assertEqual(r.returncode, 2)  # needs --solve
        self.assertIn("--solve", r.stderr)
        r = run(src, "--anon-map", self.path("m"))
        self.assertEqual(r.returncode, 2)  # needs --write-anon
        r = run(src, "--write-anon", self.path("a.mps"), "--solve", "--profile", self.path("p.json"),
                "--threads", "1")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("Status: OPTIMAL", r.stdout)
        self.assertTrue(os.path.exists(self.path("a.mps")))
        json.loads(read(self.path("p.json")))

    def check_profile(self, src, *extra):
        out = self.path("p.json")
        r = run(src, "--profile", out, "--threads", "1", "--verbosity", "0", *extra)
        self.assertEqual(r.returncode, 0, r.stderr)
        text = read(out)
        js = json.loads(text)
        self.assertEqual(js["profile_version"], 1)
        for key in ("solver", "model", "numerics", "presolve", "run"):
            self.assertIn(key, js)
        self.assertEqual(sum(js["model"]["row_classes"].values()), js["model"]["rows"])
        self.assertNotIn(os.path.basename(src), text)
        rows, cols = mps_names(src)
        strings = set()

        def walk(v):
            if isinstance(v, str):
                strings.add(v)
            elif isinstance(v, list):
                for x in v:
                    walk(x)
            elif isinstance(v, dict):
                for k, x in v.items():
                    strings.add(k)
                    walk(x)
        walk(js)
        self.assertFalse((rows | cols) & strings)
        return js

    def test_profile_milp_lp_qp(self):
        js = self.check_profile(instance("gap_tiny_1"))
        self.assertEqual(js["run"]["status"], "OPTIMAL")
        self.assertTrue(js["run"]["incumbents"])
        self.assertTrue(all(e["source"] for e in js["run"]["incumbents"]))
        js = self.check_profile(instance("lpsanity_ranges_0"))
        self.assertEqual(js["run"]["status"], "OPTIMAL")
        js = self.check_profile(os.path.join(SOURCE, "tests", "data", "maros_meszaros", "HS35.mps"))
        self.assertGreater(js["model"]["quadratic"]["nonzeros"], 0)

    def test_profile_on_time_limit(self):
        js = self.check_profile(instance("cflp_medium_1"), "--time-limit", "0.05")
        self.assertEqual(js["run"]["status"], "TIME_LIMIT")
        self.assertEqual(js["solver"]["time_limit"], 0.05)


if __name__ == "__main__":
    POPSO, SOURCE = sys.argv[1], sys.argv[2]
    unittest.main(argv=sys.argv[:1])
