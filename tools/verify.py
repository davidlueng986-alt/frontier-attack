"""Run every self-check in the repo. Exit 0 only if all pass."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

CHECKS = [
    ("hadwiger_nelson / exact toolkit", [sys.executable, "problems/hadwiger_nelson/hn.py"]),
    ("jacobian / verifier", [sys.executable, "problems/jacobian/jac.py"]),
    ("bsd / gatekeeper", [sys.executable, "problems/bsd/bsd.py"]),
]

SMOKE = [
    ("hadwiger_nelson / search smoke",
     [sys.executable, "problems/hadwiger_nelson/search.py",
      "--mode", "glue", "--trials", "5", "--budget", "24",
      "--out", "runs/_smoke_hn.json"]),
    ("jacobian / enumeration smoke",
     [sys.executable, "problems/jacobian/enumerate.py"]),
]


def run(name, cmd):
    print(f"\n=== {name} ===")
    p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    print(p.stdout.strip())
    if p.returncode != 0:
        print(p.stderr.strip()[-2000:])
    return p.returncode == 0


if __name__ == "__main__":
    ok = True
    for name, cmd in CHECKS + SMOKE:
        ok &= run(name, cmd)
    print("\nRESULT:", "PASS" if ok else "FAIL")
    raise SystemExit(0 if ok else 1)
