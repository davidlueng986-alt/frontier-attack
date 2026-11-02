"""Rigid gluing over all vertex pairs -- Hadwiger-Nelson track.

Start from the Moser spindle.  Repeatedly pick a random pair of vertices in the
current graph whose exact squared distance equals that of a random pair of
vertices in the base spindle (tolerance 1e-30), then glue a rotated (never
reflected) copy of the base onto that pair and merge coincident points.
Every candidate is scored with the exact DSATUR chromatic number.

Usage:
    python problems/hadwiger_nelson/glue2.py --trials 20 --budget 40
"""

from __future__ import annotations

import argparse
import json
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import mpmath as mp

from search import chromatic_number, spindle_points, unit_edges, merge

mp.mp.dps = 80
SEARCH_TOL = mp.mpf("1e-30")


def _dist2(p, q):
    dx = p[0] - q[0]
    dy = p[1] - q[1]
    return dx * dx + dy * dy


def glue_any(base, dst, pair_a, pair_b):
    """Rigid motion (rotation + translation, no reflection) sending
    base[pair_a[0]] -> dst[pair_b[0]] and base[pair_a[1]] -> dst[pair_b[1]].

    The two pairwise distances must agree to within 1e-30.
    """
    a0, a1 = base[pair_a[0]], base[pair_a[1]]
    b0, b1 = dst[pair_b[0]], dst[pair_b[1]]
    da = _dist2(a0, a1)
    db = _dist2(b0, b1)
    if abs(da - db) > SEARCH_TOL:
        raise ValueError(
            "pair distances differ: "
            f"{mp.nstr(da, 20)} vs {mp.nstr(db, 20)}"
        )
    theta = mp.atan2(b1[1] - b0[1], b1[0] - b0[0]) - mp.atan2(
        a1[1] - a0[1], a1[0] - a0[0]
    )
    c, s = mp.cos(theta), mp.sin(theta)
    ax, ay = a0
    bx, by = b0
    out = []
    for x, y in base:
        dx, dy = x - ax, y - ay
        out.append((bx + c * dx - s * dy, by + s * dx + c * dy))
    return out


def _base_pairs(base):
    return [
        (i, j, _dist2(base[i], base[j]))
        for i in range(len(base))
        for j in range(i + 1, len(base))
    ]


def _grow(base, pairs, budget: int, rng: random.Random):
    """Grow a copy-union of `base` up to `budget` vertices."""
    pts = [tuple(p) for p in base]
    guard = 0
    while len(pts) < budget and guard < 400:
        guard += 1
        i, j = sorted(rng.sample(range(len(pts)), 2))
        d2 = _dist2(pts[i], pts[j])
        matches = [bp for bp in pairs if abs(bp[2] - d2) <= SEARCH_TOL]
        if not matches:
            continue
        a, b, _ = rng.choice(matches)
        if rng.random() < 0.5:
            pair_a, pair_b = (a, b), (i, j)
        else:
            pair_a, pair_b = (b, a), (j, i)
        cand = glue_any(base, pts, pair_a, pair_b)
        merged = merge(pts, cand)
        if len(pts) < len(merged) <= budget:
            pts = merged
    return pts


def _score(pts):
    edges = unit_edges(pts)
    n = len(pts)
    adj = [0] * n
    for i, j in edges:
        adj[i] |= 1 << j
        adj[j] |= 1 << i
    return chromatic_number(n, adj), len(edges)


def _record(pts, chi, m):
    return {
        "n": len(pts),
        "m": m,
        "chi": chi,
        "points": [[mp.nstr(x, 60), mp.nstr(y, 60)] for x, y in pts],
    }


def run(trials: int, budget: int, seed: int, out: Path) -> None:
    rng = random.Random(seed)
    base = spindle_points()
    pairs = _base_pairs(base)
    best = None
    best_key = None
    hist: dict[int, int] = {}
    t0 = time.time()

    for t in range(trials):
        pts = _grow(base, pairs, budget, rng)
        chi, m = _score(pts)
        hist[chi] = hist.get(chi, 0) + 1
        key = (chi, -len(pts))
        if best_key is None or key > best_key:
            best_key = key
            best = _record(pts, chi, m)
            print(
                f"[{t:4d}] best chi={chi} n={len(pts)} m={m} "
                f"({time.time() - t0:.1f}s)",
                flush=True,
            )

    if best is None:
        base = [tuple(p) for p in base]
        chi, m = _score(base)
        best = _record(base, chi, m)

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(
            {"best": best, "histogram": hist, "mode": "glue2"},
            indent=2,
        ),
        encoding="utf-8",
    )
    print("histogram:", dict(sorted(hist.items())))
    print("best:", {k: best[k] for k in ("n", "m", "chi")})
    print("written:", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=20)
    ap.add_argument("--budget", type=int, default=40)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="runs/hn_glue2.json")
    a = ap.parse_args()
    run(a.trials, a.budget, a.seed, Path(a.out))
