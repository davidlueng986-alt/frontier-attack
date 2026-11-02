"""Portfolio search over unit-distance graphs -- Hadwiger-Nelson track.

Two generators
--------------
`circle`  Grow a graph by intersecting unit circles: every new vertex is at
          unit distance from two old ones.  This reproduces the classical
          construction method (Moser spindle, Golomb graph, de Grey graph).

`glue`    Start from the Moser spindle and repeatedly identify one of its unit
          edges with a unit edge of the current graph.  Rigid gluing along an
          edge preserves every existing distance and creates new cross edges,
          which is the cheapest way to raise density without breaking the
          unit-distance property.

Search runs at 80 digits; retained candidates are re-verified at 120 digits.
A hit at chi=5 would be the first small 5-chromatic unit-distance graph; a hit
at chi=7 would settle Hadwiger-Nelson as 7, because a finite 7-chromatic
unit-distance graph exists iff the plane is not 6-colourable (de Bruijn-Erdos).

Usage:
    python problems/hadwiger_nelson/search.py --mode glue --trials 200
"""

from __future__ import annotations

import argparse
import json
import random
import time
from pathlib import Path

import mpmath as mp

mp.mp.dps = 80
PI = mp.pi
SEARCH_TOL = mp.mpf("1e-30")
VERIFY_TOL = mp.mpf("1e-50")


def _bits(x: int) -> int:
    return bin(x).count("1")


def chromatic_number(n: int, adj: list[int]) -> int:
    """Exact DSATUR branch-and-bound on bitmask adjacency."""
    if n == 0:
        return 0
    order = sorted(range(n), key=lambda v: -_bits(adj[v]))
    col = {}
    for v in order:
        used = {col[u] for u in range(n) if ((adj[v] >> u) & 1) and u in col}
        c = 0
        while c in used:
            c += 1
        col[v] = c
    best = [max(col.values()) + 1]
    color = [-1] * n
    sat = [0] * n

    def rec(colored: int, used: int) -> None:
        if used >= best[0]:
            return
        if colored == n:
            best[0] = used
            return
        v = max(
            (i for i in range(n) if color[i] < 0),
            key=lambda i: (_bits(sat[i]), _bits(adj[i])),
        )
        for c in range(used + 1):
            if (sat[v] >> c) & 1:
                continue
            color[v] = c
            touched = []
            for u in range(n):
                if ((adj[v] >> u) & 1) and color[u] < 0 and not ((sat[u] >> c) & 1):
                    sat[u] |= 1 << c
                    touched.append(u)
            rec(colored + 1, max(used, c + 1))
            for u in touched:
                sat[u] &= ~(1 << c)
            color[v] = -1

    rec(0, 0)
    return best[0]


def edge_masks(points, tol=SEARCH_TOL):
    n = len(points)
    adj = [0] * n
    for i in range(n):
        xi, yi = points[i]
        for j in range(i + 1, n):
            xj, yj = points[j]
            if abs((xi - xj) ** 2 + (yi - yj) ** 2 - 1) <= tol:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    return adj


def unit_edges(points, tol=SEARCH_TOL):
    n = len(points)
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            xi, yi = points[i]
            xj, yj = points[j]
            if abs((xi - xj) ** 2 + (yi - yj) ** 2 - 1) <= tol:
                out.append((i, j))
    return out


# --------------------------------------------------------------------------
# generators
# --------------------------------------------------------------------------

def circle_intersections(p1, p2):
    x1, y1 = p1
    x2, y2 = p2
    dx, dy = x2 - x1, y2 - y1
    d2 = dx * dx + dy * dy
    d = mp.sqrt(d2)
    if d == 0 or d > 2:
        return []
    a = d2 / 2 / d
    h2 = 1 - a * a
    if h2 < 0:
        return []
    h = mp.sqrt(h2)
    xm, ym = x1 + a * dx / d, y1 + a * dy / d
    rx, ry = -dy / d * h, dx / d * h
    return [(xm + rx, ym + ry), (xm - rx, ym - ry)]


def grow(budget: int, rng: random.Random, leaf_p: float = 0.15):
    pts = [(mp.mpf(0), mp.mpf(0)), (mp.mpf(1), mp.mpf(0))]
    guard = 0
    while len(pts) < budget and guard < budget * 40:
        guard += 1
        if rng.random() < leaf_p:
            x, y = rng.choice(pts)
            th = rng.uniform(0, 2 * float(PI))
            cand = [(x + mp.cos(th), y + mp.sin(th))]
        else:
            i, j = rng.sample(range(len(pts)), 2)
            cand = circle_intersections(pts[i], pts[j])
            if not cand:
                continue
        p = rng.choice(cand)
        if all(abs(p[0] - q[0]) + abs(p[1] - q[1]) > mp.mpf("1e-25") for q in pts):
            pts.append(p)
    return pts


def spindle_points():
    s3, s11, s33 = mp.sqrt(3), mp.sqrt(11), mp.sqrt(33)
    return [
        (mp.mpf(0), mp.mpf(0)),
        (s3 / 2, mp.mpf(1) / 2),
        (s3 / 2, -mp.mpf(1) / 2),
        (s3, mp.mpf(0)),
        (5 * s3 / 12 - s11 / 12, mp.mpf(5) / 12 + s33 / 12),
        (5 * s3 / 12 + s11 / 12, -mp.mpf(5) / 12 + s33 / 12),
        (5 * s3 / 6, s33 / 6),
    ]


def rigid_from_two(base, i, j, t1, t2):
    """Rigid motion sending base[i] -> t1 and base[j] -> t2."""
    ax, ay = base[i]
    bx, by = base[j]
    vx, vy = bx - ax, by - ay
    wx, wy = t2[0] - t1[0], t2[1] - t1[1]
    theta = mp.atan2(wy, wx) - mp.atan2(vy, vx)
    c, s = mp.cos(theta), mp.sin(theta)
    out = []
    for x, y in base:
        dx, dy = x - ax, y - ay
        out.append((t1[0] + c * dx - s * dy, t1[1] + s * dx + c * dy))
    return out


def merge(pts, new, eps=mp.mpf("1e-25")):
    out = list(pts)
    for p in new:
        if all(abs(p[0] - q[0]) + abs(p[1] - q[1]) > eps for q in out):
            out.append(p)
    return out


def glue_grow(copies: int, rng: random.Random, max_vertices: int = 90):
    """Glue rotated copies of the spindle along existing unit edges."""
    base = spindle_points()
    base_edges = unit_edges(base)
    pts = list(base)
    for _ in range(copies):
        if len(pts) >= max_vertices:
            break
        edges = unit_edges(pts)
        if not edges:
            continue
        i, j = rng.choice(edges)
        a, b = rng.choice(base_edges)
        t1, t2 = (pts[i], pts[j]) if rng.random() < 0.5 else (pts[j], pts[i])
        cand = rigid_from_two(base, a, b, t1, t2)
        # only accept if the shared edge really is shared (guards against drift)
        if abs((cand[a][0] - t1[0]) ** 2 + (cand[a][1] - t1[1]) ** 2) > SEARCH_TOL:
            continue
        pts = merge(pts, cand)
    return pts


# --------------------------------------------------------------------------
# verification
# --------------------------------------------------------------------------

def verify(points):
    """Re-check a candidate at higher precision and report (chi, edges)."""
    old = mp.mp.dps
    mp.mp.dps = 120
    try:
        pts = [(mp.mpf(str(x)), mp.mpf(str(y))) for x, y in points]
        adj = edge_masks(pts, VERIFY_TOL)
        return chromatic_number(len(pts), adj), len(unit_edges(pts, VERIFY_TOL))
    finally:
        mp.mp.dps = old


def run(trials: int, budget: int, seed: int, out: Path, mode: str):
    rng = random.Random(seed)
    best = None
    hist: dict[int, int] = {}
    t0 = time.time()
    for t in range(trials):
        pts = glue_grow(max(1, budget // 6), rng) if mode == "glue" else grow(budget, rng)
        adj = edge_masks(pts)
        chi = chromatic_number(len(pts), adj)
        hist[chi] = hist.get(chi, 0) + 1
        if best is None or (chi, -len(pts)) > best[0]:
            chi2, m2 = verify(pts)
            best = ((chi, -len(pts)), {
                "n": len(pts), "m": m2, "chi": chi2,
                "points": [[mp.nstr(x, 40), mp.nstr(y, 40)] for x, y in pts],
            })
            print(f"[{t:4d}] best chi={chi2} n={len(pts)} m={m2} "
                  f"({time.time()-t0:.1f}s)")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"mode": mode, "histogram": hist, "best": best[1]},
                              indent=2), encoding="utf-8")
    print("histogram:", dict(sorted(hist.items())))
    print("best:", {k: best[1][k] for k in ("n", "m", "chi")})
    print("written:", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=100)
    ap.add_argument("--budget", type=int, default=30)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--mode", choices=["circle", "glue"], default="glue")
    ap.add_argument("--out", default="runs/hn_search.json")
    a = ap.parse_args()
    run(a.trials, a.budget, a.seed, Path(a.out), a.mode)
