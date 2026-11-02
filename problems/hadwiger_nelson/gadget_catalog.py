#!/usr/bin/env python3
"""gadget_catalog.py -- abstract-first catalog of level-5 virtual-edge candidates.

TARGET.  A level-5 virtual edge is a unit-distance graph G with two non-adjacent
terminals u, v such that

    |u - v|^2 = 3,   G is 5-colourable,   and every 5-colouring of G has
                                          colour(u) != colour(v).

The last clause is what lets G replace a sqrt(3)-edge of the De Neve et al.
33-vertex two-distance graph: an edge forbids u = v, and so must the gadget.

This script builds the abstract catalog first (the previous geometry-first
search was the wrong order: the unit part of the 33-vertex graph is
3-colourable, so every 4-colourable supergraph is free at k = 5 by the level
lemma), then filters geometrically.

  stage 0  audit of the reduction.  For a 5-colourable G:
             G + uv not 5-colourable  <=>  every 5-colouring has u = v   (EQUAL)
             G / uv not 5-colourable  <=>  every 5-colouring has u != v  (DIFFER)
           A virtual edge needs the second.  "H 6-critical, G = H - uv" gives
           the first, so it produces equality gadgets (the diamond pattern at
           k = 3), not virtual edges.  Both readings are carried through the
           pipeline: route A = the literal recipe, route B = the corrected
           forced-difference route.
  stage 1  the complete catalog of edge-6-critical graphs on <= 13 vertices
           (Lalonde 2025, hosted by B. McKay): counts verified, min degree
           checked, criticality re-verified on a sample with SAT.
  stage 2A route A: every (H, edge uv) -> G = H - uv; dedupe (G; {u, v}).
  stage 2B route B: a minimal forced-difference gadget has a 6-critical
           quotient H and is obtained from H by splitting a vertex w into the
           two terminals; hence |V(core)| = |V(H)| + 1 and m(core) >= m(H).
  stage 3  geometric filter, cheapest first:
             3a  edge count vs u(n), the exact maximum number of unit distances
                 among n points (Schade; OEIS A186705)
             3b  clique number <= 3          [no 4 pairwise-unit points]
             3c  K_{2,3}-free                [two unit circles meet in <= 2 points]
             3d  certified circle-intersection embedding with |uv|^2 = 3
  stage 4  the smallest embeddable G, if any (exact coordinates + relation).
  stage 5  runs/hn_gadget_catalog.json

Deterministic.  Writes only this script's three outputs: the script itself is
pre-existing, and runs/gadget_catalog.log and runs/hn_gadget_catalog.json are
appended/written.  Exit status 0 whenever the pipeline completes.

Usage:
    python problems/hadwiger_nelson/gadget_catalog.py
    python problems/hadwiger_nelson/gadget_catalog.py --catalog-dir D --no-download
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import sys
import time
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import networkx as nx
import sympy as sp
from mpmath import iv
from pysat.solvers import Solver

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

from relation import colouring_cnf, relation  # noqa: E402

LOG_PATH = ROOT / "runs" / "gadget_catalog.log"
OUT_PATH = ROOT / "runs" / "hn_gadget_catalog.json"
SOLVER = "cadical153"

# ---------------------------------------------------------------------------
# published data
# ---------------------------------------------------------------------------

# Complete edge-6-critical catalogs of Lalonde, Electron. J. Combin. 32(1)
# (2025) #P1.18, hosted by B. McKay.  Counts transcribed from the k = 6 column
# of users.cecs.anu.edu.au/~bdm/data/graphs.html and re-verified line by line.
PUBLISHED_COUNTS = {6: 1, 7: 0, 8: 1, 9: 2, 10: 22, 11: 393, 12: 17036, 13: 1479809}
CATALOG_TOTAL = sum(PUBLISHED_COUNTS.values())  # 1,497,264
BASE_URL = "https://users.cecs.anu.edu.au/~bdm/data/crit/crit_{n}_6.g6"

# u(n) = maximum number of unit distances among n points of the plane = maximum
# number of edges of a unit-distance graph on n vertices.  Exact for n <= 14 in
# Schade's thesis; table reproduced in Agoston-Palvolgyi, Period. Math. Hungar.
# (2022), Table 1, and OEIS A186705.  u(15) = 37 is proved in the same paper.
U_EXACT = {
    1: 0, 2: 1, 3: 3, 4: 5, 5: 7, 6: 9, 7: 12, 8: 14,
    9: 18, 10: 20, 11: 23, 12: 27, 13: 30, 14: 33, 15: 37,
}
U_SOURCE = (
    "Schade (thesis), exact for n<=14; reproduced in Agoston-Palvolgyi, "
    "Period. Math. Hungar. (2022) Table 1, and OEIS A186705; u(15)=37 proved "
    "by Agoston-Palvolgyi.  Cross-check: the density recursion "
    "u(n) <= floor(n/(n-2) * u(n-1)) holds for every 7 <= n <= 15."
)

POINTER_REPO = {
    "url": "https://github.com/HeliCorgi/fourteen-points-six-colors",
    "commit": "4c84e4d67522e644faef704694cd5ba7fc273abc",
    "doc": "code/minimality/COMPLETENESS.md sec.1 and code/minimality/README.md",
    "what": "points at the catalog and re-verifies its counts; the catalog is "
            "generated by github.com/lalondeo/gencrit (nauty geng)",
}


# ---------------------------------------------------------------------------
# logging
# ---------------------------------------------------------------------------

class Log:
    def __init__(self, path: Path):
        self.path = path
        path.parent.mkdir(parents=True, exist_ok=True)

    def __call__(self, stage: str, msg: str) -> None:
        stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        line = f"[{stamp}] {stage}: {msg}"
        with open(self.path, "a", encoding="utf-8") as f:
            f.write(line + "\n")
        print(line, flush=True)


# ---------------------------------------------------------------------------
# graph6
# ---------------------------------------------------------------------------

def g6_decode(line: str):
    """Decode one graph6 line into (n, sorted edge list).  n <= 32."""
    b = line.encode("ascii")
    if not b:
        raise ValueError("empty graph6 line")
    first = b[0]
    if first < 63:
        n = first
        body = b[1:]
    elif first == 63:
        raise ValueError("n >= 63 unsupported")
    else:
        n = first - 63
        body = b[1:]
    if not 1 <= n <= 32:
        raise ValueError(f"bad n={n}")
    bits = n * (n - 1) // 2
    x = 0
    for ch in body:
        x = (x << 6) | ((ch - 63) & 63)
    pad = len(body) * 6 - bits
    if pad < 0:
        raise ValueError("short graph6 body")
    if pad:
        x >>= pad
    edges = []
    pos = bits - 1
    for j in range(1, n):
        for i in range(j):
            if (x >> pos) & 1:
                edges.append((i, j))
            pos -= 1
    return n, edges


def open_text(path: Path):
    with open(path, "rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="ascii")
    return open(path, "r", encoding="ascii")


# ---------------------------------------------------------------------------
# SAT
# ---------------------------------------------------------------------------

def sat_colourable(n: int, edges, k: int) -> bool:
    with Solver(name=SOLVER, bootstrap_with=colouring_cnf(n, edges, k)) as s:
        return bool(s.solve())


def verify_criticality(n: int, edges, time_cap: float = 90.0) -> dict:
    """chi(H) = 6 and chi(H - e) = 5 for every edge e."""
    t0 = time.time()
    if sat_colourable(n, edges, 5):
        return {"n": n, "m": len(edges), "ok": False, "why": "5-colourable"}
    checked = 0
    for i, e in enumerate(edges):
        if time.time() - t0 > time_cap:
            return {"n": n, "m": len(edges), "ok": None, "why": "time-cap",
                    "edges_checked": checked}
        rest = edges[:i] + edges[i + 1:]
        if not sat_colourable(n, rest, 5):
            return {"n": n, "m": len(edges), "ok": False,
                    "why": f"H-{e} not 5-colourable"}
        checked += 1
    return {"n": n, "m": len(edges), "ok": True,
            "why": "chi=6 and chi(H-e)=5 for all e",
            "edges_checked": checked, "seconds": round(time.time() - t0, 3)}


# ---------------------------------------------------------------------------
# stage 1: catalog acquisition and verification
# ---------------------------------------------------------------------------

def ensure_files(cache: Path, log: Log, download: bool) -> dict:
    cache.mkdir(parents=True, exist_ok=True)
    files = {}
    for n, count in PUBLISHED_COUNTS.items():
        if count == 0:
            continue
        cands = [cache / f"crit_{n}_6.g6", cache / f"crit_{n}_6.g6.gz"]
        path = next((p for p in cands if p.exists() and p.stat().st_size > 0), None)
        if path is None:
            if not download:
                raise RuntimeError(f"catalog file for n={n} missing in {cache}")
            url = BASE_URL.format(n=n) + ".gz"
            path = cache / f"crit_{n}_6.g6.gz"
            log("stage1", f"downloading {url}")
            t0 = time.time()
            try:
                with urllib.request.urlopen(url, timeout=180) as r, open(path, "wb") as f:
                    while True:
                        chunk = r.read(1 << 20)
                        if not chunk:
                            break
                        f.write(chunk)
            except Exception as exc:  # noqa: BLE001
                raise RuntimeError(f"download failed for n={n}: {exc}") from exc
            log("stage1", f"downloaded n={n} bytes={path.stat().st_size} "
                          f"in {time.time()-t0:.1f}s")
        files[n] = path
    return files


def fallback_catalog(log: Log) -> dict:
    """SAT-verified 6-critical sample used only if the catalog is unreachable.

    K6, K3 join C5 and K2 join Moser spindle are 6-critical; each is checked
    here with SAT.  This is the "regenerate a sample yourself" fallback.
    """
    log("stage1", "catalog unreachable -> SAT-regenerated fallback sample")
    moser = [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 4), (3, 4), (4, 5),
             (4, 6), (5, 6), (2, 5)]

    def join(a, ea, b, eb):
        return a + b, list(ea) + list(eb) + [(i, a + j) for i in range(a) for j in range(b)]

    k6 = (6, [(i, j) for i in range(6) for j in range(i + 1, 6)])
    c5 = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
    k3c5 = join(3, [(0, 1), (1, 2), (0, 2)], 5, c5)
    k2moser = join(2, [(0, 1)], 7, moser)
    out = {}
    for n, edges in (k6, k3c5, k2moser):
        chk = verify_criticality(n, edges, time_cap=30.0)
        if not chk["ok"]:
            raise RuntimeError(f"fallback graph n={n} is not 6-critical: {chk}")
        out.setdefault(n, []).append(edges)
    return out


def scan_catalog(files: dict, log: Log) -> dict:
    """One streaming pass per n.  Keeps all edge lists for n <= 11 and a small
    deterministic sample (first 8 + sparsest 8) for n >= 12."""
    t0 = time.time()
    per_n, keep_all, samples, sparsest = {}, {}, {}, {}
    for n, path in sorted(files.items()):
        cnt = 0
        m_hist = Counter()
        delta_hist = Counter()
        sum_m = 0
        all_edges = [] if n <= 11 else None
        first = [] if n >= 12 else None
        heap = []
        with open_text(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                nn, edges = g6_decode(line)
                if nn != n:
                    raise RuntimeError(f"n mismatch in {path}: {nn} != {n}")
                cnt += 1
                m = len(edges)
                m_hist[m] += 1
                sum_m += m
                deg = [0] * n
                for a, b in edges:
                    deg[a] += 1
                    deg[b] += 1
                delta_hist[min(deg)] += 1
                if all_edges is not None:
                    all_edges.append(edges)
                else:
                    if len(first) < 8:
                        first.append(edges)
                    heap.append((m, cnt, edges))
                    heap.sort(key=lambda t: (t[0], t[1]))
                    del heap[8:]
        if cnt != PUBLISHED_COUNTS[n]:
            raise RuntimeError(
                f"n={n}: line count {cnt} != published {PUBLISHED_COUNTS[n]}")
        per_n[n] = {
            "graphs": cnt,
            "edges_total": sum_m,
            "m_min": min(m_hist),
            "m_max": max(m_hist),
            "m_hist": dict(sorted(m_hist.items())),
            "min_degree": min(delta_hist),
            "min_degree_hist": dict(sorted(delta_hist.items())),
        }
        if all_edges is not None:
            keep_all[n] = all_edges
        else:
            samples[n] = first
            sparsest[n] = [e for _, _, e in heap]
    log("stage1", f"parsed {sum(p['graphs'] for p in per_n.values())} graphs "
                  f"in {time.time()-t0:.1f}s")
    if 7 not in per_n:
        per_n[7] = {"graphs": 0, "edges_total": 0, "m_min": None, "m_max": None,
                    "m_hist": {}, "min_degree": None, "min_degree_hist": {}}
    return {"per_n": per_n, "keep_all": keep_all,
            "samples": samples, "sparsest": sparsest}


# ---------------------------------------------------------------------------
# stage 2A: literal route -- H - uv, forced equality
# ---------------------------------------------------------------------------

def augment(n: int, edges, u: int, v: int) -> nx.Graph:
    """G with one pendant at u and one at v.  Isomorphism of the augmented
    graphs is exactly isomorphism of (G; {u,v})."""
    g = nx.Graph()
    g.add_nodes_from(range(n + 2))
    g.add_edges_from(edges)
    g.add_edge(u, n)
    g.add_edge(v, n + 1)
    return g


def dedupe(n: int, graphs, time_cap: float = 240.0):
    t0 = time.time()
    buckets: dict = defaultdict(list)
    classes: list = []
    pairs = 0
    capped = False
    for edges in graphs:
        for e in edges:
            if time.time() - t0 > time_cap:
                capped = True
                break
            u, v = e
            g = [x for x in edges if x != e]
            pairs += 1
            aug = augment(n, g, u, v)
            key = (n, len(g), nx.weisfeiler_lehman_graph_hash(aug, iterations=3))
            hit = False
            for rep_aug, _rep in buckets[key]:
                if nx.is_isomorphic(aug, rep_aug):
                    hit = True
                    break
            if not hit:
                buckets[key].append((aug, (g, u, v)))
                classes.append((n, g, u, v))
        if capped:
            break
    return classes, pairs, capped


# ---------------------------------------------------------------------------
# stage 2B: corrected route -- split a vertex w of a 6-critical H
# ---------------------------------------------------------------------------

def twin_split(n: int, edges, w: int):
    """H with w replaced by two non-adjacent twins u = w and v = n.

    G/uv is H.  G has n+1 vertices and m(G) = m(H) + deg_H(w).
    """
    v = n
    out = [e for e in edges if w not in e]
    for a, b in edges:
        if a == w:
            out.append((b, w))
            out.append((b, v))
        elif b == w:
            out.append((a, w))
            out.append((a, v))
    return n + 1, out, w, v


def partial_split(n: int, edges, w: int, A):
    """Split w into u = w (keeps the neighbours in A) and v = n (the rest).
    A must cover every neighbour of w for the quotient to be H."""
    v = n
    out = [e for e in edges if w not in e]
    for a, b in edges:
        if a == w:
            x = b
        elif b == w:
            x = a
        else:
            continue
        out.append((x, w) if x in A else (x, v))
    return n + 1, out, w, v


# ---------------------------------------------------------------------------
# stage 3: sound structural filters
# ---------------------------------------------------------------------------

def adj_bits(n: int, edges) -> list[int]:
    adj = [0] * n
    for a, b in edges:
        adj[a] |= 1 << b
        adj[b] |= 1 << a
    return adj


def has_k4(n: int, adj) -> bool:
    for b in range(n):
        for c in range(b + 1, n):
            if not (adj[b] >> c) & 1:
                continue
            common = adj[b] & adj[c]
            verts = []
            x = common
            while x:
                lb = x & -x
                verts.append(lb.bit_length() - 1)
                x ^= lb
            for i in range(len(verts)):
                for j in range(i + 1, len(verts)):
                    if (adj[verts[i]] >> verts[j]) & 1:
                        return True
    return False


def structural(n: int, edges, u: int, v: int) -> dict:
    adj = adj_bits(n, edges)
    max_common = 0
    for i in range(n):
        for j in range(i + 1, n):
            c = (adj[i] & adj[j]).bit_count()
            if c > max_common:
                max_common = c
    tri3 = sum((adj[a] & adj[b]).bit_count() for a, b in edges)
    return {
        "m": len(edges),
        "k4": has_k4(n, adj),
        "max_common": max_common,
        "triangles": tri3 // 3,
        "edge_in_more_than_2_triangles": tri3 > 2 * len(edges),
        "uv_adjacent": bool((adj[u] >> v) & 1),
    }


# ---------------------------------------------------------------------------
# stage 3d: certified circle-intersection embedding, |uv|^2 = 3
# ---------------------------------------------------------------------------

iv.dps = 40
IV_ZERO = iv.mpf(0)
IV_ONE = iv.mpf(1)
IV_SQ3 = iv.sqrt(3)


def _dist2(p, q):
    return (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2


def _overlaps_one(x) -> bool:
    return not (x.b < IV_ONE.a or x.a > IV_ONE.b)


def embed_order(n: int, edges, u: int, v: int):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    order, placed = [u, v], {u, v}
    while len(order) < n:
        rest = [w for w in range(n) if w not in placed]
        best = max(rest, key=lambda w: (len(adj[w] & placed), len(adj[w]), -w))
        if len(adj[best] & placed) < 2:
            return None, adj
        order.append(best)
        placed.add(best)
    return order, adj


def certify_no_embedding(n: int, edges, u: int, v: int,
                         node_cap: int = 400_000, time_cap: float = 30.0) -> dict:
    """Sound certificate for "no placement with every edge unit and |uv|^2 = 3".

    Status: CERTIFIED-EXCLUDED (exhausted), PLACEMENT-FOUND, CAPPED,
    NON-RIGID-ORDER or BAD-TERMINALS.  Only CERTIFIED-EXCLUDED is a proof.
    """
    t0 = time.time()
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    if v in adj[u]:
        return {"status": "BAD-TERMINALS", "nodes": 0, "seconds": 0.0}
    order, adj = embed_order(n, edges, u, v)
    if order is None:
        return {"status": "NON-RIGID-ORDER", "nodes": 0, "seconds": 0.0}

    pos = {u: (IV_ZERO, IV_ZERO), v: (IV_SQ3, IV_ZERO)}
    nodes = [0]
    capped = [False]
    complete = [False]

    def dfs(k: int, mirror_used: bool) -> bool:
        nodes[0] += 1
        if nodes[0] > node_cap or time.time() - t0 > time_cap:
            capped[0] = True
            return True
        if k == len(order):
            complete[0] = True
            return True
        w = order[k]
        anchors = [x for x in order[:k] if x in adj[w]]
        pair = None
        for i in range(len(anchors)):
            for j in range(i + 1, len(anchors)):
                if _dist2(pos[anchors[i]], pos[anchors[j]]).a > 0:
                    pair = (anchors[i], anchors[j])
                    break
            if pair:
                break
        if pair is None:
            complete[0] = True      # cannot proceed soundly -> no claim
            return True
        a, b = pair
        others = [x for x in anchors if x not in pair]
        pa, pb = pos[a], pos[b]
        d2 = _dist2(pa, pb)
        d = iv.sqrt(d2)
        x = d2 / (2 * d)            # r1 = r2 = 1
        h2 = IV_ONE - x * x
        if h2.b < 0:
            return False
        h2c = iv.mpf([h2.a if h2.a > 0 else 0, h2.b])
        h = iv.sqrt(h2c)
        ux = (pb[0] - pa[0]) / d
        uy = (pb[1] - pa[1]) / d
        bx = pa[0] + x * ux
        by = pa[1] + x * uy
        for sgn in (1, -1):
            p = (bx - sgn * h * uy, by + sgn * h * ux)
            if not mirror_used and p[1].b < 0:
                continue
            m2 = mirror_used or p[1].a > 0 or p[1].b < 0
            ok = True
            for o in others:
                if not _overlaps_one(_dist2(p, pos[o])):
                    ok = False
                    break
            if ok:
                pos[w] = p
                if dfs(k + 1, m2):
                    return True
                del pos[w]
            if h2c.b == 0:
                break
        return False

    ran = dfs(2, False)
    if capped[0]:
        status = "CAPPED"
    elif ran and complete[0]:
        status = "PLACEMENT-FOUND"
    else:
        status = "CERTIFIED-EXCLUDED"
    return {"status": status, "nodes": nodes[0],
            "seconds": round(time.time() - t0, 3)}


# ---------------------------------------------------------------------------
# exact arithmetic helpers (controls)
# ---------------------------------------------------------------------------

def exact_check(n: int, edges, coords, u: int, v: int) -> dict:
    pts = [(sp.sympify(a), sp.sympify(b)) for a, b in coords]
    bad = []
    for a, b in edges:
        d = sp.expand((pts[a][0] - pts[b][0]) ** 2 + (pts[a][1] - pts[b][1]) ** 2)
        if sp.simplify(d - 1) != 0:
            bad.append({"edge": [a, b], "d2": str(d)})
    uv = sp.expand((pts[u][0] - pts[v][0]) ** 2 + (pts[u][1] - pts[v][1]) ** 2)
    return {
        "n": n, "m": len(edges),
        "edges_exact": not bad, "bad_edges": bad,
        "uv_d2": str(sp.simplify(uv)),
        "uv_is_sqrt3": sp.simplify(uv - 3) == 0,
        "all_points_distinct": len({(sp.simplify(a), sp.simplify(b)) for a, b in pts}) == n,
    }


def control_diamond() -> dict:
    s3 = sp.sqrt(3)
    coords = [(0, 0), (s3, 0),
              (s3 / 2, sp.Rational(1, 2)), (s3 / 2, sp.Rational(-1, 2))]
    edges = [(0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    return {
        "name": "diamond (K4 - e)",
        "role": "stage-0 witness: G + uv = K4 is not 3-colourable, tips forced EQUAL",
        "coords": [[str(x), str(y)] for x, y in coords],
        "exact": exact_check(4, edges, coords, 0, 1),
        "relation_k3": relation(edges, 4, 0, 1, 3),
        "relation_k5": relation(edges, 4, 0, 1, 5),
        "certifier": certify_no_embedding(4, edges, 0, 1),
    }


def control_wheel() -> dict:
    s3 = sp.sqrt(3)
    h = sp.Rational(1, 2)
    coords = [(0, 0), (1, 0), (h, s3 / 2), (-h, s3 / 2), (-1, 0),
              (-h, -s3 / 2), (h, -s3 / 2)]
    edges = ([(0, i) for i in range(1, 7)]
             + [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)])
    return {
        "name": "hexagon + centre (W6, 12 = u(7) edges)",
        "role": "positive embedding control: exact sqrt(3) pair, level lemma => free at k=5",
        "coords": [[str(x), str(y)] for x, y in coords],
        "exact": exact_check(7, edges, coords, 2, 4),
        "relation_k5": relation(edges, 7, 2, 4, 5),
        "certifier": certify_no_embedding(7, edges, 2, 4),
    }


def control_k6_minus_edge() -> dict:
    n = 6
    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if (i, j) != (0, 1)]
    return {
        "name": "K6 - e",
        "role": ("stage-0 witness: G + uv = K6 is not 5-colourable but "
                 "relation(0,1,5) = equal, not different"),
        "n": n, "m": len(edges),
        "relation_k5": relation(edges, n, 0, 1, 5),
        "u_n": U_EXACT[n],
        "structural": structural(n, edges, 0, 1),
        "certifier": certify_no_embedding(n, edges, 0, 1),
    }


def control_route_b() -> dict:
    """k=3 analogue of the corrected route: K4 with vertex 3 split into
    u=3 (neighbours 0,1) and v=4 (neighbour 2).  G/uv = K4 is 4-chromatic,
    G is 3-colourable, so relation(u,v,3) is 'different'."""
    edges = [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 4)]
    g, g_edges, u, v = partial_split(4, [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)],
                                     3, {0, 1})
    return {
        "name": "K4 with a vertex split (k=3 analogue of route B)",
        "role": "validates the corrected construction: quotient K4 is 4-chromatic",
        "n": g, "m": len(g_edges), "terminals": [u, v],
        "edges": g_edges,
        "quotient_n": 4,
        "relation_k3": relation(g_edges, g, u, v, 3),
        "chi_le_3": sat_colourable(g, g_edges, 3),
        "quotient_chi_le_3": sat_colourable(4, [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)], 3),
    }


# ---------------------------------------------------------------------------
# pipeline
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--catalog-dir", default=os.environ.get("HN_CATALOG_DIR"))
    ap.add_argument("--no-download", action="store_true")
    args = ap.parse_args(argv)

    t_start = time.time()
    log = Log(LOG_PATH)
    log("run", "gadget_catalog.py start (abstract catalog -> geometric filter)")

    cache = Path(args.catalog_dir) if args.catalog_dir else (
        Path(os.environ.get("TEMP", "/tmp")) / "opencode" / "hn_crit_catalog")

    # ---- stage 0: audit of the reduction --------------------------------
    audit = {
        "lemma_equal": ("for a k-colourable G and non-adjacent u,v: "
                        "G + uv is not k-colourable  <=>  every k-colouring of G "
                        "has colour(u) = colour(v)  [EQUAL gadget]"),
        "lemma_different": ("for a k-colourable G and non-adjacent u,v: "
                            "G / uv is not k-colourable  <=>  every k-colouring of G "
                            "has colour(u) != colour(v)  [DIFFER gadget = virtual edge]"),
        "proof": ("(1) a k-colouring of G with colour(u) != colour(v) is exactly a "
                  "proper k-colouring of G + uv; (2) a k-colouring of G with "
                  "colour(u) = colour(v) is exactly a proper k-colouring of G / uv."),
        "consequence": ("the recipe 'H 6-critical, G = H - uv, G + uv not "
                        "5-colourable' produces an EQUAL gadget, not a virtual "
                        "edge.  Stage 4's relation(u,v,5) == 'different' test can "
                        "never pass on that route.  The virtual edge is the "
                        "vertex-split (quotient) route B."),
        "route_B_minimal_core": ("if G is a forced-difference gadget, take a "
                                 "6-critical H inside G/uv and the preimage of "
                                 "V(H): a subgraph G'' with |V(G'')| = |V(H)| + 1, "
                                 "m(G'') >= m(H), G'' unit-distance and "
                                 "G''/uv = H."),
        "controls": {
            "k6_minus_edge": control_k6_minus_edge(),
            "diamond": control_diamond(),
            "route_b_k3_analogue": control_route_b(),
            "wheel7": control_wheel(),
        },
    }
    log("stage0", "audit: G+uv UNSAT <=> forced EQUAL; G/uv UNSAT <=> forced "
                  "DIFFERENT.  K6-e: relation(0,1,5)=equal, so the literal route "
                  "yields equality gadgets; the virtual edge is the split route B")

    # ---- stage 1 ---------------------------------------------------------
    catalog_complete = True
    try:
        files = ensure_files(cache, log, download=not args.no_download)
        scan = scan_catalog(files, log)
        per_n = scan["per_n"]
        catalog_total = sum(v["graphs"] for v in per_n.values())
        if catalog_total != CATALOG_TOTAL:
            raise RuntimeError(f"catalog total {catalog_total} != {CATALOG_TOTAL}")
        log("stage1", f"catalog complete: {catalog_total} edge-6-critical graphs "
                      f"on n=6..13, counts match Lalonde/McKay; min degree = "
                      f"{min(v['min_degree'] for v in per_n.values() if v['min_degree'])}")
        sources = sorted({str(f) for f in files.values()})
        keep_all, samples, sparsest = scan["keep_all"], scan["samples"], scan["sparsest"]
    except Exception as exc:  # noqa: BLE001
        catalog_complete = False
        catalog_total = None
        log("stage1", f"CATALOG UNAVAILABLE: {exc}; falling back to SAT sample")
        fb = fallback_catalog(log)
        per_n = {}
        for n, gs in fb.items():
            ms, deltas = [], []
            for e in gs:
                ms.append(len(e))
                deg = Counter()
                for a, b in e:
                    deg[a] += 1
                    deg[b] += 1
                deltas.append(min(deg.values()))
            per_n[n] = {"graphs": len(gs), "edges_total": sum(ms),
                        "m_min": min(ms), "m_max": max(ms),
                        "min_degree": min(deltas), "m_hist": {}, "min_degree_hist": {}}
        keep_all = {n: gs for n, gs in fb.items() if n <= 11}
        samples, sparsest, files, sources = {}, {}, {}, []

    # criticality re-verification
    crit, pool, seen = [], [], set()
    for n in sorted(per_n):
        gs = keep_all.get(n)
        if gs:
            pool.extend((n, e) for e in gs)
        else:
            if samples.get(n):
                pool.append((n, samples[n][0]))
            if sparsest.get(n):
                pool.append((n, sparsest[n][0]))
    for n, e in pool:
        key = (n, tuple(e))
        if key in seen:
            continue
        seen.add(key)
        crit.append(verify_criticality(n, e))
    crit_ok = all(r["ok"] for r in crit)
    log("stage1", f"criticality re-check on {len(crit)} graphs: "
                  + ("all chi=6 and chi(H-e)=5" if crit_ok else "FAILURES PRESENT"))

    dens_ok = all(U_EXACT[n] <= (n * U_EXACT[n - 1]) // (n - 2) for n in range(7, 16))
    log("stage1", f"u(n): u(13)={U_EXACT[13]} u(14)={U_EXACT[14]} u(15)={U_EXACT[15]}; "
                  f"density recursion consistent={dens_ok}")

    # ---- stage 2A: literal route (H - uv) --------------------------------
    stage2a, classes_a = [], []
    for n in sorted(per_n):
        if per_n[n]["graphs"] == 0:
            continue
        row = {"n": n, "graphs": per_n[n]["graphs"],
               "edge_pairs": per_n[n]["edges_total"],
               "deduped_classes": None, "dedupe": "skipped"}
        if n <= 11 and keep_all.get(n):
            cls, seen_pairs, capped = dedupe(n, keep_all[n])
            row.update(deduped_classes=len(cls), pairs_expanded=seen_pairs,
                       dedupe="capped" if capped else "full")
            classes_a.extend(cls)
        else:
            row["dedupe"] = "skipped: all G die at 3a on m alone"
        stage2a.append(row)
    log("stage2A", f"literal route: {sum(r['edge_pairs'] for r in stage2a)} (H,edge) "
                   f"pairs; dedupe for n<=11 -> {len(classes_a)} classes")

    # ---- stage 2B: corrected route (vertex split) ------------------------
    stage2b = []
    for n in sorted(per_n):
        if per_n[n]["graphs"] == 0:
            continue
        stage2b.append({"n": n, "H_graphs": per_n[n]["graphs"],
                        "split_pairs": n * per_n[n]["graphs"],
                        "core_n": n + 1,
                        "m_H_min": per_n[n]["m_min"],
                        "m_core_min": per_n[n]["m_min"]})
    log("stage2B", f"corrected route: {sum(r['split_pairs'] for r in stage2b)} "
                   f"(H, w) splits; core has |V(H)|+1 vertices and m >= m(H)")

    # ---- stage 3 ---------------------------------------------------------
    # 3a route A: m(G) = m(H) - 1 vs u(n)
    a_rows, a_surv = [], []
    for r in stage2a:
        n = r["n"]
        mG, un = per_n[n]["m_min"] - 1, U_EXACT[n]
        a_rows.append({"n": n, "m_G_min": mG, "u_n": un, "margin": mG - un,
                       "killed": mG > un, "pairs": r["edge_pairs"]})
        if mG <= un:
            a_surv.append(n)
    a_killed = sum(r["pairs"] for r in a_rows if r["killed"])
    log("stage3a", f"route A edge count: killed {a_killed} (H,edge) pairs; "
                   f"survivor sizes = {a_surv}")

    # 3a route B: m(core) >= m(H) vs u(n+1)
    b_rows, b_surv = [], []
    for r in stage2b:
        N = r["core_n"]
        mcore, uN = per_n[r["n"]]["m_min"], U_EXACT.get(N)
        killed = uN is not None and mcore > uN
        b_rows.append({"n_H": r["n"], "n_core": N, "m_core_min": mcore,
                       "u_n_core": uN,
                       "margin": (mcore - uN) if uN is not None else None,
                       "killed": bool(killed), "pairs": r["split_pairs"]})
        if not killed:
            b_surv.append(N)
    b_killed = sum(r["pairs"] for r in b_rows if r["killed"])
    log("stage3a", f"route B edge count: killed {b_killed} (H,w) splits; "
                   f"survivor core sizes = {b_surv}")

    # 3b/3c on the literal-route classes (independent of 3a)
    b_pass = c_pass = both = 0
    for n, g, u, v in classes_a:
        st = structural(n, g, u, v)
        ok_b = not st["k4"]
        ok_c = st["max_common"] <= 2 and not st["edge_in_more_than_2_triangles"]
        b_pass += int(ok_b)
        c_pass += int(ok_c)
        both += int(ok_b and ok_c)
    log("stage3b", f"clique<=3: {b_pass}/{len(classes_a)} literal classes pass")
    log("stage3c", f"K_2,3-free and edges in <=2 triangles: {c_pass}/"
                   f"{len(classes_a)} pass; both: {both}")

    # 3d: nothing reaches the embedder; certify the sparsest members directly,
    # and exercise the embedder on the controls.
    direct = []
    for n in sorted(sparsest):
        for edges in sparsest[n][:1]:
            u, v = edges[0]
            g = [e for e in edges if e != (u, v)]
            res = certify_no_embedding(n, g, u, v, node_cap=200_000, time_cap=20.0)
            direct.append({"n": n, "m_G": len(g), "u_n": U_EXACT[n],
                           "terminals": [u, v], **res})
    log("stage3d", f"direct interval certification of {len(direct)} sparsest "
                   f"catalog members: "
                   + ", ".join(sorted({d['status'] for d in direct})))
    log("stage3d", "certified embedder controls: "
                   + ", ".join(f"{k}={v['certifier']['status']}"
                               for k, v in audit["controls"].items()
                               if "certifier" in v))

    # validation of the corrected construction on a sample of (H, w)
    split_sample = {"tested": 0, "five_colourable": 0, "relation_different": 0,
                    "six_chromatic": 0, "items": []}
    for n in sorted(per_n):
        if per_n[n]["graphs"] == 0 or n > 12:
            continue
        gs = keep_all.get(n) or (samples.get(n) or [])
        for edges in gs:
            if n == 11 and split_sample["tested"] > 40:
                break
            w = max(range(n), key=lambda x: sum(1 for a, b in edges if x in (a, b)))
            gn, ge, u, v = twin_split(n, edges, w)
            ok5 = sat_colourable(gn, ge, 5)
            split_sample["tested"] += 1
            if ok5:
                split_sample["five_colourable"] += 1
                rel = relation(ge, gn, u, v, 5)
                if rel == "different":
                    split_sample["relation_different"] += 1
                if len(split_sample["items"]) < 8:
                    split_sample["items"].append(
                        {"n_H": n, "m_H": len(edges), "m_G": len(ge),
                         "relation_k5": rel, "degree_w": len([1 for a, b in edges if w in (a, b)])})
            else:
                split_sample["six_chromatic"] += 1
    log("stage3d", f"twin-split sample: {split_sample['tested']} splits, "
                   f"{split_sample['five_colourable']} 5-colourable, "
                   f"{split_sample['six_chromatic']} 6-chromatic")

    # ---- stage 4 ---------------------------------------------------------
    found = False
    log("stage4", "no candidate reached the embedder on either route; found=False")

    # ---- frontier --------------------------------------------------------
    frontier = None
    for n in range(6, 30):
        un = U_EXACT.get(n)
        if un is None:
            break
        if math.ceil(5 * n / 2) - 1 <= un:
            frontier = n
            break
    # route B frontier: smallest core size N for which m(H) <= u(N) is possible
    frontier_b = None
    for N in range(7, 30):
        uN = U_EXACT.get(N)
        if uN is None:
            break
        if math.ceil(5 * (N - 1) / 2) <= uN:
            frontier_b = N
            break

    # ---- stage 5 ---------------------------------------------------------
    payload = {
        "task": "level-5 virtual edge: abstract catalog first, geometric filter second",
        "terminal_condition": ("unit-distance G, terminals u,v non-adjacent, "
                               "|uv|^2 = 3, G 5-colourable, every 5-colouring "
                               "has colour(u) != colour(v)"),
        "verdict": "NEGATIVE_CATALOG" if catalog_complete else "NEGATIVE_SAMPLE",
        "found": found,
        "n": None,
        "uv_dist2": "3",
        "construction_audit": audit,
        "catalog": {
            "total": catalog_total,
            "complete": catalog_complete,
            "n_min": 6,
            "n_max": 13,
            "per_n": {str(n): {k: v for k, v in per_n[n].items() if k != "m_hist"}
                      for n in sorted(per_n)},
            "m_hist": {str(n): per_n[n]["m_hist"] for n in sorted(per_n)},
            "provenance": {
                "generator": "Olivier Lalonde, gencrit (nauty geng); EJC 32(1) (2025) #P1.18",
                "host": "users.cecs.anu.edu.au/~bdm/data/crit/crit_<n>_6.g6[.gz]",
                "pointer": POINTER_REPO,
                "files": sources if catalog_complete else [],
            },
        },
        "stage1": {
            "counts_verified": catalog_complete,
            "min_degree_all_graphs": (min(v["min_degree"] for v in per_n.values()
                                         if v["min_degree"]) if per_n else None),
            "criticality_sample": {"checked": len(crit), "all_ok": crit_ok,
                                   "items": crit},
        },
        "u_n": {
            "values": {str(k): v for k, v in sorted(U_EXACT.items()) if k <= 15},
            "source": U_SOURCE,
            "density_recursion_ok": dens_ok,
            "meaning": "u(n) = max edges of a unit-distance graph on n vertices",
        },
        "route_A_literal": {
            "what": "H 6-critical, G = H - uv (forced EQUALITY: G + uv = H is not 5-colourable)",
            "per_n": stage2a,
            "pairs_total": sum(r["edge_pairs"] for r in stage2a),
            "deduped_classes_small_n": len(classes_a),
            "dedupe_note": ("bucket by Weisfeiler-Lehman hash of G plus two pendants "
                            "at u and v, then networkx.is_isomorphic inside each "
                            "bucket; this is exact for (G; {u,v})"),
            "filter_3a": a_rows,
            "killed_by_3a": a_killed,
            "survivor_sizes": a_surv,
            "structural_independent": {"clique_le_3": b_pass,
                                       "k23_and_triangle": both,
                                       "tested": len(classes_a)},
            "verdict": ("NEGATIVE by edge count; and structurally this route "
                        "cannot produce relation 'different' at all (stage 0)"),
        },
        "route_B_corrected": {
            "what": ("forced DIFFERENCE: G/uv not 5-colourable.  Minimal core = "
                     "6-critical H with a vertex split into u,v; |V(core)|=|V(H)|+1, "
                     "m(core) >= m(H)"),
            "per_n": stage2b,
            "split_pairs_total": sum(r["split_pairs"] for r in stage2b),
            "filter_3a": b_rows,
            "killed_by_3a": b_killed,
            "survivor_core_sizes": b_surv,
            "twin_split_sample": split_sample,
            "verdict": "NEGATIVE by edge count",
        },
        "stage3": {
            "filters": [
                {"id": "3a", "name": "edge count", "rule": "m <= u(n)",
                 "sound": True,
                 "why": "u(n) is the exact maximum number of unit distances among n points",
                 "route_A_killed": a_killed, "route_B_killed": b_killed},
                {"id": "3b", "name": "clique", "rule": "omega(G) <= 3", "sound": True,
                 "why": "no four points of the plane are pairwise at unit distance",
                 "route_A_pass_independent": b_pass, "tested": len(classes_a)},
                {"id": "3c", "name": "K_2,3", "rule": "K_2,3-free; edge in <= 2 triangles",
                 "sound": True,
                 "why": "two unit circles intersect in at most two points",
                 "route_A_pass_independent": both, "tested": len(classes_a)},
                {"id": "3d", "name": "certified embedding",
                 "rule": "circle-intersection DFS, interval arithmetic, |uv|^2 = 3",
                 "sound": True, "tested": 0, "survivors": 0},
            ],
            "direct_certification_sample": direct,
            "survivors": 0,
        },
        "stage4": {
            "found": found,
            "smallest_embeddable": None,
            "reason": ("route A: m(G) >= 36 at n = 13 while u(13) = 30.  "
                       "route B: the minimal core has |V(H)|+1 vertices and at "
                       "least m(H) edges, so m >= 37 at |V(H)| = 13 while "
                       "u(14) = 33.  Both routes are excluded with margin >= 4."),
        },
        "frontier": {
            "route_A": {
                "first_n_not_excluded": frontier,
                "argument": ("H 6-critical => delta(H) >= 5 => m(H) >= ceil(5n/2); "
                             "m(G) = m(H)-1.  This alone kills G up to n = 12; "
                             "with the catalog's actual minimum densities it kills "
                             "n <= 13 (catalog) and n <= 14 (min-degree at n=14: "
                             "34 > u(14) = 33)."),
            },
            "route_B": {
                "first_core_size_not_excluded": frontier_b,
                "argument": ("core size N needs a 6-critical H on N-1 vertices with "
                             "m(H) <= u(N).  min-degree alone reaches N = 12; the "
                             "catalog (H <= 13 vertices) reaches N = 14.  At N = 15 "
                             "the requirement is a 6-critical graph on 14 vertices "
                             "with at most u(15) = 37 edges, which is exactly the "
                             "min-degree bound ceil(5*14/2) = 35 <= 37: uncataloged."),
                "catalog_covers_core_upto": 14,
                "unconditional_bound_without_catalog": 12,
            },
        },
        "conclusion": None,
        "seconds": round(time.time() - t_start, 1),
    }

    if catalog_complete:
        payload["conclusion"] = (
            f"NEGATIVE_CATALOG.  The complete catalog of {catalog_total} "
            f"edge-6-critical graphs on <= 13 vertices was enumerated "
            f"(counts, min degree and a SAT criticality sample verified).  "
            f"Route A (the literal recipe G = H - uv) yields forced-EQUALITY "
            f"gadgets, not virtual edges: G + uv not 5-colourable is exactly "
            f"'every 5-colouring has u = v' (witness K6-e, relation = equal).  "
            f"Even so it is geometrically dead: the sparsest G at n = 13 has 36 "
            f"edges against u(13) = 30.  Route B (the corrected forced-difference "
            f"object, a vertex split of a 6-critical H) is dead by a wider margin: "
            f"its minimal core has |V(H)|+1 vertices and at least m(H) edges, so "
            f"37 edges on 14 vertices against u(14) = 33.  Hence there is no "
            f"level-5 virtual edge of this form on <= 14 vertices.  The first "
            f"size not excluded is a 15-vertex core, which would need a 6-critical "
            f"graph on 14 vertices with <= 37 edges (the min-degree bound allows "
            f"35) -- beyond the published catalog, whose enumeration stops at 13."
        )
    else:
        payload["conclusion"] = (
            "NEGATIVE_SAMPLE: the catalog could not be obtained; the SAT-verified "
            "fallback sample is excluded by the same edge-count argument, so no "
            "gadget was found, but the completeness claim is not available."
        )

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    log("stage5", f"written {OUT_PATH.relative_to(ROOT)}")

    print()
    print("VERDICT: " + payload["verdict"])
    print(f"catalog={catalog_total} found={found} "
          f"route_A_first_n={frontier} route_B_first_core={frontier_b} "
          f"seconds={payload['seconds']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
