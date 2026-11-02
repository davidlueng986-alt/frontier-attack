"""Parse Jaan Parts' 509-vertex .vtx file (Mathematica coordinate syntax) and
build the unit-distance graph by exact/high-precision distance tests."""
import re
import sys
import json
import mpmath as mp

mp.mp.dps = 80

PATH = sys.argv[1] if len(sys.argv) > 1 else "research/scratch/509_parts.vtx"


def parse_coord(s: str) -> mp.mpf:
    """Evaluate a Mathematica-style coordinate expression at 80 digits.

    Handles: integer/rational literals, Sqrt[...] (possibly nested),
    products with *, sums, unary minus and division.
    We translate to Python and evaluate in a namespace where sqrt is mp.sqrt.
    """
    t = s.strip()
    t = t.replace("Sqrt[", "sqrt(").replace("]", ")")
    # Mathematica implicit multiplication is not used in this file; all `*`
    # are explicit.  Guard anyway: number followed by '(' -> multiply.
    t = re.sub(r"(?<=[0-9)])\s*(?=sqrt\()", "*", t)
    t = t.replace("^", "**")
    val = eval(t, {"__builtins__": {}}, {"sqrt": mp.sqrt})
    return mp.mpf(val)


def main():
    lines = [ln.strip() for ln in open(PATH, encoding="utf-8") if ln.strip()]
    pts = []
    for i, ln in enumerate(lines, 1):
        assert ln[0] == "{" and ln[-1] == "}", (i, ln)
        body = ln[1:-1]
        # split top-level ", " (no commas nested inside in this file)
        parts = body.split(", ")
        assert len(parts) == 2, (i, ln)
        pts.append((parse_coord(parts[0]), parse_coord(parts[1])))
    n = len(pts)

    # exact squared distance, 80-digit arithmetic
    d2 = [[None] * n for _ in range(n)]
    edges = []
    max_dev_edge = mp.mpf(0)
    min_nonedge_gap = mp.inf
    for i in range(n):
        xi, yi = pts[i]
        for j in range(i + 1, n):
            xj, yj = pts[j]
            v = (xi - xj) ** 2 + (yi - yj) ** 2
            dev = abs(v - 1)
            if dev < mp.mpf("1e-30"):
                edges.append((i, j))
                if dev > max_dev_edge:
                    max_dev_edge = dev
            elif dev < min_nonedge_gap:
                min_nonedge_gap = dev

    deg = [0] * n
    for u, v in edges:
        deg[u] += 1
        deg[v] += 1

    # connectivity
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    seen = [False] * n
    comps = []
    for s in range(n):
        if not seen[s]:
            stack, comp = [s], []
            seen[s] = True
            while stack:
                a = stack.pop()
                comp.append(a)
                for b in adj[a]:
                    if not seen[b]:
                        seen[b] = True
                        stack.append(b)
            comps.append(comp)

    out = {
        "n": n,
        "m": len(edges),
        "max_dev_edge": float(max_dev_edge),
        "min_nonedge_gap": float(min_nonedge_gap),
        "deg_min": min(deg),
        "deg_max": max(deg),
        "deg_hist": {str(k): deg.count(k) for k in sorted(set(deg))},
        "components": sorted(len(c) for c in comps),
    }
    print(json.dumps(out, indent=2))
    with open("research/scratch/edges.txt", "w") as f:
        for u, v in edges:
            f.write(f"{u} {v}\n")
    with open("research/scratch/points.txt", "w") as f:
        for x, y in pts:
            f.write(f"{mp.nstr(x, 30)} {mp.nstr(y, 30)}\n")
    # sanity: first few coordinates
    print("first 3 points:", [(mp.nstr(x, 17), mp.nstr(y, 17)) for x, y in pts[:3]])


if __name__ == "__main__":
    main()
