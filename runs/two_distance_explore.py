"""Scratch: parse De Neve et al. Table 1 and inspect the 33-point set."""
import itertools, re
import sympy as sp

TXT = r"runs/papers/deNeve2025.txt"

lines = open(TXT, encoding="utf-8").read().splitlines()
# Table 1 rows are four-int bracket rows between the Table 1 caption and Table 2 caption
i1 = next(i for i, l in enumerate(lines) if l.startswith("Table 1."))
i2 = next(i for i, l in enumerate(lines) if l.startswith("Table 2."))
rows = []
for l in lines[i1:i2]:
    m = re.fullmatch(r"\[(-?\d+) (-?\d+) (-?\d+) (-?\d+)\]", l.strip())
    if m:
        rows.append(tuple(int(g) for g in m.groups()))
print("rows:", len(rows))

s5, s3, s15 = sp.sqrt(5), sp.sqrt(3), sp.sqrt(15)
pts = []
for a0, a1, a2, a3 in rows:
    x = (a0 + a1 * s5) / 8
    y = (a2 * s3 + a3 * s15) / 8
    pts.append((sp.nsimplify(x), sp.nsimplify(y)))

print("distinct points:", len({sp.simplify(p[0]).__str__() + "," + sp.simplify(p[1]).__str__() for p in pts}),
      "of", len(pts))

vals = {}
for i, j in itertools.combinations(range(len(pts)), 2):
    d2 = sp.simplify(sum((pts[i][k] - pts[j][k]) ** 2 for k in range(2)))
    vals.setdefault(sp.nsimplify(d2), []).append((i, j))

print("\ndistinct squared distances:", len(vals))
for d2, prs in sorted(vals.items(), key=lambda kv: float(kv[0])):
    print(f"  d2={d2}  d={sp.sqrt(d2)}  count={len(prs)}")

unit = [p for p in vals.get(sp.Integer(1), [])]
d3 = [p for p in vals.get(sp.Integer(3), [])]
print("\nunit pairs:", len(unit), " sqrt3 pairs:", len(d3), " total edges:", len(unit) + len(d3))
