"""Replay LMFDB rank-4 curves through the BSD gatekeeper.

Queries the LMFDB elliptic-curve API for curves of rank 4 and feeds each
published pair (rank, analytic_rank) to ``bsd.check_claim``.  A mismatch would
not prove anything by itself -- the point is to confirm that the published
data is internally coherent and that every discriminant is nonzero.

The API ignores ``_limit`` and can return a full page of records, so the
payload is read with a hard byte cap and only the first ``limit`` records are
used.  If the request fails, a small hardcoded sample is replayed instead and
the output is marked ``"offline": true``.

Run `python problems/bsd/lmfdb_scan.py --limit 50`.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from bsd import check_claim, invariants  # noqa: E402

API = "https://www.lmfdb.org/api/ec_curvedata/"
# `discriminant` is not an ec_curvedata column (absD/signD are); the primary
# URL below is kept verbatim as specified, with a corrected retry.
FIELDS = "lmfdb_label,ainvs,rank,analytic_rank,conductor,discriminant"
FIELDS_FALLBACK = "lmfdb_label,ainvs,rank,analytic_rank,conductor,absD,signD"
MAX_BYTES = 32 * 1024 * 1024
TIMEOUT = 30

# 37a1 plus four rank-4 curves (234446.a1, 501029.a1, 545723.a1, 602461.a1)
# taken from an ec_curvedata response.  Used only when the network is down.
OFFLINE_CURVES = [
    {"lmfdb_label": "37a1", "ainvs": [0, 0, 1, -1, 0], "rank": 1, "analytic_rank": 1},
    {"lmfdb_label": "234446.a1", "ainvs": [1, -1, 0, -79, 289], "rank": 4, "analytic_rank": 4},
    {"lmfdb_label": "501029.a1", "ainvs": [0, 1, 1, -72, 210], "rank": 4, "analytic_rank": 4},
    {"lmfdb_label": "545723.a1", "ainvs": [0, 0, 1, -7, 36], "rank": 4, "analytic_rank": 4},
    {"lmfdb_label": "602461.a1", "ainvs": [1, 0, 0, -202, 1089], "rank": 4, "analytic_rank": 4},
]


def _api_url(fields: str) -> str:
    return f"{API}?_format=json&rank=4&_fields={fields}&_limit=50"


def _read_json(url: str, timeout: float) -> list:
    """Fetch and parse one API reply, rejecting HTML/captcha/oversized bodies."""
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "frontier-attack-bsd-scan/1.0 (research script)",
            "Accept": "application/json",
        },
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        ctype = resp.headers.get("Content-Type", "")
        raw = resp.read(MAX_BYTES + 1)
    if len(raw) > MAX_BYTES:
        raise ValueError("response larger than byte cap")
    text = raw.decode("utf-8", "replace")
    if "json" not in ctype.lower() and not text.lstrip().startswith("{"):
        raise ValueError(f"non-JSON response ({ctype or 'no content type'})")
    payload = json.loads(text)
    data = payload.get("data") if isinstance(payload, dict) else None
    if not isinstance(data, list):
        raise ValueError("JSON response has no data list")
    return data


def _curl_json(url: str, timeout: float) -> list | None:
    """Secondary transport: some environments get a bot challenge over urllib."""
    exe = shutil.which("curl")
    if not exe:
        return None
    try:
        proc = subprocess.run(
            [exe, "-sS", "-L", "--max-time", str(int(timeout)),
             "-H", "Accept: application/json", url],
            capture_output=True,
            timeout=timeout + 5,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if proc.returncode != 0 or not proc.stdout:
        return None
    try:
        payload = json.loads(proc.stdout.decode("utf-8", "replace"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        return None
    data = payload.get("data") if isinstance(payload, dict) else None
    return data if isinstance(data, list) else None


def fetch(limit: int) -> tuple[list, bool]:
    """Return (records capped at limit, offline flag)."""
    started = time.monotonic()

    def remaining(cap: float) -> float:
        return max(1.0, min(cap, 80.0 - (time.monotonic() - started)))

    for fields, budget in ((FIELDS, TIMEOUT), (FIELDS_FALLBACK, 20.0)):
        url = _api_url(fields)
        try:
            return _read_json(url, remaining(budget))[:limit], False
        except urllib.error.HTTPError as exc:
            exc.close()
            if exc.code == 404:
                continue  # unknown field, retry with a supported field list
            break
        except Exception:
            break

    data = _curl_json(_api_url(FIELDS_FALLBACK), remaining(20.0))
    if data is not None:
        return data[:limit], False

    return [dict(rec) for rec in OFFLINE_CURVES[:limit]], True


def _claim_record(rec: dict) -> dict | None:
    ainvs = rec.get("ainvs")
    rank = rec.get("rank")
    analytic_rank = rec.get("analytic_rank")
    if not isinstance(ainvs, list) or len(ainvs) != 5:
        return None
    if not isinstance(rank, int) or not isinstance(analytic_rank, int):
        return None
    try:
        curve = dict(zip(("a1", "a2", "a3", "a4", "a6"), (int(v) for v in ainvs)))
    except (TypeError, ValueError):
        return None
    return {
        "label": rec.get("lmfdb_label"),
        "curve": curve,
        "rank": rank,
        "analytic_rank": analytic_rank,
    }


def scan(limit: int) -> dict:
    records, offline = fetch(limit)

    verdicts: dict[str, int] = {}
    mismatches: list[dict] = []
    nonzero = 0
    fetched = 0

    for raw in records:
        rec = _claim_record(raw)
        if rec is None:
            continue
        fetched += 1

        if invariants(rec["curve"])["discriminant"] != 0:
            nonzero += 1

        claim = {
            "curve": rec["curve"],
            "rank": rec["rank"],
            "analytic_rank": rec["analytic_rank"],
            "rank_method": "published",
            "analytic_rank_method": "published",
        }
        verdict = check_claim(claim)["verdict"]
        verdicts[verdict] = verdicts.get(verdict, 0) + 1

        if rec["rank"] != rec["analytic_rank"]:
            mismatches.append({
                "lmfdb_label": rec["label"],
                "ainvs": [rec["curve"][k] for k in ("a1", "a2", "a3", "a4", "a6")],
                "rank": rec["rank"],
                "analytic_rank": rec["analytic_rank"],
            })

    return {
        "fetched": fetched,
        "verdicts": verdicts,
        "mismatches": mismatches,
        "offline": offline,
        "nonzero_discriminants": nonzero,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--limit", type=int, default=50,
                        help="maximum records to scan (default: 50)")
    parser.add_argument("--out", default="runs/bsd_lmfdb.json",
                        help="output JSON path (default: runs/bsd_lmfdb.json)")
    args = parser.parse_args(argv)

    limit = max(0, args.limit)
    result = scan(limit)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
