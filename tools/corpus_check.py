#!/usr/bin/env python3
"""Validate research/corpus/index.json.

Exit 0 only if every check passes. Prints counts per problem and per role.
"""
import json
import os
import sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INDEX = os.path.join(ROOT, "research", "corpus", "index.json")

FIELDS = ["id", "title", "authors", "year", "url", "kind", "problem", "role", "why", "facts", "license", "verified"]
PROBLEMS = {"hn", "jac", "bsd", "methods"}
ROLES = {"core", "method", "analogue-solved", "analogue-open"}
KINDS = {"paper", "repo", "dataset", "note"}


def fail(msg):
    print("FAIL: " + msg)
    return 1


def main():
    if not os.path.exists(INDEX):
        return fail("missing " + INDEX)

    try:
        data = json.load(open(INDEX, encoding="utf-8"))
    except Exception as e:
        return fail("index.json is not valid JSON: %s" % e)

    if not isinstance(data, list):
        return fail("index.json must be a JSON array")

    errors = []
    if len(data) < 110:
        errors.append("only %d items; need at least 110" % len(data))

    seen = set()
    for n, item in enumerate(data):
        where = "item %d (%r)" % (n, item.get("id", "?"))
        if not isinstance(item, dict):
            errors.append("%s is not an object" % where)
            continue
        missing = [f for f in FIELDS if f not in item]
        if missing:
            errors.append("%s missing fields: %s" % (where, ", ".join(missing)))
            continue

        iid = item["id"]
        if not isinstance(iid, str) or not iid.strip():
            errors.append("%s has an empty id" % where)
        elif iid in seen:
            errors.append("duplicate id %r" % iid)
        else:
            seen.add(iid)

        if not isinstance(item["title"], str) or not item["title"].strip():
            errors.append("%s has an empty title" % where)
        if not isinstance(item["authors"], list) or not item["authors"] or not all(
                isinstance(a, str) and a.strip() for a in item["authors"]):
            errors.append("%s has a bad authors list" % where)
        if not isinstance(item["year"], int) or not (1000 <= item["year"] <= 2100):
            errors.append("%s has a bad year" % where)
        if not isinstance(item["url"], str) or not item["url"].startswith("http"):
            errors.append("%s has a non-URL url" % where)
        if item["kind"] not in KINDS:
            errors.append("%s has kind %r" % (where, item["kind"]))
        pl = item["problem"]
        if not isinstance(pl, list) or not pl or any(p not in PROBLEMS for p in pl):
            errors.append("%s has bad problem %r" % (where, pl))
        if item["role"] not in ROLES:
            errors.append("%s has bad role %r" % (where, item["role"]))
        if not isinstance(item["why"], str) or not item["why"].strip():
            errors.append("%s has an empty why" % where)
        if not isinstance(item["facts"], list) or not all(isinstance(f, str) and f.strip() for f in item["facts"]):
            errors.append("%s has a bad facts list" % where)
        if not item["facts"]:
            errors.append("%s has no facts and no NO FACT AVAILABLE marker" % where)
        if not isinstance(item["verified"], bool):
            errors.append("%s has a non-boolean verified" % where)
        elif item["verified"] and not item["url"].startswith("https://"):
            errors.append("%s is verified but its url is not https" % where)

    if errors:
        for e in errors:
            print("FAIL: " + e)
        print("\n%d error(s); %d item(s) present" % (len(errors), len(data)))
        return 1

    print("OK: %d items" % len(data))
    print("\nby problem:")
    pc = Counter(p for item in data for p in item["problem"])
    for k in sorted(pc):
        print("  %-8s %d" % (k, pc[k]))
    print("\nby role:")
    rc = Counter(item["role"] for item in data)
    for k in sorted(rc):
        print("  %-16s %d" % (k, rc[k]))
    print("\nunverified: %d" % sum(1 for item in data if not item["verified"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
