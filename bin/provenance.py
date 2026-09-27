#!/usr/bin/env python3
"""Provenance record: tool versions (from every task's versions.yml), reference build, database
versions, pipeline git commit and dirty state, Nextflow run metadata, and timestamp.

A tool reported at more than one version across tasks lists all of them (e.g. samtools in the
lane-3 exception documented in docs/01), so nothing is silently collapsed.

Usage: provenance.py --versions software_versions.yml --meta run_meta.json [--pharmcat-report report.json]
                     [--ot-summary pgx_ot.summary.json] --out provenance.json
"""
import argparse
import json
import os
import re
from collections import defaultdict
from datetime import datetime, timezone

ap = argparse.ArgumentParser()
ap.add_argument("--versions", required=True)
ap.add_argument("--meta", required=True)
ap.add_argument("--pharmcat-report")
ap.add_argument("--ot-summary")
ap.add_argument("--out", required=True)
a = ap.parse_args()

tools, per_process, proc = defaultdict(set), defaultdict(dict), None
for line in open(a.versions):
    if not line.strip():
        continue
    m = re.match(r'^\s*"?([^":]+(?::[^":]+)*)"?:\s*$', line)
    if m and not line.startswith(" " * 8) and line.rstrip().endswith(":"):
        proc = m.group(1).strip()
        continue
    kv = re.match(r"^\s+([A-Za-z0-9_.\-]+):\s*(.*)$", line)
    if kv and proc:
        k, v = kv.group(1), kv.group(2).strip().strip('"')
        per_process[proc][k] = v
        if not k.endswith(("_args", "_source", "_md5", "_sha256", "_b64")) and v and v != "stub":
            tools[k].add(v)

meta = json.load(open(a.meta))
dbs = dict(meta.get("databases", {}))
if a.pharmcat_report and os.path.exists(a.pharmcat_report):
    r = json.load(open(a.pharmcat_report))
    dbs["pharmcat_data_version (CPIC/DPWG/FDA via PharmGKB)"] = r.get("dataVersion")
if a.ot_summary and os.path.exists(a.ot_summary):
    ot = json.load(open(a.ot_summary)).get("open_targets") or {}
    dbs["open_targets"] = f"release {ot.get('data_release')} (API {ot.get('api_version')}, queried {ot.get('queried_utc')})"
else:
    dbs["open_targets"] = "stage 1b not run in this output directory"

record = {
    "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    "pipeline": meta.get("pipeline", {}),
    "reference": meta.get("reference", {}),
    "databases": dbs,
    "params": meta.get("params", {}),
    "tools": {k: " / ".join(sorted(v)) for k, v in sorted(tools.items())},
    "per_process": per_process,
}
json.dump(record, open(a.out, "w"), indent=2)
print(json.dumps({"tools": len(record["tools"]), "multi_version_tools": [k for k, v in tools.items() if len(v) > 1]}))
