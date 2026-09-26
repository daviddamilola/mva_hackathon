#!/usr/bin/env python3
"""Validate that a CRAM was encoded against the given reference.

Compares every @SQ line in the CRAM header (SN, LN, M5) with the reference sequence dictionary
(samtools dict output, which carries M5 = MD5 of the uppercase sequence). Any missing M5, length
difference or MD5 difference fails the run, because CRAM decoding against the wrong reference
silently corrupts bases.

Usage: validate_cram_m5.py <cram_header.sam> <reference.dict> <out.json>
"""
import json
import sys


def sq_records(path):
    out = {}
    order = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("@SQ"):
                continue
            tags = dict(f.split(":", 1) for f in line.rstrip("\n").split("\t")[1:] if ":" in f)
            out[tags["SN"]] = {"LN": int(tags["LN"]), "M5": tags.get("M5")}
            order.append(tags["SN"])
    return out, order


header_path, dict_path, out_path = sys.argv[1:4]
cram, cram_order = sq_records(header_path)
ref, ref_order = sq_records(dict_path)

missing_m5 = [sn for sn, r in cram.items() if not r["M5"]]
only_cram = sorted(set(cram) - set(ref))
only_ref = sorted(set(ref) - set(cram))
len_mismatch = [sn for sn in set(cram) & set(ref) if cram[sn]["LN"] != ref[sn]["LN"]]
md5_mismatch = [sn for sn in set(cram) & set(ref)
                if cram[sn]["M5"] and ref[sn]["M5"] and cram[sn]["M5"].lower() != ref[sn]["M5"].lower()]

result = {
    "n_sq_cram": len(cram),
    "n_sq_reference": len(ref),
    "n_m5_checked": sum(1 for sn in set(cram) & set(ref) if cram[sn]["M5"] and ref[sn]["M5"]),
    "same_order": cram_order == ref_order,
    "cram_sq_missing_m5": missing_m5[:50],
    "only_in_cram": only_cram[:50],
    "only_in_reference": only_ref[:50],
    "length_mismatch": sorted(len_mismatch)[:50],
    "md5_mismatch": sorted(md5_mismatch)[:50],
}
result["pass"] = not (missing_m5 or only_cram or only_ref or len_mismatch or md5_mismatch)

with open(out_path, "w") as out:
    json.dump(result, out, indent=2)
print(json.dumps({k: v for k, v in result.items() if not isinstance(v, list) or v}))
if not result["pass"]:
    sys.exit(f"CRAM/reference M5 validation failed; see {out_path}")
