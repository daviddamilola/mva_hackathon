#!/usr/bin/env python3
"""Verify a file's checksum. Usage: verify_checksum.py <md5|sha256> <expected> <file> [record.tsv]"""
import hashlib
import sys

algo, expected, path = sys.argv[1:4]
h = hashlib.new(algo)
with open(path, "rb") as f:
    for chunk in iter(lambda: f.read(1 << 22), b""):
        h.update(chunk)
got = h.hexdigest()
if len(sys.argv) > 4:
    with open(sys.argv[4], "a") as out:
        out.write(f"{path}\t{algo}\t{got}\t{'OK' if got == expected else 'MISMATCH'}\n")
if got != expected:
    sys.exit(f"{algo} mismatch for {path}: got {got}, expected {expected}")
print(f"{path}: {algo} OK")
