#!/usr/bin/env python3
"""Synthetic tests for reference-allele naming across callers (no subject data).

PharmVar renamed CYP2C19*1 to *38. Aldy ships the new name, PharmCAT reports the old one, so the
same genotype would be scored as a caller disagreement unless the two names are equated. Since the
report presents caller agreement as a quality signal, a false disagreement in the most frequently
called gene would discredit the whole table.

The equivalence is deliberately gene-scoped and narrow, so these tests pin both directions: that
CYP2C19 *1 and *38 ARE equated, and that lookalike cases in other genes are NOT.

Both bin scripts carry their own copy of the table (they are standalone CLI tools staged
individually by Nextflow), so one test asserts the copies are identical.

Run:  python3 tests/unit/test_star_allele_normalisation.py
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RECONCILE = os.path.join(ROOT, "bin", "pgx_reconcile.py")
COMPARE = os.path.join(ROOT, "bin", "compare_star_alleles.py")


def load_region(path, start_marker, end_marker, extra=""):
    """Exec the pure-helper region of a script that otherwise runs as a CLI at import time."""
    src = open(path).read()
    body = src[src.index(start_marker):src.index(end_marker)]
    ns = {"re": re}
    exec(compile(extra + body, os.path.basename(path), "exec"), ns)  # noqa: S102 - test harness
    return ns


rec = load_region(RECONCILE, "# Reference-allele synonyms", "rows = []")
cmp_src = open(COMPARE).read()
cmp_ns = load_region(
    COMPARE, "# Reference-allele synonyms", "ap = argparse.ArgumentParser",
    extra=cmp_src[cmp_src.index("def core("):cmp_src.index("# Reference-allele synonyms")],
)

norm_pharmcat, norm_aldy = rec["norm_pharmcat"], rec["norm_aldy"]
semantic_pharmcat, semantic_aldy = rec["semantic_pharmcat"], rec["semantic_aldy"]
normalise = cmp_ns["normalise"]


def test_cyp2c19_reference_allele_is_not_a_disagreement():
    """PharmCAT *1/*2 and Aldy *38/*2 are the same genotype under two naming conventions."""
    assert norm_pharmcat("*1/*2", "CYP2C19") == norm_aldy("*38/*2", "CYP2C19")


def test_cyp2c19_both_alleles_reference():
    assert norm_pharmcat("*1/*1", "CYP2C19") == norm_aldy("*38/*38", "CYP2C19")


def test_cyp2c19_without_a_reference_allele_still_agrees():
    """The control: this project's own sample, where neither allele is reference."""
    assert norm_pharmcat("*2/*17", "CYP2C19") == norm_aldy("*2/*17", "CYP2C19")


def test_a_real_difference_is_still_a_difference():
    """The mapping must not make genuinely different CYP2C19 genotypes compare equal."""
    assert norm_pharmcat("*1/*2", "CYP2C19") != norm_aldy("*38/*17", "CYP2C19")


def test_equivalence_does_not_leak_to_other_genes():
    """*38 means something else outside CYP2C19; the table is keyed by gene for that reason."""
    assert norm_pharmcat("*1/*2", "CYP2D6") != norm_aldy("*38/*2", "CYP2D6")


def test_nat2_reference_defaults_are_not_equated():
    """Aldy's no-mutation NAT2 allele is *1 while PyPGx defaults to *4/*12. Those are distinct
    haplotypes, not synonyms, so equating them would manufacture agreement."""
    assert norm_pharmcat("*1/*1", "NAT2") != norm_aldy("*4/*4", "NAT2")
    assert norm_pharmcat("*1/*1", "NAT2") != norm_aldy("*12/*12", "NAT2")


def test_semantic_path_also_canonicalises():
    """A gene routed through the non-star semantic path must get the same treatment."""
    assert semantic_pharmcat("*1/*2", "CYP2C19") == semantic_aldy("*38/*2", "CYP2C19")


def test_compare_script_agrees_with_reconcile():
    """The Aldy-vs-PyPGx comparison uses a second normaliser; it must behave the same."""
    assert normalise("*38/*2", "CYP2C19") == normalise("*1/*2", "CYP2C19")
    assert normalise("*38/*2", "CYP2D6") != normalise("*1/*2", "CYP2D6")


def test_suballele_is_reduced_before_mapping():
    """Aldy emits suballeles such as *38.001; the core allele is what gets mapped."""
    assert norm_aldy("*38.001/*2", "CYP2C19") == norm_pharmcat("*1/*2", "CYP2C19")


def test_both_scripts_carry_the_same_table():
    """The table is duplicated across two standalone scripts; it must not drift."""
    assert rec["REF_ALLELE_SYNONYMS"] == cmp_ns["REF_ALLELE_SYNONYMS"], (
        rec["REF_ALLELE_SYNONYMS"], cmp_ns["REF_ALLELE_SYNONYMS"])


def test_unknown_gene_is_a_no_op():
    assert norm_aldy("*2/*17", "") == norm_aldy("*2/*17", "SOMETHING_ELSE")


if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn(); print(f"PASS {name}")
            except AssertionError as e:
                fails += 1; print(f"FAIL {name}: {e}")
    print(f"\n{fails} failure(s)")
    sys.exit(1 if fails else 0)
