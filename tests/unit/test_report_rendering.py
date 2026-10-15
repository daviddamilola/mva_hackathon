#!/usr/bin/env python3
"""Synthetic tests for the report's display helpers (no subject data).

Four display defects were found in a rendered report, and each is fixed in bin/render_report.py
rather than in the data, because the pipeline's own tables are the record and must stay raw:

  collapse_rows   The Open Targets join emits one row per (allele tested, evidence level,
                  annotation category). A drug whose guideline text is identical for every allele
                  therefore printed the same paragraph three or four times. Collapsing is display
                  only, so the test pins that nothing is dropped: every allele, level and matched
                  allele of the raw rows survives in the merged cell.
  keyed_notes     A long guideline paragraph moves out of the table cell to a keyed note quoted
                  once; a short basis stays in the cell. The test pins the threshold behaviour and
                  that a note records any OTHER drug its text names, which is how the simvastatin
                  row carrying DPWG's atorvastatin wording is caught instead of reading as a join
                  error.
  starname        Aldy truncates the reference-allele name to '*ference'. Only the printed name is
                  repaired; no other allele name may be rewritten.
  tex             Long single tokens (md5, accession, URL) have no break opportunity and stood past
                  the margin. The escaper inserts break points. The test pins that escaping is still
                  correct and that no break point lands inside a LaTeX control sequence.

Run:  python3 tests/unit/test_report_rendering.py
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RENDER = os.path.join(ROOT, "bin", "render_report.py")

src = open(RENDER).read()


def region(start, end, ns=None):
    """Exec one helper region of a script that otherwise parses argv at import time."""
    body = src[src.index(start):src.index(end)]
    ns = dict(ns or {}, re=re)
    exec(compile(body, "render_report.py", "exec"), ns)  # noqa: S102 - test harness
    return ns


rows_ns = region("# ---- stage 1b verdict rows, as printed", "ot_by_axis = []",
                 {"pc_report": {"drugs": {"DPWG Guideline Annotation": {"atorvastatin": {}, "simvastatin": {}}}}})
collapse_rows, keyed_notes, note_key = rows_ns["collapse_rows"], rows_ns["keyed_notes"], rows_ns["note_key"]
tex = region("def tex(s):", "def fmt(")["tex"]
starname = region("# Aldy truncates the reference-allele name", "def mdcell(")["starname"]

fails = []


def check(name, got, want):
    if got != want:
        fails.append(f"{name}: got {got!r}, want {want!r}")


LONG = ("[CPIC; Strong] Prescribe the desired starting dose and adjust doses based on "
        "disease-specific guidelines, which is a sentence long enough to belong in a note.")
SHORT = "no PGx rows in Open Targets"

# ---- collapse_rows: three alleles of one gene, one guideline paragraph -> one row, nothing lost
raw = [{"gene": "MT-RNR1", "caller": "MT_RNR1_CALL", "flag": "no flag", "guideline_text": LONG,
        "tier": "2", "patient_support": "mechanism_only", "called": c, "evidence_level": "1A",
        "matched_ot_allele": m}
       for c, m in (("T (m.1095T>C not detected)", "rs267606618 = T"),
                    ("C (m.1494C>T not detected)", "rs267606619 = C"),
                    ("A (m.1555A>G not detected)", "rs267606617 = A"))]
got = collapse_rows(raw)
check("three alleles collapse to one row", len(got), 1)
for r in raw:
    if r["called"] not in got[0]["called_merged"]:
        fails.append(f"collapse dropped a call: {r['called']}")
    if r["matched_ot_allele"] not in got[0]["matched_merged"]:
        fails.append(f"collapse dropped a matched allele: {r['matched_ot_allele']}")
check("one level, printed once", got[0]["levels_merged"], "1A")

# differing evidence levels merge into one cell rather than into separate rows
lv = collapse_rows([dict(raw[0], evidence_level=x) for x in ("1B", "1A", "1B")])
check("levels merge", (len(lv), lv[0]["levels_merged"]), (1, "1A, 1B"))

# a different caller is a different row: the two are not the same measurement
two = collapse_rows([dict(raw[0], caller="PharmCAT"), dict(raw[0], caller="force-call (CRAM)")])
check("different callers stay apart", len(two), 2)

# ---- keyed_notes: long text keyed and quoted once, short basis left in the cell
drugs = [{"name": "GENTAMICIN", "rows_collapsed": collapse_rows(raw)},
         {"name": "AMIKACIN", "rows_collapsed": collapse_rows([dict(raw[0], guideline_text=LONG)])},
         {"name": "ATALUREN", "rows_collapsed": collapse_rows([dict(raw[0], guideline_text=SHORT)])}]
notes = keyed_notes(drugs)
check("one note for one shared paragraph", len(notes), 1)
check("shared paragraph keyed on both drugs",
      [d["rows_collapsed"][0]["note"] for d in drugs[:2]], ["a", "a"])
check("short basis is not keyed", drugs[2]["rows_collapsed"][0]["note"], "")
check("short basis stays in the cell", drugs[2]["rows_collapsed"][0]["basis_short"], SHORT)
check("note lists the drugs it serves", notes[0]["for"], ["amikacin", "gentamicin"])

# a paragraph naming a drug other than its own row's is reported, not silently printed
simva = [{"name": "SIMVASTATIN", "rows_collapsed": collapse_rows(
    [dict(raw[0], guideline_text="[DPWG; No recommendation] The guideline does not provide a "
                                 "recommendation for atorvastatin in patients with normal function.")])}]
n = keyed_notes(simva)[0]
check("other drug named in the text is flagged", n["names_other_drug"], ["atorvastatin"])
check("its own drug is not flagged as foreign", "simvastatin" in n["names_other_drug"], False)
own = keyed_notes([{"name": "atorvastatin", "rows_collapsed": collapse_rows(
    [dict(raw[0], guideline_text="[DPWG] no recommendation for atorvastatin in normal function, a "
                                 "sentence long enough to be keyed as a note.")])}])[0]
check("a drug's own name is never foreign", own["names_other_drug"], [])

# statements concatenated by the join are printed apart
st = keyed_notes([{"name": "SIMVASTATIN", "rows_collapsed": collapse_rows(
    [dict(raw[0], guideline_text="[CPIC; Strong] First statement, long enough to be keyed. || "
                                 "[DPWG; No recommendation] Second statement.")])}])[0]
check("concatenated statements split", len(st["statements"]), 2)
check("no separator survives", any("||" in x for x in st["statements"]), False)

# keys keep going past 26 notes
check("note keys past z", [note_key(0), note_key(25), note_key(26)], ["a", "z", "aa"])

# ---- starname: repair the truncated reference name, touch nothing else
check("truncated reference repaired", starname("*6178GT/*ference"), "*6178GT/Reference")
check("real allele names untouched", starname("*4/*35x2"), "*4/*35x2")
check("already-correct name untouched", starname("Reference/Reference"), "Reference/Reference")
check("empty is passed through", starname(""), "")

# ---- tex: escaping stays correct, and break points never land inside a control sequence
out = tex("a_b & 50% $x #1 {y} ~ ^ < > |")
for frag in ("\\_", "\\&", "\\%", "\\$", "\\#", "\\{", "\\}",
             "\\textasciitilde{}", "\\textasciicircum{}", "\\textless{}", "\\textgreater{}", "\\textbar{}"):
    if frag not in out:
        fails.append(f"tex() lost the escape {frag}")
md5 = tex("a056c57649f3c9964c68aead3849bbf8")
check("a long unbreakable run gains break points", "\\allowbreak{}" in md5, True)
check("its characters are unchanged", re.sub(r"\\allowbreak\{\}", "", md5),
      "a056c57649f3c9964c68aead3849bbf8")
for m in re.finditer(r"\\[a-zA-Z]+", tex("GCA_000001405.15_GRCh38_no_alt/sub;x")):
    if m.group(0) not in ("\\allowbreak", "\\_"):
        fails.append(f"tex() produced an unexpected control sequence {m.group(0)}")
check("none is passed through as empty", tex(None), "")

if fails:
    print("FAIL")
    for f in fails:
        print("  -", f)
    sys.exit(1)
print("ok: report display helpers (collapse, keyed notes, allele names, LaTeX escaping)")
