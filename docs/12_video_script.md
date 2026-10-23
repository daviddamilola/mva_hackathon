# 12 — Three-minute video script

Spoken script for the Track 2 submission video. Square brackets are delivery notes, not spoken.
Timing note at the end.

---

**[0:00 — thanks, ~30s]**

Thank you to Sage Bionetworks, the MVA Society, Hugging Face and BEACON for running this hackathon,
and to AWS and Anthropic for supporting it.

Most of all, thank you to the family. They shared their child's genome and their story so that other
families might get answers sooner. Everything I am about to show rests on that decision.

I'm David Oluwasusi. This is my Track 2 submission.

**[0:30 — the problem, ~30s]**

This child has mosaic variegated aneuploidy, caused by variants in BUB1B. The honest starting point
is that BUB1B cannot be drugged — an essential gene, with no approved drug and no clinical candidate
in any modality.

So I asked whether anything near it could be drugged instead, and checked what those drugs actually
do. Of seventy mechanism records against the mitotic machinery, sixty-nine are inhibitors —
developed to *reduce* mitotic fidelity in tumour cells.

This disease *is* a loss of mitotic fidelity, so that pharmacology would make it worse. The strategy
cannot be to fix the broken protein. It has to be compensation.

**[1:00 — method, ~40s]**

I built a Nextflow pipeline. It realigns the genome against a reference then annotates every variant in a thirty-two-gene panel — the genes reported in
this syndrome, and the checkpoint machinery around them — for population frequency, clinical
significance and predicted effect.

It finds the causal variants: two rare variants in BUB1B, and both hit the same part of the
protein. The nonsense allele truncates it before the kinase domain. The missense sits inside that
domain. Between them, very little working BubR1 is left — and that is the checkpoint failure the
disease is built on.

And it encodes one rule as code: any drug that inhibits a mitotic-fidelity gene is excluded, and the
pipeline reports what it rejected.

**[1:40 — the key result, ~65s]**

The strongest axis is readthrough of the nonsense allele. Strongest, because it is the only one
computed from this child's own genome rather than reasoned from the syndrome.

The pipeline doesn't assume a nonsense allele is druggable. It checks. The stop codon is UGA, the
most permissive of the three. The three proteins a readthrough agent would produce were each scored,
and all three come back tolerated.

But two findings argue against it, and I report both. The base immediately after the stop is the
least permissive in the published hierarchy. And the transcript is destroyed by nonsense-mediated
decay before it is translated — so readthrough must be paired with NMD inhibition, and no drug has
ever been developed for that.

That pairing is where the real answer sits. Gene-specific antisense inhibition silences decay of one
transcript instead of the whole pathway — which matters enormously in a cancer-predisposition
syndrome. It isn't approved, so I didn't propose it as a candidate. I reported it as the direction
the mechanism points.

**[2:45 — close, ~15s]**

Eleven approved medicines, grouped by how well each is supported *for this child*, every one tiered
with the evidence that would raise it.

Nothing here is a treatment recommendation. These are research hypotheses, stated with their
evidence and their assumptions, for a clinical team to take further.

Thank you.

---

## Notes

- **Length.** 502 spoken words: about 200 seconds at a normal 150 words-per-minute pace,
  about 215 seconds if you take it slowly. Time one read-through before recording. If you are
  over, cut the "And it encodes one rule as code…" sentence rather than rushing the readthrough
  section — that section carries the argument.
- **The one line to land.** "It doesn't assume a nonsense allele is druggable. It checks." That is
  the difference between this submission and a literature search, and it is worth a beat either side.
- **Reporting the two findings against the axis is deliberate** — say it as a strength, not an
  apology. A judge who knows the field will notice that the weaknesses are reported at all.
- **If you have time in hand**, add this after the variant-finding sentence: "It calls the
  pharmacogenes three separate ways and measures the downstream aneuploidy with a stated detection
  limit." It costs about seven seconds and shows the breadth of the pipeline.
- **If you have slides**, the natural four are: the tractability table (nothing to target), the
  domain map with both variants, the readthrough product scores, and the candidate axes grouped by
  support level.
