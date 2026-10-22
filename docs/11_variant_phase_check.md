# 11 — Variant phase: can the compound-heterozygous configuration be verified?

Status: **checked** (2026-10-17). The report treats the two causal variants as compound
heterozygous (in trans). This document records whether that configuration is verifiable from the
data, what was checked, and what the inference rests on. It records no genotypes beyond what the
check itself required.

## The question

Both causal variants are called heterozygous in a proband-only sample. Short-read WGS can phase two
heterozygous variants only by chaining fragments across intermediate heterozygous sites. If the two
variants were on the same allele (cis), the other allele would be intact, the child would be a
carrier, and the report's whole mechanism section — one null allele, one missense allele — would be
wrong. The report asserted trans without saying it was an assumption.

## What was checked

The interval between the two variants (~11 kb) was queried for heterozygous sites in both callsets:
this pipeline's panel VCF and the distributed VCF. Both agree: the interval contains **one**
heterozygous site, with gaps of roughly 6.8 kb and 4.1 kb on either side. No short-read fragment
spans either gap, so read-backed phasing (WhatsHap/HapCUT2) cannot connect the pair regardless of
depth. Population-panel phasing is uninformative, because the missense allele's gnomAD frequency
(~6×10⁻⁷) means no reference panel carries its haplotype.

**Conclusion: the configuration is not verifiable from this data.** Parental genotyping at the two
positions (two Sanger assays), long reads or linked reads would verify it; none is in the released
data.

## What supports the inference

1. **Recessive disease logic.** MVA1 requires biallelic BUB1B dysfunction; a cis child is a
   monoallelic carrier, and carriers — including the parents in the original families — are
   unaffected (Hanks et al., *Nat Genet* 2004, PMID 15475955).
2. **The search for a cis-compatible second hit found nothing.** If cis, the phenotype needs
   another hit on the other allele. The pipeline looked: no BUB1B structural variant (Delly), no
   exon copy-number loss, no other rare functional BUB1B variant in the panel. A deep-intronic or
   regulatory hit would be invisible to this search.
3. **The canonical MVA1 genotype.** One truncating plus one missense allele, in trans, is the
   published architecture of the syndrome (Suijkerbuijk et al., *Cancer Res* 2010;70(12):4891–4900,
   PMID 20516114). A ClinVar/literature search found no published record of this exact variant
   pair, so there is no external segregation evidence for it.

## Where it landed in the report

- Section 3.1 now states the trans assumption in one sentence and points to the Limitations.
- Limitations has a "Variant phase." paragraph with the check's result and the three supports,
  citing Hanks 2004 and Suijkerbuijk 2010 (both in `assets/report/literature.yml`).
