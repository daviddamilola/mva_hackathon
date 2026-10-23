Hi, I’m David Oluwasusi, and I’m thankful to the organisers of the hackathon, and also thankful to the parents for releasing the data to work with.

I’ll be making a submission for Track 2, which is the drug repurposing track. The child has mosaic variegated aneuploidy, which was found to be caused by variants in the BUB1B gene. But the honest starting point is that the BUB1B gene cannot be drugged directly. It’s an essential gene with no approved drugs and no clinical candidates in any modality.

So I asked myself the question: is there anything near it that could be drugged instead? And I checked what those drugs actually do. Of the 70 mechanism records against the mitotic machinery, 69 are inhibitors, and they are developed to reduce mitotic fidelity in tumour cells. Since this disease is a loss of mitotic fidelity, that pharmacology would make it worse. So the strategy for me was not to fix the broken protein; it was to compensate for it.

So what I have done is to build a Nextflow pipeline that realigns the genome against the reference and annotates every variant in a 32-gene panel. The genes reported in this syndrome and the checkpoint machinery around them are what make up this 32-gene panel. I’ve also looked at the population frequency, clinical significance and predicted effect.

What it finds is that the causal variants are two rare variants in the BUB1B gene, and both hit the same part of the protein. There is one that is a nonsense allele, which truncates before the kinase domain, and there is a second, which is a missense allele that sits inside that domain. So between these two alleles, there is very little working BUBR1 left, and that is the checkpoint failure the disease is built on.

The strongest axis is readthrough of the nonsense allele. It’s strongest because it is the only one computed from this child’s genome rather than reasoned from the syndrome. The pipeline doesn’t assume that the nonsense allele is druggable. It checks, and the stop codon here is UGA, which is also the most permissive of the three. The three proteins a readthrough agent would produce were each scored, and all three came back as tolerated.

Two findings argue against this approach. The first is that the base immediately after the stop is the least permissive in the published hierarchy. So the transcript is destroyed by nonsense-mediated decay before it is translated. So readthrough must be paired with NMD inhibition. But there are no direct or specific NMD inhibitors, although some drugs perform NMD inhibition as a secondary mechanism.

The pairing of NMD inhibitors and readthrough agents is a potential answer to this problem. Gene-specific antisense inhibition silences decay of one transcript instead of the whole pathway, which matters enormously in a cancer predisposition syndrome. It isn’t approved at the moment, so I’ve not proposed it as a candidate, but I’ve reported this as a direction the mechanism points to.

So there are 11 approved medicines, which are grouped by how each is supported for this child, and every one is tiered with the evidence that would raise it. Thank you.