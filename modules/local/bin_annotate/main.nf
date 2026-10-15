/*
 * Fixed-width genome bins (primary chromosomes) and panel MANE exons, annotated with GC fraction
 * (bedtools nuc on the explicit reference), uniquely-mappable fraction (Umap k100, the longest k
 * published; conservative for 150 bp reads) and blacklist overlap (ENCODE v2).
 */
process BIN_ANNOTATE {
    tag "${params.bin_size} bp"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--h13024bc_3' :
        'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3' }"

    input:
    tuple val(ref_meta), path(fasta), path(fai)
    tuple val(umap_name), path(umap)
    tuple val(bl_name), path(blacklist)
    path spliceai_annotation      // MANE Select exons of the panel genes

    output:
    path "bins.bed"          , emit: bins_bed
    path "bins.annot.tsv"    , emit: bins
    path "exons.bed"         , emit: exons_bed
    path "exons.annot.tsv"   , emit: exons
    path "versions.yml"      , emit: versions

    script:
    """
    set -o pipefail
    awk -v OFS='\\t' '\$1 ~ /^([0-9]+|X|Y)\$/ {print \$1, \$2}' ${fai} > genome.txt
    bedtools makewindows -g genome.txt -w ${params.bin_size} > bins.bed
    cut -f1-3 ${umap} | sort -k1,1 -k2,2n > umap.bed
    cut -f1-3 ${blacklist} | sort -k1,1 -k2,2n > blacklist.bed

    bedtools nuc -fi ${fasta} -bed bins.bed | tail -n +2 | cut -f1-3,5 > gc.tsv
    bedtools coverage -a bins.bed -b umap.bed | awk -v OFS='\\t' '{print \$NF}' > map.tsv
    bedtools coverage -a bins.bed -b blacklist.bed | awk -v OFS='\\t' '{print \$(NF-2)}' > bl.tsv
    ( printf 'chrom\\tstart\\tend\\tgc\\tmap_frac\\tblacklist_bp\\n'; paste gc.tsv map.tsv bl.tsv ) > bins.annot.tsv

    # MANE exons: gene, exon number in transcript order
    awk -F'\\t' -v OFS='\\t' 'NR>1 {n=split(\$6,s,","); split(\$7,e,","); for(i=1;i<n;i++){ex=(\$3=="+")?i:n-i; print \$2, s[i], e[i], \$1"|exon"ex}}' ${spliceai_annotation} \\
        | sort -k1,1V -k2,2n > exons.bed
    bedtools nuc -fi ${fasta} -bed exons.bed | tail -n +2 | cut -f1-4,6 > exons.gc.tsv
    bedtools coverage -a exons.bed -b umap.bed | awk '{print \$NF}' > exons.map.tsv
    ( printf 'chrom\\tstart\\tend\\tname\\tgc\\tmap_frac\\n'; paste exons.gc.tsv exons.map.tsv ) > exons.annot.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
        bin_size: "${params.bin_size}"
    END_VERSIONS
    """

    stub:
    """
    printf '1\\t0\\t100000\\n' > bins.bed
    printf 'chrom\\tstart\\tend\\tgc\\tmap_frac\\tblacklist_bp\\n1\\t0\\t100000\\t0.4\\t1\\t0\\n' > bins.annot.tsv
    printf '1\\t100\\t200\\tG|exon1\\n' > exons.bed
    printf 'chrom\\tstart\\tend\\tname\\tgc\\tmap_frac\\n1\\t100\\t200\\tG|exon1\\t0.4\\t1\\n' > exons.annot.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.31.1
        bin_size: "${params.bin_size}"
    END_VERSIONS
    """
}
