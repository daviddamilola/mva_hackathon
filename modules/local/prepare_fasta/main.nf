process PREPARE_FASTA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--h13024bc_3' :
        'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3' }"

    input:
    val meta        // [ id, url, md5 ]
    path mask_bed

    output:
    tuple val(meta), path("${meta.id}.fa"), emit: fasta
    path "${meta.id}.masked_regions.bed"  , emit: mask
    path "versions.yml"                   , emit: versions

    script:
    """
    curl -fsSL --retry 5 --retry-delay 10 -o source.fna.gz '${meta.url}'
    verify_checksum.py md5 ${meta.md5} source.fna.gz

    # UCSC names -> no-chr names used by the subject VCF (chrM -> M, chrEBV -> EBV, chrUn_* -> Un_*)
    gzip -dc source.fna.gz | sed -E 's/^>chr/>/; s/^>([^ \\t]+).*/>\\1/' > unmasked.fa
    rm source.fna.gz

    grep -v '^#' ${mask_bed} > ${meta.id}.masked_regions.bed
    bedtools maskfasta -fi unmasked.fa -bed ${meta.id}.masked_regions.bed -fo ${meta.id}.fa
    rm -f unmasked.fa unmasked.fa.fai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
        curl: \$(curl --version | head -1 | cut -d' ' -f2)
        reference_source: "${meta.url}"
        reference_md5: "${meta.md5}"
    END_VERSIONS
    """

    stub:
    """
    printf '>1\\nACGTACGTNN\\n' > ${meta.id}.fa
    grep -v '^#' ${mask_bed} > ${meta.id}.masked_regions.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.31.1
        curl: 8.14.1
        reference_source: "${meta.url}"
        reference_md5: "${meta.md5}"
    END_VERSIONS
    """
}
