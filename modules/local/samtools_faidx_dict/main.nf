process SAMTOOLS_FAIDX_DICT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'quay.io/biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta}.fai")         , emit: fai
    tuple val(meta), path("${fasta.baseName}.dict"), emit: dict
    path "versions.yml"                            , emit: versions

    script:
    """
    samtools faidx ${fasta}
    samtools dict -a ${params.genome} -s 'Homo sapiens' -o ${fasta.baseName}.dict ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    printf '1\\t1000\\t3\\t60\\t61\\n' > ${fasta}.fai
    touch ${fasta.baseName}.dict

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.22.1
    END_VERSIONS
    """
}
