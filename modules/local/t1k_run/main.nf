process T1K_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/t1k:1.0.10--h5814d7d_0' :
        'quay.io/biocontainers/t1k:1.0.10--h5814d7d_0' }"

    input:
    tuple val(meta), path(r1), path(r2)
    path t1k_ref

    output:
    tuple val(meta), path("${meta.id}.t1k_genotype.tsv"), emit: genotype
    path "${meta.id}.t1k_allele.tsv"                    , emit: alleles
    path "versions.yml"                                 , emit: versions

    script:
    """
    mkdir -p t1k_out   # run-t1k segfaults in fastq-extractor if the output directory is missing
    run-t1k -1 ${r1} -2 ${r2} --preset hla-wgs -f ${t1k_ref}/hla_dna_seq.fa -t ${task.cpus} -o t1k_out/${meta.id}
    cp t1k_out/${meta.id}_genotype.tsv ${meta.id}.t1k_genotype.tsv
    cp t1k_out/${meta.id}_allele.tsv ${meta.id}.t1k_allele.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: "1.0.10"
        t1k_reference: "\$(head -1 ${t1k_ref}/SOURCE.txt)"
        t1k_args: "--preset hla-wgs"
    END_VERSIONS
    """

    stub:
    """
    printf 'HLA-A\\t2\\tHLA-A*01:01:01\\t10\\t20\\tHLA-A*02:01:01\\t10\\t20\\n' > ${meta.id}.t1k_genotype.tsv
    touch ${meta.id}.t1k_allele.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: "1.0.10"
    END_VERSIONS
    """
}
