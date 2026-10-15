/*
 * Runs of homozygosity (bcftools roh, HMM on genotype likelihoods -G30, default AF 0.4, autosomes)
 * plus a model-free cross-check: heterozygous-SNV density in 1 Mb windows. F_ROH = total ROH length
 * (>= 1.5 Mb) / autosomal length. Informs consanguinity, i.e. the prior between homozygous and
 * compound-heterozygous recessive models.
 */
process BCFTOOLS_ROH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.22--h3a4d415_1' :
        'quay.io/biocontainers/bcftools:1.22--h3a4d415_1' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(ref_meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.roh.summary.json"), emit: summary
    path "${meta.id}.roh.txt"                           , emit: roh
    path "versions.yml"                                 , emit: versions

    script:
    """
    bcftools roh -G30 --AF-dflt 0.4 -O r -o ${meta.id}.roh.txt --threads ${task.cpus} \\
        -r \$(seq -s, 1 22) ${vcf}
    roh_summary.py --roh ${meta.id}.roh.txt --vcf ${vcf} --fai ${fai} --out ${meta.id}.roh.summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | cut -d' ' -f2)
        roh_args: "roh -G30 --AF-dflt 0.4 (autosomes)"
    END_VERSIONS
    """

    stub:
    """
    echo '{"stub": true}' > ${meta.id}.roh.summary.json; touch ${meta.id}.roh.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.22
    END_VERSIONS
    """
}
