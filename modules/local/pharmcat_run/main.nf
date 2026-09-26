/*
 * PharmCAT on the force-called VCF: rename to chr, preprocess against our reference, then match,
 * phenotype and report with Aldy outside calls (e.g. CYP2D6).
 * Missing-to-ref options are deliberately NOT used: every position is genotyped from the CRAM.
 */
process PHARMCAT_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/pharmcat3:3.4.0--py313h106432d_0' :
        'quay.io/biocontainers/pharmcat3:3.4.0--py313h106432d_0' }"

    input:
    tuple val(meta), path(vcf), path(csi), path(outside_calls)
    path chr_map
    tuple path(ref_chr), path(ref_chr_fai), path(ref_chr_gzi)

    output:
    tuple val(meta), path("${meta.id}.phenotype.json"), emit: phenotype
    tuple val(meta), path("${meta.id}.report.json")   , emit: report
    path "${meta.id}.match.json"                      , emit: match
    path "${meta.id}.report.html"                     , emit: html
    path "${meta.id}.preprocessed.vcf.bgz"            , emit: preprocessed
    path "versions.yml"                               , emit: versions

    script:
    def outside = "-po ${outside_calls}"
    """
    set -o pipefail
    awk '{print \$2, \$1}' ${chr_map} > nochr_to_chr.txt
    CONTIGS=\$(cut -f1 ${ref_chr_fai} | paste -sd, -)
    bcftools annotate --rename-chrs nochr_to_chr.txt ${vcf} -Ou | bcftools view -t \$CONTIGS -Oz -o input.chr.vcf.gz
    bcftools index input.chr.vcf.gz

    pharmcat_vcf_preprocessor -vcf input.chr.vcf.gz -refFna ${ref_chr} -o pre -bf ${meta.id}
    mv pre/${meta.id}.preprocessed.vcf.bgz .

    # an empty outside-call file is valid input for PharmCAT
    pharmcat -vcf ${meta.id}.preprocessed.vcf.bgz ${outside} -o . -bf ${meta.id} \\
        -reporterJson -reporterHtml -rs ${params.pharmcat_reporter_sources}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: \$(pharmcat -version 2>&1 | head -1 | sed 's/PharmCAT //')
        pharmcat_data_version: "\$(python3 -c 'import json; print(json.load(open("${meta.id}.report.json"))["dataVersion"])')"
        bcftools: \$(bcftools --version | head -1 | cut -d' ' -f2)
        pharmcat_reporter_sources: "${params.pharmcat_reporter_sources}"
        missing_to_ref: "not used"
    END_VERSIONS
    """

    stub:
    """
    echo '{"geneReports": {}}' > ${meta.id}.phenotype.json
    echo '{"drugs": {}, "dataVersion": "stub"}' > ${meta.id}.report.json
    touch ${meta.id}.match.json ${meta.id}.report.html ${meta.id}.preprocessed.vcf.bgz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: 3.4.0
        pharmcat_data_version: "stub"
        bcftools: 1.22
        pharmcat_reporter_sources: "${params.pharmcat_reporter_sources}"
        missing_to_ref: "not used"
    END_VERSIONS
    """
}
