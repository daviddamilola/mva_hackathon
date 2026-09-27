/*
 * VerifyBamID2 2.0.3: DNA contamination (FREEMIX) with the bundled hg38 SVD resource (1000 Genomes
 * phase 3, 100k sites), CRAM decoded with the explicit reference. The resource's .bed/.mu use chr
 * names and are renamed to no-chr in-task; .UD/.V carry no names.
 */
process VERIFYBAMID2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/verifybamid2:2.0.3--hc004090_0' :
        'quay.io/biocontainers/verifybamid2:2.0.3--hc004090_0' }"

    input:
    tuple val(meta), path(cram), path(crai), path(ref_check)
    tuple val(ref_meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.selfSM"), emit: selfsm
    path "${meta.id}.Ancestry"                , emit: ancestry
    path "versions.yml"                       , emit: versions

    script:
    def res = params.verifybamid2_resource
    """
    RES=\$(ls -d \$(dirname \$(command -v verifybamid2))/../share/verifybamid2*/resource | head -1)
    for x in UD V; do cp \$RES/${res}.vcf.gz.dat.\$x svd.dat.\$x; done
    for x in bed mu; do sed 's/^chr//' \$RES/${res}.vcf.gz.dat.\$x > svd.dat.\$x; done
    verifybamid2 --SVDPrefix svd.dat --Reference ${fasta} --BamFile ${cram} --NumThread ${task.cpus} --Output ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verifybamid2: "2.0.3"
        svd_resource: "${res}"
    END_VERSIONS
    """

    stub:
    """
    printf '#SEQ_ID\\tRG\\tCHIP_ID\\t#SNPS\\t#READS\\tAVG_DP\\tFREEMIX\\n${meta.id}\\tNA\\tNA\\t100000\\t1\\t40\\t0.0\\n' > ${meta.id}.selfSM
    touch ${meta.id}.Ancestry
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verifybamid2: "2.0.3"
    END_VERSIONS
    """
}
