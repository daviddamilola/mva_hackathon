/*
 * MT-RNR1 (12S rRNA) for the CPIC aminoglycoside guideline. Not in the Aldy panel or PharmCAT's
 * matcher, so it is genotyped here and passed to PharmCAT as an outside call, as HLA is.
 * mpileup over the whole gene on the mitochondrial contig (haploid; no diploid genotype model):
 * allele fractions from FORMAT/AD, so heteroplasmy is reported as a fraction. MAPQ >= 20 limits
 * NUMT-derived reads; max depth raised because mtDNA depth is in the thousands.
 */
process MT_RNR1_CALL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(cram), path(crai), path(ref_check)
    path cpic_alleles
    tuple val(ref_meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.mt_rnr1.outside.tsv"), emit: outside
    path "${meta.id}.mt_rnr1.json"                         , emit: json
    path "${meta.id}.mt_rnr1.pileup.vcf.gz"                , emit: pileup
    path "versions.yml"                                    , emit: versions

    script:
    // digest of the bin script: part of the task hash (Nextflow does not hash bin/), and recorded in provenance
    def script_md5 = file("${projectDir}/bin/mt_rnr1_call.py").text.md5()
    def region = "${params.mt_contig}:${params.mt_rnr1_start}-${params.mt_rnr1_end}"
    """
    bcftools mpileup --fasta-ref ${fasta} -r ${region} -d 100000 -q 20 -Q 20 -a FORMAT/AD,FORMAT/DP \\
        -Oz -o ${meta.id}.mt_rnr1.pileup.vcf.gz ${cram}
    mt_rnr1_call.py --pileup ${meta.id}.mt_rnr1.pileup.vcf.gz --alleles ${cpic_alleles} --fasta ${fasta} \\
        --contig ${params.mt_contig} --gene-start ${params.mt_rnr1_start} --gene-end ${params.mt_rnr1_end} \\
        --sample ${meta.id} --out-json ${meta.id}.mt_rnr1.json --out-outside ${meta.id}.mt_rnr1.outside.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | cut -d' ' -f2)
        mpileup_args: "-r ${region} -d 100000 -q 20 -Q 20 -a FORMAT/AD,FORMAT/DP"
        cpic_mt_rnr1_alleles: "\$(grep -m1 sha256 ${cpic_alleles} | sed 's/.*sha256 of the JSON response: //')"
        mt_rnr1_call_py_md5: "${script_md5}"
    END_VERSIONS
    """

    stub:
    """
    printf 'MT-RNR1\\tReference\\n' > ${meta.id}.mt_rnr1.outside.tsv
    echo '{"stub": true}' > ${meta.id}.mt_rnr1.json
    touch ${meta.id}.mt_rnr1.pileup.vcf.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.22
    END_VERSIONS
    """
}
