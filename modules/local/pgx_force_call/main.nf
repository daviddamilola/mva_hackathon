/*
 * Force-call PharmCAT positions from the CRAM so reference genotypes are observed, not assumed.
 *
 *   constrained   : bcftools call -m -A -C alleles -T alleles.tsv.gz
 *                   every PharmCAT definition allele is genotyped, including hom-ref at indel
 *                   sites, and records keep PharmCAT's REF/ALT (needed by the preprocessor)
 *   unconstrained : bcftools call -m (sees alleles PharmCAT does not define)
 *   merged        : constrained, except where the unconstrained call carries an undefined
 *                   non-reference allele (merge_forcecalls.py). Those sites are listed in
 *                   off_definition.tsv.
 * mpileup: explicit --fasta-ref, -R merged windows (index jump) + -T exact positions, -a AD,DP.
 * Ploidy: GRCh38 preset (no-chr names) with the sample's sex, so chrX is haploid in XY samples.
 */
process PGX_FORCE_CALL {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.22--h3a4d415_1' :
        'quay.io/biocontainers/bcftools:1.22--h3a4d415_1' }"

    input:
    tuple val(meta), path(cram), path(crai), path(ref_check)
    tuple path(positions), path(positions_csi), path(alleles), path(alleles_tbi), path(regions)
    tuple val(ref_meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.pgx_forcecall.vcf.gz"), path("${meta.id}.pgx_forcecall.vcf.gz.csi"), emit: vcf
    path "${meta.id}.pgx_forcecall.constrained.vcf.gz"  , emit: constrained
    path "${meta.id}.pgx_forcecall.unconstrained.vcf.gz", emit: unconstrained
    tuple val(meta), path("${meta.id}.off_definition.tsv"), emit: off_definition
    path "${meta.id}.pgx_forcecall.evidence.tsv"        , emit: evidence
    path "${meta.id}.pgx_forcecall.evidence.summary.json", emit: evidence_summary
    path "${meta.id}.pgx_forcecall.summary.json"        , emit: summary
    path "versions.yml"                                 , emit: versions

    script:
    def sex = meta.sex ?: 'F'
    """
    printf '%s\\t%s\\n' '${meta.id}' '${sex}' > sample_sex.tsv

    bcftools mpileup --threads ${task.cpus} --fasta-ref ${fasta} -R ${regions} -T ${positions} \\
        -a FORMAT/AD,FORMAT/DP -Ou ${cram} > pileup.bcf
    bcftools call -m -A -C alleles -T ${alleles} --ploidy GRCh38 -S sample_sex.tsv \\
        -Oz -o ${meta.id}.pgx_forcecall.constrained.vcf.gz pileup.bcf
    bcftools call -m --ploidy GRCh38 -S sample_sex.tsv \\
        -Oz -o ${meta.id}.pgx_forcecall.unconstrained.vcf.gz pileup.bcf
    rm pileup.bcf

    merge_forcecalls.py ${meta.id}.pgx_forcecall.constrained.vcf.gz ${meta.id}.pgx_forcecall.unconstrained.vcf.gz \\
        merged.vcf ${meta.id}.off_definition.tsv ${meta.id}.pgx_forcecall.evidence.tsv \\
        ${params.pgx_min_depth} ${params.pgx_min_gq}
    bgzip -c merged.vcf > ${meta.id}.pgx_forcecall.vcf.gz
    bcftools index ${meta.id}.pgx_forcecall.vcf.gz

    python3 - <<'PYEOF' > ${meta.id}.pgx_forcecall.summary.json
    import gzip, json, collections
    gts, dps, n_pos = collections.Counter(), [], 0
    for l in gzip.open("${positions}", "rt"):
        if not l.startswith("#"): n_pos += 1
    for l in gzip.open("${meta.id}.pgx_forcecall.vcf.gz", "rt"):
        if l.startswith("#"): continue
        f = l.rstrip().split("\\t"); d = dict(zip(f[8].split(":"), f[9].split(":")))
        gts[d.get("GT", ".")] += 1
        if d.get("DP", ".") not in (".", ""): dps.append(int(d["DP"]))
    dps.sort()
    n_off = sum(1 for _ in open("${meta.id}.off_definition.tsv")) - 1
    print(json.dumps({"pharmcat_positions": n_pos, "records": sum(gts.values()), "genotypes": dict(gts),
                      "missing_gt": sum(v for k, v in gts.items() if "." in k),
                      "dp_min": dps[0] if dps else None, "dp_median": dps[len(dps)//2] if dps else None,
                      "dp_max": dps[-1] if dps else None, "off_definition_sites": n_off,
                      "sex_used_for_ploidy": "${sex}"}, indent=2))
    PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | cut -d' ' -f2)
        mpileup_args: "--fasta-ref -R regions -T positions -a FORMAT/AD,FORMAT/DP"
        call_args_constrained: "-m -A -C alleles -T alleles --ploidy GRCh38 -S sex"
        call_args_unconstrained: "-m --ploidy GRCh38 -S sex"
        min_evidence: "depth >= ${params.pgx_min_depth}, GQ >= ${params.pgx_min_gq} (GQ derived from PL); failing positions become no-calls"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.pgx_forcecall.vcf.gz ${meta.id}.pgx_forcecall.vcf.gz.csi
    touch ${meta.id}.pgx_forcecall.constrained.vcf.gz ${meta.id}.pgx_forcecall.unconstrained.vcf.gz
    printf 'chrom\\tpos\\tid\\tdefinition_ref\\tdefinition_alts\\tobserved_ref\\tobserved_undefined_alt\\tobserved_gt\\tobserved_ad\\tobserved_dp\\tconstrained_gt\\n' > ${meta.id}.off_definition.tsv
    echo '{"stub": true}' > ${meta.id}.pgx_forcecall.summary.json
    touch ${meta.id}.pgx_forcecall.evidence.tsv
    echo '{"stub": true}' > ${meta.id}.pgx_forcecall.evidence.summary.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.22
    END_VERSIONS
    """
}
