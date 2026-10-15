/*
 * somalier 0.3.5: extract genotypes at its 17k hg38 sites from the CRAM (explicit reference), then
 *  - ancestry: SVM on PCA of the 1000 Genomes profiles (5 superpopulations, pinned labels)
 *  - relate: per-sample stats incl. X/Y depth and X heterozygosity (sex check)
 * Resources (sha256-pinned) are fetched in-task; the sites VCF is renamed to no-chr. Site order is
 * unchanged, so the chr-named 1000 Genomes profiles stay comparable.
 */
process SOMALIER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/somalier:0.3.5--h5205c93_0' :
        'quay.io/biocontainers/somalier:0.3.5--h5205c93_0' }"

    input:
    tuple val(meta), path(cram), path(crai), path(ref_check)
    tuple val(ref_meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.somalier-ancestry.tsv"), path("${meta.id}.relate.samples.tsv"), emit: results
    path "${meta.id}.somalier-ancestry.html", emit: html
    path "*.somalier"                       , emit: profile
    path "versions.yml"                     , emit: versions

    script:
    """
    set -o pipefail
    fetch() { python3 -c "import sys,urllib.request; urllib.request.urlretrieve(sys.argv[1], sys.argv[2])" "\$1" "\$2" && verify_checksum.py sha256 "\$3" "\$2"; }
    fetch '${params.somalier_sites_url}'  sites.hg38.vcf.gz   ${params.somalier_sites_sha256}
    fetch '${params.somalier_1kg_url}'    1kg.somalier.tar.gz ${params.somalier_1kg_sha256}
    fetch '${params.somalier_labels_url}' labels.tsv          ${params.somalier_labels_sha256}
    gzip -dc sites.hg38.vcf.gz | sed -E 's/^chr//; s/##contig=<ID=chr/##contig=<ID=/' | bgzip -c > sites.nochr.vcf.gz
    tar -xzf 1kg.somalier.tar.gz

    somalier extract -d extracted --sites sites.nochr.vcf.gz -f ${fasta} ${cram}
    cp extracted/*.somalier .
    somalier ancestry --labels labels.tsv 1kg-somalier/*.somalier ++ extracted/*.somalier -o ${meta.id}
    somalier relate -o ${meta.id}.relate extracted/*.somalier
    rm -rf 1kg-somalier 1kg.somalier.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(somalier 2>&1 | grep -oE 'version: [0-9.]+' | head -1 | sed 's/version: //')
        ancestry_reference: "1000 Genomes phase 3 (2504 samples, somalier profiles, zenodo 3479773)"
    END_VERSIONS
    """

    stub:
    """
    printf '#sample_id\\tpredicted_ancestry\\tgiven_ancestry\\tAFR_prob\\tEUR_prob\\tEAS_prob\\tAMR_prob\\tSAS_prob\\n${meta.id}\\tEUR\\t\\t0\\t1\\t0\\t0\\t0\\n' > ${meta.id}.somalier-ancestry.tsv
    printf '#family_id\\tsample_id\\tdepth_mean\\tX_depth_mean\\tX_n\\tX_het\\tY_depth_mean\\tp_middling_ab\\n${meta.id}\\t${meta.id}\\t45\\t22\\t300\\t0\\t20\\t0.01\\n' > ${meta.id}.relate.samples.tsv
    touch ${meta.id}.somalier-ancestry.html ${meta.id}.somalier
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: 0.3.5
    END_VERSIONS
    """
}
