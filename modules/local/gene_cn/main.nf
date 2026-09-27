process GENE_CN {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), val(label), path(exon_regions)
    path exons_annot
    path bins_with_depth
    tuple val(meta2), path(delly_vcf), path(delly_tbi)
    path panel_bed
    path arms_json

    output:
    tuple val(meta), path("${meta.id}.panel.gene_cn.tsv"), emit: genes
    path "${meta.id}.panel.exon_cn.tsv"                   , emit: exons
    path "${meta.id}.panel.sv_panel.tsv"                  , emit: sv
    path "${meta.id}.panel.cnv_sv.summary.json"           , emit: summary
    path "versions.yml"                                   , emit: versions

    script:
    """
    gene_cn.py --exon-depth ${exon_regions} --exons ${exons_annot} --bins ${bins_with_depth} \\
        --delly ${delly_vcf} --panel ${panel_bed} --arms ${arms_json} --prefix ${meta.id}.panel \\
        --min-map ${params.aneuploidy_min_map}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.panel.gene_cn.tsv ${meta.id}.panel.exon_cn.tsv ${meta.id}.panel.sv_panel.tsv
    echo '{"stub": true}' > ${meta.id}.panel.cnv_sv.summary.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.11
    END_VERSIONS
    """
}
