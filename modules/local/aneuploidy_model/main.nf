process ANEUPLOIDY_MODEL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // no single biocontainer with the pinned numpy/scipy/matplotlib set; Wave builds from environment.yml

    input:
    tuple val(meta), val(label), path(regions)
    path bins_annot
    tuple val(meta2), path(vcf), path(tbi)
    tuple val(cyto_name), path(cytoband)
    tuple val(sd_name), path(segdup)

    output:
    tuple val(meta), path("${meta.id}.aneuploidy.summary.json"), emit: summary
    path "${meta.id}.aneuploidy.arms.json"                     , emit: arms_json
    path "${meta.id}.aneuploidy.arms.tsv"                      , emit: arms
    path "${meta.id}.aneuploidy.bins.tsv.gz"                   , emit: bins
    path "${meta.id}.aneuploidy.*.png"                         , emit: figures
    path "bins_with_depth.tsv"                                 , emit: bins_depth
    path "versions.yml"                                        , emit: versions

    script:
    """
    # join mosdepth bin means onto the annotated bins (same makewindows order, checked)
    python3 - <<'PYEOF'
    import gzip
    depth = {}
    for l in gzip.open("${regions}", "rt"):
        c, s, e, m = l.rstrip("\\n").split("\\t")[:4]
        depth[(c, s)] = m
    n_miss = 0
    with open("${bins_annot}") as fi, open("bins_with_depth.tsv", "w") as fo:
        fo.write(next(fi).rstrip("\\n") + "\\tdepth\\n")
        for l in fi:
            f = l.rstrip("\\n").split("\\t")
            d = depth.get((f[0], f[1]))
            if d is None:
                n_miss += 1
                continue
            fo.write("\\t".join(f + [d]) + "\\n")
    if n_miss:
        raise SystemExit(f"{n_miss} bins lack a mosdepth value")
    PYEOF

    aneuploidy_model.py --bins bins_with_depth.tsv --vcf ${vcf} --cytoband ${cytoband} --segdup ${segdup} --sex ${meta.sex} \\
        --prefix ${meta.id}.aneuploidy --min-map ${params.aneuploidy_min_map}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
        scipy: \$(python3 -c 'import scipy; print(scipy.__version__)')
        model: "direct implementation (bin/aneuploidy_model.py), bins ${params.bin_size} bp, min mappability ${params.aneuploidy_min_map}"
    END_VERSIONS
    """

    stub:
    """
    echo '{"stub": true}' > ${meta.id}.aneuploidy.summary.json
    echo '[]' > ${meta.id}.aneuploidy.arms.json
    touch ${meta.id}.aneuploidy.arms.tsv ${meta.id}.aneuploidy.bins.tsv.gz ${meta.id}.aneuploidy.genome_lrr_baf.png bins_with_depth.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.11
        model: "stub"
    END_VERSIONS
    """
}
