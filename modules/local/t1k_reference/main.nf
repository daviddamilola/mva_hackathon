/*
 * T1K HLA reference from a pinned IPD-IMGT/HLA release: hla.dat.zip fetched by commit URL, verified by
 * sha256, then built with t1k-build.pl. storeDir: built once.
 */
process T1K_REFERENCE {
    tag "IPD-IMGT/HLA ${params.imgt_hla_version}"
    label 'process_single'
    storeDir "${params.reference_dir}/t1k"

    conda "${moduleDir}/environment.yml"

    output:
    path "t1k_hla_${params.imgt_hla_version}", emit: ref

    script:
    """
    python3 - <<'PYEOF'
    import urllib.request, zipfile
    urllib.request.urlretrieve("${params.imgt_hla_dat_url}", "hla.dat.zip")
    PYEOF
    verify_checksum.py sha256 ${params.imgt_hla_dat_sha256} hla.dat.zip
    python3 -c "import zipfile; zipfile.ZipFile('hla.dat.zip').extractall('.')"
    mkdir -p t1k_hla_${params.imgt_hla_version}
    t1k-build.pl -d hla.dat -o t1k_hla_${params.imgt_hla_version} --prefix hla
    printf 'IPD-IMGT/HLA ${params.imgt_hla_version}\\nsource=${params.imgt_hla_dat_url}\\nsha256=${params.imgt_hla_dat_sha256}\\n' > t1k_hla_${params.imgt_hla_version}/SOURCE.txt
    rm -f hla.dat hla.dat.zip
    """

    stub:
    """
    mkdir -p t1k_hla_${params.imgt_hla_version} && touch t1k_hla_${params.imgt_hla_version}/hla_dna_seq.fa
    """
}
