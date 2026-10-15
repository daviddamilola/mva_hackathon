/*
 * One lane end to end, to bound disk use: download R1/R2 from the gated HF dataset into the task
 * directory, verify sha256 against the upstream LFS oid, align with bwa mem, fixmate, then
 * coordinate-sort to CRAM. The FASTQ is deleted when the task exits, on success or failure.
 * Only the lane CRAM (~5 GB) is kept, so -resume never needs the FASTQ again.
 * HF auth: the task inherits HF_TOKEN or ~/.cache/huggingface/token. The token is never written.
 */
process FETCH_ALIGN_LANE {
    tag "${meta.id}.L${meta.lane}"
    label 'process_high'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    // no single biocontainer ships hf + bwa + samtools; docker/singularity profiles build this from
    // environment.yml with Wave

    input:
    tuple val(meta), val(r1), val(r1_sha256), val(r2), val(r2_sha256)
    tuple val(ref_meta), path(fasta), path(fai), path(bwa_index)

    output:
    tuple val(meta), path("${prefix}.cram"), emit: cram
    path "${prefix}.fastq.sha256.tsv"      , emit: checksums
    path "versions.yml"                    , emit: versions

    script:
    prefix = "${meta.id}.${meta.flowcell}.L${meta.lane}"
    def rg = "@RG\\tID:${meta.flowcell}.${meta.lane}\\tSM:${meta.id}\\tLB:${meta.id}\\tPL:ILLUMINA\\tPU:${meta.flowcell}.${meta.lane}.${meta.sample_index}"
    def sort_threads = Math.max(1, (task.cpus as int).intdiv(4))
    def sort_mem_mb = 768
    """
    trap 'rm -rf fastq' EXIT

    hf download ${params.hf_repo} ${r1} ${r2} --repo-type ${params.hf_repo_type} --local-dir fastq
    verify_checksum.py sha256 ${r1_sha256} fastq/${r1} ${prefix}.fastq.sha256.tsv
    verify_checksum.py sha256 ${r2_sha256} fastq/${r2} ${prefix}.fastq.sha256.tsv

    INDEX=\$(find -L ${bwa_index} -name '*.bwt' | sed 's/\\.bwt\$//')

    # -K fixes the batch size so output is independent of thread count; -Y soft-clips supplementary
    bwa mem -t ${task.cpus} -K 100000000 -Y -R '${rg}' \$INDEX fastq/${r1} fastq/${r2} \\
        | samtools fixmate -@ 2 -m -u - - \\
        | samtools sort -@ ${sort_threads} -m ${sort_mem_mb}M -T sort_tmp \\
            --reference ${fasta} -O cram -o ${prefix}.cram -

    rm -rf fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | sed -n 's/^Version: //p')
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
        huggingface_hub: \$(python3 -c 'import huggingface_hub as h; print(h.__version__)')
        bwa_args: "mem -K 100000000 -Y"
    END_VERSIONS
    """

    stub:
    prefix = "${meta.id}.${meta.flowcell}.L${meta.lane}"
    """
    touch ${prefix}.cram
    printf 'fastq/${r1}\\tsha256\\t${r1_sha256}\\tSTUB\\n' > ${prefix}.fastq.sha256.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: 0.7.19
        samtools: 1.22.1
        huggingface_hub: 0.35.3
        bwa_args: "mem -K 100000000 -Y"
    END_VERSIONS
    """
}
