#!/usr/bin/env nextflow
/*
 * mva-track2: MVA Hackathon 2026, Track 2 (mechanism and drug repurposing).
 * Not a variant-discovery pipeline. Research use only; outputs are hypotheses, not
 * treatment recommendations.
 */
nextflow.enable.dsl = 2

include { MVA } from './workflows/mva'

workflow {
    if (workflow.profile.tokenize(',').contains('stub') && !workflow.stubRun) {
        error "The 'stub' profile uses synthetic inputs and must be run with -stub-run"
    }
    // fail fast on resource parameters that modules read but config never defined (an undefined
    // param is null in Nextflow, which once surfaced only as a download error deep in a task)
    def required = ['ref_fasta_url', 'ref_fasta_md5', 'imgt_hla_dat_url', 'imgt_hla_dat_sha256', 'umap_url', 'umap_sha256',
                    'blacklist_url', 'blacklist_sha256', 'cytoband_url', 'cytoband_sha256', 'segdup_url', 'segdup_sha256',
                    'delly_exclude_url', 'delly_exclude_sha256', 'somalier_sites_url', 'somalier_sites_sha256',
                    'somalier_1kg_url', 'somalier_1kg_sha256', 'somalier_labels_url', 'somalier_labels_sha256',
                    'alphamissense_url', 'alphamissense_md5_b64', 'clinvar_release', 'pypgx_bundle_sha256',
                    'verifybamid2_resource', 'skip_qc', 'skip_report', 'skip_aneuploidy', 'skip_modifier_panel', 'skip_pgx', 'skip_hla']
    def undefined = required.findAll { params[it] == null }
    if (undefined) error "Parameters not defined in config: ${undefined.join(', ')}"

    def on_apple_silicon = System.getProperty('os.name') == 'Mac OS X' && System.getProperty('os.arch') == 'aarch64'
    if (workflow.profile.tokenize(',').contains('conda') && on_apple_silicon && System.getenv('CONDA_SUBDIR') != 'osx-arm64') {
        log.warn "Apple Silicon without CONDA_SUBDIR=osx-arm64: an x86_64 conda will build Intel envs that run under Rosetta"
    }
    log.info """
    mva-track2 ${workflow.manifest.version}
      sample      : ${params.sample_id}
      genome      : ${params.genome} (${params.ref_name})
      outdir      : ${params.outdir}
      stub run    : ${workflow.stubRun}
    """.stripIndent()

    MVA()
}
