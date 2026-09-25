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
