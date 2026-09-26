/*
 * Analysis 2: modifier genotypes in the mechanistic panel (BUB1B axis; gene +/-10 kb).
 *   PANEL_VARIANTS (subject VCF, panel regions, split/normalised)
 *   GNOMAD_PANEL (v4.1 joint, remote region fetch) | CLINVAR_FETCH (dated, md5) |
 *   ALPHAMISSENSE_PANEL (md5, panel subset) | SPLICEAI (local, panel variants only)
 *   -> VEP_ANNOTATE (offline cache 116, MANE Select) -> MODIFIER_REPORT (all + rare & functional)
 * Hypothesis generation, not diagnostic classification: everything found is reported, VUS included.
 */
include { PANEL_VARIANTS      } from '../../../modules/local/panel_variants'
include { GNOMAD_PANEL        } from '../../../modules/local/gnomad_panel'
include { CLINVAR_FETCH       } from '../../../modules/local/clinvar_fetch'
include { ALPHAMISSENSE_PANEL } from '../../../modules/local/alphamissense_panel'
include { SPLICEAI            } from '../../../modules/local/spliceai'
include { VEP_ANNOTATE        } from '../../../modules/local/vep_annotate'
include { MODIFIER_REPORT     } from '../../../modules/local/modifier_report'

workflow MODIFIER_PANEL {
    take:
    ch_vcf      // [ meta, vcf, tbi ]
    ch_fasta    // [ ref_meta, fa ]
    ch_fai      // [ ref_meta, fai ]

    main:
    panel_bed = file(params.gene_panel_bed, checkIfExists: true)
    vep_cache = file(params.vep_cache_dir, checkIfExists: true)
    def marker = file("${params.vep_cache_dir}/homo_sapiens/${params.vep_cache_version}_GRCh38/.download_ok")
    if (!workflow.stubRun && !marker.exists()) {
        error "VEP cache ${params.vep_cache_version} not verified (${marker} missing). Run: bin/fetch_vep_cache.sh homo_sapiens GRCh38 ${params.vep_cache_version} ${params.vep_cache_dir}"
    }
    ch_ref = ch_fasta.join(ch_fai).first()

    PANEL_VARIANTS(ch_vcf, panel_bed, ch_ref)
    GNOMAD_PANEL(panel_bed)
    CLINVAR_FETCH(params.clinvar_release)
    ALPHAMISSENSE_PANEL(panel_bed)
    SPLICEAI(PANEL_VARIANTS.out.vcf, ch_ref)

    VEP_ANNOTATE(
        PANEL_VARIANTS.out.vcf.join(SPLICEAI.out.vcf),
        ch_ref,
        vep_cache,
        GNOMAD_PANEL.out.vcf,
        CLINVAR_FETCH.out.vcf,
        ALPHAMISSENSE_PANEL.out.vcf
    )
    MODIFIER_REPORT(VEP_ANNOTATE.out.vcf, panel_bed)

    ch_versions = Channel.empty().mix(PANEL_VARIANTS.out.versions, GNOMAD_PANEL.out.versions, ALPHAMISSENSE_PANEL.out.versions,
                                      SPLICEAI.out.versions, VEP_ANNOTATE.out.versions, MODIFIER_REPORT.out.versions)

    emit:
    rare_functional = MODIFIER_REPORT.out.rare_functional
    all             = MODIFIER_REPORT.out.all
    summary         = MODIFIER_REPORT.out.summary
    versions        = ch_versions
}
