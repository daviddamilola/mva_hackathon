/*
 * Stage 1a: fixed-panel pharmacogenomic calling from the CRAM.
 *   VALIDATE_CRAM_REF  (first node: CRAM @SQ M5 == reference dict; nothing below starts without it)
 *   PHARMCAT_SITES -> PGX_FORCE_CALL (constrained + unconstrained, merged) -> PHARMCAT_RUN
 *   ALDY_GENOTYPE (per gene) -> ALDY_PARSE (outside calls for PharmCAT, e.g. CYP2D6)
 *   PGX_RECONCILE (PharmCAT vs Aldy per gene; disagreements reported, neither caller preferred)
 *   PYPGX_* (Stargazer successor, in parallel with Aldy) -> STAR_ALLELE_COMPARE (Aldy vs PyPGx)
 */
include { VALIDATE_CRAM_REF  } from '../../../modules/local/validate_cram_ref'
include { PHARMCAT_SITES     } from '../../../modules/local/pharmcat_sites'
include { PGX_FORCE_CALL     } from '../../../modules/local/pgx_force_call'
include { PHARMCAT_REFERENCE } from '../../../modules/local/pharmcat_reference'
include { ALDY_GENOTYPE      } from '../../../modules/local/aldy_genotype'
include { ALDY_PARSE         } from '../../../modules/local/aldy_parse'
include { PHARMCAT_RUN       } from '../../../modules/local/pharmcat_run'
include { PGX_RECONCILE      } from '../../../modules/local/pgx_reconcile'
include { PYPGX_BUNDLE        } from '../../../modules/local/pypgx_bundle'
include { PYPGX_PREPARE       } from '../../../modules/local/pypgx_prepare'
include { PYPGX_GENOTYPE      } from '../../../modules/local/pypgx_genotype'
include { STAR_ALLELE_COMPARE } from '../../../modules/local/star_allele_compare'

workflow PGX_CALL {
    take:
    ch_cram     // [ meta, cram, crai ]   meta: [ id, sex ]
    ch_fasta    // [ ref_meta, fa ]
    ch_fai      // [ ref_meta, fai ]
    ch_dict     // [ ref_meta, dict ]

    main:
    ch_ref      = ch_fasta.join(ch_fai).first()
    ch_ref_dict = ch_fasta.join(ch_fai).join(ch_dict).first()

    VALIDATE_CRAM_REF(ch_cram, ch_ref_dict)
    ch_valid_cram = VALIDATE_CRAM_REF.out.cram      // [ meta, cram, crai, ref_check ]

    PHARMCAT_SITES()
    PGX_FORCE_CALL(ch_valid_cram, PHARMCAT_SITES.out.sites, ch_ref)

    ch_genes = Channel.fromList(params.aldy_genes.tokenize(',').collect { it.trim() })
    ALDY_GENOTYPE(ch_valid_cram.combine(ch_genes), ch_ref)
    ALDY_PARSE(ALDY_GENOTYPE.out.aldy.groupTuple())

    PHARMCAT_REFERENCE(ch_ref)
    PHARMCAT_RUN(
        PGX_FORCE_CALL.out.vcf
            .join(ALDY_PARSE.out.calls)
            .map { meta, vcf, csi, aldy_json, outside -> [ meta, vcf, csi, outside ] },
        PHARMCAT_SITES.out.chr_map,
        PHARMCAT_REFERENCE.out.fasta
    )

    PGX_RECONCILE(
        PHARMCAT_RUN.out.phenotype
            .join(ALDY_PARSE.out.calls.map { meta, aldy_json, outside -> [ meta, aldy_json ] })
            .join(PGX_FORCE_CALL.out.off_definition)
    )

    // second star-allele caller, in parallel with Aldy
    ch_pypgx_versions = Channel.empty()
    ch_star_compare   = Channel.empty()
    if (!params.skip_pypgx) {
        PYPGX_BUNDLE(Channel.value([ params.pypgx_bundle_version, params.pypgx_bundle_sha256 ]))
        PYPGX_PREPARE(ch_valid_cram, ch_ref, PYPGX_BUNDLE.out.bundle)
        ch_pgenes = Channel.fromList(params.pypgx_genes.tokenize(',').collect { it.trim() })
        PYPGX_GENOTYPE(PYPGX_PREPARE.out.inputs.combine(ch_pgenes), PYPGX_BUNDLE.out.bundle)
        STAR_ALLELE_COMPARE(
            ALDY_PARSE.out.calls
                .map { meta, aldy_json, outside -> [ meta, aldy_json ] }
                .join(PYPGX_GENOTYPE.out.calls.groupTuple())
        )
        ch_star_compare   = STAR_ALLELE_COMPARE.out.tsv
        ch_pypgx_versions = ch_pypgx_versions.mix(PYPGX_PREPARE.out.versions, PYPGX_GENOTYPE.out.versions.first(),
                                                  STAR_ALLELE_COMPARE.out.versions)
    }

    ch_versions = ch_pypgx_versions.mix(
        VALIDATE_CRAM_REF.out.versions,
        PHARMCAT_SITES.out.versions,
        PGX_FORCE_CALL.out.versions,
        ALDY_GENOTYPE.out.versions.first(),
        ALDY_PARSE.out.versions,
        PHARMCAT_REFERENCE.out.versions,
        PHARMCAT_RUN.out.versions,
        PGX_RECONCILE.out.versions
    )

    emit:
    valid_cram = ch_valid_cram                  // [ meta, cram, crai, ref_check ] for other branches
    calls      = PGX_RECONCILE.out.calls        // [ meta, tsv, json ]
    report     = PHARMCAT_RUN.out.report        // [ meta, report.json ]
    forcecall  = PGX_FORCE_CALL.out.vcf
    star_compare = ch_star_compare               // [ meta, aldy_vs_pypgx.tsv ]
    versions   = ch_versions
}
