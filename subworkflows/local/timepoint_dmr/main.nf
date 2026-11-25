// IMPORT MODULES
include { DMR_CHECK_MODIFICATIONS } from '../../../modules/local/dmr/check_modifications/main'
include { HTSLIB_BGZIPTABIX       } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { DMR_HAPLOTYPE_REGIONS   } from '../../../modules/local/dmr/haplotype_regions/main'
include { MODKIT_DMR              } from '../../../modules/nf-core/modkit/dmr/main'
include { DMR_NEAREST_GENE        } from '../../../modules/local/dmr/nearest_gene/main'

// Differential methylation between two timepoints of one patient. Reuses the haplotype DMR modules
// unchanged: they take any two bedMethyl files, so the baseline goes where the haplotype pair would.
workflow TIMEPOINT_DMR {

    take:
    comparisons      // [meta, treatment_bed, control_bed]  -- meta = [id: comparison]; treatment = later sample, control = baseline
    fasta            // [[:], fasta]  -- value channel
    fai              // [[:], fai]    -- value channel
    cpg_islands_bed  // path          -- value channel
    gencode_gene_bed // path          -- value channel

    main:
    ch_versions = channel.empty()

    // One item per file, so every file is checked and indexed on its own
    comparisons
        .flatMap { meta, treatment, control ->
            [
                [ meta + [ id: "${meta.id}_treatment", comparison: meta.id, role: 'treatment' ], treatment ],
                [ meta + [ id: "${meta.id}_control",   comparison: meta.id, role: 'control'   ], control ]
            ]
        }
        .set { pileups }
    // pileups: [meta + comparison + role, bedmethyl]

    //
    // MODULE: DMR_CHECK_MODIFICATIONS (label: process_single)
    // Stops the run early when modkit dmr pair would reject the pileup (other modifications, no 5hmC)
    //
    DMR_CHECK_MODIFICATIONS ( pileups )

    //
    // MODULE: HTSLIB_BGZIPTABIX (label: process_low)
    // Tabix-index each pileup (already bgzipped by modkit --bgzf)
    //
    DMR_CHECK_MODIFICATIONS.out.checked
        .map { meta, bedmethyl -> [ meta, bedmethyl, [], [] ] }
        .set { bgziptabix_input }

    HTSLIB_BGZIPTABIX (
        bgziptabix_input,
        'compress',
        true,
        'bed'
    )

    HTSLIB_BGZIPTABIX.out.output
        .join(HTSLIB_BGZIPTABIX.out.index)
        .map { meta, bedgz, tbi -> [ meta.comparison, meta.role, bedgz, tbi ] }
        .branch { _comparison, role, _bedgz, _tbi ->
            treatment: role == 'treatment'
            control:   role == 'control'
        }
        .set { indexed }

    indexed.treatment
        .map { comparison, _role, bedgz, tbi -> [ comparison, bedgz, tbi ] }
        .set { treatment_indexed }
    indexed.control
        .map { comparison, _role, bedgz, tbi -> [ comparison, bedgz, tbi ] }
        .set { control_indexed }

    control_indexed
        .join(treatment_indexed)
        .map { comparison, control_bedgz, control_tbi, treatment_bedgz, treatment_tbi ->
            [ [ id: comparison ], control_bedgz, control_tbi, treatment_bedgz, treatment_tbi ]
        }
        .set { dmr_input }
    // dmr_input: [meta, control_bedgz, control_tbi, treatment_bedgz, treatment_tbi]

    //
    // MODULE: DMR_HAPLOTYPE_REGIONS (label: process_single)
    // Restrict the CpG islands to those covered in both samples (the baseline takes the first slot).
    //
    DMR_HAPLOTYPE_REGIONS (
        dmr_input,
        cpg_islands_bed,
        fai.first()
    )
    ch_versions = ch_versions.mix(DMR_HAPLOTYPE_REGIONS.out.versions)

    // modkit dmr pair fails on an empty regions file: skip those comparisons
    DMR_HAPLOTYPE_REGIONS.out.regions_bed
        .filter { meta, regions_bed ->
            if (regions_bed.size() == 0) {
                log.warn("DMR: no CpG island is covered in both samples of ${meta.id}; skipping modkit dmr pair")
                return false
            }
            return true
        }
        .set { nonempty_regions_bed }

    dmr_input
        .join(nonempty_regions_bed)
        .multiMap { meta, control_bedgz, control_tbi, treatment_bedgz, treatment_tbi, regions_bed ->
            treatment: [ meta, treatment_bedgz, treatment_tbi ]
            control:   [ meta, control_bedgz, control_tbi ]
            regions:   [ meta, regions_bed ]
        }
        .set { modkit_dmr_input }

    // .first(): as a queue channel MODKIT_DMR would run for one comparison only
    fasta
        .combine(fai)
        .map { meta, fasta_file, _meta2, fai_file -> [ meta, fasta_file, fai_file ] }
        .first()
        .set { fasta_with_fai }

    //
    // MODULE: MODKIT_DMR (label: process_medium)
    // Later timepoint (treatment) as -a, baseline (control) as -b
    //
    MODKIT_DMR (
        modkit_dmr_input.treatment,
        modkit_dmr_input.control,
        modkit_dmr_input.regions,
        fasta_with_fai
    )

    //
    // MODULE: DMR_NEAREST_GENE (label: process_single)
    //
    DMR_NEAREST_GENE (
        MODKIT_DMR.out.bed,
        gencode_gene_bed
    )
    ch_versions = ch_versions.mix(DMR_NEAREST_GENE.out.versions)

    emit:
    nearest_gene_bed = DMR_NEAREST_GENE.out.nearest_gene_bed  // [meta, bed]  -- DMRs annotated with the nearest gene
    dmr_log          = MODKIT_DMR.out.log                      // [meta, log]
    versions         = ch_versions
}
