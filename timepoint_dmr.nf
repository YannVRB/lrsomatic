#!/usr/bin/env nextflow
/*
    Timepoint DMR: differentially methylated regions of each tumor against a baseline tumor of the same patient.
    Entry script of the Cirro process "Timepoint DMR (lrsomatic dev, test)", run as `nextflow run . -main-script timepoint_dmr.nf`.

    --comparisons  CSV with columns comparison,treatment_name,treatment_bed,control_name,control_bed (one row per
                   comparison; treatment = later sample, control = baseline), written by the Cirro preprocess step
    --outdir       output directory; --genome is GRCh38 (the DMR reference BEDs and fasta come from conf/igenomes.config)
*/

include { TIMEPOINT_DMR  } from './subworkflows/local/timepoint_dmr/main'
include { SAMTOOLS_FAIDX as TIMEPOINT_FAIDX } from './modules/nf-core/samtools/faidx/main'

workflow {

    main:
    if (!params.comparisons) {
        error("--comparisons is required: a CSV with columns comparison,treatment_name,treatment_bed,control_name,control_bed")
    }
    // GRCh38 only; the Cirro process also sets it. dev leaves --genome null by default.
    def genome_name = params.genome ?: 'GRCh38'
    def genome_attrs = params.genomes[genome_name]
    if (!genome_attrs) {
        error("--genome '${genome_name}' is not a known genome")
    }

    channel
        .fromPath(params.comparisons, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            [ [ id: row.comparison ], file(row.treatment_bed, checkIfExists: true), file(row.control_bed, checkIfExists: true) ]
        }
        .set { ch_comparisons }

    ch_fasta = channel.value([ [:], file(genome_attrs.fasta, checkIfExists: true) ])
    TIMEPOINT_FAIDX ( ch_fasta.map { meta, fasta -> [ meta, fasta, [] ] }, false )
    ch_fai = TIMEPOINT_FAIDX.out.fai

    TIMEPOINT_DMR (
        ch_comparisons,
        ch_fasta,
        ch_fai,
        file(genome_attrs.dmr_cpg_islands_bed, checkIfExists: true),
        file(genome_attrs.dmr_gencode_gene_bed, checkIfExists: true)
    )
}
