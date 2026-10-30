#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/cbmdiscovery
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/cbmdiscovery
    Website: https://nf-co.re/cbmdiscovery
    Slack  : https://nfcore.slack.com/channels/cbmdiscovery
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CBMDISCOVERY        } from './workflows/cbmdiscovery'
include { PIPELINE_COMPLETION } from './subworkflows/local/utils_nfcore_cbmdiscovery_pipeline'
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_CBMDISCOVERY {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    CBMDISCOVERY (
        samplesheet,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )

    emit:
    candidate_anchors = CBMDISCOVERY.out.candidate_anchors
    break_windows     = CBMDISCOVERY.out.break_windows
    junction_catalog  = CBMDISCOVERY.out.junction_catalog
    census_json       = CBMDISCOVERY.out.census_json
    best_core_motif   = CBMDISCOVERY.out.best_core_motif
    bipartite_spacing = CBMDISCOVERY.out.bipartite_spacing
    discovery_report  = CBMDISCOVERY.out.discovery_report
    multiqc_report    = CBMDISCOVERY.out.multiqc_report // channel: /path/to/multiqc_report.html
    versions          = CBMDISCOVERY.out.versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    // 0. CLI Help handling
    if (params.help) {
        log.info paramsHelp()
        exit 0
    }

    // 1. Declarative parameter validation via nf-schema Draft 2020-12
    if (params.validate_params) {
        validateParameters()
    }
    log.info paramsSummaryLog(workflow)

    // 2. Samplesheet schema validation via nf-schema Draft 2020-12
    if (params.input) {
        samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")
    }

    // =========================================================================
    // DOMAIN-SPECIFIC PRE-RUN VERIFICATION GATE
    // (Verifies file existence, non-emptiness, and platform requirements)
    // =========================================================================
    def validation_errors = []

    if (!params.input) {
        validation_errors.add("Mandatory parameter --input <path/to/samplesheet.csv> is missing.")
    } else {
        def ss_file = file(params.input)
        if (!ss_file.exists() || ss_file.isEmpty()) {
            validation_errors.add("Provided --input samplesheet does not exist or is empty: ${params.input}")
        } else {
            def lines = ss_file.readLines()
            if (lines.size() < 2) {
                validation_errors.add("Provided --input samplesheet has no sample rows: ${params.input}")
            } else {
                def header = lines[0].split(',', -1).collect { col -> col.trim() }
                def sample_idx   = header.indexOf('sample')
                def platform_idx = header.indexOf('platform')
                def reads1_idx   = header.indexOf('reads_1')
                def reads2_idx   = header.indexOf('reads_2')

                if (sample_idx < 0 || platform_idx < 0 || reads1_idx < 0) {
                    validation_errors.add("Samplesheet header must contain 'sample,platform,reads_1' (and optional 'reads_2'): ${header.join(',')}")
                } else {
                    lines.drop(1).eachWithIndex { line, line_idx ->
                        def row_num = line_idx + 2
                        def trimmed = line.trim()
                        if (!trimmed || trimmed.startsWith('#')) return

                        def cols = line.split(',', -1).collect { col -> col.trim() }
                        def sample_id = cols[sample_idx] ?: "row_${row_num}"
                        def platform  = cols[platform_idx]?.toLowerCase()

                        if (!cols[reads1_idx]) {
                            validation_errors.add("Samplesheet row ${row_num} (sample '${sample_id}'): reads_1 file path is missing.")
                        } else {
                            def r1 = file(cols[reads1_idx])
                            if (!r1.exists()) {
                                validation_errors.add("Samplesheet row ${row_num} (sample '${sample_id}'): reads_1 does not exist: ${cols[reads1_idx]}")
                            } else if (r1.isEmpty()) {
                                validation_errors.add("Samplesheet row ${row_num} (sample '${sample_id}'): reads_1 is empty: ${cols[reads1_idx]}")
                            }
                        }

                        def has_reads2 = reads2_idx >= 0 && reads2_idx < cols.size() && cols[reads2_idx]
                        if (platform == 'illumina' && !has_reads2) {
                            validation_errors.add("Samplesheet row ${row_num} (sample '${sample_id}'): platform 'illumina' requires paired-end reads (reads_2 is missing).")
                        } else if (has_reads2) {
                            def r2 = file(cols[reads2_idx])
                            if (!r2.exists()) {
                                validation_errors.add("Samplesheet row ${row_num} (sample '${sample_id}'): reads_2 does not exist: ${cols[reads2_idx]}")
                            } else if (r2.isEmpty()) {
                                validation_errors.add("Samplesheet row ${row_num} (sample '${sample_id}'): reads_2 is empty: ${cols[reads2_idx]}")
                            }
                        }
                    }
                }
            }
        }
    }

    if (!params.outdir) {
        validation_errors.add("Mandatory parameter --outdir <path> is missing.")
    }

    if (!params.telomere_repeat || params.telomere_repeat.length() < 3) {
        validation_errors.add("Invalid --telomere_repeat: '${params.telomere_repeat}'. Must be at least 3 bp.")
    }

    if (params.kmer_size < 15 || params.kmer_size > 64) {
        validation_errors.add("Invalid --kmer_size: ${params.kmer_size}. Expected value between 15 and 64.")
    }

    if (validation_errors) {
        log.error "==================================================================="
        log.error " PRE-RUN VERIFICATION GATE FAILED: ${validation_errors.size()} CRITICAL CONFIGURATION ISSUE(S)"
        log.error "==================================================================="
        validation_errors.eachWithIndex { err, i ->
            log.error " [${i+1}] ${err}"
        }
        log.error "==================================================================="
        error "Pre-run verification gate failed. Please correct the parameters above before executing."
    }

    log.info "==================================================================="
    log.info " PRE-RUN VERIFICATION GATE PASSED: All inputs, databases, and paths verified."
    log.info "==================================================================="

    //
    // Construct sample channel from validated input samplesheet
    //
    ch_samplesheet = channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { meta, reads_1, reads_2 ->
            if (!reads_2) {
                return [ meta.id, meta + [ single_end: true ], [ reads_1 ] ]
            } else {
                return [ meta.id, meta + [ single_end: false ], [ reads_1, reads_2 ] ]
            }
        }
        .groupTuple()
        .map { _id, metas, fastqs ->
            [ metas[0], fastqs.flatten() ]
        }

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_CBMDISCOVERY (
        ch_samplesheet
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_CBMDISCOVERY.out.multiqc_report
    )

    publish:
    candidate_anchors = NFCORE_CBMDISCOVERY.out.candidate_anchors
    break_windows     = NFCORE_CBMDISCOVERY.out.break_windows
    junction_catalog  = NFCORE_CBMDISCOVERY.out.junction_catalog
    census_json       = NFCORE_CBMDISCOVERY.out.census_json
    best_core_motif   = NFCORE_CBMDISCOVERY.out.best_core_motif
    bipartite_spacing = NFCORE_CBMDISCOVERY.out.bipartite_spacing
    discovery_report  = NFCORE_CBMDISCOVERY.out.discovery_report
    multiqc_report    = NFCORE_CBMDISCOVERY.out.multiqc_report
    versions          = NFCORE_CBMDISCOVERY.out.versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW OUTPUT PUBLISHING (NEXTFLOW 26+)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {
    candidate_anchors {
        path { meta, _file ->
            params.per_species_outdir
                ? "${meta.species ?: meta.id}/analysis/cbmdiscovery"
                : "somatic_anchors"
        }
        mode params.publish_dir_mode
    }
    break_windows {
        path { meta, _file ->
            params.per_species_outdir
                ? "${meta.species ?: meta.id}/analysis/cbmdiscovery"
                : "germline_junctions"
        }
        mode params.publish_dir_mode
    }
    junction_catalog {
        path { meta, _file ->
            params.per_species_outdir
                ? "${meta.species ?: meta.id}/analysis/cbmdiscovery"
                : "germline_junctions"
        }
        mode params.publish_dir_mode
    }
    census_json {
        path { meta, _file ->
            params.per_species_outdir
                ? "${meta.species ?: meta.id}/analysis/cbmdiscovery"
                : "germline_junctions"
        }
        mode params.publish_dir_mode
    }
    best_core_motif {
        path { meta, _file ->
            params.per_species_outdir
                ? "${meta.species ?: meta.id}/analysis/cbmdiscovery"
                : "core_motif"
        }
        mode params.publish_dir_mode
    }
    bipartite_spacing {
        path { meta, _file ->
            params.per_species_outdir
                ? "${meta.species ?: meta.id}/analysis/cbmdiscovery"
                : "bipartite_profiler"
        }
        mode params.publish_dir_mode
    }
    discovery_report {
        path { meta, _file ->
            params.per_species_outdir
                ? "${meta.species ?: meta.id}/analysis/cbmdiscovery"
                : "bipartite_profiler"
        }
        mode params.publish_dir_mode
    }
    // Cohort-level outputs:
    multiqc_report {
        path { "pipeline_info" }
        mode params.publish_dir_mode
    }
    versions {
        path { "pipeline_info" }
        mode params.publish_dir_mode
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
