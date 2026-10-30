/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { SAMTOOLS_BAM2FASTQ     } from '../modules/local/samtools_bam2fastq/main'
include { FISHTELOMERICANCHORS   } from '../modules/local/fishtelomericanchors/main'
include { BAITGERMLINEJUNCTIONS  } from '../modules/local/baitgermlinejunctions/main'
include { MEMECOREGRID           } from '../modules/local/memecoregrid/main'
include { BIPARTITEPROFILER      } from '../modules/local/bipartiteprofiler/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_cbmdiscovery_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CBMDISCOVERY {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // STAGE 0: On-the-fly unaligned BAM conversion
    //
    ch_reads_branch = ch_samplesheet.branch { meta, reads ->
        def r1 = meta.single_end ? reads : reads[0]
        bam: r1.name.endsWith('.bam')
        fastq: true
    }

    SAMTOOLS_BAM2FASTQ(ch_reads_branch.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_BAM2FASTQ.out.versions_samtools.first())

    ch_normalized_samplesheet = ch_reads_branch.fastq
        .mix(SAMTOOLS_BAM2FASTQ.out.fastq)

    //
    // MODULE: Run FastQC
    //
    FASTQC(ch_normalized_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })

    //
    // STAGE 1: Fish telomeric reads and construct candidate somatic anchors
    //
    FISHTELOMERICANCHORS (
        ch_normalized_samplesheet,
        params.telomere_repeat
    )

    //
    // STAGE 2 & 3: K-mer baiting, germline reconstruction, microhomology & junction census
    //
    ch_bait_input = FISHTELOMERICANCHORS.out.candidate_anchors
        .join(ch_normalized_samplesheet)

    BAITGERMLINEJUNCTIONS (
        ch_bait_input,
        params.telomere_repeat,
        params.kmer_size,
        params.min_germ_reads,
        params.min_telo_reads
    )

    //
    // STAGE 4: Tier-1 Core motif parameter grid search
    //
    MEMECOREGRID (
        BAITGERMLINEJUNCTIONS.out.break_windows_fasta
    )

    //
    // STAGE 5: Core-masking, Tier-2 auxiliary grid search & bipartite profiler
    //
    ch_bipartite_input = BAITGERMLINEJUNCTIONS.out.break_windows_fasta
        .join(MEMECOREGRID.out.best_core_motif)
        .join(BAITGERMLINEJUNCTIONS.out.census_json)
        .join(MEMECOREGRID.out.motif_metrics_json)
        .join(MEMECOREGRID.out.best_background)

    BIPARTITEPROFILER (
        ch_bipartite_input
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'cbmdiscovery_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'cbmdiscovery'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    telomeric_reads_fastq = FISHTELOMERICANCHORS.out.telomeric_reads_fastq
    baited_reads_fastq    = BAITGERMLINEJUNCTIONS.out.baited_reads_fastq
    candidate_anchors     = FISHTELOMERICANCHORS.out.candidate_anchors
    break_windows_fasta   = BAITGERMLINEJUNCTIONS.out.break_windows_fasta
    junction_catalog      = BAITGERMLINEJUNCTIONS.out.junction_catalog
    census_json           = BAITGERMLINEJUNCTIONS.out.census_json
    best_core_motif       = MEMECOREGRID.out.best_core_motif
    grid_summary_tsv      = MEMECOREGRID.out.grid_summary_tsv
    best_auxiliary_motif  = BIPARTITEPROFILER.out.best_auxiliary_motif
    bipartite_spacing_tsv = BIPARTITEPROFILER.out.bipartite_spacing_tsv
    discovery_report_md   = BIPARTITEPROFILER.out.discovery_report_md
    multiqc_report        = MULTIQC.out.report.map { _meta, report -> [report] }.toList()
    versions              = ch_collated_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
