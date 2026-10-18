process FISHTELOMERICANCHORS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community.wave.seqera.io/library/bbmap_vsearch:6fca92c6c19c5c7d':
        'community.wave.seqera.io/library/bbmap_vsearch:6fca92c6c19c5c7d' }"

    input:
    tuple val(meta), path(reads)
    val telomere_repeat

    output:
    tuple val(meta), path("*.candidate_anchors.fasta")              , emit: candidate_anchors
    tuple val(meta), path("*.pass1_summary.json")                   , emit: summary_json
    tuple val(meta), path("*.pass1_telomeric_reads.fastq.gz")       , emit: telomeric_reads_fastq
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def platform = meta.platform ?: "hifi"
    def is_single = meta.single_end
    def r1 = is_single ? "${reads}" : "${reads[0]}"
    def r2 = is_single ? "" : "${reads[1]}"

    """
    # Tool: cbm_fish_telomeric_anchors.py | Script-SHA: 56ad31ef629cddb5

    # 1. Fishing telomeric reads using bbduk.sh
    REPEAT_3X="${telomere_repeat}${telomere_repeat}${telomere_repeat}"
    bbduk.sh \\
        in=${r1} \\
        ${is_single ? '' : "in2=${r2}"} \\
        outm=${prefix}.pass1_telomeric_reads.fastq.gz \\
        literal="\$REPEAT_3X" \\
        k=18 \\
        rcomp=t \\
        hdist=0 \\
        threads=${task.cpus} \\
        overwrite=true

    # 2. Extract somatic anchors from fished reads
    cbm_fish_telomeric_anchors.py \\
        --reads-1 ${prefix}.pass1_telomeric_reads.fastq.gz \\
        --platform ${platform} \\
        --telomere-repeat ${telomere_repeat} \\
        --output-fasta "${prefix}.candidate_anchors.fasta" \\
        --output-summary "${prefix}.pass1_summary.json" \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.candidate_anchors.fasta"
    touch "${prefix}.pass1_telomeric_reads.fastq.gz"
    echo '{"platform": "hifi", "candidate_anchors_count": 0}' > "${prefix}.pass1_summary.json"
    """
}
