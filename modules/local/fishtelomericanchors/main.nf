process FISHTELOMERICANCHORS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'quay.io/biocontainers/python:3.10.4' }"

    input:
    tuple val(meta), path(reads)
    val telomere_repeat

    output:
    tuple val(meta), path("*.candidate_anchors.fasta"), emit: candidate_anchors
    tuple val(meta), path("*.pass1_summary.json")     , emit: summary_json
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def platform = meta.platform ?: "hifi"
    def reads_arg = meta.single_end ? "--reads-1 ${reads}" : "--reads-1 ${reads[0]} --reads-2 ${reads[1]}"

    """
    # Tool: cbm_fish_telomeric_anchors.py | Script-SHA: d3909feccd93f96f
    cbm_fish_telomeric_anchors.py \\
        ${reads_arg} \\
        --platform ${platform} \\
        --telomere-repeat ${telomere_repeat} \\
        --output-fasta "${prefix}.candidate_anchors.fasta" \\
        --output-summary "${prefix}.pass1_summary.json" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.candidate_anchors.fasta"
    echo '{"platform": "hifi", "candidate_anchors_count": 0}' > "${prefix}.pass1_summary.json"
    """
}
