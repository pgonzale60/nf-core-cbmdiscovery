process SAMTOOLS_BAM2FASTQ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/samtools_seqkit_gawk:5f1679612f236815':
        'community.wave.seqera.io/library/samtools_seqkit_gawk:5f1679612f236815' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | head -1 | sed 's/samtools //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools fastq -@ ${task.cpus} "${bam}" | gzip -c > "${prefix}.fastq.gz"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.fastq.gz"
    """
}
