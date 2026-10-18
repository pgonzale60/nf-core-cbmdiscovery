process FISHTELOMERICANCHORS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/kmc_python:d875681ce2658826'

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
    def is_single = meta.single_end
    def r1 = is_single ? "${reads}" : "${reads[0]}"
    def r2 = is_single ? "" : "${reads[1]}"
    def mem_gb = task.memory ? task.memory.toGiga() : 4

    """
    # Tool: cbm_fish_telomeric_anchors.py | Script-SHA: 5b3a02e5b09d5e59

    # 1. Unconditionally build KMC 18-mer database from 4 tandem repeats
    mkdir -p kmc_tmp
    python3 -c "
repeat = '${telomere_repeat}'.upper()
telo_4x = repeat * 4
with open('telo_seed.fasta', 'w') as f:
    f.write(f'>telo_seed\\n{telo_4x}\\n')
"
    kmc -k18 -m${mem_gb} -t${task.cpus} -ci1 -fa telo_seed.fasta telo_kmc kmc_tmp
    kmc_tools filter -t${task.cpus} telo_kmc -ci1 "${r1}" -fq "fished_1.fastq"
    if [ -n "${r2}" ]; then
        kmc_tools filter -t${task.cpus} telo_kmc -ci1 "${r2}" -fq "fished_2.fastq"
        READS_ARGS="--reads-1 fished_1.fastq --reads-2 fished_2.fastq"
    else
        READS_ARGS="--reads-1 fished_1.fastq"
    fi

    # 2. Extract somatic anchors from fished reads
    cbm_fish_telomeric_anchors.py \\
        \${READS_ARGS} \\
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
