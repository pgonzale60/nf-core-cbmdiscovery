process BAITGERMLINEJUNCTIONS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/kmc_python:d875681ce2658826'

    input:
    tuple val(meta), path(candidate_anchors), path(reads)
    val telomere_repeat
    val kmer_size
    val min_germ_reads

    output:
    tuple val(meta), path("*.break_windows_300bp.fasta")            , emit: break_windows_fasta
    tuple val(meta), path("*.junction_catalog.tsv")                 , emit: junction_catalog
    tuple val(meta), path("*.coverage_and_chromosome_census.json")  , emit: census_json
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
    # Tool: cbm_bait_germline_junctions.py | Script-SHA: 7ab4ab88eb3d87dc

    # 1. Extract terminal 60 bp of candidate anchors into bait FASTA
    python3 -c "
with open('${candidate_anchors}') as f_in, open('anchors_terminal_60bp.fasta', 'w') as f_out:
    cur_id = None
    cur_seq = []
    for line in f_in:
        line = line.strip()
        if not line: continue
        if line.startswith('>'):
            if cur_id:
                s = ''.join(cur_seq)
                f_out.write(f'>{cur_id}\\n{s[-60:]}\\n')
            cur_id = line[1:].split()[0]
            cur_seq = []
        else:
            cur_seq.append(line.upper())
    if cur_id:
        s = ''.join(cur_seq)
        f_out.write(f'>{cur_id}\\n{s[-60:]}\\n')
"

    # 2. High-speed C++ k-mer baiting via KMC if available, else native Python streaming
    mkdir -p kmc_tmp
    if command -v kmc >/dev/null 2>&1 && command -v kmc_tools >/dev/null 2>&1; then
        kmc -k${kmer_size} -m${mem_gb} -t${task.cpus} -ci1 -fa anchors_terminal_60bp.fasta anchors_kmc kmc_tmp
        kmc_tools filter -t${task.cpus} anchors_kmc -ci1 "${r1}" -fq "baited_reads_1.fastq"
        if [ -n "${r2}" ]; then
            kmc_tools filter -t${task.cpus} anchors_kmc -ci1 "${r2}" -fq "baited_reads_2.fastq"
            READS_ARGS="--reads-1 baited_reads_1.fastq --reads-2 baited_reads_2.fastq"
        else
            READS_ARGS="--reads-1 baited_reads_1.fastq"
        fi
    else
        if [ -n "${r2}" ]; then
            READS_ARGS="--reads-1 ${r1} --reads-2 ${r2}"
        else
            READS_ARGS="--reads-1 ${r1}"
        fi
    fi

    # 3. Exact bifurcation partitioning, microhomology annotation, and junction census
    cbm_bait_germline_junctions.py \\
        --candidate-anchors "${candidate_anchors}" \\
        \${READS_ARGS} \\
        --platform ${platform} \\
        --telomere-repeat ${telomere_repeat} \\
        --kmer-size ${kmer_size} \\
        --min-germ-reads ${min_germ_reads} \\
        --output-catalog "${prefix}.junction_catalog.tsv" \\
        --output-fasta "${prefix}.break_windows_300bp.fasta" \\
        --output-census "${prefix}.coverage_and_chromosome_census.json" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.break_windows_300bp.fasta"
    touch "${prefix}.junction_catalog.tsv"
    echo '{"n_breaks": 0, "n_germ_telo": 0}' > "${prefix}.coverage_and_chromosome_census.json"
    """
}
