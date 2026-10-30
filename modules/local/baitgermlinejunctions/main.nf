process BAITGERMLINEJUNCTIONS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/baitgermlinejunctions:b348bbd6f96216b6':
        'community.wave.seqera.io/library/baitgermlinejunctions:b348bbd6f96216b6' }"

    input:
    tuple val(meta), path(candidate_anchors), path(reads)
    val telomere_repeat
    val kmer_size
    val min_germ_reads
    val min_telo_reads

    output:
    tuple val(meta), path("*.break_windows_300bp.fasta")            , emit: break_windows_fasta
    tuple val(meta), path("*.junction_catalog.tsv")                 , emit: junction_catalog
    tuple val(meta), path("*.coverage_and_chromosome_census.json")  , emit: census_json
    tuple val(meta), path("*.pass2_baited_reads.fastq.gz")          , emit: baited_reads_fastq
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
    # Tool: cbm_bait_germline_junctions.py | Script-SHA: 2e1b56c8181b80e3

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

    # 2. Mate-synchronized read baiting using bbduk.sh
    if [ "${is_single}" = "true" ]; then
        bbduk.sh \\
            in=${r1} \\
            ref=anchors_terminal_60bp.fasta \\
            outm=baited_R1.fastq \\
            k=27 \\
            hdist=1 \\
            rcomp=t \\
            minkmerhits=1 \\
            entropy=0.35 \\
            entropywindow=50 \\
            threads=${task.cpus} \\
            overwrite=true

        gzip -c baited_R1.fastq > ${prefix}.pass2_baited_reads.fastq.gz
        READS_ARGS="--reads-1 baited_R1.fastq"
    else
        bbduk.sh \\
            in1=${r1} \\
            in2=${r2} \\
            ref=anchors_terminal_60bp.fasta \\
            outm1=baited_R1.fastq \\
            outm2=baited_R2.fastq \\
            k=27 \\
            hdist=1 \\
            rcomp=t \\
            minkmerhits=1 \\
            entropy=0.35 \\
            entropywindow=50 \\
            threads=${task.cpus} \\
            overwrite=true

        cat baited_R1.fastq baited_R2.fastq | gzip -c > ${prefix}.pass2_baited_reads.fastq.gz
        READS_ARGS="--reads-1 baited_R1.fastq --reads-2 baited_R2.fastq"
    fi

    # 3. Exact bifurcation partitioning, microhomology annotation, and junction census
    cbm_bait_germline_junctions.py \\
        --candidate-anchors "${candidate_anchors}" \\
        \${READS_ARGS} \\
        --platform ${platform} \\
        --telomere-repeat ${telomere_repeat} \\
        --kmer-size ${kmer_size} \\
        --min-germ-reads ${min_germ_reads} \\
        --min-telo-reads ${min_telo_reads} \\
        --threads ${task.cpus} \\
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
    touch "${prefix}.pass2_baited_reads.fastq.gz"
    echo '{"n_breaks": 0, "n_germ_telo": 0}' > "${prefix}.coverage_and_chromosome_census.json"
    """
}
