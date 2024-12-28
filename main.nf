#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process TRIMMOMATIC {
    container "biocontainers/trimmomatic:0.39--hdfd78af_2"

    publishDir "${params.outdir}/trimmomatic/", mode: 'copy'

    input:
    path(reads)

    output:
    path("*.paired.trim*.fastq.gz")   , emit: trimmed_reads
    path("*.unpaired.trim_*.fastq.gz"), emit: unpaired_reads, optional:true
    path("*_trim.log")                , emit: trim_log
    path("*_out.log")                 , emit: out_log
    path("*.summary")                 , emit: summary
    path "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trimmed = meta.single_end ? "SE" : "PE"
    def output = meta.single_end ?
        "${prefix}.SE.paired.trim.fastq.gz" // HACK to avoid unpaired and paired in the trimmed_reads output
        : "${prefix}.paired.trim_1.fastq.gz ${prefix}.unpaired.trim_1.fastq.gz ${prefix}.paired.trim_2.fastq.gz ${prefix}.unpaired.trim_2.fastq.gz"
    def qual_trim = task.ext.args2 ?: ''
    """
    trimmomatic \\
        $trimmed \\
        -threads $task.cpus \\
        -trimlog ${prefix}_trim.log \\
        -summary ${prefix}.summary \\
        $reads \\
        $output \\
        $qual_trim \\
        $args 2> >(tee ${prefix}_out.log >&2)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$(trimmomatic -version)
    END_VERSIONS
    """
}

process METASPADES {
    container "biocontainers/spades:4.0.0--h5fb382e_1"

    publishDir "${params.outdir}/metaspades/", mode: 'copy'

    input:
    path(reads)

    output:
    path('*.scaffolds.fa.gz')    , optional:true, emit: scaffolds
    path('*.contigs.fa.gz')      , optional:true, emit: contigs
    path('*.transcripts.fa.gz')  , optional:true, emit: transcripts
    path('*.gene_clusters.fa.gz'), optional:true, emit: gene_clusters
    path('*.assembly.gfa.gz')    , optional:true, emit: gfa
    path('*.warnings.log')       , optional:true, emit: warnings
    path('*.spades.log')         , emit: log
    path  "versions.yml"         , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maxmem = task.memory.toGiga()
    def illumina_reads = illumina ? ( meta.single_end ? "-s $illumina" : "-1 ${illumina[0]} -2 ${illumina[1]}" ) : ""
    """
    spades.py \\
        $args \\
        --threads $task.cpus \\
        --memory $maxmem \\
        $custom_hmms \\
        $reads \\
        -o ./
    mv spades.log ${prefix}.spades.log

    if [ -f scaffolds.fasta ]; then
        mv scaffolds.fasta ${prefix}.scaffolds.fa
        gzip -n ${prefix}.scaffolds.fa
    fi
    if [ -f contigs.fasta ]; then
        mv contigs.fasta ${prefix}.contigs.fa
        gzip -n ${prefix}.contigs.fa
    fi
    if [ -f transcripts.fasta ]; then
        mv transcripts.fasta ${prefix}.transcripts.fa
        gzip -n ${prefix}.transcripts.fa
    fi
    if [ -f assembly_graph_with_scaffolds.gfa ]; then
        mv assembly_graph_with_scaffolds.gfa ${prefix}.assembly.gfa
        gzip -n ${prefix}.assembly.gfa
    fi

    if [ -f gene_clusters.fasta ]; then
        mv gene_clusters.fasta ${prefix}.gene_clusters.fa
        gzip -n ${prefix}.gene_clusters.fa
    fi

    if [ -f warnings.log ]; then
        mv warnings.log ${prefix}.warnings.log
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed -n 's/^.*SPAdes genome assembler v//p')
    END_VERSIONS
    """
}

workflow {
    reads_ch    = Channel.fromPath(params.reads)

    TRIMMOMATIC(reads_ch)
    METASPADES(TRIMMOMATIC.out.trimmed_reads)
}