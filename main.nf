#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process TRIMMOMATIC {
    container 'community.wave.seqera.io/library/trimmomatic:0.39--a688969e471089d7'

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
    def prefix = reads[0].name.replace('_R1.fastq.gz', '')
    def qual_trim = task.ext.args2 ?: ''
    """
    trimmomatic PE \\
        -threads $task.cpus \\
        $reads \\
        ${prefix}.paired.trim_1.fastq.gz ${prefix}.unpaired.trim_1.fastq.gz ${prefix}.paired.trim_2.fastq.gz ${prefix}.unpaired.trim_2.fastq.gz \\
        -trimlog ${prefix}_trim.log \\
        -summary ${prefix}.summary \\
        LEADING:30 TRAILING:30 SLIDINGWINDOW:4:20 MINLEN:35 \\
        2> >(tee ${prefix}_out.log >&2)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$(trimmomatic -version)
    END_VERSIONS
    """
}

process METASPADES {
    container 'community.wave.seqera.io/library/spades:4.0.0--dc56d3b41f13769d'

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
    def prefix = reads[0].name.replace('.paired.trim_1.fastq.gz', '')
    def maxmem = task.memory.toGiga()
    def reads = "-1 ${reads[0]} -2 ${reads[1]}"
    """
    spades.py \\
        --meta \\
        --threads $task.cpus \\
        --memory $maxmem \\
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
    reads_ch = Channel.fromPath(params.reads)

    TRIMMOMATIC(reads_ch.collect())
    METASPADES(TRIMMOMATIC.out.trimmed_reads)
}
