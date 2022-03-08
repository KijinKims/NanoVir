nextflow.enable.dsl=2

workflow {
    main:
    Channel.fromPath(params.input).set { fastq }
    Channel.fromPath(params.hmmdb).set { hmmdb }
    Canu(fastq)
    Translate(Canu.out.contigs)
    HMMscan(Translate.out.translated_contigs, hmmdb)
    Correct(Canu.out.contigs, Canu.out.graph, hmmdb, HMMscan.out.domtbl)
    Virsorter2(Correct.out.corrected_contigs)
    CheckV(Virsorter2.out.viral_contigs)
}

process Canu {
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path fastq
    output:
        path "canu.${params.prefix}/*"
        path "canu.${params.prefix}/${params.prefix}.contigs.fasta", emit: contigs
        path "canu.${params.prefix}/${params.prefix}.graph.dot", emit: graph
    """
    ${params.canu_binary_path} -p ${params.prefix} -d canu.${params.prefix} -nanopore-raw $fastq ${params.canu_high_sens_opts}
    """

}

process Translate {
    input:
        path contigs
    output:
        path "${params.prefix}.pep.fasta", emit: translated_contigs
    """
    python ${params.nanovir_dir}/modules/dna2pep.py -i $contigs -o ${params.prefix}.pep.fasta
    """
}

process HMMscan {
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path translated_contigs
        path hmmdb
    output:
        path "hmmscan.${params.prefix}/${params.prefix}.domtbl", emit: domtbl
    """
    mkdir hmmscan.${params.prefix}
    hmmscan --domtblout hmmscan.${params.prefix}/${params.prefix}.domtbl -E ${params.hmmscan_evalue} $hmmdb $translated_contigs 
    """
}

process Correct {
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path contigs
        path graphs
        path hmmdb
        path hmmscan_result    
    output:
        path "correct.${params.prefix}/${params.prefix}.corrected_contigs.fasta", emit: corrected_contigs
        path "correct.${params.prefix}/${params.prefix}.corrected_graphs.dot", emit: corrected_graphs
    """
    python ${params.nanovir_dir}/modules/correct.py --contigs $contigs --DAG $graphs --hmmdb $hmmdb --hmmscan_result $hmmscan_result --minimum_edge_weight ${params.min_weight} -o correct.${params.prefix}/${params.prefix}.corrected_contigs.fasta -og correct.${params.prefix}/${params.prefix}.corrected_graphs.dot
    """
}

process Virsorter2 {
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path contigs
    output:
        path "virsorter.${params.prefix}/*"
        path "virsorter.${params.prefix}/final-viral-combined.fa", emit: viral_contigs
    """
    virsorter run -d ${params.virsorterdb} -i $viral_contigs -j ${params.threads} ${params.virsorter_viruses} -w virsorter.${params.prefix}
    """
}

process CheckV {
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path viral_contigs
    output:
        path "checkv.${params.prefix}/*"
    """
    checkv ${params.checkv_mode} -d ${params.checkvdb} $viral_contigs checkv.${params.prefix} -t ${params.threads}
    """
}