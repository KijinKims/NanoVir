nextflow.enable.dsl=2

workflow {
    main:
    Channel.fromPath(params.input).set { fastq }
    Canu(fastq)
    Translate(Canu.out())
    HMMscan(Translate.out())
    Correct(HMMscan.out())
    Virsorter2(Correct.out())
    CheckV(Virsorter2.out())
}

process Canu {
    input:
        path fastq
    output:
        path "canu/${params.prefix}.contigs.fasta"
        path graphs
    """
    canu $fastq -d canu -p ${params.prefix} ${params.canu_high_sens_opts}
    """

}

process Translate {
    input:
        path contigs
        path graphs
    output:
        path contigs
        path "${contigs.baseName}.pep.fasta"
        path graphs
    """
    python lab_scripts/dna2pep.py -r all --fasta ${contigs.baseName}.pep.fasta $contigs
    """
}

process HMMscan {
    input:
        path contigs
        path translated_contigs
        path graphs
    output:
        path contigs
        path graphs
        path "${params.prefix}_${params.hmmdb}.domtbl"
    """
    hmmscan $translated_contigs $graphs ${params.hmmdb} --domtbl ${params.prefix}_${params.hmmdb}.domtbl
    """
}

process Correct {
    input:
        path contigs
        path graphs
        path hmmscan_result    
    output:
        path corrected_contigs
    """
    python correct.py --contigs $contigs --DAG $graphs --hmmscan_result $hmmscan_result --DAG_coef ${params.DAG_coef} --HMM_coef ${params.HMM_coef}
    """
}

process Virsorter2 {
    input:
        path contigs
    output:
        path viral_contigs
    """
    virsorter run -d ${params.virsorterdb} -i $viral_contigs -j ${params.threads} all
    """
}

process Checkv {
    input:
        path viral_contigs
    output:
        path result
    """
    checkv end_to_end -d ${params.checkvdb} $viral_contigs -t ${params.threads}
    """
}