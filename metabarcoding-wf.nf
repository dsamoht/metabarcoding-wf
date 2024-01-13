include { FIGARO } from './modules/figaro'
include { DADA2  } from './modules/dada2'

workflow METABARCODING_WF {

    raw_fastq_dir = Channel.fromPath(params.input, type: 'dir')
    FIGARO(raw_fastq_dir)
    DADA2(raw_fastq_dir, FIGARO.out)

}

workflow {

    METABARCODING_WF()

}