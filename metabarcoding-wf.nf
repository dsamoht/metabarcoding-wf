include { FIGARO } from './modules/figaro'
include { DADA2  } from './modules/dada2'

workflow METABARCODING_WF {

    if (params.input == ""){
        exit 1, "Enter a valid fastq directory (--input /path/to/fastq/)."
    }
    if (params.output == ""){
        exit 1, "Enter a valid output directory (--output /path/to/output/)."
    }
    if (params.amplicon_length == ""){
        exit 1, "Enter a valid amplicon length (--amplicon_length [AMPLICONLENGTH])."
    }
    if (params.fwd_primer_length == ""){
        exit 1, "Enter a valid forward primer length (--fwd_primer_length [FWDPRIMERLENGTH])."
    }
    if (params.rev_primer_length == ""){
        exit 1, "Enter a valid reverse primer length (--rev_primer_length [REVPRIMERLENGTH])."
    }

    raw_fastq_dir = Channel.fromPath(params.input, type: 'dir')
    FIGARO(raw_fastq_dir)
    DADA2(raw_fastq_dir, FIGARO.out)

}

workflow {

    METABARCODING_WF()

}