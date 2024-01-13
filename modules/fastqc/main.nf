process FASTQC {

    if (workflow.containerEngine == 'singularity') {
        container = params.fastqc_singularity
    } else {
        container = params.fastqc_docker
    }

    publishDir "${params.output}/fastqc", mode: 'copy'

    input:
    path inputDir

    output:
    path "R1_combined_fastqc.html", emit: fastqc_fwd_report
    path "R2_combined_fastqc.html", emit: fastqc_rev_report

    script:
    """
    cat ${inputDir}/*R1_001.fastq.gz > R1_combined.fastq.gz
    cat ${inputDir}/*R2_001.fastq.gz > R2_combined.fastq.gz
    fastqc R*_combined.fastq.gz
    """

}