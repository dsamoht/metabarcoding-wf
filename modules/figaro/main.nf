process FIGARO {

    if (workflow.containerEngine == 'singularity') {
        container = params.figaro_singularity
    } else {
        container = params.figaro_docker
    }

    publishDir "${params.outdir}/figaro", mode: 'copy'

    input:
    path inputDir

    output:
    path "filterAndTrimParameters.txt", emit: figaro_out

    script:
    def python_cmd = """import json
    with open('${inputDir}_figaro_out/trimParameters.json') as json_file:
        data = json.load(json_file)
    print(data[0]['trimPosition'][0])
    print(data[0]['trimPosition'][1])
    print(data[0]['maxExpectedError'][0])
    print(data[0]['maxExpectedError'][1])
    """
    """
    figaro.py -a ${params.amplicon_length} -f ${params.fwd_primer_length} -r ${params.rev_primer_length} -i ${inputDir} -o ${inputDir}_figaro_out
    python -c "${python_cmd}" > filterAndTrimParameters.txt
    """

}