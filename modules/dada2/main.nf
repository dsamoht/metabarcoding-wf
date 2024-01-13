process DADA2 {
  
    if (workflow.containerEngine == 'singularity') {
        container = params.dada2_singularity
    } else {
        container = params.dada2_docker
    }

    publishDir "${params.output}/dada2", mode: 'copy'

    input:
    path rawReads
    path figaroOutput

    output:
    path "reads_tracker_figaro_dada2.rds", emit: track
    path "seqtab_nochim_figaro_dada2.rds", emit: seqtab_nochim

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    figaro_out = readLines("${figaroOutput}")
    parameters = c("truncLen_fwd", "truncLen_rev", "ee_fwd", "ee_rev")
    for (i in seq_along(figaro_out)) {
        assign(parameters[i], as.integer(figaro_out[i]))
    }

    path = "${rawReads}"

    fnFs = sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
    fnRs = sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

    sample_names = sapply(strsplit(basename(fnFs), "_L001_"), `[`, 1)

    filtFs = file.path(".", "dada2_filtered", paste0(sample_names, "_F_filt.fastq.gz"))
    filtRs = file.path(".", "dada2_filtered", paste0(sample_names, "_R_filt.fastq.gz"))
    names(filtFs) = sample_names
    names(filtRs) = sample_names

    out = filterAndTrim(fnFs, filtFs, fnRs, filtRs,
              truncLen=c(truncLen_fwd, truncLen_rev),
              maxN=0, trimLeft=c(${params.fwd_primer_length}, ${params.rev_primer_length}),
              maxEE=c(ee_fwd, ee_rev),
              truncQ=2,
              rm.phix=TRUE,
              compress=TRUE,
              multithread=TRUE)
    head(out)

    errF = learnErrors(filtFs, multithread=TRUE)
    errR = learnErrors(filtRs, multithread=TRUE)

    dadaFs = dada(filtFs, err=errF, multithread=TRUE)
    dadaRs = dada(filtRs, err=errR, multithread=TRUE)

    mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
    seqtab = makeSequenceTable(mergers)
    seqtab_nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

    getN = function(x) sum(getUniques(x))
    track = cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab_nochim))
    
    colnames(track) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
    rownames(track) = sample_names
    head(track)

    saveRDS(track, "reads_tracker_figaro_dada2.rds")
    saveRDS(seqtab_nochim, "seqtab_nochim_figaro_dada2.rds")

    """
}