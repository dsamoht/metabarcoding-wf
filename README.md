# Automation of Amplicon Data Treatment for Microbial Ecology
## Introduction
A Nextflow pipeline automating metabarcoding data treatment from __raw__ Illumina paired reads to __ASVs__
> *__Note__* : Input reads must have these characteristics:
> - Demultiplexed
> - Typical Illumina file name (i.e. \*_L001_R{1,2}_001.fastq.gz)  
> - Reads of the same orientation (forward / reverse) must all have the same length

## Dependencies
- __Software :__  
  [Nextflow](https://www.nextflow.io/)  
  [Docker](https://www.docker.com/) and/or [Apptainer/Singularity](https://apptainer.org/)  

- (TODO) __A DADA2-formatted taxonomy database :__  
[SILVA]()  
[UNITE]()  
[PR2]()  
...TODO
- __Edit__ *nextflow.config* :  
  ```
  TODO    
  ```
## How to run the pipeline
__This command will test the setup and download the containers for off-line use__:  
```
nextflow run metabarcoding-wf.nf -profile {docker,singularity},local,test
```
__Run on your data__:  
```
nextflow run metabarcoding-wf.nf -profile {docker,singularity},{local,hpc} --input /path/to/fastq/files/ --output /path/to/output/ --amplicon_length [AMPLICONLENGTH] --fwd_primer_length [FWDPRIMERLENGTH] --rev_primer_length [REVPRIMERLENGTH]
```
