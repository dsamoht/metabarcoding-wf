# How to use the script
1) Navigate to the directory containing the .ab1 files
2) Copy the script `paired_reads_analysis.py` in the directory
3) Launch the following script

```
#!/bin/bash
#SBATCH --job-name=16s_sanger_analysis_job
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=9
#SBATCH --account=def-yourpi
#SBATCH --out=16s_sanger_analysis_job.stdout
#SBATCH --error=16s_sanger_analysis_job.stderr

module load StdEnv/2020 gcc/9.3.0 blast+/2.13.0 emboss fastqc python
pip install biopython
pip install pandas
git clone https://github.com/lh3/seqtk.git
cd seqtk; make; cd ..                                          
export PATH=$PATH:$PWD/seqtk

python paired_reads_analysis.py
```
