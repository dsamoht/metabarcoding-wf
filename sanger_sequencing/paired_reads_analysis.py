"""Sanger Analysis - 16S"""
import glob
import os
import sys
from pathlib import Path
import subprocess
from datetime import datetime
import logging

from Bio import SeqIO
import pandas as pd

__author__ = "Thomas Deschênes"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("16s_sanger_analysis_{:%Y_%m_%d}.log".format(datetime.now()), mode="w"),
        logging.StreamHandler()
    ]
)
INFO = """
================== Sanger Analysis - 16S =============================

This script is made to be executed on the Alliance clusters.
Execute to following commands before executing the script.
You can ignore this message if it's already done.

Commands to execute (copy-paste the following lines in your terminal):

module load StdEnv/2020  gcc/9.3.0 blast+/2.13.0 emboss fastqc python
pip install biopython
pip install pandas
git clone https://github.com/lh3/seqtk.git
cd seqtk; make; cd ..                                          
export PATH=$PATH:$PWD/seqtk                                                      


What this script does:

1) Converts .ab1 files to .fastq files;
2) Runs fastqc on forward and reverse files to log quality statistics; 
3) Trims the 5' and 3' ends to remove low quality bases
4) Merges forward and reverse reads into a consensus sequence;
5) Runs blast locally for taxonomy assignment;

To run the script:

Be in the directory where you have your .ab1 files.
Have a copy of this script in this directory.

python analysis_16s_sanger.py

(CTRL+C to quit)
======================================================================
"""

print(INFO)

def quality_assessment_warning(fastq_report):
    """
    Return the mean quality score as reported
    in `fastqc_data.txt`
    """
    lines = []
    with open(fastq_report, "r", encoding="UTF-8") as file:
        found_start = False
        for line in file:
            if found_start:
                if line.startswith(">>END_MODULE"):
                    break
                lines.append(line.strip())
            elif line.startswith(">>Per base sequence quality"):
                found_start = True

    quality_df = pd.DataFrame([line.split("\t") for line in lines],
                              columns=["Base", "Mean", "Median", "Lower Quartile", "Upper Quartile", "10th Percentile", "90th Percentile"]).drop(0, axis=0)
    quality_df["Mean"] = pd.to_numeric(quality_df["Mean"])
    return quality_df["Mean"].mean()

FWD_FILE_ID = "For"
REV_FILE_ID = "Rev"

TRIM_LEFT = "40"
TRIM_RIGHT = "100"

logging.info("Searching for FWD files with expression 'For'")
logging.info("Searching for REV files with expression 'Rev'")

fwd_ab1_files = glob.glob(f"*{FWD_FILE_ID}*ab1")
rev_ab1_files = glob.glob(f"*{REV_FILE_ID}*ab1")

logging.info(f"Found {len(fwd_ab1_files)} FWD .ab1 files")
logging.info(f"Found {len(rev_ab1_files)} REV .ab1 files")

logging.info("Converting .ab1 files to fastq")
Path("./fastq/").mkdir(parents=True, exist_ok=True)
for _file in sorted(fwd_ab1_files+rev_ab1_files):
    logging.info(f"Converting file `{_file}` from .ab1 to .fastq")
    if Path(f"./fastq/{_file.split('ab1')[0]}fastq").is_file():
        continue
    else:
        try:
            rec = SeqIO.parse(_file, "abi")
            count = SeqIO.write(rec, Path(f"./fastq/{_file.split('ab1')[0]}fastq"), "fastq")
        except Exception as err:
            logging.error(f"A problem occured when trying to convert file `{_file}`")
            sys.exit()

logging.info("Done. Files are in `./fastq/` directory")

fwd_fq_files = glob.glob(f"./fastq/*{FWD_FILE_ID}*fastq")
rev_fq_files = glob.glob(f"./fastq/*{REV_FILE_ID}*fastq")

logging.info(f"{len(fwd_fq_files)} FWD .fastq files were produced")
logging.info(f"{len(rev_fq_files)} REV .fastq files were produced")

logging.info("Running `fastqc` on .fastq files")
Path("./fastqc/").mkdir(parents=True, exist_ok=True)
for _file in sorted(fwd_fq_files+rev_fq_files):
    logging.info(f"Running `fastqc` on file `{_file}`")
    if (Path(f"./fastqc/{Path(_file).name.split('.fastq')[0]}_fastqc/fastqc_data.txt").is_file() and
        Path(f"./fastqc/{Path(_file).name.split('.fastq')[0]}_fastqc.html").is_file()):
        pass
    else:
        try:
            subprocess.run(["fastqc", _file, "--extract", 
                            "--outdir", f"{Path(__file__).parent.joinpath('./fastqc/')}"],
                            check=False,
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL)
        except Exception as err:
            logging.error(f"A problem occured when trying to run `fastqc` on file `{_file}`")
            sys.exit()

    try:
        mean_quality = quality_assessment_warning(f"./fastqc/{Path(_file).name.split('.fastq')[0]}_fastqc/fastqc_data.txt")
        if mean_quality < 20:
            logging.warning(f"Sequence of `{_file}` has an average quality score of {round(mean_quality, 2)} (< 20). Inspect the `fastqc` report for more details.")
    except Exception as err:
        logging.error(f"A problem occured when trying to read the `fastqc` report of `{_file}`: {err}")
        sys.exit()

logging.info("Done. Profiles are in `./fastqc/` directory")
logging.info("Running `seqtk` on .fastq files with default parameters")
logging.info(f"Trimming the left end of sequences at {TRIM_LEFT} bp and right end at {TRIM_RIGHT} bp")

Path("./fastq/trimmed/").mkdir(parents=True, exist_ok=True)
for _file in sorted(fwd_fq_files+rev_fq_files):
    file_name = Path(_file).name.split(".fastq")[0] + "_trimmed.fastq"
    logging.info(f"Trimming file `{_file}` with `seqtk`")
    if Path("./fastq/trimmed/"+ file_name).is_file():
        continue
    else:
        try:
            with open(f"{Path(__file__).parent.joinpath('./fastq/trimmed/'+ file_name)}", "w") as trimmed_output:
                subprocess.run(["seqtk", "trimfq", "-b", TRIM_LEFT, "-e", TRIM_RIGHT, _file],
                                stdout=trimmed_output,
                                check=False)
        except Exception as err:
            logging.error(f"A problem occured when trying to run `seqtk` on file {_file}")
            sys.exit()

logging.info("Done. Files are in `./fastqc/trimmed/` directory")
logging.info("Merging the FWD and the reverse complement of the REV files")

fwd_trimmed_files = glob.glob(f"./fastq/trimmed/*{FWD_FILE_ID}*trimmed*fastq")
sample_names = [Path(i.split(f"_{FWD_FILE_ID}_")[0]).name for i in fwd_trimmed_files]

if not Path("./Merge_Sanger_v2.py").exists():
    logging.info("Downloading `Merge_Sanger_v2.py` from https://peerj.com/articles/11354/")
    try:
        subprocess.run(["wget", "https://dfzljdn9uc3pi.cloudfront.net/2021/11354/1/Merge_Sanger_v2.py"],
                       check=True)
    except Exception as err:
        logging.error("A problem occured when downloading `https://dfzljdn9uc3pi.cloudfront.net/2021/11354/1/Merge_Sanger_v2.py`")
        sys.exit()

rev_trimmed_files = glob.glob(f"./fastq/trimmed/*{REV_FILE_ID}*trimmed*fastq")

logging.info("Executing `Merge_Sanger_v2.py` on all samples")
with open("all_merged_sequences.fasta", "w"):
    pass

for sample_name in sorted(sample_names):
    logging.info(f"Executing `Merge_Sanger_v2.py` on sample `{sample_name}`")
    try:
        Path(f"./{sample_name}_merger/").mkdir(parents=True, exist_ok=True)
        fwd_file = glob.glob(f"./fastq/trimmed/{sample_name}*{FWD_FILE_ID}*trimmed*fastq")[-1]
        rev_file = glob.glob(f"./fastq/trimmed/{sample_name}*{REV_FILE_ID}*trimmed*fastq")[-1]
        for rec_fwd in SeqIO.parse(fwd_file, "fastq"):
            if len(rec_fwd.seq) < 20:
                logging.warning(f"`{fwd_file}` has a length of {len(rec_fwd.seq)}. Skipping this file")
            else:
                with open(f"./{sample_name}_merger/000F-{sample_name}.seq", "w", encoding="UTF-8") as fastq_fwd_output:
                    fastq_fwd_output.write(str(rec_fwd.seq))
        for rec_rev in SeqIO.parse(rev_file, "fastq"):
            if len(rec_rev.seq) < 20:
                logging.warning(f"`{rev_file}` has a length of {len(rec_rev.seq)}. Skipping this file")
            else:
                with open(f"./{sample_name}_merger/001R-{sample_name}.seq", "w", encoding="UTF-8") as fastq_rev_output:
                    fastq_rev_output.write(str(rec_rev.seq))
    
        os.chdir(f"./{sample_name}_merger/")
        with open(f"{sample_name}_merger.log", "w") as log:
            subprocess.run(["python", "../Merge_Sanger_v2.py", f"000F-{sample_name}.seq", f"001R-{sample_name}.seq"],
                            check=True,
                            stdout=log,
                            stderr=subprocess.DEVNULL)
        seq = ""
        with open("./merged_sequence/merged.seq", "r") as seq_input:
            for _line in seq_input:
                seq += _line.strip()
            with open("../all_merged_sequences.fasta", "a") as fasta_output:
                fasta_output.write(f">{sample_name}\n{seq}\n")
        os.chdir("..")
    except Exception as err:
        logging.error(f"A problem occured when executing the script `Merge_Sanger_v2.py` on the sample `{sample_name}`. Error: '{err}'")
        os.chdir("..")
        #sys.exit()

logging.info("Done. All merged sequences are in the file `all_merged_sequences.fasta`")
logging.info("Running blastn on nt database... ")

try:
    subprocess.run(["blastn", "-task", "megablast", "-query", "all_merged_sequences.fasta", "-db", "/cvmfs/ref.mugqic/genomes/blast_db/LATEST/nt", "-out", "all_seq_blast_output.txt", "-outfmt",
                    "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore", "-num_threads", "4", "-max_target_seqs", "10"])
except Exception as err:
    logging.error(f"A problem occured when running blastn on nt database with the query file `all_merged_sequences.fasta`")
    sys.exit()

logging.info("Done. Results are reported in the file `all_seq_blast_output.txt`")
logging.info("Run completed!")
