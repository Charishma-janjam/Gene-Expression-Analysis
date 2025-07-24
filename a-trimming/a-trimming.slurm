#!/bin/bash
#SBATCH --job-name=trimming
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16

### SCRIPT TRIMS RAW SEQUENCE DATA
### SOFTWARE USED: TrimGalore!, cutadapt, FastQC

##JOB LOG HEADER - Puts a copy of this script &  other debugging info in "slurm.out"
  perl -E 'say"="x80'; echo "JOB STARTED: $(date)"; echo "NODE: $(hostname)"; echo "SCRIPT ${0}:"; cat $0; perl -E 'say"="x80'

## SOFTWARE REQUIREMENTS
  export PATH=$PATH:/usr/local/FastQC

## PRIMARY VARIABLES - commonly need to change                                                    # TIP: ${HOME} inserts your home directory
  FASTQ_IN="/work/binf-courses/binf694/project_1/raw_data/placenta/placenta_6b.s.fastq.gz"     # FASTQ FILE LOCATION
  OUTDIR="${HOME}/exercises/gene_expression/a-trimming"                                               # OUTPUT DIRECTORY

## SECONDARY VARIABLES - rarely need to change
  MIN_LEN=50            # MINIMUM READ LENGTH TO KEEP
  MIN_QUAL=28           # MIN QUAL SCORE (PHRED) TO TRIM
  CPU=4                 # CPU CORES - TG oddly uses 4x this # - make sure #SBATCH requests enough (e.g. 16 for 4)

## COMMAND(s) TO RUN
  mkdir -p ${OUTDIR}    # MAKE SURE OUTPUT DIRECTORY EXISTS
  cd ${OUTDIR}          # CHANGE TO OUTPUT DIRECTORY

## Use _TrimGalore!_ to remove low quality sequences and adapters using _cutadapt_
 # It will also rerun _FastQC_ so you can see how the trimming worked...
    trim_galore \
      --length ${MIN_LEN} \
      --quality ${MIN_QUAL} \
      --cores ${CPU} \
      --fastqc \
      --fastqc_args "--threads ${SLURM_CPUS_PER_TASK}" \
      ${FASTQ_IN}

## JOB CLEANUP
  perl -E 'say"="x80'; echo "JOB COMPLETED: $(date)"; perl -E 'say"="x80'
