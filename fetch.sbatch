#!/bin/bash
#SBATCH --job-name=fetchSRA
#SBATCH -p long
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=32gb
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/Users/jewe1055/exp_files/exp30/eofiles/%x%j.out
#SBATCH --error=/scratch/Users/jewe1055/exp_files/exp30/eofiles/%x%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jewe1055@colorado.edu

### Displays the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Script directory is `pwd`
echo Tfit bedfile directory is ${INDIR}
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

### Load module
module load sra/2.8.0

### Parameters
OUTDIR=/scratch/Users/jewe1055/exp_files/exp12/input

### Fetch sra
# SRR from ChIP-seq in K562 anti-RUNX1 (Pencovich et al 2011)
#fastq-dump -O ${OUTDIR} -split-3 SRR069898 #Input control for K562-RUNX1
#fastq-dump -O ${OUTDIR} -split-3 SRR069899 #Input control for K562-RUNX1
#fastq-dump -O ${OUTDIR} -split-3 SRR069902 #Input control for K562-TPA-RUNX1
#fastq-dump -O ${OUTDIR} -split-3 SRR069903 #Input control for K562-TPA-RUNX1

# SRR from ChIP-seq in Jurkat anti-RUNX1 (Hollenhorst et al 2009)
fastq-dump -O ${OUTDIR} -split-3 SRR034158 #Input control for Jurkat
fastq-dump -O ${OUTDIR} -split-3 SRR034159 #Input control for Jurkat
fastq-dump -O ${OUTDIR} -split-3 SRR034154 #anti-RUNX1 for Jurkat-RUNX1
fastq-dump -O ${OUTDIR} -split-3 SRR034155 #anti-RUNX1 for Jurkat-RUNX1

date
echo Job $SLURM_JOB_ID finished at `date +"%T %a %d %b %Y"`
