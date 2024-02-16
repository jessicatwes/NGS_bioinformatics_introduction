#!/bin/bash
#SBATCH --job-name=bedtofasta
#SBATCH -p long
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=32gb
#SBATCH --time=00:30:00
#SBATCH --output=/Users/jewe1055/experiments/exp12/eofiles/%x%j.out
#SBATCH --error=/Users/jewe1055/experiments/exp12/eofiles/%x%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jewe1055@colorado.edu

### Displays the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Script directory is `pwd`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

### Load module
module load bedtools/2.25.0
module load python/2.7.14

### Assign variables
INDIR=/Users/jewe1055/experiments/exp12/output/macs_peak
OUTDIR=/Users/jewe1055/experiments/exp12/output/macs_peak/fasta
hg38_fasta=/scratch/Shares/dowell/genomes/hg38/hg38.fa
ROOTNAME=SRR069900

### Get fasta file for bed regions
bedtools getfasta \
 -fi ${hg38_fasta} \
 -bed ${INDIR}/${ROOTNAME}/${ROOTNAME}_peaks.narrowPeak \
 -fo ${OUTDIR}/${ROOTNAME}.fasta

date
echo Job $SLURM_JOB_ID finished at `date +"%T %a %d %b %Y"`
