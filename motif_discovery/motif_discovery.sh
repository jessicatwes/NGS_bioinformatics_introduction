#!/bin/bash
#SBATCH --job-name=motif_discovery
#SBATCH -p long
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=32gb
#SBATCH --time=04:00:00
#SBATCH --output=/Users/jewe1055/experiments/exp12/eofiles/%x%j.out
#SBATCH --error=/Users/jewe1055/experiments/exp12/eofiles/%x%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jewe1055@colorado.edu

### Displays the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Script directory is `pwd`
echo Tfit bedfile directory is ${INDIR}
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

### Load Modules
module load meme/5.1.1

### Set Parameters
INDIR=/Users/jewe1055/experiments/exp12/output/macs_peak/fasta
OUTDIR=/Users/jewe1055/experiments/exp12/output/macs_peak/fasta
ROOTNAME=SRR069904

### Run MEME
/opt/meme/5.1.1/bin/meme \
 ${INDIR}/${ROOTNAME}.fasta \
 -dna \
 -oc ${OUTDIR}/meme/${ROOTNAME}

### Run DREME
#/opt/meme/4.12.0/bin/dreme \
# -p ${INDIR}/${ROOTNAME}.fasta \
# -oc ${OUTDIR}/dreme/${ROOTNAME}


### Run TOMTOM
/opt/meme/5.1.1/bin/tomtom \
 ${OUTDIR}/meme/${ROOTNAME}/meme.html \
 -no-ssc \
 -oc ${OUTDIR}/tomtom/${ROOTNAME} \
 -verbosity 1 \
 -min-overlap 5 \
 -dist pearson \
 -evalue \
 -thresh 10.0 \
 /opt/meme/4.12.0/db/motif_databases/JASPAR/JASPAR_CORE_2014_vertebrates.meme \
 /opt/meme/4.12.0/db/motif_databases/MOUSE/uniprobe_mouse.meme \
 /opt/meme/4.12.0/db/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme
