#!/bin/bash
#SBATCH --job-name=macs
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jewe1055@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --time=10:00:00
#SBATCH --output=/Users/jewe1055/experiments/exp12/eofiles/%x%j.out
#SBATCH --error=/Users/jewe1055/experiments/exp12/eofiles/%x%j.err

#### Display job content in output file
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo "Running on host `hostname`"
echo "Job started at `date +"%T %a %d %b %Y"`"
echo Directory is `pwd`
echo Using $SLURM_NTASKS processors across $SLURM_NODES nodes

#### Set parameters
INDIR=/Users/jewe1055/experiments/exp12/output/mapped/bams
OUTDIR=/Users/jewe1055/experiments/exp12/output/macs_peak
ROOTNAME=SRR034154
CONTROL=SRR034158

#### Load module
module load python/2.7.14/MACS/2.1.1

mkdir -p ${OUTDIR}
mkdir -p ${OUTDIR}/${ROOTNAME}

#### Call peaks without controls
#macs2 callpeak -t ${INDIR}${PROJECT}/${ROOTNAME}_trimmed.sorted.bam.tri.BedGraph -f BED -g hs -n ${PROJECT}_${ROOTNAME} --nomodel --shift -100 --extsize 200 -B --broad --outdir ${OUTDIR}${PROJECT}/
#macs2 callpeak -t ${INDIR}/${ROOTNAME}.sorted.bam -f BAM -g hs -n ${ROOTNAME} --nomodel --shift -100 --extsize 200 -B -q 0.01 --outdir ${OUTDIR}

#### Call peaks with control
#macs2 callpeak -t ChIP.bam -c Control.bam -f BAM -g hs -n name -B -q 0.01 for regular peaks
macs2 callpeak -t ${INDIR}/${ROOTNAME}.sorted.bam -c ${INDIR}/${CONTROL}.sorted.bam -f BAM -g  hs -n ${ROOTNAME} --nomodel --shift -100 --extsize 200 -B -q 0.01 --outdir ${OUTDIR}/${ROOTNAME}
