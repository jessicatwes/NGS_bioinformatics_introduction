#!/bin/bash
#SBATCH --job-name=SS_map 				# Job name
#SBATCH --mail-type=FAIL,END 				# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=margaret.gruca@colorado.edu 	# Where to send mail
#SBATCH --nodes=1					# Number of nodes
#SBATCH --ntasks=16 					# Number of CPU (processer cores i.e. tasks)
#SBATCH --time=08:00:00 				# Time limit hrs:min:sec
#SBATCH --partition short				# Job queue
#SBATCH --mem=100gb 					# Memory limit
#SBATCH --output=/scratch/Users/magr0763/Andrysik2017/RNA-seq/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/magr0763/Andrysik2017/RNA-seq/e_and_o/%x_%j.err

#This script will take you from a raw, unedited fastq file to tdf

ROOTNAME=$1
PROJECT='/scratch/Users/magr0763/Andrysik2017/RNA-seq'
GENOME='/scratch/Users/magr0763/genomes/mapping/hg38/'
DEEP_CONTAINER='singularity exec -H /scratch/Users/magr0763 /scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img'

pwd; hostname; date
date


#Quality check of fastqc file before trimming

module load fastqc/0.11.5
fastqc ${PROJECT}/fastq/${ROOTNAME}.fastq -o ${PROJECT}/fastq/FastQC/
echo pre-trim qual
date
date

###Trimming: there are a number of different tools available (Trimmomatic, cudadapt, TrimGalore). My favorite currently for most purposes in bbduk due to its rapid runtime and large number of trimming options. "literal" is commented out below. The current sequence I have there is for rrna contaminants, however you can specify any sequence (similar to -a flag in cutadapt). This is useful if you have to add an additional sequence/contaminant to trim rather than adding it to the fa file. Because of the large number of options, bbduk is particularly useful if you really know your data. For a first pass if you are not sure of the protocol, you will more than likely miss adapters as it has no autodetect feature.
###Go here for full list of options : https://github.com/BioInfoTools/BBMap/blob/master/sh/bbduk2.sh

module load bbmap/38.05

bbduk.sh -Xmx20g \
t=16 \
overwrite= t \
in=${PROJECT}/fastq/${ROOTNAME}.fastq \
out=${PROJECT}/trimmed/${ROOTNAME}.trim.fastq \
ref=/opt/bbmap/38.05/adapters.fa \
ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
maq=10 minlen=20 \
stats=${PROJECT}/trimmed/${ROOTNAME}.trimstats.txt \
refstats=${PROJECT}/trimmed/${ROOTNAME}.refstats.txt \
ehist=${PROJECT}/trimmed/${ROOTNAME}.ehist.txt

echo bbduk trim
date
date

#Quality check trimmed reads --> if first QC bad, check these before moving forward and adjust parameters if needed

module load fastqc/0.11.5
fastqc ${PROJECT}/trimmed/${ROOTNAME}.trim.fastq -o ${PROJECT}/trimmed/FastQC/
echo post-trim qual
date
date

#Map trimmed fastq files to reference genome --> sam files include all reads

module load samtools/1.3.1
module load hisat2/2.1.0

hisat2 -p 16 \
 --very-sensitive \
 -x ${GENOME}/HISAT2_indexes/genome \
 --pen-noncansplice 14 \
 --mp 1,0 \
 -U ${PROJECT}/trimmed/${ROOTNAME}.trim.fastq \
 > ${PROJECT}/mapped/sams/${ROOTNAME}.trim.sam \
 2> ${PROJECT}/mapped/sams/${ROOTNAME}.trim.stderr

##### These are the flags that will singificantly affect mapping percentages -- adjust as necessary

# --pen-noncansplice 14        
# --mp 1,0 
# --sp 3,1 \

echo mapped, sam
date
date

#Word count sam files, convert sam files to bam files (compressed, binary sam files)

wc \
 -l ${PROJECT}/mapped/sams/${ROOTNAME}.trim.sam \
 > ${PROJECT}/mapped/sams/${ROOTNAME}.trim.sam.wc

samtools view \
 -S -b -o \
 ${PROJECT}/mapped/bams/${ROOTNAME}.trim.bam \
 ${PROJECT}/mapped/sams/${ROOTNAME}.trim.sam \
 2> ${PROJECT}/mapped/bams/${ROOTNAME}.trim.bam.err

samtools flagstat \
 ${PROJECT}/mapped/bams/${ROOTNAME}.trim.bam \
 > ${PROJECT}/mapped/bams/${ROOTNAME}.trim.bam.flagstat \
 2> ${PROJECT}/mapped/bams/${ROOTNAME}.trim.bam.flagstat.err

echo bam
date
date

#Sort bam files, flagstat useful for rcc (read count correction) and QC

samtools sort \
 -m 100G ${PROJECT}/mapped/bams/${ROOTNAME}.trim.bam \
 > ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam

samtools flagstat \
 ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam \
 > ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam.flagstat \
 2> ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam.flagstat.err

echo sorted.bam
date
date

#Index sorted bam files for use in mutlicov OR if you decide to sort out a specific region (e.g. use samtools view to select only chromosome 1)

samtools index \
 ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam \
 ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam.bai

echo indexed
date
date

module load preseq/2.0.3

#c_curve : can use inputs from file types [bed/bam] and will plot the estimated complexity of a sample. Specify -B for sorted bam files, -P for paired end read files

preseq c_curve -B \
 -o ${PROJECT}/qc/preseq/${ROOTNAME}.c_curve.txt \
 ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam

echo c_curve
echo complexity
date

#lc_extrap : can use inputs from from tyles [bed/bam] and will estimate future yields for complexity if the sample is sequenced at higher read depths. Specify -B for sorted bam files, -P for paired end reads

preseq lc_extrap -B \
 -o ${PROJECT}/qc/preseq/${ROOTNAME}.lc_extrap.txt \
 ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam

echo lc_extrap
echo future yield
date

#analyze read distributions using RSeQC -- will give you number of reads over different regions of genome dependent on the input annotation file (typically will be over genes, but you can imagine the applicability to eRNA calls eventually for quick comparison of relative eRNA coverage between samples)

module load python/2.7.14/rseqc

read_distribution.py \
 -i ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam \
 -r /scratch/Users/magr0763/genomes/GeneAnnotations/hg38_refseq.bed \
 > ${PROJECT}/qc/rseqc/${ROOTNAME}.read_dist.txt

echo rseqc
date
date

#BedGraph generator -- generating positive and negative strand mapping; awk coverts values to negative values for negative strand mapping

#NOTE: For the current input of FStitch, it is VERY important to include "pos" somewhere in your output BedGraph file name, otherwise you will run into problems with your training data

$DEEP_CONTAINER bamCoverage \
 --numberOfProcessors 16 \
 -b ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam \
 --filterRNAstrand forward \
 --normalizeUsing RPKM \
 --effectiveGenomeSize 2913022398 \
 -of bedgraph \
 -o ${PROJECT}/mapped/bedgraph/${ROOTNAME}.pos.trim.rpkm.bedGraph

echo forward bedgraph
date
date

$DEEP_CONTAINER bamCoverage \
 --numberOfProcessors 16 \
 -b ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam \
 --filterRNAstrand reverse \
 --normalizeUsing RPKM \
 --effectiveGenomeSize 2913022398 \
 -of bedgraph \
 -o ${PROJECT}/mapped/bedgraph/${ROOTNAME}.tmp.neg.trim.rpkm.bedGraph

echo reverse bedgraph
date
date

awk 'BEGIN{FS=OFS="\t"} {$4=-$4}1' ${PROJECT}/mapped/bedgraph/${ROOTNAME}.tmp.neg.trim.rpkm.bedGraph \
 > ${PROJECT}/mapped/bedgraph/${ROOTNAME}.neg.trim.rpkm.bedGraph

rm ${PROJECT}/mapped/bedgraph/${ROOTNAME}.tmp.neg.trim.rpkm.bedGraph


cat \
 ${PROJECT}/mapped/bedgraph/${ROOTNAME}.pos.trim.rpkm.bedGraph \
 <(grep -v '^@' ${PROJECT}/mapped/bedgraph/${ROOTNAME}.neg.trim.rpkm.bedGraph) \
 | sortBed \
 > ${PROJECT}/mapped/bedgraph/${ROOTNAME}.trim.rpkm.bedGraph


echo BEDGRAPH
date
date


#Generate tdfs (compressed bedgraphs) -- mapped reads easily viewable in IGV

/opt/igvtools/2.3.75/igvtools toTDF \
 ${PROJECT}/mapped/bedgraph/${ROOTNAME}.trim.rpkm.bedGraph \
 ${PROJECT}/mapped/tdfs/${ROOTNAME}.trim.tdf \
 /Users/magr0763/igv/genomes/sizes/hg38.chrom.sizes

echo tdf
date
date

