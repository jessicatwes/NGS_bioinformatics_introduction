#!/bin/bash
#SBATCH --job-name=Nascent_map # Job name
#SBATCH --mail-type=FAIL,END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=margaret.gruca@colorado.edu # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=16 # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=16:00:00 # Time limit hrs:min:sec
#SBATCH --partition short
#SBATCH --mem=100gb # Memory limit
#SBATCH --output=/scratch/Users/magr0763/Andrysik2017/GRO-seq/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/magr0763/Andrysik2017/GRO-seq/e_and_o/%x_%j.err

#This script will take you from a raw, unedited fastq file to tdf


################### SET VARIABLES ######################################

ROOTNAME=$1
PROJECT='/scratch/Users/magr0763/Christa/RNA-seq'
GENOME='/scratch/Shares/dowell/rn6'
DEEP_CONTAINER='singularity exec -H /scratch/Users/magr0763 /scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img'

########################################################################
################### LOAD NECESSARY MODULES #############################

module load singularity/2.4
module load fastx-toolkit/0.0.13
module load fastqc/0.11.5
module load bbmap/38.05
module load samtools/1.8
module load hisat2/2.1.0
module load preseq/2.0.3
module load python/2.7.14/rseqc
module load bedtools/2.25.0
module load igvtools/2.3.75

########################################################################

printf "Sample ID: $ROOTNAME"
printf "Directory: $PROJECT"
printf "Genome annotation files: $GENOME"
printf "Run on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

#Produces reverse complement of each short-read seqeuence -- the need for this is protocol dependent

fastx_reverse_complement \
 -Q33 -i ${PROJECT}/fastq/${ROOTNAME}.fastq \
 -o ${PROJECT}/fastq/${ROOTNAME}.flip.fastq 

echo flipped
date
date

#Quality check of fastqc file before trimming

fastqc ${PROJECT}/fastq/${ROOTNAME}.fastq -o ${PROJECT}/fastq/FastQC/

echo pre-trim qual
date
date

###Trimming: there are a number of different tools available (Trimmomatic, cudadapt, TrimGalore). My favorite currently for most purposes in bbduk due to its rapid runtime and large number of trimming options. "literal" is commented out below. The current sequence I have there is for rrna contaminants, however you can specify any sequence (similar to -a flag in cutadapt). This is useful if you have to add an additional sequence/contaminant to trim rather than adding it to the fa file. Because of the large number of options, bbduk is particularly useful if you really know your data. For a first pass if you are not sure of the protocol, you will more than likely miss adapters as it has no autodetect feature.
###Go here for full list of options : https://github.com/BioInfoTools/BBMap/blob/master/sh/bbduk2.sh


bbduk.sh -Xmx20g \
t=16 \
overwrite= t \
in=${PROJECT}/fastq/${ROOTNAME}.fastq \
out=${PROJECT}/trimmed/${ROOTNAME}.trim.fastq \
ref=/opt/bbmap/38.05/resources/adapters.fa \
ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
maq=10 minlen=20 \
literal=AAAAAAAAAAAAAAAAAAAAAAA \
stats=${PROJECT}/trimmed/${ROOTNAME}.trimstats.txt \
refstats=${PROJECT}/trimmed/${ROOTNAME}.refstats.txt \
ehist=${PROJECT}/trimmed/${ROOTNAME}.ehist.txt

echo bbduk trim
date
date

#Quality check trimmed reads --> if first QC bad, check these before moving forward and adjust parameters if needed

fastqc ${PROJECT}/trimmed/${ROOTNAME}.trim.fastq -o ${PROJECT}/trimmed/FastQC/

echo post-trim qual
date
date

#Map trimmed fastq files to reference genome --> sam files include all reads


hisat2 -p 16 \
 --very-sensitive \
 --no-spliced-alignment \
 -x ${GENOME}/HISAT2_indexes/genome \
 -U ${PROJECT}/trimmed/${ROOTNAME}.trim.fastq \
 > ${PROJECT}/mapped/sams/${ROOTNAME}.trim.sam \
 2> ${PROJECT}/mapped/sams/${ROOTNAME}.trim.stderr

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
 > ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam

samtools flagstat \
 ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam \
 > ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam.flagstat \
 2> ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam.flagstat.err

rm ${PROJECT}/mapped/bams/${ROOTNAME}.trim.bam

echo sorted.bam
date
date

#Index sorted bam files for use in mutlicov OR if you decide to sort out a specific region (e.g. use samtools view to select only chromosome 1)

samtools index \
 ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam \
 ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam.bai

echo indexed
date
date


#c_curve : can use inputs from file types [bed/bam] and will plot the estimated complexity of a sample. Specify -B for sorted bam files, -P for paired end read files

preseq c_curve -B \
 -o ${PROJECT}/qc/preseq/${ROOTNAME}.c_curve.txt \
 ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam

echo c_curve
echo complexity
date

#lc_extrap : can use inputs from from tyles [bed/bam] and will estimate future yields for complexity if the sample is sequenced at higher read depths. Specify -B for sorted bam files, -P for paired end reads

preseq lc_extrap -B \
 -o ${PROJECT}/qc/preseq/${ROOTNAME}.lc_extrap.txt \
 ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam

echo lc_extrap
echo future yield
date

#analyze read distributions using RSeQC -- will give you number of reads over different regions of genome dependent on the input annotation file (typically will be over genes, but you can imagine the applicability to eRNA calls eventually for quick comparison of relative eRNA coverage between samples)


read_distribution.py \
 -i ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam \
 -r ${GENOME}/hg38_refseq.bed \
 > ${PROJECT}/qc/rseqc/${ROOTNAME}.read_dist.txt

echo rseqc
date
date

#### Analyze GC content over sample, per-chromosome coverage stats, etc.

pileup.sh \
 in=${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam \
 out=${PROJECT}/qc/coverage/${ROOTNAME}.coverage.stats.txt \
 hist=${PROJECT}/qc/coverage/${ROOTNAME}.coverage.hist.txt

echo pileup coverage
date
date

#BedGraph generator -- generating positive and negative strand mapping; awk coverts values to negative values for negative strand mapping

#NOTE: For the current input of FStitch, it is VERY important to include "pos" somewhere in your output BedGraph file name, otherwise you will run into problems with your training data

#DEEP_CONTAINER is a singularity package for deeptools

# .bt extension on bedgraphs generated from bedtools

genomeCoverageBed \
 -bg \
 -strand + \
 -g hg38 \
 -ibam ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam \
  > ${PROJECT}/mapped/bedgraph/${ROOTNAME}.bt.trim.pos.bedGraph

genomeCoverageBed \
 -bg \
 -strand - \
 -g hg38 \
 -ibam ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam \
 | awk -F '	' -v OFS='	' '{ $4 = - $4 ; print $0 }' \
 > ${PROJECT}/mapped/bedgraph/${ROOTNAME}.bt.trim.neg.bedGraph

cat ${PROJECT}/mapped/bedgraph/${ROOTNAME}.bt.trim.pos.bedGraph \
 ${PROJECT}/mapped/bedgraph/${ROOTNAME}.bt.trim.neg.bedGraph \
 > ${PROJECT}/BedGraph/${ROOTNAME}.bt.trim.unsorted.bedGraph

sortBed \
 -i ${PROJECT}/mapped/bedgraph/${ROOTNAME}.bt.trim.unsorted.bedGraph \
 > ${PROJECT}/mapped/bedgraph/${ROOTNAME}.bt.trim.bedGraph

echo BedGraph.pos.neg
date
date

################# NORMALIZED BEDGRAPHS AND BIGWIGS FOR VISUALIZATION ###################

$DEEP_CONTAINER bamCoverage \
 --numberOfProcessors 16 \
 -b ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam \
 --filterRNAstrand forward \
 --normalizeUsing RPKM \
 --effectiveGenomeSize 2913022398 \
 -of bigwig \
 -o ${PROJECT}/mapped/bigwig/${ROOTNAME}.pos.trim.rpkm.bw

echo forward bw
date
date

$DEEP_CONTAINER bamCoverage \
 --numberOfProcessors 16 \
 -b ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam \
 --filterRNAstrand reverse \
 --normalizeUsing RPKM \
 --effectiveGenomeSize 2913022398 \
 -of bigwig \
 -o ${PROJECT}/mapped/bigwig/${ROOTNAME}.neg.trim.rpkm.bw

echo reverse bigwig
date
date

$DEEP_CONTAINER bamCoverage \
 --numberOfProcessors 16 \
 -b ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam \
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
 -b ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam \
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
 ${PROJECT}/mapped/tdfs/${ROOTNAME}.trim.rpkm.tdf \
 /Users/magr0763/igv/genomes/sizes/hg38.chrom.sizes

echo tdf
date
date

