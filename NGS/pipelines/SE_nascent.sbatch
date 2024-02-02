#!/bin/bash
#SBATCH --job-name=nascent		# Job name
#SBATCH --mail-type=FAIL,END            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremail@host.com  # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=16                     # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=16:00:00                 # Time limit hrs:min:sec
#SBATCH --partition short
#SBATCH --mem=100gb                     # Memory limit
#SBATCH --output=/path/to/directory/%x_%j.out   # Standard output
#SBATCH --error=/path/to/directory/%x_%j.err    # Standard error

########################################################################
########################### Nassent transcription ###########################
### Nascent transcription sequencing assays such as GRO-seq and PRO-seq pipeline. Workflow takes raw, unedited SE fastq file to TDF file format to use in the IGV browser.

########################################################################
################### Set up variable and environment ######################################

##### Set variables #####
ROOTNAME=$1
PROJECT='/path/to/directory/'
ADAPTER='/path/to/opt/bbmap/38.05/adapter.fa'   #adapter from bbmap/38.05

##GENOME='/scratch/Shares/dowell/rn6'
##DEEP_CONTAINER='singularity exec -H /scratch/Users/magr0763 /scratch/Shares/public/singularity/deeptools-3.0.1-py35_1.img'

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
########################### QC and Pre-processing ###########################

##### FASTX-toolkit to  producee  reverse complement #####
### Depending on protocol, your short reads may need to be  flipped.

fastx_reverse_complement \
 -Q33 -i ${PROJECT}/fastq/${ROOTNAME}.fastq \
 -o ${PROJECT}/fastq/${ROOTNAME}.flip.fastq 

echo flipped
date
date

##### Pre-trim quality check of fastqc file before trimming #####

fastqc ${PROJECT}/fastq/${ROOTNAME}.fastq -o ${PROJECT}/fastq/FastQC/

echo pre-trim qual
date
date

### Trimming of additional sequences including adapters, rRNA contaminants, or any specified sequence. 'literal' parameter to add additional sequence to trim. There are multiple trimming tools that can be used for this step (Trimmomatic, Cutadapt, TrimGalore, BBDuk).
## See https://github.com/BioInfoTools/BBMap/blob/master/sh/bbduk2.sh for parameters for bbduk.

bbduk.sh -Xmx20g \
t=16 \
overwrite= t \
in=${PROJECT}/fastq/${ROOTNAME}.fastq \
out=${PROJECT}/trimmed/${ROOTNAME}.trim.fastq \
ref=${ADAPTER} \
ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
maq=10 minlen=20 \
literal=AAAAAAAAAAAAAAAAAAAAAAA \
stats=${PROJECT}/trimmed/${ROOTNAME}.trimstats.txt \
refstats=${PROJECT}/trimmed/${ROOTNAME}.refstats.txt \
ehist=${PROJECT}/trimmed/${ROOTNAME}.ehist.txt

echo bbduk trim
date
date

##### Quality check of reads after trim  #####
# Run QC on reads to compare to pre-trim reads and adjust parameters as needed

fastqc ${PROJECT}/trimmed/${ROOTNAME}.trim.fastq -o ${PROJECT}/trimmed/FastQC/

echo post-trim qual
date
date

########################################################################
########################### Alignment to reference genome ###########################

##### Map trimmed fastq files to reference genome #####

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

########################### File formatting ###########################

##### Word count SAM files, convert sam files to bam files (compressed, binary sam files) #####
wc \
 -l ${PROJECT}/mapped/sams/${ROOTNAME}.trim.sam \
 > ${PROJECT}/mapped/sams/${ROOTNAME}.trim.sam.wc

## SAMtools view converts file from SAM format into binary  BAM format. Can also conver to CRAM which is compressed, index version of BAM
samtools view \
 -S -b -o \
 ${PROJECT}/mapped/bams/${ROOTNAME}.trim.bam \
 ${PROJECT}/mapped/sams/${ROOTNAME}.trim.sam \
 2> ${PROJECT}/mapped/bams/${ROOTNAME}.trim.bam.err

## SAMtools flagstat calculate statistics based on the information in the FLAG field of each alignment. Useful for read count correction (rcc) and QC.
samtools flagstat \
 ${PROJECT}/mapped/bams/${ROOTNAME}.trim.bam \
 > ${PROJECT}/mapped/bams/${ROOTNAME}.trim.bam.flagstat \
 2> ${PROJECT}/mapped/bams/${ROOTNAME}.trim.bam.flagstat.err

echo bam
date
date

## SAMtools sort sequence alignments in file based on leftmost coordinates
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

## SAMtools index creates a file to enable rapid retrieval of alignment of genomic regions. Used in multicov or to select specific genomic regions
samtools index \
 ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam \
 ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam.bai

echo indexed
date
date


########################################################################
########################### QC for report summary ###########################
##### Preseq is tool designed to estimate and predict complexity  of sequencing libary. #####

## preseq c_curve takes bam file and plot the estimated complexity of a sample. Parameters -B for sorted bam files and -P for paired end read files.
preseq c_curve -B \
 -o ${PROJECT}/qc/preseq/${ROOTNAME}.c_curve.txt \
 ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam

echo c_curve
echo complexity
date

## preseq lc_extrap takes bam file and estimate predicted yield for complexity if the sample were to be sequenced at higher read depths. Parameters -B for sorted bam files and -P for paired end reads
preseq lc_extrap -B \
 -o ${PROJECT}/qc/preseq/${ROOTNAME}.lc_extrap.txt \
 ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam

echo lc_extrap
echo future yield
date

##### RSeQC evaluates sequence quality including nucleotide composition bias, PCR and GC bias, mapped reads distribution, strand specificity, and other quality output. Read distributions gives the number of reads over different region of genome dependent on input annotation file per samples. #####

#analyze read distributions using RSeQC -- will give you number of reads over different regions of genome dependent on the input annotation file (typically will be over genes, but you can imagine the applicability to eRNA calls eventually for quick comparison of relative eRNA coverage between samples)


read_distribution.py \
 -i ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam \
 -r ${GENOME}/hg38_refseq.bed \
 > ${PROJECT}/qc/rseqc/${ROOTNAME}.read_dist.txt

echo rseqc
date
date

##### Coverage stats. Use tool in the bbmap suite to get the GC content over sample per chromosome and the number of bases covered by reads (coverage, not depth) #####

pileup.sh \
 in=${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam \
 out=${PROJECT}/qc/coverage/${ROOTNAME}.coverage.stats.txt \
 hist=${PROJECT}/qc/coverage/${ROOTNAME}.coverage.hist.txt

echo pileup coverage
date
date

##### BedTools generates bedGraph as output; This is a coverage track of the positive and negative strand mapping #####

#NOTE: If you want to use libraries such as FStitch(https://github.com/Dowell-Lab/FStitch), it is important to include "pos" in your output BedGraph file name in order to run your training data.


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
##### deepTools bamCoverage used to generate bedGraph (or bigWig), a coverage track calculated as number of reads per bin of a defined size. Use bamCoverage to normalize by scaling factor, Reads Per Kilobase per Million mapped reads (RPKM) #####.
## Install deeptools. Information on software found https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html

## RPKM forward strand bigWig
bamCoverage \
 --numberOfProcessors 16 \
 -b ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam \
 --filterRNAstrand forward \
 --normalizeUsing RPKM \
 --effectiveGenomeSize 2913022398 \
 -of bigwig \
 -o ${PROJECT}/mapped/bigwig/${ROOTNAME}.pos.trim.rpkm.bw

echo forward bigWig bw
date
date

## RPKM reverse strand bigWig
bamCoverage \
 --numberOfProcessors 16 \
 -b ${PROJECT}/mapped/bams/${ROOTNAME}.trim.sorted.bam \
 --filterRNAstrand reverse \
 --normalizeUsing RPKM \
 --effectiveGenomeSize 2913022398 \
 -of bigwig \
 -o ${PROJECT}/mapped/bigwig/${ROOTNAME}.neg.trim.rpkm.bw

echo reverse bigWig bw
date
date

## RPKM forward strand bedGraph
bamCoverage \
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

## RPKM reverse strand bedGraph
Coverage \
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

## awk coverts values to negative values for negative strand mapping
awk 'BEGIN{FS=OFS="\t"} {$4=-$4}1' ${PROJECT}/mapped/bedgraph/${ROOTNAME}.tmp.neg.trim.rpkm.bedGraph \
 > ${PROJECT}/mapped/bedgraph/${ROOTNAME}.neg.trim.rpkm.bedGraph

rm ${PROJECT}/mapped/bedgraph/${ROOTNAME}.tmp.neg.trim.rpkm.bedGraph

## concatenate positive and negative strands
cat \
 ${PROJECT}/mapped/bedgraph/${ROOTNAME}.pos.trim.rpkm.bedGraph \
 <(grep -v '^@' ${PROJECT}/mapped/bedgraph/${ROOTNAME}.neg.trim.rpkm.bedGraph) \
 | sortBed \
 > ${PROJECT}/mapped/bedgraph/${ROOTNAME}.trim.rpkm.bedGraph

echo BEDGRAPH
date
date

##### Generate compressed bedgraphs, TDF files, for easier use in IGV #####
## Need to install IGVtools to use script below

/path/to/opt/igvtools/2.3.75/igvtools toTDF \
 ${PROJECT}/mapped/bedgraph/${ROOTNAME}.trim.rpkm.bedGraph \
 ${PROJECT}/mapped/tdfs/${ROOTNAME}.trim.rpkm.tdf \
 /Users/magr0763/igv/genomes/sizes/hg38.chrom.sizes

echo tdf
date
date

