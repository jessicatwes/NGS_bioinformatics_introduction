#!/bin/bash
#SBATCH --job-name=RNAseq		# Job name
#SBATCH --mail-type=FAIL,END            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremail@host.com  # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=16                     # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=08:00:00                 # Time limit hrs:min:sec
#SBATCH --partition short
#SBATCH --mem=100gb                     # Memory limit
#SBATCH --output=/path/to/directory/%x_%j.out   # Standard output
#SBATCH --error=/path/to/directory/%x_%j.err    # Standard error

########################### RNA-seq ########################### 
### Single strand RNA-sequencing pipeline. This script will take the SE raw, unedited fastq file to TDF file format to use in the IGV browser.

#This script will take you from a raw, unedited fastq file to tdf

ROOTNAME=$1
PROJECT='/path/to/RNAseq/directory'
GENOME='/scratch/Users/magr0763/genomes/mapping/hg38/'
ADAPTER='/path/to/opt/bbmap/38.05/adapter.fa'   #adapter from bbmap/38.05
GENOMEINDEX='resources/HISAT2/genome'   #see resource for detail on where to find index

pwd; hostname; date
date

########################### QC and Pre-processing ###########################

module load fastqc/0.11.5
fastqc ${PROJECT}/fastq/${ROOTNAME}.fastq -o ${PROJECT}/fastq/FastQC/
echo pre-trim qual
date
date

### Trimming of additional sequences including adapters, rRNA contaminants, or any specified sequence. The sequence from adapter.fa is for rRNA contaminants, but you can specify any sequence (similar to -a flag in Cutadapt). 

## There are multiple trimming tools that can be used for this step (Trimmomatic, Cutadapt, TrimGalore, BBDuk). BBDuk has a rapid runtime and useful if you know your data. BBDuk can not autodetect adapters so you may miss the sequence.


module load bbmap/38.05

bbduk.sh -Xmx20g \
t=16 \
overwrite= t \
in=${PROJECT}/fastq/${ROOTNAME}.fastq \
out=${PROJECT}/trimmed/${ROOTNAME}.trim.fastq \
ref=${ADAPTER}
ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
maq=10 minlen=20 \
stats=${PROJECT}/trimmed/${ROOTNAME}.trimstats.txt \
refstats=${PROJECT}/trimmed/${ROOTNAME}.refstats.txt \
ehist=${PROJECT}/trimmed/${ROOTNAME}.ehist.txt

echo bbduk trim
date
date

### Quality check of reads after trim. 
# Run QC on reads to compare to pre-trim reads and adjust parameters as needed

module load fastqc/0.11.5
fastqc ${PROJECT}/trimmed/${ROOTNAME}.trim.fastq -o ${PROJECT}/trimmed/FastQC/
echo post-trim qual
date
date

########################### Alignment to reference genome ###########################

module load samtools/1.3.1
module load hisat2/2.1.0

### HISAT2 is a good aligner for RNA-seq as it is splice-aware. For RNA-seq data, you can adjust parameters as they are related to how the aligner handles certain aspects of the alignment.
# --pen-noncansplice 14	#penalty for intron-containing alignment where splice site not at canonical GT-AG boundaries. Higher penalty of 14 means alignment with non-canonical splice site less favored in mapping
# --mp 1,0	#penalty for maximum and minimum mismatches in alignment. Higher penalty value makes aligner more sensitive to mismatch but more stringent
# --sp 3,1	#penalty for short and long insertions and deletions. Higher penalty for short indels may make aligner more strigent and reduce false positive alignments

hisat2 -p 16 \
 --very-sensitive \
 -x ${GENOMEINDEX}
 --pen-noncansplice 14 \
 --mp 1,0 \
 -U ${PROJECT}/trimmed/${ROOTNAME}.trim.fastq \
 > ${PROJECT}/mapped/sams/${ROOTNAME}.trim.sam \
 2> ${PROJECT}/mapped/sams/${ROOTNAME}.trim.stderr

echo mapped, sam
date
date

########################### File formatting ###########################

### Word count SAM files, convert sam files to bam files (compressed, binary sam files)i
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
 > ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam

samtools flagstat \
 ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam \
 > ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam.flagstat \
 2> ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam.flagstat.err

echo sorted.bam
date
date

## SAMtools index creates a file to enable rapid retrieval of alignment of genomic regions. Used in multicov or to select specific genomic regions
samtools index \
 ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam \
 ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam.bai

echo indexed
date
date

########################### QC for report summary ###########################
### Preseq is tool designed to estimate and predict complexity  of sequencing libary.
module load preseq/2.0.3

## preseq c_curve takes bam file and plot the estimated complexity of a sample. Parameters -B for sorted bam files and -P for paired end read files.
preseq c_curve -B \
 -o ${PROJECT}/qc/preseq/${ROOTNAME}.c_curve.txt \
 ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam

echo c_curve
echo complexity
date

## preseq lc_extrap takes bam file and estimate predicted yield for complexity if the sample were to be sequenced at higher read depths. Parameters -B for sorted bam files and -P for paired end reads
preseq lc_extrap -B \
 -o ${PROJECT}/qc/preseq/${ROOTNAME}.lc_extrap.txt \
 ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam

echo lc_extrap
echo future yield
date

### RSeQC evaluates sequence quality including nucleotide composition bias, PCR and GC bias, mapped reads distribution, strand specificity, and other quality output. Read distributions gives the number of reads over different region of genome dependent on input annotation file per samples.

module load python/2.7.14/rseqc

read_distribution.py \
 -i ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam \
 -r /scratch/Users/magr0763/genomes/GeneAnnotations/hg38_refseq.bed \
 > ${PROJECT}/qc/rseqc/${ROOTNAME}.read_dist.txt

echo rseqc
date
date

### deepTools bamCoverage used to generate bedGraph (or bigWig), a coverage track calculated as number of reads per bin of a defined size. Use bamCoverage to normalize by scaling factor, Reads Per Kilobase per Million mapped reads (RPKM).
## Install deeptools. Information on software found https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html

## RPKM forward strand bedGraph
bamCoverage \
 --numberOfProcessors 16 \
 -b ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam \
 --filterRNAstrand forward \
 --normalizeUsing RPKM \
 --effectiveGenomeSize 2913022398 \
 -of bedgraph \
 -o ${PROJECT}/mapped/bedgraph/${ROOTNAME}.pos.trim.rpkm.bedGraph

echo forward bedGraph
date
date

## RPKM reverse strand bedGraph
bamCoverage \
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

### Generate compressed bedgraphs, TDF files, for easier use in IGV
## Need to install igvtools to use script below

/path/to/opt/igvtools/2.3.75/igvtools toTDF \
 ${PROJECT}/mapped/bedgraph/${ROOTNAME}.trim.rpkm.bedGraph \
 ${PROJECT}/mapped/tdfs/${ROOTNAME}.trim.tdf \
 /Users/magr0763/igv/genomes/sizes/hg38.chrom.sizes

echo tdf
date
date

