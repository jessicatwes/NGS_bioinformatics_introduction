#!/bin/bash
#SBATCH --job-name=ChIP_map 		# Job name
#SBATCH --mail-type=FAIL,END 		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremail@host.com	# Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=32 			# Number of CPU (processer cores i.e. tasks)
#SBATCH --time=08:00:00 		# Time limit hrs:min:sec
#SBATCH --partition short
#SBATCH --mem=100gb 			# Memory limit
#SBATCH --output=/path/to/directory/%x_%j.out	# Standard output
#SBATCH --error=/path/to/directory/%x_%j.err	# Standard error

########################### SE ChIP-seq ########################### 
### Single-end (SE) ChIP-sequencing pipeline. This script will take the SE raw, unedited fastq file to TDF file format to use in the IGV browser.

rootname=				#name of file root
project=				#/path/to/project/directory/

pwd; hostname; date
date

########################### QC and Pre-processing ###########################
### Pre-trim quality check of fastqc file before trimming

module load fastqc/0.11.5
fastqc ${project}/fastq/${rootname}.fastq -o ${project}/fastq/FastQC/
echo pre-trim qual
date
date

### Trimming of additional sequences including adapters, rRNA contaminants, or any specified sequence. There are multiple trimming tools that can be used for this step (Trimmomatic, Cutadapt, TrimGalore, BBDuk). 
## See https://github.com/BioInfoTools/BBMap/blob/master/sh/bbduk2.sh for parameters for bbduk.

module load bbmap/38.05

bbduk.sh -Xmx20g \
 t=16 \
 overwrite= t \
 in=${project}/fastq/${rootname}.fastq \
 out=${project}/trimmed/${rootname}.trim.fastq \
 ref=/scratch/Users/magr0763/bbmap/resources/adapters.fa \
 ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
 maq=10 minlen=20 \
 stats=${project}/trimmed/${rootname}.trimstats.txt \
 refstats=${project}/trimmed/${rootname}.refstats.txt \
 ehist=${project}/trimmed/${rootname}.ehist.txt

echo bbduk trim
date
date

### Quality check of reads after trim. 
# Run QC on reads to compare to pre-trim reads and adjust parameters as needed

module load fastqc/0.11.5
fastqc ${project}/trimmed/${rootname}.bbduk.trim.fastq -o ${project}/trimmed/FastQC/
echo post-trim qual
date
date

########################### Alignment to reference genome ###########################

module load samtools/1.3.1
module load hisat2/2.1.0

### Map trimmed fastq files to reference genome
hisat2 -p 16 \
 --very-sensitive \
 --no-spliced-alignment \
 -x ${ref}/HISAT2_indexes/genome \
 -U ${project}/trimmed/${rootname}.trim.fastq \
 > ${project}/mapped/sams/${rootname}.trim.sam \
 2> ${project}/mapped/sams/${rootname}.trim.stderr

echo mapped, sam
date
date

########################### File formatting ###########################

### Word count SAM files, convert sam files to bam files (compressed, binary sam files)

wc \
 -l ${project}/mapped/sams/${rootname}.trim.sam \
 > ${project}/mapped/sams/${rootname}.trim.sam.wc

## samtools view converts file from SAM format into binary  BAM format. Can also conver to CRAM which is compressed, index version of BAM
samtools view \
 -S -b -o \
 ${project}/mapped/bams/${rootname}.trim.bam \
 ${project}/mapped/sams/${rootname}.trim.sam \
 2> ${project}/mapped/bams/${rootname}.trim.bam.err

## samtools flagstat calculate statistics based on the information in the FLAG field of each alignment. Useful for read count correction (rcc) and QC.
samtools flagstat \
 ${project}/mapped/bams/${rootname}.trim.bam \
 > ${project}/mapped/bams/${rootname}.trim.bam.flagstat \
 2> ${project}/mapped/bams/${rootname}.trim.bam.flagstat.err

echo bam
date
date

## samtools sort sequence alignments in file based on leftmost coordinates
samtools sort \
 -m 100G ${project}/mapped/bams/${rootname}.trim.bam \
 > ${project}/mapped/bams/sorted/${rootname}.trim.sorted.bam

samtools flagstat \
 ${project}/mapped/bams/sorted/${rootname}.trim.sorted.bam \
 > ${project}/mapped/bams/sorted/${rootname}.trim.sorted.bam.flagstat \
 2> ${project}/mapped/bams/sorted/${rootname}.trim.sorted.bam.flagstat.err

echo sorted.bam
date
date

## samtools index creates a file to enable rapid retrieval of alignment of genomic regions. Used in multicov or to select specific genomic regions
samtools index \
 ${project}/mapped/bams/sorted/${rootname}.trim.sorted.bam \
 ${project}/mapped/bams/sorted/${rootname}.trim.sorted.bam.bai

echo indexed
date
date

########################### QC for report summary ###########################
### Preseq is tool designed to estimate and predict complexity  of sequencing libary.
module load preseq/2.0.3

## preseq c_curve takes bam file and plot the estimated complexity of a sample. Parameters -B for sorted bam files and -P for paired end read files. 

preseq c_curve -B \
 -o ${project}/preseq/${rootname}.c_curve.txt \
 ${project}/mapped/bams/sorted/${rootname}.trim.sorted.bam

echo c_curve
echo complexity
date

## preseq lc_extrap takes bam file and estimate predicted yield for complexity if the sample were to be sequenced at higher read depths. Parameters -B for sorted bam files and -P for paired end reads

preseq lc_extrap -B \
 -o ${project}/preseq/${rootname}.lc_extrap.txt \
 ${project}/mapped/bams/sorted/${rootname}.trim.sorted.bam

echo lc_extrap
echo future yield
date

### RSeQC evaluates sequence quality including nucleotide composition bias, PCR and GC bias, mapped reads distribution, strand specificity, and other quality output. Read distributions gives the number of reads over different region of genome dependent on input annotation file per samples.

module load python/2.7.14/rseqc

read_distribution.py \
 -i ${project}/mapped/bams/sorted/${rootname}.trim.sorted.bam \
 -r /scratch/Users/magr0763/genomes/GeneAnnotations/hg38_refseq.bed \
 > ${project}/rseqc/${rootname}.read_dist.txt

echo rseqc
date
date

### Coverage stats. Use tool in the bbmap suite to get the GC content over sample per chromosome and the number of bases covered by reads (coverage, not depth)

echo generating coverage stats
date
date

module load java/1.8.0_101

pileup.sh \
 in=${project}/mapped/bams/sorted/${rootname}.trim.sorted.bam \
 out=${project}/qc/coverage/${rootname}.coverage.stats.txt \
 hist=${project}/qc/coverage/${rootname}.coverage.hist.txt

echo coverage stats done
date
date

### BedTools generates BedGraph; positive and negative strand mapping

module load bedtools/2.25.0

genomeCoverageBed \
 -bg -strand + -g hg38 \
 -ibam ${project}/mapped/bams/sorted/${rootname}.trim.sorted.bam \
 > ${project}/mapped/bedgraph/${rootname}.trim.unsorted.BedGraph 

sortBed \
 -i ${project}/mapped/bedgraph/${rootname}.trim.unsorted.BedGraph \
 > ${project}/mapped/bedgraph/${rootname}.trim.BedGraph

echo BedGraph.pos.neg
date
date

### Read count correcting (rcc) reduces sequencing bias based on depth for visualization in IGV

module load python/2.7.14

python \
 /Users/magr0763/mappingandvis/readcountcorrectBG.py \
 ${project}/mapped/bedgraph/${rootname}.trim.BedGraph \
 ${project}/mapped/bams/sorted/${rootname}.trim.sorted.bam.flagstat \
 ${project}/mapped/bedgraph/${rootname}.trim.rcc.BedGraph

echo readcountcorrectedbedgraph
date
date

### Generate compressed bedgraphs, TDF files, for easier use in IGV

/opt/igvtools/2.3.75/igvtools toTDF \
 ${project}/mapped/bedgraph/${rootname}.trim.rcc.BedGraph \
 ${project}/mapped/tdfs/${rootname}.trim.tdf \
 /Users/magr0763/igv/genomes/sizes/hg38.chrom.sizes

echo tdf
date
date
