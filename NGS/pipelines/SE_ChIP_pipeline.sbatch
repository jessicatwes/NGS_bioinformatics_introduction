#!/bin/bash
#SBATCH --job-name=ChIPseq 		# Job name
#SBATCH --mail-type=FAIL,END 		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremail@host.com	# Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=32 			# Number of CPU (processer cores i.e. tasks)
#SBATCH --time=08:00:00 		# Time limit hrs:min:sec
#SBATCH --partition short
#SBATCH --mem=100gb 			# Memory limit
#SBATCH --output=/path/to/directory/%x_%j.out	# Standard output
#SBATCH --error=/path/to/directory/%x_%j.err	# Standard error

########################################################################
########################### SE ChIP-seq ########################### 
### Single-end (SE) ChIP-sequencing pipeline. This script will take the SE raw, unedited fastq file to TDF file format to use in the IGV browser.

########################################################################
################### Set up variable and environment ######################################

##### Set variables #####
ROOTNAME='rootname'			#root name for your project/files
PROJECT='/path/to/directory/'
ADAPTER='/path/to/opt/bbmap/38.05/adapter.fa'	#adapter from bbmap/38.05	
GENOMEINDEX='resources/HISAT2/genome'	#see resource for detail on where to find index
REFGENOME='resources/hg38_refseq.bed'	
CHROMESIZE='resources/hg38.chrom.sizes'

##### Load modules #####
module load fastqc/0.11.5	#QC
module load bbmap/38.05  	#trimming
module load samtools/1.3.1	
module load hisat2/2.1.0	#alignment
module load preseq/2.0.3	
module load python/2.7.14/rseqc	#read distribution stats
module load java/1.8.0_101	#coverage stats
module load bedtools/2.25.0
module load python/2.7.14       #Python readcountcorrectBG.py script

pwd; hostname; date
date

########################################################################
########################### QC and Pre-processing ###########################
##### Pre-trim quality check of fastqc file before trimming #####

fastqc ${tPROJECT}/fastq/${ROOTNAME}.fastq -o ${PROJECT}/fastq/FastQC/
echo pre-trim qual
date
date

### Trimming of additional sequences including adapters, rRNA contaminants, or any specified sequence. There are multiple trimming tools that can be used for this step (Trimmomatic, Cutadapt, TrimGalore, BBDuk). 
## See https://github.com/BioInfoTools/BBMap/blob/master/sh/bbduk2.sh for parameters for bbduk.

bbduk.sh -Xmx20g \
 t=16 \
 overwrite= t \
 in=${PROJECT}/fastq/${rootname}.fastq \
 out=${PROJECT}/trimmed/${rootname}.trim.fastq \
 ref=${ADAPTER}  \
 ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
 maq=10 minlen=20 \
 stats=${PROJECT}/trimmed/${rootname}.trimstats.txt \
 refstats=${PROJECT}/trimmed/${rootname}.refstats.txt \
 ehist=${PROJECT}/trimmed/${rootname}.ehist.txt

echo bbduk trim
date
date

##### Quality check of reads after trim  #####
# Run QC on reads to compare to pre-trim reads and adjust parameters as needed

fastqc ${PROJECT}/trimmed/${ROOTNAME}.bbduk.trim.fastq -o ${PROJECT}/trimmed/FastQC/
echo post-trim qual
date
date

########################################################################
########################### Alignment to reference genome ###########################

##### Map trimmed fastq files to reference genome #####
hisat2 -p 16 \
 --very-sensitive \
 --no-spliced-alignment \
 -x ${GENOMEINDEX} \
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

########################################################################
########################### QC for report summary ###########################
##### Preseq is tool designed to estimate and predict complexity  of sequencing libary. #####

## preseq c_curve takes bam file and plot the estimated complexity of a sample. Parameters -B for sorted bam files and -P for paired end read files. 
preseq c_curve -B \
 -o ${PROJECT}/preseq/${ROOTNAME}.c_curve.txt \
 ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam

echo c_curve
echo complexity
date

## preseq lc_extrap takes bam file and estimate predicted yield for complexity if the sample were to be sequenced at higher read depths. Parameters -B for sorted bam files and -P for paired end reads
preseq lc_extrap -B \
 -o ${PROJECT}/preseq/${ROOTNAME}.lc_extrap.txt \
 ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam

echo lc_extrap
echo future yield
date

##### RSeQC evaluates sequence quality including nucleotide composition bias, PCR and GC bias, mapped reads distribution, strand specificity, and other quality output. Read distributions gives the number of reads over different region of genome dependent on input annotation file per samples. #####

read_distribution.py \
 -i ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam \
 -r ${REFGENOME} \
 > ${PROJECT}/rseqc/${ROOTNAME}.read_dist.txt

echo rseqc
date
date

##### Coverage stats. Use tool in the bbmap suite to get the GC content over sample per chromosome and the number of bases covered by reads (coverage, not depth) #####

echo generating coverage stats
date
date

pileup.sh \
 in=${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam \
 out=${PROJECT}/qc/coverage/${ROOTNAME}.coverage.stats.txt \
 hist=${PROJECT}/qc/coverage/${ROOTNAME}.coverage.hist.txt

echo coverage stats done
date
date

##### BedTools generates bedGraph as output; This is a coverage track of the positive and negative strand mapping #####

genomeCoverageBed \
B -bg -strand + -g hg38 \
 -ibam ${PROJECT}/mapped/bams/sorted/${ROOTNAME}.trim.sorted.bam \
 > ${PROJECT}/mapped/bedgraph/${ROOTNAME}.trim.unsorted.BedGraph 

sortBed \
 -i ${PROJECT}/mapped/bedgraph/${ROOTNAME}.trim.unsorted.BedGraph \
 > ${PROJECT}/mapped/bedgraph/${ROOTNAME}.trim.BedGraph

echo BedGraph.pos.neg
date
date

##### Read count correcting (rcc) reduces sequencing bias based on depth for visualization in IGV. To run you need to git clone https://github.com/maallen3/mapandvischainonfijii. #####

python \
 /path/to/readcountcorrectBG.py \
 ${PROJECT}/mapped/bedgraph/${ROOTNAME}.trim.BedGraph \
 ${PROJECT}maapped/bams/sorted/${ROOTNAME}.trim.sorted.bam.flagstat \
 ${PROJECT}/mapped/bedgraph/${ROOTNAME}.trim.rcc.BedGraph

echo read count corrected bedgraph
date
date

##### Generate compressed bedgraphs, TDF files, for easier use in IGV  #####
## Need to install IGVtools to use script below

/path/to/opt/igvtools/2.3.75/igvtools toTDF \
 ${PROJECT}/mapped/bedgraph/${ROOTNAME}.trim.rcc.BedGraph \
 ${PROJECT}/mapped/tdfs/${ROOTNAME}.trim.tdf \
 ${CHROMESIZE}

echo tdf
date
date
