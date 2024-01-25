#!/bin/bash
#SBATCH --job-name=ChIP_map # Job name
#SBATCH --mail-type=FAIL,END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=margaret.gruca@colorado.edu # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=32 # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=08:00:00 # Time limit hrs:min:sec
#SBATCH --partition short
#SBATCH --mem=100gb # Memory limit
#SBATCH --output=/scratch/Users/magr0763/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/magr0763/e_and_o/%x_%j.err

#This script will take you from a raw, unedited fastq file to tdf

rootname=
project=
ref=

pwd; hostname; date
date

#Quality check of fastqc file before trimming

module load fastqc/0.11.5
fastqc ${project}/fastq/${rootname}.fastq -o ${project}/fastq/FastQC/
echo pre-trim qual
date
date

###Trimming: there are a number of different tools available (Trimmomatic, cudadapt, TrimGalore). My favorite currently for most purposes in bbduk due to its rapid runtime and large number of trimming options. "literal" is commented out below. The current sequence I have there is for rrna contaminants, however you can specify any sequence (similar to -a flag in cutadapt). This is useful if you have to add an additional sequence/contaminant to trim rather than adding it to the fa file. Because of the large number of options, bbduk is particularly useful if you really know your data. For a first pass if you are not sure of the protocol, you will more than likely miss adapters as it has no autodetect feature.
###Go here for full list of options : https://github.com/BioInfoTools/BBMap/blob/master/sh/bbduk2.sh

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

#Quality check trimmed reads --> if first QC bad, check these before moving forward and adjust parameters if needed

module load fastqc/0.11.5
fastqc ${project}/trimmed/${rootname}.bbduk.trim.fastq -o ${project}/trimmed/FastQC/
echo post-trim qual
date
date

#Map trimmed fastq files to reference genome --> sam files include all reads

module load samtools/1.3.1
module load hisat2/2.1.0

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

#Word count sam files, convert sam files to bam files (compressed, binary sam files)

wc \
 -l ${project}/mapped/sams/${rootname}.trim.sam \
 > ${project}/mapped/sams/${rootname}.trim.sam.wc

samtools view \
 -S -b -o \
 ${project}/mapped/bams/${rootname}.trim.bam \
 ${project}/mapped/sams/${rootname}.trim.sam \
 2> ${project}/mapped/bams/${rootname}.trim.bam.err

samtools flagstat \
 ${project}/mapped/bams/${rootname}.trim.bam \
 > ${project}/mapped/bams/${rootname}.trim.bam.flagstat \
 2> ${project}/mapped/bams/${rootname}.trim.bam.flagstat.err

echo bam
date
date

#Sort bam files, flagstat useful for rcc (read count correction) and QC

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

#Index sorted bam files for use in mutlicov OR if you decide to sort out a specific region (e.g. use samtools view to select only chromosome 1)

samtools index \
 ${project}/mapped/bams/sorted/${rootname}.trim.sorted.bam \
 ${project}/mapped/bams/sorted/${rootname}.trim.sorted.bam.bai

echo indexed
date
date

module load preseq/2.0.3

#c_curve : can use inputs from file types [bed/bam] and will plot the estimated complexity of a sample. Specify -B for sorted bam files, -P for paired end read files

preseq c_curve -B \
 -o ${project}/preseq/${rootname}.c_curve.txt \
 ${project}/mapped/bams/sorted/${rootname}.trim.sorted.bam

echo c_curve
echo complexity
date

#lc_extrap : can use inputs from from tyles [bed/bam] and will estimate future yields for complexity if the sample is sequenced at higher read depths. Specify -B for sorted bam files, -P for paired end reads

preseq lc_extrap -B \
 -o ${project}/preseq/${rootname}.lc_extrap.txt \
 ${project}/mapped/bams/sorted/${rootname}.trim.sorted.bam

echo lc_extrap
echo future yield
date

#analyze read distributions using RSeQC -- will give you number of reads over different regions of genome dependent on the input annotation file (typically will be over genes, but you can imagine the applicability to eRNA calls eventually for quick comparison of relative eRNA coverage between samples)

module load python/2.7.14/rseqc

read_distribution.py \
 -i ${project}/mapped/bams/sorted/${rootname}.trim.sorted.bam \
 -r /scratch/Users/magr0763/genomes/GeneAnnotations/hg38_refseq.bed \
 > ${project}/rseqc/${rootname}.read_dist.txt

echo rseqc
date
date

### Tool in the bbmap suite which will give you GC content info over the sample per chromosome (required bai file), and amount of coverage which is define as total number of bases covered by reads (NOT DEPTH)

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

#BedGraph generator -- generating positive and negative strand mapping; awk coverts values to negative values for negative strand mapping

#NOTE: For the current input of FStitch, it is VERY important to include "pos" somewhere in your output BedGraph file name, otherwise you will run into problems with your training data

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

#Read count correcting (rcc) -- reducing sequencing bias based on depth for visualization in IGV

module load python/2.7.14

python \
 /Users/magr0763/mappingandvis/readcountcorrectBG.py \
 ${project}/mapped/bedgraph/${rootname}.trim.BedGraph \
 ${project}/mapped/bams/sorted/${rootname}.trim.sorted.bam.flagstat \
 ${project}/mapped/bedgraph/${rootname}.trim.rcc.BedGraph

echo readcountcorrectedbedgraph
date
date

#Generate tdfs (compressed bedgraphs) -- mapped reads easily viewable in IGV

/opt/igvtools/2.3.75/igvtools toTDF \
 ${project}/mapped/bedgraph/${rootname}.trim.rcc.BedGraph \
 ${project}/mapped/tdfs/${rootname}.trim.tdf \
 /Users/magr0763/igv/genomes/sizes/hg38.chrom.sizes

echo tdf
date
date

