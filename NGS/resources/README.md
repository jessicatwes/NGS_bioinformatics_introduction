# NGS Resources

## Details on resources and files needed for NGS  alignment  workflow
* BBMap adapter: Adapter used with  BBMap. File 'adapters.fa' in directory is for BBMap. File was downloadad from [BBMap resource page](https://github.com/BioInfoTools/BBMap/blob/master/resources/adapters.fa).
* HISAT2 index: Used with HISAT2 during alignment to map reads to reference genome. Resource to download [HISAT2 index genome/](https://daehwankimlab.github.io/hisat2/download/).
* RefSeq reference genome: NCBI RefSeq h38 reference genome. File 'hg38_refseq.bed' is human GRCh38 reference genome. Additional reference can be found from [NCBI refseq](https://www.ncbi.nlm.nih.gov/refseq/).
* Chromosome size: Human chrome size file, 'hg38.chrom.sizes' downloaded from IGV. Additional chrom.size file can be found through [IGV](https://github.com/igvteam/igv/blob/master/genomes/sizes/hg38.chrom.sizes).
* IGVtools: Details on [igvtools toTDF](https://github.com/philres/IGV/blob/master/docs/igvtools_readme.txt).
* read correction count (rcc): readcountcorrectBG.py program available on [mapandvischainonfiji github](https://github.com/maallen3/mapandvischainonfiji). Use rcc to reduce sequencing bias based on depth for visualization in IGV.
* bamCoverage: Part of [deepTools suite](https://github.com/deeptools/deepTools). Used for generating coverage tracks (bedGraph or bigWig). Details on [bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) found on deepTools site.
