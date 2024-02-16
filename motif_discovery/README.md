# Motif discovery

Motif discovery allows us to identify sequences (motifs) that are found in our data and have a biological significance within the larger DNA, RNA , or protein sequences. The significance can be indicative of regulatory elements such as transcription binding sites that play a crucial role in gene expression and other cellular processes.

Using MEME (Multiple EM for Motif Elicitation) and TOMTOM, we can discover motifs in two steps
<br>Step 1: MEME is used to discover motifs in a set of sequences. It uses the expectation maximization (EM) algorithm to identify motifs that are statistically significant and have potential biological relevance.
<br>Step 2: TOMTOM takes the output of motifs discovered by MEME (or other motif discovery tools) and compare it against a database of known motifs (JASPAR, TRANSFAC) to identify potential matches.  

## Resources
[MEME](https://meme-suite.org/meme/): Suite of tools for motif-based sequence analysis
<br>[TOMTOM](https://meme-suite.org/meme/tools/tomtom): A motif comparison tool; part of the MEME suite
<br>[HOMER](http://homer.ucsd.edu/homer/motif/): Alternative to motif discovery for DNA sequences
