.. _backbone-label:

Small RNA-seq
==============================
Small noncoding RNAs act in gene silencing and post-transcriptional regulation of gene expression. Small RNA sequencing (RNA-Seq) is a technique to isolate and sequence small RNA species, such as microRNAs (miRNAs). Small RNA-Seq can query thousands of small RNA and miRNA sequences with unprecedented sensitivity and dynamic range. 
With small RNA-Seq you can discover novel miRNAs and other small noncoding RNAs, and examine the differential expression of all small RNAs in any sample. You can characterize variations such as isomiRs with single-base resolution, as well as analyze any small RNA or miRNA without prior sequence or secondary structure information.

Sequencing steps
~~~~~~~~~~~~~~~~~~~
- Total RNA isolation
- Small RNA library construction
- Deep sequencing

Bioinformatics analysis
~~~~~~~~~~~~~~~~~~~~~~~~~
https://blog.csdn.net/Candle_light/article/details/104320026
``mirDeep2`` provides a perl script ``extract_miRNA.pl`` to extract target organism's microRNA fasta sequences::

	extract_miRNA.pl mature.fa hsa > mature_hsa.fa
	extract_miRNA.pl hairpin.fa hsa > hairpin_hsa.fa

you only need to change the tree letter code, in this case ``hsa`` stand for human.

You can also convert the file to DNA with::

	perl rna2dna.pl mature_hsa.fa > mature.hsa.dna.fa

It´s the same for the hairpin sequences::

	perl rna2dna.pl hairpin_hsa.fa > hairpin.hsa.dna.fa

