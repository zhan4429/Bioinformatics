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
https://blog.csdn.net/weixin_43569478/article/details/108079336?utm_medium=distribute.pc_relevant.none-task-blog-2~default~baidujs_baidulandingword~default-1-108079336-blog-104320026.235^v29^pc_relevant_default_base3&spm=1001.2101.3001.4242.2&utm_relevant_index=4


Remove adaptors
+++++++++++++++++
Usually you will have gotten a small RNA sequencing data file from a collaborator that wants you to analyze the data file. Before you can start with any kind of analysis you should either already know the small RNA sequencing adapter that was used for the sequencing of the sample or ask your collaborator to sent it to you. If you don't clip the adapter then the majority of the reads having an adapter are likely not to be aligned to anywhere. 

Once you know the adapter sequence you should do a simple check to see how many of your sequences contain the adapter. This you can do by typing::

	grep -c TGGAATTC example_small_rna_file.fastq

where 'TGGAATTC' are the first 8 nucleotides of the adapter that has been used for this sample. Replace it with your own sequencing adapter. MicroRNAs have a mean length of 22 nucleotides in animals so if you have sequenced one of those it will likely have the sequencing adapter attached to it. If the resulting number of sequences with an adapter is around 70% of the number of your input sequences the data set can be considered as reasonably good. Note: In case that only adapters have been sequenced predominantly you will also get a high number which is obviously not good. If you only get 10% of sequences with an adpater then likely something went wrong during the sequnecing library preparation or your sample doesn't contain too many small RNAs. For novel miRNA prediction we need to map the reads against a reference database which has to be indexed by bowtie 1. For this we take a reference database file, lets call it refdb.fa (This can be a genome file or simply a file with scaffolds) and build a bowtie index by typing::

	bowtie-build refdb.fa refdb.fa

The first argument is here the actual file to index, the second argument is the prefix for the bowtie index files. You can name it differently but for ease of use I use the same name as my reference database file. Depending on the input file size it can take several hours (for the human genome for example) to be indexed. However, the bowtie website has already some prebuild index files for download. If you decide to download index files you will also need to download the fasta file with which the index was build. Otherwise the results in the miRDeep2 prediction will be not reliable.

Data preprocessing for novel miRNA prediction
+++++++++++++++++++++++++++++++++++++++++++++++++++
Since the miRDeep2 package was designed as a complete solution for miRNA prediction and quantification it also contains data preprocessing routines that will also clip the sequencing adapter. The main function of the mapping module is the mapping of the preprocessed reads file to reference database. The reference database is typically an annotated genome sequence but can also be simply a scaffold assembly if no genome is available. The scaffolds itself however should be at least 200 nucleotides long so that a sane miRNA precursor plus some flanking region fits into it. Apart from clipping adapters the module does sanity checks on your sequencing reads and also collapse read sequences to reduce the file size which will save computing time::

	mapper.pl example_small_rna_file.fastq -e -h -i -j -k TGGAATTC -l 18 -m -p refdb.fa -s /reads_collapsed.fa -t reads_vs_refdb.arf -v -o 4


What does this command do? The first argument needs to be your sequencing file. Typically, this will be a fastq file. The format of the fastq file is designated by specifying option '-e'. If your file is in fasta format already you specify option '-c' instead. If your reads file is not in fasta format you need to specify option '-h' which advises the mapper module to parse your file to fasta format. Option -i will convert RNA to DNA and option '-j' will remove sequences that contain characters other than ACGTN. 

Now comes the actual adapter clipping which is only done if a adapter sequence is given by option -k. Only the first 6 nucleotides of this sequence will be used to search for an exact match in the sequencing reads. Option '-m' will collapse the reads to remove redundancy and decrease the file size. A sequnecing read seen 10 times in your raw file will occur only once in the collapsed file and have a _x10 in its identifier. 

After that the reads will be mapped to the given refence genome which index file was specified by option '-p'. Option -s indicates the preprocessed read file name which is output by the mapper module and option -t is the file name of the read mappings to the reference database ('refdb.fa') in miRDeep2's arf format. A mapping file in arf format can be easily obtained from a standard bowtie 1 output file (This is NOT in 'sam' format but a proprietary bowtie text file format) by typing::

	convert_bowtie_output.pl reads_vs_refdb.bwt > reads_vs_refdb.arf

However, if you used the mapping module then the mapped output file is already in arf format.

Identification of known and novel miRNAs
++++++++++++++++++++++++++++++++++++++++++
For predicting novel miRNAs the miRDeep2 module from the package is called with a collapsed reads file and a reference genome file in fasta format. For better prediction results reference files of miRNAs and related miRNAs should be given since miRDeep2 considers predicted miRNAs with conserved seeds in other species more reliable that miRNAs with non-conserved seeds::

	miRDeep2.pl reads_collapsed.fa refdb.fa reads_vs_refdb.arf mature_ref.fa mature_other.fa /hairpin_ref.fa -t hsa 2>report.log


``mirDeep2`` provides a perl script ``extract_miRNA.pl`` to extract target organism's microRNA fasta sequences::

	extract_miRNA.pl mature.fa hsa > mature_hsa.fa
	extract_miRNA.pl hairpin.fa hsa > hairpin_hsa.fa

you only need to change the tree letter code, in this case ``hsa`` stand for human.

You can also convert the file to DNA with::

	perl rna2dna.pl mature_hsa.fa > mature.hsa.dna.fa

It´s the same for the hairpin sequences::

	perl rna2dna.pl hairpin_hsa.fa > hairpin.hsa.dna.fa

``miRdeep2`` can also extract sequcences from multiple organisms::

	extract_miRNA.pl mature.fa mmu,ptr > mature_mmu_ptr.fa

