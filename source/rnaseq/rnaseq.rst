.. _backbone-label:

RNAseq
==============================



edgeR
~~~~~~~~~~~~~
https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
Normalizing the data
+++++++++++++++++++++++
Read this short blog entry about normalizing RNA Seq data: http://www.rna-seqblog.com/data-analysis/which-method-should-you-use-for-normalization-of-rna-seq-data/ . edgeR normalizes by total count.

edgeR is concerned with differential expression analysis rather than with the quantification of expression levels. It is concerned with relative changes in expression levels between conditions, but not directly with estimating absolute expression levels.

The calcNormFactors() function normalizes for RNA composition by finding a set of scaling factors for the library sizes that minimize the log-fold changes between the samples for most genes. The default method for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples. We call the product of the original library size and the scaling factor the effective library size. The effective library size replaces the original library size in all downsteam analyses.


Gene set enrichment analysis 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_02_gsea.html