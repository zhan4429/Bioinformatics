.. _backbone-label:

Ensembl genomes
==============================
http://genomespot.blogspot.com/2015/06/mapping-ngs-data-which-genome-version.html

When you look at all the genome files available from Ensembl. You are presented with a bunch of options. Which one is the best to use/download?

You have a combination of choices.

First part options:

- dna_sm - Repeats soft-masked (converts repeat nucleotides to lowercase)
- dna_rm - Repeats masked (converts repeats to to N's)
- dna - No masking

Second part options:

- toplevel - Includes haplotype information (not sure how aligners deal with this)
- primary_assembly - Single reference base per position

Right now I usually use a non-masked primary assembly for analysis, so in the case of humans: Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

Does this make sense for standard RNA-Seq, ChIP-Seq, ATAC-Seq, CLIP-Seq, scRNA-Seq, etc... ?

In what cases would I prefer other genomes? Which tools/aligners take into account softmasked repeat regions?

Generally, you should use the soft-masked or unmasked primary assembly. Cross-species whole-genome aligners, especially older ones, do need to know soft-masked regions; otherwise they can be impractically slow for mammalian genomes. Modern read aligners are designed to work with repeats efficiently and therefore they don't need to see the soft mask.

For GRCh38, though, I would recommend to use the official build at GRC FTP. Most people will probably choose "no_alt_analysis_set". Using the Ensembl version is discouraged due to its chromosome naming. We more often use "chr1" instead of "1" for GRCh38. At one point, Ensembl actually agreed to use "chr1" as well, but didn't make that happen due to technical issues, I guess.

As to alternate haplotypes, most aligners can't work with them; no variant callers can take the advantage of these sequences, either. When you align to a reference genome containing haplotypes with an aligner not supporting these extra sequences, you will get poor mapping results.

There's rarely a good reason to use a hard-masked genome (sometimes for blast, but that's it). For that reason, we use soft-masked genomes, which only have the benefit of showing roughly where repeats are (we never make use of this for our *-seq experiments, but it's there in case we ever want to).

For primary vs. toplevel, very few aligners can properly handle additional haplotypes. If you happen to be using BWA, then the toplevel assembly would benefit you. For STAR/hisat2/bowtie2/BBmap/etc. the haplotypes will just cause you problems due to increasing multimapper rates incorrectly. Note that none of these actually use soft-masking.

