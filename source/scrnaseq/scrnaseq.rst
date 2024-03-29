scRNA-seq
==============================

Analyze single-cell RNA-Seq data in R 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
- Count matrix
- QC and filtering
- Normalization
- Identify highly variable genes
- Scale data
- Linear dimensinality reduction(PCA)
- Clustering
- Non-linear dimensinality reduction (UMAP/t-SNE)


Downstream Analysis
~~~~~~~~~~~~~~~~~~~~~
- Cluster identification
- Differential gene expression between clusters/differential chromatiin accessibility analysis
- Inferring trajectories/lineage

PCA
+++++++
I'm working on a project involving analyzing scRNA-seq data. A large part of the project involves clustering cells, identifying DE genes between clusters, pathway analysis of the DE genes, etc. To do the analysis, I am planning to use Seurat. By default, Seurat uses the graph-based Louvain algorithm to cluster cells. So that would seem to indicate that it is important that the 2D embedding generated by t-SNE or UMAP is as accurate as possible so that the clusters are also maximally accurate.

Prior to doing t-SNE or UMAP, Seurat's vignettes recommend doing PCA to perform an initial reduction in the dimensionality of the input dataset while still preserving most of the important data structure. Seurat is definitely not the only pipeline to do this; it seems to me that most analysis pipelines use PCA prior to t-SNE / UMAP basically like Seurat does. However, it also seems to me that ICA is generally better at dividing cells based on the activation of gene modules than PCA. This seems to me to make sense in principle - i.e. gene modules behave more like independent gene combinations (as modeled by ICA) than orthogonal gene combinations (as modeled by PCA) - and also in practice - i.e. I've read a few papers presenting empirical evidence that ICA is better than PCA for differentiating cells based on gene module activation. Assuming this is correct, would it make more sense to use ICA rather than PCA to do the pre-t-SNE / UMAP dimensionality reduction? Or is there a compelling reason that most people seem to use PCA for this that I am simply unaware of?

Answer:
t-SNE scales poorly to datasets in tens of thousands data points and more than 30-50 features, and struggles mightily with millions of data points and 50+ features. I have done UMAP easily with 2-5 million data points and 200+ features, so you may not need any initial dimensionality reduction with UMAP.
I don't think it matters much whether you use ICA or PCA for initial reduction, as long as it is followed by t-SNE. If PCA with 30 or so PCs explains >80-90% of variance, that should be good enough. If you provide exact data dimensions and the variance that is explained by each PC up to 50, I may be able able to offer better advice.

I think the main reason dimensionality reduction is performed before t-SNE is because of the poor performance of t-SNE with high dimensional data (this could be due to the difficulty in finding the right parameters in such situation). UMAP seems better in this respect but this is anecdotal. This paper compares t-SNE and UMAP on single cell data. I would also add that UMAP can do metric learning: it can be used to learn a projection that best separates annotated samples then used to project unannotated samples in this space. This can be quite useful if one has annotated samples.

The choice is typically guided by some assumptions about the data and what the goal of the transformation is. PCA assumes that the only relevant components to explain the variability in the data are the uncorrelated ones. This generally works well for multivariate Gaussian distributions because in this case, uncorrelated also means statistically independent. ICA assumes the data to be generated by statistically independent non-Gaussian sources (it's a form of blind source separation like NMF) and tries to identify them by minimizing statistical dependence of the components. Why use ICA for non-Gaussian sources? Because in this case, uncorrelated doesn't imply statistical independence so PCA wouldn't necessarily recover the desired components. The downside of ICA is that there's no ranking of the components (i.e. there's no relationship between an ICA with k-1 components and one with k components unlike in PCA).


