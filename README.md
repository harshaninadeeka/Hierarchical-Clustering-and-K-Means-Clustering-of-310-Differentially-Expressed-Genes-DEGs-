# Hierarchical-Clustering-and-K-Means-Clustering-of-310-Differentially-Expressed-Genes-DEGs-
This analysis compares two widely used unsupervised machine learning methods—Hierarchical Clustering and K-Means Clustering (K = 4)—to investigate expression patterns among 310 Differentially Expressed Genes (DEGs) identified from the RNA-Seq analysis of the airway dataset.
Both methods grouped the DEGs into four clusters, allowing the identification of genes with similar expression profiles across untreated and dexamethasone-treated samples. Although the clustering results were similar, each method provides different biological insights.

##Work flow
Import RNA seq count data,sample meta data,  Filtering low count genes, Differential gene expression analysis through DESeq2, Identification of differentially expressed significant genes based on (logFC) > 1.5 & FDR < 0.05),Hierarchical-Clustering, Heatmap, k-means clustering, Visulization of matplots, saving cluster excel sheets for downstream analysis (Enrichment Analysis).

##Visulization 
#Hierarchical-Clustering (Dendogram-cluster 01)

<img width="507" height="407" alt="Dendogram-Cluster 01" src="https://github.com/user-attachments/assets/2fd30184-238a-4eca-9353-5e2029115cc4" />




##Biological Interpretation

Both clustering approaches successfully grouped the 310 DEGs into four distinct expression patterns, indicating that the genes exhibit similar transcriptional responses following dexamethasone treatment.

Hierarchical clustering provides detailed information about the relationships among genes by organizing them into a dendrogram based on expression similarity. This makes it valuable for identifying co-expressed genes and exploring potential functional relationships.
K-Means clustering groups genes according to their average expression profiles, producing clusters with clear upregulated and downregulated expression patterns. This method simplifies the interpretation of global transcriptional responses and is particularly useful for visualization.
