# Hierarchical-Clustering-and-K-Means-Clustering-of-310-Differentially-Expressed-Genes-DEGs-
This analysis compares two widely used unsupervised machine learning methods—Hierarchical Clustering and K-Means Clustering (K = 4)—to investigate expression patterns among 310 Differentially Expressed Genes (DEGs) identified from the RNA-Seq analysis of the airway dataset.
Both methods grouped the DEGs into four clusters, allowing the identification of genes with similar expression profiles across untreated and dexamethasone-treated samples. Although the clustering results were similar, each method provides different biological insights.

##Work flow
Import RNA seq count data,sample meta data,  Filtering low count genes, Differential gene expression analysis through DESeq2, Identification of differentially expressed significant genes based on (logFC) > 1.5 & FDR < 0.05),Hierarchical-Clustering, Heatmap, k-means clustering, Visulization of matplots, saving cluster excel sheets for downstream analysis (Enrichment Analysis).

##Visulization 
#Hierarchical-Clustering (Dendogram-cluster 01)
<img width="631" height="507" alt="Dendogram-Cluster 01" src="https://github.com/user-attachments/assets/bccbafff-e297-4df5-8005-1f6e90f57b73" />
#Hierarchical-Clustering (Dendogram-cluster 02)
<img width="631" height="507" alt="Dendogram-Cluster 02" src="https://github.com/user-attachments/assets/fd9ab8f9-b0db-4f16-ab88-bcab8e49316a" />

#Hierarchical-Clustering (Dendogram-cluster 03)
<img width="631" height="507" alt="Dendogram-Cluster 03" src="https://github.com/user-attachments/assets/0382822f-bc48-4f3d-bf0a-3610b8a79fa7" />



#K means clustering (cluster 01)
<img width="631" height="507" alt="K means clustering-Cluster 01" src="https://github.com/user-attachments/assets/112d232b-f23e-4364-8a63-5520666fe77a" />

#K means clustering (cluster 02)
<img width="631" height="507" alt="K means clustering-Cluster 02" src="https://github.com/user-attachments/assets/a54fa76a-5ea5-4066-931e-c8caab2415f3" />

#K means clustering (cluster 03)
<img width="631" height="507" alt="K means clustering-Cluster 03" src="https://github.com/user-attachments/assets/cc6b0a6b-7dec-497b-b37b-f494989d90ff" />

#K means clustering (cluster 04)
<img width="631" height="507" alt="K means clustering-Cluster 04" src="https://github.com/user-attachments/assets/a7063840-2c99-4553-9594-673350ca8e54" />

##Biological Interpretation

Both clustering approaches successfully grouped the 310 DEGs into four distinct expression patterns, indicating that the genes exhibit similar transcriptional responses following dexamethasone treatment.

Hierarchical clustering provides detailed information about the relationships among genes by organizing them into a dendrogram based on expression similarity. This makes it valuable for identifying co-expressed genes and exploring potential functional relationships.
K-Means clustering groups genes according to their average expression profiles, producing clusters with clear upregulated and downregulated expression patterns. This method simplifies the interpretation of global transcriptional responses and is particularly useful for visualization.
