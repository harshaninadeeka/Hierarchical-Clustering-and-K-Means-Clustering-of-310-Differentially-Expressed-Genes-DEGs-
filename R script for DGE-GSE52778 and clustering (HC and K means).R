
#Analysis differential gene expression for data set GSE52778

# script to get data from airway package

#load packages

library(airway)
library(dplyr)
library(tidyverse)
library(GEOquery)
library(tidyverse)


data(airway)
airway

sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)

# Step 1: preparing count data ----------------

# read in counts data
countTable <- read.csv('counts_data.csv')
head(countTable)


# read in sample info
sampleData <- read.csv('sample_info.csv')
head(sampleData)

# making sure the row names in sampleData matches to column names in counts_data
all(colnames(countTable) %in% rownames(sampleData))

# are they in the same order?
all(colnames(countTable) == rownames(sampleData))


#Define sample groups
treatment<-factor (sampleData$dexamethasone)
treatment
pairs<-factor(sampleData$cellLine)
pairs


#Filter low counts
library(edgeR)
meanLog2CPM <- rowMeans(log2(cpm(countTable) + 1))
sum(meanLog2CPM <= 1)

#Histogram
hist(meanLog2CPM)

#save histogram before filtering
ggsave(file="Histogram-GSE52778-Before filtering.pdf",width = 10, height = 8)

#keep log2cpm values greater than 1 for downstream analysis
countsTable <- countTable[meanLog2CPM > 1, ]
dim(countsTable)

#Check distribution
library(DESeq2)
dds <- DESeqDataSetFromMatrix(
  countData = countsTable,
  design = ~ treatment + pairs,
  colData = data.frame(
    treatment = treatment,
    pairs = pairs
  ))
normCounts <- rlog(dds)

#Histogram (To see negative binomial distribution)
hist(assay(normCounts), main = "GSE52778")
ggsave(file="Histogram-GSE52778-Negative Binomial.pdf",width = 10, height = 8)


#Heatmap (to visualize the correlations between samples)
library(pheatmap)
sampleDist <- cor(assay(normCounts), method = "spearman")
p<-pheatmap(sampleDist, clustering_distance_rows = as.dist(1 - sampleDist),
         clustering_distance_cols = as.dist(1 - sampleDist))

#PCA
library(DESeq2)
pcaPlot <- plotPCA(normCounts, intgroup = c("treatment"))
pcaPlot

#Prepare count data for analysis
dge <- DGEList(countsTable)

#Calculate normalization factors
dge <- calcNormFactors(dge)

#Define design matrix
designDF <- data.frame(Treatment = treatment, Pair = pairs)
design <- model.matrix(~ 0 + Treatment + Pair, data = designDF) 
design


#Fit GLM
dge <- estimateDisp(dge, design, robust = TRUE)
fit <- glmQLFit(dge, design, robust = TRUE)
fit

#Define contrast matrix
contrastMat <- makeContrasts(TreatmentVscontrol = Treatmenttreated - Treatmentuntreated,
                             levels = design)


#Statistical testing
qlfRes <- glmQLFTest(fit, contrast = contrastMat)
topRes <- topTags(qlfRes, n = nrow(fit$counts))
topRes <- subset(topRes$table, abs(logFC) > 1.5 & FDR < 0.05)
print(head(topRes))
print(nrow(topRes))
topRes


# Get the list of 310 DEG gene names (assuming rownames are Entrez/Symbols)
sig_gene_ids <- rownames(topRes)

# Calculate log2 CPM (using your edgeR DGEList object 'dge' or 'y')

logCPM <- cpm(dge, log = TRUE, prior.count = 2)

# Subset the log2 CPM matrix to ONLY your 310 significant DEGs
clust_data <- logCPM[sig_gene_ids, ]


#delete missing values
clust_data <- na.omit(clust_data)
#if the dataset contains any missing values

#plot a histogram
hist(clust_data) 
ggsave(file="Histogram-GSE52778.pdf",width = 10, height = 8)

# standardize the dataset
clust_data <- scale(clust_data)
hist(clust_data)
dim(clust_data)
ggsave(file="Histogram-GSE52778-After scaling.pdf",width = 10, height = 8)

#Hierarchical clustering Euclidian distance
#calculate distance matrix

d <- dist(clust_data, method="euclidean") #by genes

#  Hierarchical Clustering
#first set the margins of for plot by using the par() function
par(mar=c(7,7,7,7))
hc <- hclust(d, method="average") 

##convert the hc to a dendrogram object using dendextend package
library(dendextend) 

hc<-as.dendrogram(hc)



#plot the clustering tree
plot(hc, main="Airway data-GSE52778 (Differentially Expressed Genes =310", cex=0.2,cex.lab=0.2 ) # display dendrogram


#cut the tree and plot the cluster profiles
clusters <- cutree(hc, k=4)

# draw dendogram with red borders around the 4 clusters 
rect.dendrogram(hc , k=4, border="red")



#Plot 4 clusters
cl1 <- clust_data[which(clusters==1),]
matplot(t(cl1),type="l",xaxt = "n", main="Dendrogram- Cluster 01")
axis(1,at=1:8,labels=colnames(clust_data),las=2)

cl2 <- clust_data[which(clusters==2),]
matplot(t(cl2),type="l",xaxt = "n", main="Dendrogram- Cluster 02")
axis(1,at=1:8,labels=colnames(clust_data),las=2)


cl3 <- clust_data[which(clusters==3),]
matplot(t(cl3),type="l",xaxt = "n", main="Dendrogram- Cluter 03")
axis(1,at=1:8,labels=colnames(clust_data),las=2)

cl4 <- clust_data[which(clusters==4),]
matplot(t(cl4),type="l",xaxt = "n", main="Dendrogram- Cluster 04")
axis(1,at=1:8,labels=colnames(clust_data),las=2)




#check how many clusters got?
max(clusters)
#save content of one cluster to file
library(openxlsx)
write.xlsx(cl1,file="cl1-Dendrogram.xlsx")
write.xlsx(cl2,file="cl2-Dendrogram.xlsx")
write.xlsx(cl3,file="cl3-Dendrogram.xlsx")
write.xlsx(cl4,file="cl4-Dendrogram.xlsx")




#visualization using heatmap

heatmap(
  as.matrix(clust_data),
  Rowv = TRUE,        # Explicitly clusters rows (Genes)
  Colv = TRUE,        # Clusters columns (Samples) 
  scale = "row",      # CRITICAL: Standardizes expression per gene
  margins = c(5, 5),
  cexCol = 0.7)


#K-means clustering of genes

#To produce a K-means clustering with 4 clusters type:
k <- 4
km <- kmeans(clust_data, k, iter.max=10000)


#To extract the genes in cluster nr 1 you type:
clk1 <- clust_data[which(km$cluster==1),]

#Visualize the genes in a cluster

#Use matplot() to plot gene profiles in one cluster
#Transpose the matrix before plotting
clk1 <- t(clust_data[which(km$cluster==1),])
par(mar=c(4,4,4,4))#set the margins
matplot(clk1,type="l") 
#delete and type in ””

#Adding more informative labels on X axis
matplot(clk1,type="l",xaxt = "n")
#add more informative axis
axis(1,at=1:8,labels=colnames(clust_data),las=2)

#For wider bottom margins
par(mar=c(7,4,4,2))
matplot(clk1,type="l",xaxt = "n", main="K-means clustering (K=4), Cluster 1")
axis(1,at=1:8,labels=colnames(clust_data),las=2)


clk2 <- t(clust_data[which(km$cluster==2),])
par(mar=c(7,4,4,2))
matplot(clk2,type="l",xaxt = "n", main="K-means clustering (K=4), Cluster 2")
axis(1,at=1:8,labels=colnames(clust_data),las=2)

clk3 <- t(clust_data[which(km$cluster==3),])
par(mar=c(7,4,4,2))
matplot(clk3,type="l",xaxt = "n", main="K-means clustering (K=4), Cluster 3")
axis(1,at=1:8,labels=colnames(clust_data),las=2)

clk4 <- t(clust_data[which(km$cluster==4),])
par(mar=c(7,4,4,2))
matplot(clk4,type="l",xaxt = "n", main="K-means clustering (K=4), Cluster 4")
axis(1,at=1:8,labels=colnames(clust_data),las=2)




#save cluster files
library(openxlsx)
write.xlsx(clk1,file="clk1-K means clustering-cluster1.xlsx")
write.xlsx(clk2,file="clk1-K means clustering-cluster2.xlsx")
write.xlsx(clk3,file="clk1-K means clustering-cluster3.xlsx")
write.xlsx(clk4,file="clk1-K means clustering-cluster4.xlsx")





















