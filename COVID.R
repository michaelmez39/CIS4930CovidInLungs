#Subset data.
library(data.table)
library(dplyr)
library(purrr)
library(DESeq2)
infected_raw <- fread("data/slam_inf_params.txt", check.names=TRUE)
uninfected_raw <- fread("data/slam_uninf_params.txt", check.names=TRUE)
infected_1 = infected_raw[,c(1, 3:8)]
uninfected_1 = uninfected_raw[,c(1, 3:8)]

infected_dedup <- infected_1 %>%
  group_by(gene.symbol) %>%
  summarize(
    one.hour.inf = as.integer(median(X1h.4sU.Readcount, na.rm=TRUE)),
    two.hour.inf = as.integer(median(X2h.4sU.rep1.Readcount)),
    three.hour.inf = as.integer(median(X3h.4sU.Readcount)),
    four.hour.inf = as.integer(median(X4h.4sU.Readcount))
  )

uninfected_dedup <- uninfected_1 %>%
  group_by(gene.symbol) %>%
  summarize(
    one.hour.uninf = as.integer(median(X1h.4sU.Readcount, na.rm=TRUE)),
    two.hour.uninf = as.integer(median(X2h.4sU.rep2.Readcount)),
    three.hour.uninf = as.integer(median(X3h.4sU.Readcount)),
    four.hour.uninf = as.integer(median(X4h.4sU.Readcount))
  )

cts_merged <- merge(x=infected_dedup, y=uninfected_dedup, by="gene.symbol")
cts <- cts_merged[,-1]
rownames(cts) <- cts_merged[,1]
coldata <- as.data.frame(names(cts))
coldata$condition <- c("infected", "infected", "infected", "infected", "uninfected", "uninfected", "uninfected", "uninfected")

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds) #this is doing a lot of number crunching
res <- results(dds)
vsd <- vst(dds, blind=FALSE)
subset <- head(cts[order(res$log2FoldChange),], 5000)

#Algorithm of K nearest neighbors.
kmeans(as.matrix(subset), 2, iter.max = 10, nstart = 1) #When k = 2.
#The ratio of BSS/TSS is 76.3%.
kmeans(as.matrix(subset), 3, iter.max = 10, nstart = 1) #when k = 3.
#The ratio of BSS/TSS is 86.7%.
kmeans(as.matrix(subset), 4, iter.max = 10, nstart = 1) #When k = 4.
#The ratio of BSS/TSS is 88.6%.
kmeans(as.matrix(subset), 5, iter.max = 10, nstart = 1) #When k = 5.
#The ratio of BSS/TSS is 89.4%.
kmeans(as.matrix(subset), 6, iter.max = 10, nstart = 1) #When k = 6.
##The ratio of BSS/TSS is 96.4%
#K-Means is a method that requires me to select the number of clusters (k).
#By changing the value of clusters (k), the ratio of BSS/TSS increases.

#Rerun the clustering method using 10, 100, 1000, and 10000.
subset_rev1 <- head(cts[order(res$log2FoldChange),], 10)
kmeans(as.matrix(subset_rev1), 2, iter.max = 10, nstart = 1)
#The BSS/TSS ratio is 82.2%.
subset_rev2 <- head(cts[order(res$log2FoldChange),], 100)
kmeans(as.matrix(subset_rev2), 2, iter.max = 10, nstart = 1)
#The BSS/TSS ratio is 97.8%.
subset_rev3 <- head(cts[order(res$log2FoldChange),], 1000)
kmeans(as.matrix(subset_rev3), 2, iter.max = 10, nstart = 1)
#The BSS/TSS ratio is 72.2%.
subset_rev4 <- head(cts[order(res$log2FoldChange),], 10000)
kmeans(as.matrix(subset_rev4), 2, iter.max = 10, nstart = 1)
#The BSS/TSS ratio is 71.9%.
#By increasing the number of genes, the BSS/TSS ratio increases from 82.2% to 97.8% and then decreases from 72.2% to 71.9%.

#The alluvial diagram.

library(ggplot2)
ggplot(alluvial.table, aes(y="Freq", x="Cluster") + geom_alluvium(aes(fill="Remain"), width=1/12))

#The heat-map of 5000 genes.
heatmap(as.matrix(subset), name=subset, column_km = 2)

#The chi-squared test.
chisq.test(cbind(c(4,4), c(3,5)))

#Compare different clustering results.
chisq.test(cbind(c(2,6), c(3,5)))

#Adjust all statistical test results for multiple hypothesis testing.
p.adjust()
