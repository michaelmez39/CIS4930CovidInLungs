library(data.table)
library(dplyr)
library(purrr)
library(DESeq2)
# Data Preparation
prepare_data <- function(infected_filepath, uninfected_filepath) {
    infected_raw <- fread(file=infected_filepath, check.names=TRUE)
    uninfected_raw <- fread(uninfected_filepath, check.names=TRUE)
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
    return(cts)
}
cts <- prepare_data("data/GSE162323_slam_inf_params.txt", "data/GSE162323_slam_uninf_params.txt")
coldata <- as.data.frame(names(cts))
coldata$condition <- c("infected", "infected", "infected", "infected", "uninfected", "uninfected", "uninfected", "uninfected")
coldata$batches <-  c("one", "two", "three", "four", "one", "two", "three", "four");

# GSEA Analysis
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds) #this is doing a lot of number crunching
res <- results(dds)
vsd <- vst(dds, blind=FALSE)

# Most Differentially Expressed
most.variable <- head(cts[order(res$log2FoldChange),], 10000)

# # KMeans
# c2 <- kmeans(as.matrix(subset), 3, iter.max = 10, nstart=1)

# Consensus Cluster Analysis
library(ConsensusClusterPlus)
gene.samples = c(10, 100, 1000, 10000)
for (num.genes in gene.samples) {
    ConsensusClusterPlus(as.matrix(most.variable[1:num.genes,]),
                           maxK=6,
                           reps=50,
                           pItem=0.8,
                           pFeature=1,
                           title=paste("./plots/ConsensusClusterPlus", num.genes, sep=""),
                           clusterAlg="hc",
                           distance="pearson",
                           seed=1223618388.7149,
                           plot="png"
                        )
}
library(ComplexHeatmap)
results = ConsensusClusterPlus(as.matrix(most.variable[1:5000,]),
                     maxK=6,
                     reps=50,
                     pItem=0.8,
                     pFeature=1,
                     title="./plots/ConsensusClusterPlus",
                     clusterAlg="hc",
                     distance="pearson",
                     seed=1223618388.7149,
                     plot="png"
)
Heatmap(as.matrix(most.variable), name = "Counts", 
    column_split = results[[4]]$consensusClass)