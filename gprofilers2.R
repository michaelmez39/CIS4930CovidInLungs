packages.install("dplyr")
BiocManager::install("DESeq2",update=FALSE)
install.packages("gprofilers2")
library(gprofilers2)
library(DESeq2)
multi_raw<-gost(query=list("infected"=c("gene symbol","ucsc / ensembl ID","1h 4sU Readcount","2h 4sU rep1 Readcount",
                                        "2h 4sU rep2 Readcount","3h 4sU Readcount",
                                        "4h 4sU Readcount","no 4sU Readcount","1h 4sU 0.05 quantile"),
                           "uninfected"=c("gene symbol","ucsc / ensembl ID","1h 4sU Readcount","2h 4sU rep1 Readcount",
                                          "2h 4sU rep2 Readcount","3h 4sU Readcount",
                                          "4h 4sU Readcount","no 4sU Readcount","1h 4sU 0.05 quantile"),evcodes = TRUE, multi_query = FALSE, )
                gostplot(multi_raw, capped = TRUE, interactive = TRUE)
                