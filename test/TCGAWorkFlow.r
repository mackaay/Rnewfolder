
library(TCGAbiolinks)
# Obs: The data in the legacy database has been aligned to hg19
query.met.gbm <- GDCquery(project = "TCGA-GBM", 
                          legacy = TRUE,
                          data.category = "DNA methylation",
                          platform = "Illumina Human Methylation 450", 
                          barcode = c("TCGA-76-4926-01B-01D-1481-05", "TCGA-28-5211-01C-11D-1844-05"))
GDCdownload(query.met.gbm)

met.gbm.450 <- GDCprepare(query = query.met.gbm,
                          save = TRUE, 
                          save.filename = "gbmDNAmet450k.rda",
                          summarizedExperiment = TRUE)

query.met.lgg <- GDCquery(project = "TCGA-LGG", 
                          legacy = TRUE,
                          data.category = "DNA methylation",
                          platform = "Illumina Human Methylation 450",
                          barcode = c("TCGA-HT-7879-01A-11D-2399-05", "TCGA-HT-8113-01A-11D-2399-05"))
GDCdownload(query.met.lgg)
met.lgg.450 <- GDCprepare(query = query.met.lgg,
                          save = TRUE, 
                          save.filename = "lggDNAmet450k.rda",
                          summarizedExperiment = TRUE)
met.gbm.lgg <- SummarizedExperiment::cbind(met.lgg.450, met.gbm.450)


query.exp.lgg <- GDCquery(project = "TCGA-LGG", 
                          legacy = TRUE,
                          data.category = "Gene expression",
                          data.type = "Gene expression quantification",
                          platform = "Illumina HiSeq", 
                          file.type = "results",
                          sample.type = "Primary solid Tumor")
GDCdownload(query.exp.lgg)
exp.lgg <- GDCprepare(query = query.exp.lgg, save = TRUE, save.filename = "lggExp.rda")

query.exp.gbm <- GDCquery(project = "TCGA-GBM", 
                          legacy = TRUE,
                          data.category = "Gene expression",
                          data.type = "Gene expression quantification",
                          platform = "Illumina HiSeq", 
                          file.type = "results",
                          sample.type = "Primary solid Tumor")
GDCdownload(query.exp.gbm)
exp.gbm <- GDCprepare(query = query.exp.gbm, save = TRUE, save.filename = "gbmExp.rda")
exp.gbm.lgg <- SummarizedExperiment::cbind(exp.lgg, exp.gbm)


#-----------------------------------------------------------------------------
#                   Data.category: Copy number variation aligned to hg38
#-----------------------------------------------------------------------------
query <- GDCquery(project = "TCGA-ACC",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment",
                  barcode = c( "TCGA-OR-A5KU-01A-11D-A29H-01", "TCGA-OR-A5JK-01A-11D-A29H-01"))
GDCdownload(query)
data <- GDCprepare(query)

query <- GDCquery("TCGA-ACC",
                  "Copy Number Variation",
                  data.type = "Masked Copy Number Segment",
                  sample.type = c("Primary solid Tumor")) # see the barcodes with getResults(query)$cases
GDCdownload(query)
data <- GDCprepare(query)
