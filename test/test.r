source("https://bioconductor.org/biocLite.R")
biocLite("methylationArrayAnalysis")

source("https://bioconductor.org/biocLite.R")
biocLite("TCGAWorkflow")


library(ggplot2)

library(limma)
install.packages("ggplot2")
source("https://bioconductor.org/biocLite.R")
biocLite("limma")
install.packages("colorspace")


#source("http://bioconductor.org/workflows.R")
#workflowInstall("methylationArrayAnalysis")
#source("http://bioconductor.org/workflows.R")
#workflowInstall("TCGAWorkflow")

browseVignettes("TCGAWorkflow")
library(TCGAWorkflowData)
library(DT)

version


remove.packages("BiocInstaller", lib=.libPaths())
source("https://bioconductor.org/biocLite.R")
biocValid()
biocLite("limma")

library(BiocInstaller)
biocValid()  
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("GenomicRanges")
library(GenomicRanges)
biocLite("SummarizedExperiment")
library(SummarizedExperiment)




ip <- as.data.frame(installed.packages())
head(ip)
dim(ip)
# if you use MRO, make sure that no packages in this library will be removed
ip <- subset(ip, !grepl("MRO", ip$LibPath))
# we don't want to remove base or recommended packages either\
ip <- ip[!(ip[,"Priority"] %in% c("base", "recommended")),]
# determine the library where the packages are installed
path.lib <- unique(ip$LibPath)
# create a vector with all the names of the packages you want to remove
pkgs.to.remove <- ip[,1]
head(pkgs.to.remove)
# remove the packages
sapply(pkgs.to.remove, remove.packages, lib = path.lib)


library(readr)

library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)

#ggfortify#
library(ggplot2)
library(rlang)
library(devtools)
library(digest)
if (!require("devtools")) install.packages("devtools")
install_github('sinhrks/ggfortify')
library(ggfortify)

source("https://bioconductor.org/biocLite.R")
biocLite("DOSE") 
biocLite("topGO")
biocLite("clusterProfiler")
biocLite("pathview")
biocLite("org.Mm.eg.db")
library(DOSE)
library(org.Mm.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")
library(clusterProfiler)
biocLite("mygene")
library(mygene)

source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
library(GenomicRanges)
biocLite("SummarizedExperiment")
