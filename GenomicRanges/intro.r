### R code from vignette source 'GenomicRangesIntroduction.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: biocLite (eval = FALSE)
###################################################
## source("https://bioconductor.org/biocLite.R")
## biocLite("GenomicRanges")


###################################################
### code chunk number 3: initialize
###################################################
library(GenomicRanges)


###################################################
### code chunk number 4: example-GRanges
###################################################
gr <- GRanges(
  seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
  strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  score = 1:10,
  GC = seq(1, 0, length=10))
gr
options(warn=2)

#The components of the genomic coordinates within a GRanges object can be extracted using
#the seqnames, ranges, and strand accessor functions.


###################################################
### code chunk number 5: GRanges-location-accessors
###################################################
seqnames(gr)
ranges(gr)
strand(gr)


###################################################
### code chunk number 6: granges-accessor
###################################################
#The genomic ranges can be extracted without corresponding metadata with granges
granges(gr)


###################################################
### code chunk number 7: metadataAccess
###################################################
#Annotations for these coordinates can be extracted as a DataFrame object using the mcols accessor.
mcols(gr)
mcols(gr)$score


###################################################
### code chunk number 8: setSeqLengths
###################################################
#Information about the lengths of the various sequences that the ranges are aligned to can
#also be stored in the GRanges object. So if this is data from Homo sapiens, we can set the
#values as:

seqlengths(gr) <- c(249250621, 243199373, 198022430)


###################################################
### code chunk number 9: setSeqLengths2
###################################################
seqlengths(gr)


###################################################
### code chunk number 10: names
###################################################
names(gr)
length(gr)


###################################################
### code chunk number 11: splitAppendGRanges
###################################################
sp <- split(gr, rep(1:2, each=5))
sp


###################################################
### code chunk number 12: combine
###################################################
c(sp[[1]], sp[[2]])


###################################################
### code chunk number 13: subset1
###################################################
gr[2:3]


###################################################
### code chunk number 14: subset2
###################################################
gr[2:3, "GC"]


###################################################
### code chunk number 15: assign1
###################################################
#Elements can also be assigned to the GRanges object. Here is an example where the second
#row of a GRanges object is replaced with the first row of gr.
#singles <- split(gr, names(gr))
grMod <- gr
grMod[2] <- singles[[1]]
head(grMod, n=3)


###################################################
### code chunk number 16: other
###################################################
rep(singles[[2]], times = 3)
rev(gr)
head(gr,n=2)
tail(gr,n=2)
window(gr, start=2,end=4)
gr[IRanges(start=c(2,7), end=c(3,9))]


###################################################
### code chunk number 17: IRangesStuff
###################################################
#Basic interval characteristics of GRanges objects can be extracted using the start, end,
#width, and range methods.
g <- gr[1:3]
g <- append(g, singles[[10]])
start(g)
end(g)
width(g)
range(g)


###################################################
### code chunk number 18: flank
###################################################
flank(g, 10)


###################################################
### code chunk number 19: flank2
###################################################
flank(g, 10, start=FALSE)


###################################################
### code chunk number 20: shiftAndResize
###################################################
shift(g, 5)
resize(g, 30)


###################################################
### code chunk number 21: reduce
###################################################
reduce(g)


###################################################
### code chunk number 22: gaps
###################################################
gaps(g)


###################################################
### code chunk number 23: disjoin
###################################################
disjoin(g)


###################################################
### code chunk number 24: coverage
###################################################
coverage(g)


###################################################
### code chunk number 25: intervals1
###################################################
g2 <- head(gr, n=2)
union(g, g2)
intersect(g, g2)
setdiff(g, g2)


###################################################
### code chunk number 26: intervals2
###################################################
g3 <- g[1:2]
ranges(g3[1]) <- IRanges(start=105, end=112)
punion(g2, g3)
pintersect(g2, g3)
psetdiff(g2, g3)


###################################################
### code chunk number 27: manPage (eval = FALSE)
###################################################
## ?GRanges


###################################################
### code chunk number 28: granges-methods (eval = FALSE)
###################################################
## methods(class="GRanges")


###################################################
### code chunk number 29: example-GRangesList
###################################################
gr1 <- GRanges(
  seqnames = "chr2", 
  ranges = IRanges(103, 106),
  strand = "+", 
  score = 5L, GC = 0.45)
gr2 <- GRanges(
  seqnames = c("chr1", "chr1"),
  ranges = IRanges(c(107, 113), width = 3),
  strand = c("+", "-"), 
  score = 3:4, GC = c(0.3, 0.5))
grl <- GRangesList("txA" = gr1, "txB" = gr2)
grl


###################################################
### code chunk number 30: basicGRLAccessors
###################################################
seqnames(grl)
ranges(grl)
strand(grl)


###################################################
### code chunk number 31: exceptions
###################################################
length(grl)
names(grl)
seqlengths(grl)


###################################################
### code chunk number 32: elementNROWS
###################################################
elementNROWS(grl)


###################################################
### code chunk number 33: isEmpty
###################################################
isEmpty(grl)


###################################################
### code chunk number 34: mcolsGRL
###################################################
mcols(grl) <- c("Transcript A","Transcript B")
mcols(grl)


###################################################
### code chunk number 35: mcolsGRL-unlist
###################################################
mcols(unlist(grl))


###################################################
### code chunk number 36: unlistGRL
###################################################
ul <- unlist(grl)
ul


###################################################
### code chunk number 37: pc-grl
###################################################
grl1 <- GRangesList(
  gr1 = GRanges("chr2", IRanges(3, 6)),
  gr2 = GRanges("chr1", IRanges(c(7,13), width = 3)))
grl2 <- GRangesList(
  gr1 = GRanges("chr2", IRanges(9, 12)),
  gr2 = GRanges("chr1", IRanges(c(25,38), width = 3)))

pc(grl1, grl2)

grl3 <- c(grl1, grl2)
regroup(grl3, names(grl3))


###################################################
### code chunk number 38: intOpsGRL
###################################################
start(grl)
end(grl)
width(grl)


###################################################
### code chunk number 39: List-ops
###################################################
sum(width(grl))  # sum of widths of each grl element


###################################################
### code chunk number 40: coverageGRL
###################################################
shift(grl, 20)
coverage(grl)


###################################################
### code chunk number 41: subsetGRL (eval = FALSE)
###################################################
## grl[1]
## grl[[1]]
## grl["txA"]
## grl$txB


###################################################
### code chunk number 42: subsetGRL2
###################################################
grl[1, "score"]
grl["txB", "GC"]


###################################################
### code chunk number 43: otherSubsetGRL
###################################################
rep(grl[[1]], times = 3)
rev(grl)
head(grl, n=1)
tail(grl, n=1)
window(grl, start=1, end=1)
grl[IRanges(start=2, end=2)]


###################################################
### code chunk number 44: lapply
###################################################
lapply(grl, length)
sapply(grl, length)


###################################################
### code chunk number 45: mapply
###################################################
grl2 <- shift(grl, 10)
names(grl2) <- c("shiftTxA", "shiftTxB")

mapply(c, grl, grl2)
Map(c, grl, grl2)


###################################################
### code chunk number 46: endoapply
###################################################
endoapply(grl, rev)
mendoapply(c, grl, grl2)


###################################################
### code chunk number 47: ReduceGRL
###################################################
Reduce(c, grl)


###################################################
### code chunk number 48: unlist-relist
###################################################
gr <- unlist(grl)
gr$log_score <- log(gr$score)
grl <- relist(gr, grl)
grl


###################################################
### code chunk number 49: manPage2 (eval = FALSE)
###################################################
## ?GRangesList
## methods(class="GRangesList")   # _partial_ list


###################################################
### code chunk number 50: findOverlaps
###################################################
mtch <- findOverlaps(gr, grl)
as.matrix(mtch)


###################################################
### code chunk number 51: countOL
###################################################
countOverlaps(gr, grl)


###################################################
### code chunk number 52: subsetByOverlaps
###################################################
subsetByOverlaps(gr,grl)


###################################################
### code chunk number 53: select-first
###################################################
findOverlaps(gr, grl, select="first")
findOverlaps(grl, gr, select="first")


###################################################
### code chunk number 54: SessionInfo
###################################################
sessionInfo()

