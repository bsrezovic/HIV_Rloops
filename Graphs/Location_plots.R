library(data.table)
#library(xlsx)    # failing to connect to java lol
library(rtracklayer)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(AnnotationHub)
library(plotgardener)


### Import data:

# use import.bed and import.bedGraph (for chromatin marks and Drip): s177fw, s177rv, k36, k27ac, k27me3
# consensus peaks are data.table in format seqnames, start, end, strand: cons
# IS same as peaks: is 


### Set plotting parameters:

# locations of interest for plotting:
# (gene name, chromosome, start, end)
# this is the main location thingy with gene and stuff
location <- c("GRIN2D", "chr19", 48835613,49017192)

# colors, in order of plotting (Drip, consensus, H3K36me3, H3K27ac, H3K27me3, IS):
colorscale <- c("#C43E2E", "#4C74B1", "#137F17", "#E67017", "#87187C", "#525252", "#42A6AD")

# reasonable cutoffs, in order of plotting (Drip+, Drip-, H3K36me3, H3K27ac, H3K27me3):
# (vector of max y values for each track - needs adjustment depending on data)
lims <- c(300, 300, 40, 250, 30)  # changed to 300 from 500 upon request
# dimensions and genomic location for tracks:
params_d <- 
  pgParams(chrom = location[2], chromstart = as.numeric(location[3]), chromend = as.numeric(location[4]), assembly = "hg19", 
           x = 0.3, width = 5.6, default.units = "inches")

# signal track data ranges:
s177fw_range <- pgParams(range = c(0, lims[1]), assembly = "hg19")
s177rv_range <- pgParams(range = c(0, lims[2]), assembly = "hg19")
k36_range <- pgParams(range = c(0, lims[3]), assembly = "hg19")
k27ac_range <- pgParams(range = c(0, lims[4]), assembly = "hg19")
k27me3_range <- pgParams(range = c(0, lims[5]), assembly = "hg19")

## load consensus peaks
cons <-readRDS("C:/Users/38598/Desktop/hivint/MajasCode/Rloops/R-loops-code/Consensus-peaks/consensusPeaks_withMeanSignal.RDS")
# plot forward and reverse DRIPc sample, both with label and y-axis:
s177fw <- import.bedGraph("C:/Users/38598/Desktop/hivint/plotgardnerfw19.bedgraph")
s177rw <- import.bedGraph("C:/Users/38598/Desktop/hivint/plotgardnerrw19.bedgraph")
# load chromatin data
k36 <- import.bedGraph("C:/Users/38598/Desktop/hivint/H3K36me19.bdg")
k27ac <- import.bedGraph("C:/Users/38598/Desktop/hivint/H3K27aC19.bdg")
k27me3 <- import.bedGraph("C:/Users/38598/Desktop/hivint/H3K27me19.bdg")

#making sure the table is formatted just for this
cons$name <- NULL
cons$score <- NULL
cons <- cons[, c(1:4)]
cons <- data.table(cons)
colnames(cons) <- c("chrom","start","end","strand")

##load insertion sites (there should be one insertion among these)
is <- combined # this is from NewGraphs.R, containing AQR and NTC insertion sites from both runs
is <- is[,c(3,5,6,10)]
colnames(is) <- c("chrom","start","end","strand")
is$end <- is$start + is$end  # this was isl_widt orginally, just fixing after renaming
#fixing strand
is$strand <- as.character(is$strand)
is[is$strand=="-1",]$strand <- "-"
is[is$strand=="1",]$strand <- "+"

# load consensus insertion sites from all the paperinos
is2 <- prevint   # prevint is also from New_graphs.R
colnames(is2) <- c("chrom","start","end","source","score")
is2$end <- is2$end +500  # thickening them up just so we have something to look at?
# left margin for labels (position on x, actually):
labmarg <- 0.4

# data track heights (plusd 0.1 internal margin for each):
tr_h <- c(rep(0.9, 2), 0.1, rep(0.2, 2), rep(0.9, 3), 0.1, 0.2,0.1, 0.2, 0.9,0.9)
names(tr_h) <- c(rep("DRIP", 2), "peaks label", rep("consensus peaks", 2), rep("chromatin marks", 3),
                 "IS label", "IS","IS label2", "IS2", "genes","genome")
# you may need to adjust, depending on data

# note: if you don't need IS label to be in its own track, remove the last 0.1 from tr_h (from names as well) and comment out the
# counter immediately after label definition

# final plot dimensions in inches (add a margin to plotgardener page) and resolution for jpeg:
w <- 6 + 0.5
h <- ceiling(sum(tr_h) + length(tr_h)*0.1) + 0.5
res <- 200



### Plot everything:





# create a page and hide guide lines immediately:
pageCreate(width = 6, height = ceiling(sum(tr_h) + length(tr_h)*0.1), default.units = "inches")

i <- 1   # counter for positions - used to help place each track on appropriate spot

s177fw_tr <- plotSignal(data = s177fw, params = c(params_d, s177fw_range),
                        fill = colorscale[1], linecolor = colorscale[1],
                        y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1),
                        height = tr_h[i], just = c("left","bottom"))
plotText(label = "DRIPc +", fontcolor = colorscale[1], fontsize = 8, fontface = "bold",
         x = labmarg, y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1 - tr_h[i]),
         just = c("left","top"), default.units = "inches")   # plot label
annoYaxis(plot = s177fw_tr, at = seq(0, lims[1], length.out=5), axisLine = TRUE, fontsize = 8)   # plot y-axis
i <- i+1

s177rv_tr <- plotSignal(data = s177rw, params = c(params_d, s177rv_range),
                        fill = colorscale[1], linecolor = colorscale[1],
                        y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1),
                        height = tr_h[i], just = c("left","bottom"))
plotText(label = "DRIPc -", fontcolor = colorscale[1], fontsize = 8, fontface = "bold",
         x = labmarg, y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1 - tr_h[i]),
         just = c("left","top"), default.units = "inches")
annoYaxis(plot = s177rv_tr, at = seq(0, lims[2], length.out=5), axisLine = TRUE, fontsize = 8)
i <- i+1


# plot consensus peaks (each strand separately because the coloring is stupid otherwise):
plotText(label = "consensus peaks", fontcolor = colorscale[6], fontsize = 8, fontface = "bold",
         x = labmarg, y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1 - tr_h[i]),
         just = c("left","top"), default.units = "inches")
i <- i+1
#
cons_pl <- plotRanges(data = cons[cons$strand == "+", c(1:4)], params = params_d, collapse = TRUE,
                      fill = colorscale[2], linecolor = "fill",
                      y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1),
                      height = tr_h[i], just = c("left","bottom"))
i <- i+1
cons_mi <- plotRanges(data = cons[cons$strand == "-", c(1:4)] , params = params_d, collapse = TRUE,
                      fill = colorscale[7], linecolor = "fill",
                      y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1),
                      height = tr_h[i], just = c("left","bottom"))
i <- i+1

# plot chromatin marks:
k36_tr <- plotSignal(data = k36, params = c(params_d, k36_range),
                     fill = colorscale[3], linecolor = colorscale[3],
                     y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1),
                     height = tr_h[i], just = c("left","bottom"))
plotText(label = "H3K36me3", fontcolor = colorscale[3], fontsize = 8, fontface = "bold",
         x = labmarg, y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1 - tr_h[i]),
         just = c("left","top"), default.units = "inches")
annoYaxis(plot = k36_tr, at = seq(0, lims[3], length.out=5), axisLine = TRUE, fontsize = 8)
i <- i+1
k27ac_tr <- plotSignal(data = k27ac, params = c(params_d, k27ac_range),
                       fill = colorscale[4], linecolor = colorscale[4],
                       y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1),
                       height = tr_h[i], just = c("left","bottom"))
plotText(label = "H3K27ac", fontcolor = colorscale[4], fontsize = 8, fontface = "bold",
         x = labmarg, y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1 - tr_h[i]),
         just = c("left","top"), default.units = "inches")
annoYaxis(plot = k27ac_tr, at = seq(0, lims[4], length.out=5), axisLine = TRUE, fontsize = 8)
i <- i+1

k27me3_tr <- plotSignal(data = k27me3, params = c(params_d, k27me3_range),
                        fill = colorscale[5], linecolor = colorscale[5],
                        y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1),
                        height = tr_h[i], just = c("left","bottom"))
plotText(label = "H3K27me3", fontcolor = colorscale[5], fontsize = 8, fontface = "bold",
         x = labmarg, y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1 - tr_h[i]),
         just = c("left","top"), default.units = "inches")
annoYaxis(plot = k27me3_tr, at = seq(0, lims[5], length.out=5), axisLine = TRUE, fontsize = 8)
i <- i+1




# insertion sites:
plotText(label = "CD4 insertion sites", fontcolor = colorscale[6], fontsize = 8, fontface = "bold",
         x = labmarg, y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1 - tr_h[i]), 
         just = c("left","top"), default.units = "inches")
i <- i+1   # comment this out if label doesn't need to be on its own line

is_tr <- plotRanges(data = is, params = params_d, # collapse = TRUE, 
                    fill = colorscale[6], linecolor = colorscale[6], lwd = 6,
                    y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1), 
                    height = tr_h[i], spaceWidth=0.01, just = c("left","bottom"))
i <- i+1


plotText(label = "Known insertion sites", fontcolor = colorscale[6], fontsize = 8, fontface = "bold",
         x = labmarg, y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1 - tr_h[i]), 
         just = c("left","top"), default.units = "inches")
i <- i+1

is_tr <- plotRanges(data = is2, params = params_d, # collapse = TRUE, 
                    fill = colorscale[6], linecolor = colorscale[6], lwd = 6,
                    y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1), 
                    height = tr_h[i], spaceWidth=0.01, just = c("left","bottom"))
i <- i+1


# gene track:
genes_tr <- plotGenes(params = params_d, stroke = 1, fontsize = 8, strandLabels = TRUE, 
                      fontcolor = c(colorscale[2], colorscale[7]), 
                      fill = c(colorscale[2], colorscale[7]),
                      y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1), 
                      height = tr_h[i], just = c("left","bottom"))
i <- i+1   # comment this out if label doesn't need to be on its own line
# genome scale:
annoGenomeLabel(plot = genes_tr, params = params_d, scale = "Mb", digits = 1, 
                fontsize = 8, y = as.numeric(cumsum(tr_h)[i] + length(tr_h[1:i])*0.1 + 0.05), 
                just = c("left", "top"))


# hide page guides:
pageGuideHide()   

# and save the file
dev.copy(pdf, "plotgardnerv2.pdf", width = w, height = h)
dev.off()
