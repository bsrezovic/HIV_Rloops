library(data.table)
library(ggplot2)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(org.Hs.eg.db)

save.plt.pdf <- function(p, prefix, nplots = c(1, 1)) {
  # save plot as pdf
  pdf(paste0("figs/", prefix, ".pdf"),
      width = 9 * nplots[2],
      height = 9 * nplots[1])
  print(p)
  dev.off()
}


###### LOAD GENIDF ###### 
genidf <- readRDS("genidf_hg19_smaller.RDS")
setnames(genidf,
         c("names", "RIGSLISTE", "ExpressionGroups"),
         c("id", "rigLists", "expGroup"))

genidf <- genidf[seqnames %in% paste0("chr", c(1:22, "X", "Y"))]
genidf[expGroup == "High10%", expGroup := "High 10%"]
genidf[expGroup == "Low10%", expGroup := "Low 10%"]
genidf[, expGroup := factor(expGroup, levels = c("High 10%", "Mid", "Low 10%", "Not expressed"))]
genidf[rigLists == "2 or more", rigLists := "2+"]
genidf <- makeGRangesFromDataFrame(genidf, keep.extra.columns = TRUE)


###### FUNCTIONS FOR READING AND ANNOTATING ###### 
readPeaks <- function(filename) {
  peaks <- fread(filename, col.names = paste0("V", 1:6))
  peaks <- peaks[V1 %in% paste0("chr", c(1:22, "X", "Y"))]
  peaks[, GRanges(V1, IRanges(V2, V3), peak_id = V4)]
}

annotateGenes <- function(genes, peaks) {
  counts <- countOverlaps(genes, peaks)
  scaled <- 1000*counts/width(genes)
  dists <- mcols(distanceToNearest(genes, peaks))$distance
  
  genes$peakCounts <- counts
  genes$peaksPerGeneKb <- scaled
  genes$distToNear <- dists
  as.data.table(genes)
}

###### ANNOTATE CONSENSUS PEAKS ###### 
annots <- annotateGenes(
  genidf,
  readPeaks("Rloops_ConsensusPeaksFromTriplicates_hg19.txt")
)


###### PLOTS BY GENE EXPRESSION FOR CONSENSUS PEAKS ###### 

a <- ggplot(data = annots,
       aes(x = expGroup, y = peakCounts, fill = expGroup)) +
  geom_boxplot(colour = "black", outlier.color = NA) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(lineend = "round"),
        text = element_text(size = 30)) +
  scale_fill_manual(values = c("#000000", "#525252", "#969696", "#cccccc")) +
  labs(x = "Expression group",
       y = "No. of R-loop peaks overlapping genes") +
  coord_cartesian(ylim = c(0, 9))

b <- ggplot(data = annots,
            aes(x = expGroup, y = peaksPerGeneKb, fill = expGroup)) +
  geom_boxplot(colour = "black", outlier.color = NA) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(lineend = "round"),
        text = element_text(size = 30)) +
  scale_fill_manual(values = c("#000000", "#525252", "#969696", "#cccccc")) +
  labs(x = "Expression group",
       y = "No. of R-loop peaks overlapping genes per kb") +
  coord_cartesian(ylim = c(0, 0.4))


###### PLOTS BY RIG/NON-RIG FOR CONSENSUS PEAKS ###### 

c <- ggplot(data = annots,
       aes(x = rigLists, y = peakCounts, fill = rigLists)) +
  geom_boxplot(colour = "black", outlier.color = NA) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(lineend = "round"),
        text = element_text(size = 30)) +
  scale_fill_manual(values = c("#cccccc", "#969696", "#525252")) +
  labs(x = "No. of datasets gene is in",
       y = "No. of R-loop peaks overlapping genes") +
  coord_cartesian(ylim = c(0, 8))

d <- ggplot(data = annots,
       aes(x = rigLists, y = peaksPerGeneKb, fill = rigLists)) +
  geom_boxplot(colour = "black", outlier.color = NA) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(lineend = "round"),
        text = element_text(size = 25)) +
  scale_fill_manual(values = c("#cccccc", "#969696", "#525252")) +
  labs(x = "No. of datasets gene is in",
       y = "No. of R-loop peaks overlapping genes per kb") +
  coord_cartesian(ylim = c(0, 0.08))

###### PLOTS RIG/NON-RIG ~ EXPRESSION GROUP FOR CONSENSUS PEAKS ###### 

e <- ggplot(data = annots,
            aes(x = rigLists, y = peakCounts, fill = rigLists)) +
  geom_boxplot(colour = "black", outlier.color = NA) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(lineend = "round"),
        text = element_text(size = 30)) +
 scale_fill_manual(values = c("#cccccc", "#969696", "#525252")) +
  labs(x = "",
       y = "No. of R-loop peaks overlapping genes") +
  coord_cartesian(ylim = c(0, 10)) +
  facet_wrap( ~ expGroup, nrow = 1)

f <- ggplot(data = annots,
            aes(x = rigLists, y = peaksPerGeneKb, fill = rigLists)) +
  geom_boxplot(colour = "black", outlier.color = NA) +
  theme_bw() +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(lineend = "round"),
        text = element_text(size = 30)) +
  scale_y_continuous(breaks = seq(0, 0.25, 0.025)) +
  scale_fill_manual(values = c("#cccccc", "#969696", "#525252")) +
  labs(x = "",
       y = "No. of R-loop peaks overlapping genes per kb") +
  coord_cartesian(ylim = c(0, 0.25)) +
  facet_wrap( ~ expGroup, nrow = 1)

###### DISTANCE TO NEAREST PEAKS PLOTS ######

g <- ggplot(data = annots,
            aes(x = expGroup, y = distToNear, fill = expGroup)) +
  geom_boxplot(colour = "black", outlier.color = NA) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(lineend = "round"),
        text = element_text(size = 30)) +
  scale_fill_manual(values = c("#000000", "#525252", "#969696", "#cccccc")) +
  labs(x = "Expression group",
       y = "Distance to nearest R-loop peak") +
  coord_cartesian(ylim = c(0, 4e+05))

h <- ggplot(data = annots,
            aes(x = rigLists, y = distToNear, fill = rigLists)) +
  geom_boxplot(colour = "black", outlier.color = NA) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(lineend = "round"),
        text = element_text(size = 30)) +
  scale_fill_manual(values = c("#cccccc", "#969696", "#525252")) +
  labs(x = "No. of datasets gene is in",
       y = "Distance to nearest R-loop peak") +
  coord_cartesian(ylim = c(0, 2e+05))

i <- ggplot(data = annots,
            aes(x = rigLists, y = distToNear, fill = rigLists)) +
  geom_boxplot(colour = "black", outlier.color = NA) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(lineend = "round"),
        text = element_text(size = 30)) +
  scale_fill_manual(values = c("#cccccc", "#969696", "#525252")) +
  labs(x = "",
       y = "Distance to nearest R-loop peak") +
  coord_cartesian(ylim = c(0, 5e+05)) +
  facet_wrap( ~ expGroup, nrow = 1)

###### SAVE CONSENSUS PEAKS PLOTS ###### 
save.plt.pdf(a, prefix = "Consensus_expression_unscaled")
save.plt.pdf(b, prefix = "Consensus_expression_scaled")
save.plt.pdf(c, prefix = "Consensus_integration_unscaled")
save.plt.pdf(d, prefix = "Consensus_integration_scaled")
save.plt.pdf(e, prefix = "Consensus_integration_expression_unscaled", nplots = c(1, 4))
save.plt.pdf(f, prefix = "Consensus_integration_expression_scaled", nplots = c(1, 2))
save.plt.pdf(g, prefix = "Consensus_expression_distToNear")
save.plt.pdf(h, prefix = "Consensus_integration_distToNear")
save.plt.pdf(i, prefix = "Consensus_integration_expression_distToNear", nplots = c(1, 4))


###### INDIVIDUAL SAMPLE PEAKS PLOTTING FUNCTION ###### 

sampleColors <- list(expression = list(c("#9b9db0", "#eb7f7f", "#de0202", "#7d0000"),
                                       c("#9b9db0", "#85b356", "#4c9401", "#274d00"),
                                       c("#9b9db0", "#f5b576", "#ff8000", "#a35202")),
                     targeted = list(c("#eb7f7f", "#de0202", "#7d0000"),
                                     c("#85b356", "#4c9401", "#274d00"),
                                     c("#f5b576", "#ff8000", "#a35202")))

plotPeaks <- function(dirpath) {
  filelist <- list.files(dirpath, pattern = "^peaks.+txt$", full.names = TRUE)
  for (i in 1:length(filelist)) {
    annots <- annotateGenes(genidf, readPeaks(filelist[i]))
    namevar <- paste0("D", i)
    
    
    ###### PLOTS BY GENE EXPRESSION ######
    
    p1 <- ggplot(data = annots,
                 aes(x = expGroup, y = peaksPerGeneKb, fill = expGroup)) +
            geom_boxplot(colour = "black", outlier.color = NA) +
            theme_classic() +
            theme(legend.position = "none",
                  axis.line = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_line(lineend = "round"),
                  text = element_text(size = 30)) +
            scale_fill_manual(values = rev(sampleColors$expression[[i]])) +
            labs(x = "Expression group",
                 y = "No. of R-loop peaks overlapping genes per kb") +
            coord_cartesian(ylim = c(0, 0.5))
    
    p2 <- ggplot(data = annots,
                 aes(x = expGroup, y = peakCounts, fill = expGroup)) +
      geom_boxplot(colour = "black", outlier.color = NA) +
      theme_classic() +
      theme(legend.position = "none",
            axis.line = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(lineend = "round"),
            text = element_text(size = 30)) +
      scale_fill_manual(values = rev(sampleColors$expression[[i]])) +
      labs(x = "Expression group",
           y = "No. of R-loop peaks overlapping genes") +
      coord_cartesian(ylim = c(0, 20))
    
    
    ###### PLOTS BY RIG/NON-RIG ######
    
    p3 <- ggplot(data = annots,
                 aes(x = rigLists, y = peakCounts, fill = rigLists)) +
            geom_boxplot(colour = "black", outlier.color = NA) +
            theme_classic() +
            theme(legend.position = "none",
                  axis.line = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_line(lineend = "round"),
                  text = element_text(size = 30)) +
            scale_fill_manual(values = sampleColors$targeted[[i]]) +
            labs(x = "No. of datasets gene is in",
                 y = "No. of R-loop peaks overlapping genes") +
            coord_cartesian(ylim = c(0, 15))
    
    p4 <- ggplot(data = annots,
                 aes(x = rigLists, y = peaksPerGeneKb, fill = rigLists)) +
            geom_boxplot(colour = "black", outlier.color = NA) +
            theme_classic() +
            theme(legend.position = "none",
                  axis.line = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_line(lineend = "round"),
                  text = element_text(size = 30)) +
            scale_fill_manual(values = sampleColors$targeted[[i]]) +
            labs(x = "No. of datasets gene is in",
                 y = "No. of R-loop peaks overlapping genes per kb") +
            coord_cartesian(ylim = c(0, 0.2))
    
    ###### PLOTS RIG/NON-RIG ~ EXPRESSION GROUP ######
    
    p5 <- ggplot(data = annots,
                aes(x = rigLists, y = peakCounts, fill = rigLists)) +
      geom_boxplot(colour = "black", outlier.color = NA) +
      theme_classic() +
      theme(legend.position = "none",
            axis.line = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(lineend = "round"),
            text = element_text(size = 30)) +
      scale_fill_manual(values = sampleColors$targeted[[i]]) +
      labs(x = "",
           y = "No. of R-loop peaks overlapping genes") +
      coord_cartesian(ylim = c(0, 25)) +
      facet_wrap( ~ expGroup, nrow = 1)
    
    p6 <- ggplot(data = annots,
                aes(x = rigLists, y = peaksPerGeneKb, fill = rigLists)) +
      geom_boxplot(colour = "black", outlier.color = NA) +
      theme_classic() +
      theme(legend.position = "none",
            axis.line = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(lineend = "round"),
            text = element_text(size = 30)) +
      scale_fill_manual(values = sampleColors$targeted[[i]]) +
      labs(x = "",
           y = "No. of R-loop peaks overlapping genes per kb") +
      coord_cartesian(ylim = c(0, 0.4)) +
      facet_wrap( ~ expGroup, nrow = 1)
    
    ###### DISTANCE TO NEAREST PEAKS PLOTS ######
    
    p7 <- ggplot(data = annots,
                 aes(x = rigLists, y = distToNear, fill = rigLists)) +
      geom_boxplot(colour = "black", outlier.color = NA) +
      theme_classic() +
      theme(legend.position = "none",
            axis.line = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(lineend = "round"),
            text = element_text(size = 30)) +
      scale_fill_manual(values = sampleColors$targeted[[i]]) +
      labs(x = "No. of datasets gene is in",
           y = "Distance to nearest R-loop peak") +
      coord_cartesian(ylim = c(0, 2e+05))
    
    p8 <- ggplot(data = annots,
                 aes(x = expGroup, y = distToNear, fill = expGroup)) +
      geom_boxplot(colour = "black", outlier.color = NA) +
      theme_classic() +
      theme(legend.position = "none",
            axis.line = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(lineend = "round"),
            text = element_text(size = 30)) +
      scale_fill_manual(values = rev(sampleColors$expression[[i]])) +
      labs(x = "Expression group",
           y = "Distance to nearest R-loop peak") +
      coord_cartesian(ylim = c(0, 4e+05))
    
    p9 <- ggplot(data = annots,
                 aes(x = rigLists, y = distToNear, fill = rigLists)) +
      geom_boxplot(colour = "black", outlier.color = NA) +
      theme_classic() +
      theme(legend.position = "none",
            axis.line = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(lineend = "round"),
            text = element_text(size = 30)) +
      scale_fill_manual(values = sampleColors$targeted[[i]]) +
      labs(x = "",
           y = "Distance to nearest R-loop peak") +
      coord_cartesian(ylim = c(0, 4e+05)) +
      facet_wrap( ~ expGroup, nrow = 1)
    
    ###### SAVE PLOTS ######
    save.plt.pdf(p1, prefix = paste0(namevar, "_expression_scaled"))
    save.plt.pdf(p2, prefix = paste0(namevar, "_expression_unscaled"))
    save.plt.pdf(p3, prefix = paste0(namevar, "_integration_unscaled"))
    save.plt.pdf(p4, prefix = paste0(namevar, "_integration_scaled"))
    save.plt.pdf(p5, prefix = paste0(namevar, "_integration_expression_unscaled"), nplots = c(1, 4))
    save.plt.pdf(p6, prefix = paste0(namevar, "_integration_expression_scaled"), nplots = c(1, 4))
    save.plt.pdf(p7, prefix = paste0(namevar, "_integration_distToNear"))
    save.plt.pdf(p8, prefix = paste0(namevar, "_expression_distToNear"))
    save.plt.pdf(p9, prefix = paste0(namevar, "_integration_expression_distToNear"), nplots = c(1, 4))
  }
}

plotPeaks(".")


################## PLOTS FOR INTRA- AND INTERGENIC ANNOTATIONS ##################


###### ANNOTATION PLOTS BY RIG/NON-RIG FOR CONSENSUS PEAKS ###### 

plotPeakAnnots <- function(genidf, peaks, prefix, title) {
  rig_peaks <- peaks[overlapsAny(peaks, genidf[genidf$rigLists == "2+"])]
  other_peaks <- peaks[overlapsAny(peaks, genidf[genidf$rigLists != "2+"])]
  peakList <- list(RIG = rig_peaks, Other = other_peaks)
  peakAnnoList <- lapply(peakList,
                         annotatePeak,
                         level = "gene",
                         tssRegion = c(-3000, 3000),
                         TxDb = txdb,
                         annoDb = "org.Hs.eg.db",
                         verbose = FALSE)
  save.plt.pdf(plotAnnoBar(peakAnnoList, title = title),
               prefix = prefix, nplots = c(0.5, 1))
}

plotPeakAnnots(genidf,
               peaks = readPeaks("Rloops_ConsensusPeaksFromTriplicates_hg19.txt"),
               prefix = "Consensus_Rloop_profile",
               title = "Consensus peaks feature distribution")

###### ANNOTATION PLOTS BY RIG/NON-RIG FOR INDIVIDUAL SAMPLES ###### 

plotProfile <- function(dirpath) {
  filelist <- list.files(dirpath, pattern = "^peaks.+txt$", full.names = TRUE)
  for (i in 1:length(filelist)) {
    peaks <- readPeaks(filelist[i])
    namevar <- paste0("D", i)
    plotPeakAnnots(genidf, peaks = peaks,
                   prefix = paste0(namevar, "_Rloop_profile"),
                   title = paste0(namevar, " peaks feature distribution"))
  }
}

plotProfile(".")
