setwd("~/Documents/WORK/rloop_paper_figs")

library(data.table)
library(readxl)
library(stringr)
library(ggplot2)
library(GenomicRanges)

peaks <- fread("Rloops_ConsensusPeaksFromTriplicates_hg19.txt")
colnames(peaks)[1] <- "V1"
peaks <- peaks[, GRanges(V1, IRanges(V2, V3))]

genidf <- readRDS("genidf_hg19_smaller.RDS")
setnames(genidf,
         c("names", "RIGSLISTE", "ExpressionGroups"),
         c("gene_id", "rigLists", "expGroup"))

exon <- read_excel("exon_AQR_KO-vs-NT.xlsx", sheet = "DEGs")
exon <- as.data.table(exon)
exon[, gene_id := gsub("\\.\\d*$", "", gene_id)]
exon <- merge(exon, genidf, by = "gene_id", all = TRUE)
exon[is.na(rigLists), rigLists := "0"]

body <- read_excel("Aquarius_2022_AQR_KO_genebody.xlsx", sheet = "DEGs")
body <- as.data.table(body)
body[, gene_id := gsub("\\.\\d*$", "", gene_id)]
body <- merge(body, genidf, by = "gene_id", all = TRUE)
body[is.na(rigLists), rigLists := "0"]

make_cat <- function(x) {
  cat <- rep("ns", nrow(x))
  cat[x$FDR < 0.05 & x$logFC > 0] <- "Up"
  cat[x$FDR < 0.05 & x$logFC < 0] <- "Down"
  x$cat <- cat
  x
}

exon <- make_cat(exon)
body <- make_cat(body)

exon <- makeGRangesFromDataFrame(exon, na.rm = T, keep.extra.columns = T)
exon$RloopOverlaps <- suppressWarnings(countOverlaps(exon, peaks))

body <- makeGRangesFromDataFrame(body, na.rm = T, keep.extra.columns = T)
body$RloopOverlaps <- suppressWarnings(countOverlaps(body, peaks))

ggplot(as.data.frame(table(cat = exon$cat, rig = exon$rigLists)),
       aes(Freq, cat, fill = rig, label = Freq)) +
  geom_col(position = "fill") +
  geom_text(position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values = c("grey80", "skyblue2", "steelblue4"),
                    name = "RIG category") +
  labs(x = "Fraction", y = "AQR KO response") +
  theme_bw() +
  theme(aspect.ratio = 2/(1 + sqrt(5))) +
  ggtitle("Exons")

ggsave("./aqr_figs/rigs_in_degs_exon.pdf", width = 6, height = 4)

ggplot(as.data.frame(table(cat = body$cat, rig = body$rigLists)),
       aes(Freq, cat, fill = rig, label = Freq)) +
  geom_col(position = "fill") +
  geom_text(position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values = c("grey80", "skyblue2", "steelblue4"),
                    name = "RIG category") +
  labs(x = "Fraction", y = "AQR KO response") +
  theme_bw() +
  theme(aspect.ratio = 2/(1 + sqrt(5))) +
  ggtitle("Gene bodies")

ggsave("./aqr_figs/rigs_in_degs_body.pdf", width = 6, height = 4)

ggplot(as.data.frame(table(cat = exon$cat, rloop = exon$RloopOverlaps > 0)),
       aes(Freq, cat, fill = rloop, label = Freq)) +
  geom_col(position = "fill") +
  geom_text(position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values = c("grey80", "skyblue2"),
                    labels = c("Not overlapping R-loop",
                               "Overlapping R-loop"),
                    name = NULL) +
  labs(x = "Fraction", y = "AQR KO response") +
  theme_bw() +
  theme(aspect.ratio = 2/(1 + sqrt(5))) +
  ggtitle("Exons")

ggsave("./aqr_figs/rloops_in_degs_exon.pdf", width = 6, height = 4)

ggplot(as.data.frame(table(cat = body$cat, rloop = body$RloopOverlaps > 0)),
       aes(Freq, cat, fill = rloop, label = Freq)) +
  geom_col(position = "fill") +
  geom_text(position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values = c("grey80", "skyblue2"),
                    labels = c("Not overlapping R-loops",
                               "Overlapping R-loop"),
                    name = NULL) +
  labs(x = "Fraction", y = "AQR KO response") +
  theme_bw() +
  theme(aspect.ratio = 2/(1 + sqrt(5))) +
  ggtitle("Gene bodies")

ggsave("./aqr_figs/rloops_in_degs_body.pdf", width = 6, height = 4)

ggplot(as.data.frame(table(cat = exon$cat, rloop = exon$RloopOverlaps > 0)),
       aes(Freq, rloop, fill = cat, label = Freq)) +
  geom_col(position = "fill") +
  geom_text(position = position_fill(vjust = 0.5)) +
  scale_y_discrete(labels = c("Not overlapping R-Loops",
                              "Overlapping R-Loops")) +
  scale_fill_manual(values = c("deepskyblue3", "grey80", "brown2"),
                    name = "AQR KO response") +
  labs(x = "Fraction", y = "Gene category") +
  theme_bw() +
  theme(aspect.ratio = 2/(1 + sqrt(5))) +
  ggtitle("Exons")

ggsave("./aqr_figs/proportion_of_exon_degs_in_rloop_genes.pdf", width = 6, height = 4)

ggplot(as.data.frame(table(cat = body$cat, rloop = body$RloopOverlaps > 0)),
       aes(Freq, rloop, fill = cat, label = Freq)) +
  geom_col(position = "fill") +
  geom_text(position = position_fill(vjust = 0.5)) +
  scale_y_discrete(labels = c("Not overlapping R-Loops",
                              "Overlapping R-Loops")) +
  scale_fill_manual(values = c("deepskyblue3", "grey80", "brown2"),
                    name = "AQR KO response") +
  labs(x = "Fraction", y = "Gene category") +
  theme_bw() +
  theme(aspect.ratio = 2/(1 + sqrt(5))) +
  ggtitle("Gene bodies")

ggsave("./aqr_figs/proportion_of_body_degs_in_rloop_genes.pdf", width = 6, height = 4)

ggplot(as.data.frame(table(genidf$RloopOverlaps > 0, genidf$rigLists)),
       aes(Freq, Var1, fill = Var2, label = Freq)) +
  geom_col(position = "fill") +
  geom_text(position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values = c("#cccccc", "#969696", "#525252"),
                    name = "RIG category") +
  labs(x = "Fraction", y = "Gene category") +
  theme_bw() +
  theme(aspect.ratio = 2/(1 + sqrt(5)))

ggsave("./aqr_figs/proportion_of_rigs_in_rloop_genes.pdf", width = 6, height = 4)


ggplot(as.data.frame(table(cat = exon$cat, rig = exon$rigLists)),
       aes(Freq, rig, fill = cat, label = Freq)) +
  geom_col(position = "fill") +
  geom_text(position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values = c("deepskyblue3", "grey80", "brown2"),
                    name = "AQR KO response") +
  labs(x = "Fraction", y = "RIG category") +
  theme_bw() +
  theme(aspect.ratio = 2/(1 + sqrt(5))) +
  ggtitle("Exons")

ggsave("./aqr_figs/proportion_of_exon_degs_in_rigs.pdf", width = 6, height = 4)

ggplot(as.data.frame(table(cat = body$cat, rig = body$rigLists)),
       aes(Freq, rig, fill = cat, label = Freq)) +
  geom_col(position = "fill") +
  geom_text(position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values = c("deepskyblue3", "grey80", "brown2"),
                    name = "AQR KO response") +
  labs(x = "Fraction", y = "RIG category") +
  theme_bw() +
  theme(aspect.ratio = 2/(1 + sqrt(5))) +
  ggtitle("Gene bodies")

ggsave("./aqr_figs/proportion_of_body_degs_in_rigs.pdf", width = 6, height = 4)

ggplot(as.data.frame(exon),
       aes(cat, RloopOverlaps, colour = cat)) +
  geom_boxplot(outlier.colour = NA) +
  scale_colour_manual(values = c("deepskyblue3", "grey80", "brown2"),
                      name = "AQR KO response") +
  theme_bw() +
  theme(aspect.ratio = 2/(1 + sqrt(5))) +
  ggtitle("Exons") +
  labs(x = "AQR KO response", y = "# R-loops") +
  facet_wrap(~rigLists) +
  coord_cartesian(ylim = c(0, 8))

ggsave("./aqr_figs/rloop_number_in_exon_degs_and_rigs.pdf", width = 7, height = 4)

ggplot(as.data.frame(body),
       aes(cat, RloopOverlaps, colour = cat)) +
  geom_boxplot(outlier.colour = NA) +
  scale_colour_manual(values = c("deepskyblue3", "grey80", "brown2"),
                      name = "AQR KO response") +
  theme_bw() +
  theme(aspect.ratio = 2/(1 + sqrt(5))) +
  ggtitle("Gene bodies") +
  labs(x = "AQR KO response", y = "# R-loops") +
  facet_wrap(~rigLists) +
  coord_cartesian(ylim = c(0, 8))

ggsave("./aqr_figs/rloop_number_in_body_degs_and_rigs.pdf", width = 7, height = 4)

ggplot(as.data.frame(exon),
       aes(cat, 1000*RloopOverlaps/width(exon), colour = cat)) +
  geom_boxplot(outlier.colour = NA) +
  scale_colour_manual(values = c("deepskyblue3", "grey80", "brown2"),
                      name = "AQR KO response") +
  theme_bw() +
  theme(aspect.ratio = 2/(1 + sqrt(5))) +
  ggtitle("Exon DEGs") +
  labs(x = "AQR KO response", y = "# R-loops per 1 kb of gene") +
  facet_wrap(~rigLists) +
  coord_cartesian(ylim = c(0, 0.2))

ggsave("./aqr_figs/rloop_density_in_exon_degs_and_rigs.pdf", width = 7, height = 4)

ggplot(as.data.frame(body),
       aes(cat, 1000*RloopOverlaps/width(body), colour = cat)) +
  geom_boxplot(outlier.colour = NA) +
  scale_colour_manual(values = c("deepskyblue3", "grey80", "brown2"),
                      name = "AQR KO response") +
  theme_bw() +
  theme(aspect.ratio = 2/(1 + sqrt(5))) +
  ggtitle("Gene body DEGs") +
  labs(x = "AQR KO response", y = "# R-loops per 1 kb of gene") +
  facet_wrap(~rigLists) +
  coord_cartesian(ylim = c(0, 0.2))

ggsave("./aqr_figs/rloop_density_in_body_degs_and_rigs.pdf", width = 7, height = 4)


df <- as.data.frame(table(
  cat = exon$cat,
  rig = exon$rigLists,
  rloop = exon$RloopOverlaps > 0
))

ggplot(df, aes(Freq, rloop, fill = cat, label = Freq)) +
  geom_col(position = "fill") +
  geom_text(position = position_fill(vjust = 0.5)) +
  scale_y_discrete(labels = c("Not overlapping R-loops",
                              "Overlapping R-loops")) +
  scale_fill_manual(values = c("deepskyblue3", "grey80", "brown2"),
                    name = "AQR KO response") +
  theme_bw() +
  theme(aspect.ratio = 1/(1 + sqrt(5))) +
  ggtitle("Exon DEGs") +
  labs(x = "Fraction", y = "Gene category") +
  facet_grid(rig~.)

ggsave("./aqr_figs/deg_exon_proportions_in_rig_rloop_categories.pdf", width = 7, height = 4)

df <- as.data.frame(table(
  cat = body$cat,
  rig = body$rigLists,
  rloop = body$RloopOverlaps > 0
))

ggplot(df, aes(Freq, rloop, fill = cat, label = Freq)) +
  geom_col(position = "fill") +
  geom_text(position = position_fill(vjust = 0.5)) +
  scale_y_discrete(labels = c("Not overlapping R-loops",
                              "Overlapping R-loops")) +
  scale_fill_manual(values = c("deepskyblue3", "grey80", "brown2"),
                      name = "AQR KO response") +
  theme_bw() +
  theme(aspect.ratio = 1/(1 + sqrt(5))) +
  ggtitle("Gene body DEGs") +
  labs(x = "Fraction", y = "Gene category") +
  facet_grid(rig~.)

ggsave("./aqr_figs/deg_body_proportions_in_rig_rloop_categories.pdf", width = 7, height = 4)


writexl::write_xlsx(as.data.frame(body), "./tables/aqr_body_degs_rigs_rloops.xlsx")
writexl::write_xlsx(as.data.frame(exon), "./tables/aqr_exon_degs_rigs_rloops.xlsx")
