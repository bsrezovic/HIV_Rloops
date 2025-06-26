library(data.table)
library(stringr)
library(Biostrings)
library(dplyr)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(biomaRt)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)  # Replace with your genome
library(MASS)  # For fitting distributions
library(ggplot2)  #
library(scales)

## Making the Random matched controls

#1. load the fragemntation data from our files


experimental <- import.bed("C:/Users/38598/Desktop/hivint/Miseq_part2/miseq_fragments.bed")
experimental_granges <- GRanges(experimental)

# Calculate fragment lengths
fragment_lengths <- width(experimental_granges)

# Fit a normal distribution to fragment lengths
fit <- fitdistr(fragment_lengths, "normal")
mean_length <- fit$estimate["mean"]
sd_length <- fit$estimate["sd"]

# Plot observed vs. fitted distribution
ggplot(data.frame(Length = fragment_lengths), aes(x = Length)) +
  geom_density() +
  stat_function(fun = dnorm, args = list(mean = mean_length, sd = sd_length), color = "red")


library(regioneR)

# Define the genome (e.g., hg3

genome <- getGenome("hg19")

# Generate random regions with matched lengths
random_controls <- createRandomRegions(
  nregions = 3*length(experimental_granges),
  length.mean = mean_length,
  length.sd = sd_length,
  genome = genome,
  mask = NULL  # Exclude regions like gaps if needed
)


# Get GC content for experimental fragments
experimental_gc <- letterFrequency(
  getSeq(BSgenome.Hsapiens.UCSC.hg19, experimental_granges),
  letters = "GC",
  as.prob = TRUE
)

# Bin GC content and sample controls to match
gc_bins <- cut(experimental_gc, breaks = 10)
random_controls <- random_controls[sample(length(random_controls)/3, prob = gc_bins)]


export.bed(random_controls, "matched_controls.bed")

# we  have 33534 sitesies
# so we select 100602 random controls

chosen_rand <- sample(random_controls,100602,replace = T)
chosen_rand <- as.data.table(chosen_rand)
## try to redo their tier one graphs (just the histograms of control vs active/resting)
write.csv(chosen_rand, "random_miseq.csv")
#1
#All integration sites and matched random controls were annotated for gene density in the 1 Mb region surrounding the integration site
# lets do just that


# Load gene annotations 
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes_gr <- genes(txdb)

# load our integration sites to go along the simulated controls
peaks1 <- fread("C:/Users/38598/Desktop/hivint/Miseq_part2/HIV1_data3.bed")
peaks2 <- fread("C:/Users/38598/Desktop/hivint/Miseq2/HIV1_data3.bed")
peaks <- rbind(peaks1,peaks2)

# add in the "matched random controls" and pursue 
peaks <- peaks[,c(1,2,3,7)]

chosen_rand <- chosen_rand[,c(1,2,3)]
colnames(chosen_rand) <- c("V1","V2","V3")
chosen_rand$V7 <- "MRC"
chosen_rand$V8 <- "MRC"
peaks$V8 <- "experimental"

peaks_all <- rbind(peaks,chosen_rand)
# Convert to GRanges and set chromosome info
coords_gr <- makeGRangesFromDataFrame(peaks_all, keep.extra.columns = TRUE, start.field = "V2", end.field = "V3", seqnames.field = "V1")
seqinfo(coords_gr) <- seqinfo(txdb)[seqlevels(coords_gr)]

# Create 1 Mb windows centered on each coordinate - this also decides that integrations are centers of simulated fragments !?
# actually create multiple windows for comparisons sake
windows_gr1 <- resize(coords_gr, width = 1e6, fix = "center")
windows_gr2 <- resize(coords_gr, width = 1e5, fix = "center")
windows_gr3 <- resize(coords_gr, width = 1e4, fix = "center")
windows_gr4 <- resize(coords_gr, width = 1e3, fix = "center")
windows_gr1 <- trim(windows_gr1)  # Ensure windows stay within chromosome bounds
windows_gr2 <- trim(windows_gr2)
windows_gr3 <- trim(windows_gr3)
windows_gr4 <- trim(windows_gr4)

# Count genes overlapping each window
coords_gr$gene_density1MB <- countOverlaps(windows_gr1, genes_gr)
coords_gr$gene_density100kB <- countOverlaps(windows_gr2, genes_gr)
coords_gr$gene_density10kB <- countOverlaps(windows_gr3, genes_gr)
coords_gr$gene_density1kB <- countOverlaps(windows_gr4, genes_gr)


# Convert back to data frame if needed
annotated_df <- as.data.frame(coords_gr)

# density of the gene density lol
annotated_df %>%
  ggplot( aes(x=gene_density100kB)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)
# okay time to bin this into 10 bins?
#  10 equal bis but how?


annotated_df <- annotated_df %>% mutate(gene_density_bin_1MB = ntile(gene_density1MB, 10),
                                        gene_density_bin_100kB = ntile(gene_density100kB, 10),
                                        gene_density_bin_10kB = ntile(gene_density10kB, 10),
                                        gene_density_bin_1kB = ntile(gene_density1kB, 10))

# Compute counts and the proportion (relative to total sites) for each bin and group

annotated_df_summary <- annotated_df %>%
  group_by(gene_density_bin_1MB, V8) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(V8) %>%
  mutate(total_in_group = sum(count),
         proportion = count / total_in_group)

# Plot side-by-side bars for the two groups in each bin
p1 <- ggplot(annotated_df_summary, aes(x = factor(gene_density_bin_1MB), y = proportion, fill = V8)) +
  geom_bar(stat = "identity",width = 0.8, position = position_dodge(width = 0.8), color = "black") +
  labs(title = "",
       x = "Gene Density, 1MB Windows",
       y = "Proportion of Sites") +
  theme_minimal()+ guides(fill=guide_legend(title=""))

p1


# splitting this all up by group
annotated_df <- as.data.table(annotated_df)
annotated_df[,group:= ifelse(V7%in% c("D183_mock","D184_mock"),"Resting mock",
                    ifelse(V7%in% c("D183A_mock","D185A_mock"),"Active mock",
                           ifelse(V7%in% c("D183_NT","D184_NT"),"Resting NT",
                                  ifelse(V7%in% c("D183A_NT","D185A_NT"),"Active NT",
                                         ifelse(V7%in% c( "D183_Nup153KO","D184_Nup153KO"),"Resting Nup153KO",
                                                ifelse(V7%in% c("D183A_Nup153KO","D185A_Nup153KO"),"Active Nup153KO",
                                                       ifelse(V7%in% c("D183A_Nup62KO","D185A_Nup62KO"),"Active Nup62KO",
                                                       ifelse(V7%in% c("D183_Nup62KO","D184_Nup62KO"),"Resting Nup62KO","MRC"
                                                       ))))))))]


# drop the mocks when dping this
annotated_df_summary <- annotated_df[!group %in% c("Resting mock","Active mock"),] %>%
  group_by(gene_density_bin_1MB, group) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(group) %>%
  mutate(total_in_group = sum(count),
         proportion = count / total_in_group)
annotated_df_summary <- as.data.table(annotated_df_summary)
# Plot side-by-side bars for the two groups in each bin
annotated_df_summary$group <- factor(annotated_df_summary$group,
                                     levels = c("MRC","Resting NT","Active NT",
                                                "Resting Nup153KO","Active Nup153KO",
                                                "Resting Nup62KO","Active Nup62KO"))

p2 <- ggplot(annotated_df_summary[!group %in% c("Active mock","Resting mock"),], aes(x = factor(gene_density_bin_1MB), y = proportion, fill = group)) +
  geom_bar(stat = "identity",width = 0.8, position = position_dodge(width = 0.8), color = "black") +
  labs(title = "",
       x = "Gene Density, 1MB Windows",
       y = "Proportion of Sites") +
  scale_fill_brewer(palette="Blues", direction = 1)+
  theme_minimal()+ guides(fill=guide_legend(title=""))
p2

p3 <- ggplot(annotated_df_summary[!group %in% c("Active NT","Active Nup153KO","Active Nup62KO"), ],
             aes(x = factor(gene_density_bin), y = proportion, fill = group)) +
  geom_bar(stat = "identity",width = 0.8, position = position_dodge(width = 0.8), color = "black") +
  labs(title = "Resting CD4",
       x = "Gene Density, 1MB Windows",
       y = "Proportion of Sites") +
  scale_fill_brewer(palette="Blues", direction = 1)+
  theme_minimal()+ guides(fill=guide_legend(title=""))
p3

p4 <- ggplot(annotated_df_summary[group %in% c("MRC","Active NT","Active Nup153KO","Active Nup62KO"), ],
             aes(x = factor(gene_density_bin), y = proportion, fill = group)) +
  geom_bar(stat = "identity",width = 0.8, position = position_dodge(width = 0.8), color = "black") +
  labs(title = "Active CD4",
       x = "Gene Density, 1MB Windows",
       y = "Proportion of Sites") +
  scale_fill_brewer(palette="Blues", direction = 1)+
  theme_minimal()+ guides(fill=guide_legend(title=""))
p4



# Continuing the heatmap procedure

# trying to calculate rank information for each integration site by comapring to three nearest MRCs

annotated_df$ROC_rank_dens_1MB <- 0  # this will be the default value, also ergo the value the MRCs will keep
annotated_df$ROC_rank_dens_100kB <- 0 
annotated_df$ROC_rank_dens_10kB <- 0 
annotated_df$ROC_rank_dens_1kB <- 0 
MRC <- annotated_df[group == "MRC",]
MRCr <- IRanges(start = MRC$start, end = MRC$end, names = MRC$seqnames)
for (i in 1:nrow(annotated_df[group!="MRC",])){
  exp_dens1 <- annotated_df[i,]$gene_density1MB# find the nearest three MRCs, making of IRanges, distancing them, remembering gene ranks
  exp_dens2 <- annotated_df[i,]$gene_density100kB
  exp_dens3 <- annotated_df[i,]$gene_density10kB
  exp_dens4 <- annotated_df[i,]$gene_density1kB
  exp <- IRanges(start = annotated_df[i,]$start, end = annotated_df[i,]$end, names = annotated_df[i,]$seqnames)
  nearest1 <- nearest(exp,MRCr)
  n1_dens1 <- MRC[nearest1,]$gene_density1MB
  n1_dens2 <- MRC[nearest1,]$gene_density100kB
  n1_dens3 <- MRC[nearest1,]$gene_density10kB
  n1_dens4 <- MRC[nearest1,]$gene_density1kB
  nearest2 <- nearest(exp,MRCr[-nearest1])
  n2_dens1 <- MRC[nearest2,]$gene_density1MB
  n2_dens2 <- MRC[nearest2,]$gene_density100kB
  n2_dens3 <- MRC[nearest2,]$gene_density10kB
  n2_dens4 <- MRC[nearest2,]$gene_density1kB
  nearest3 <- nearest(exp,MRCr[-c(nearest1,nearest2)])
  n3_dens1 <- MRC[nearest3,]$gene_density1MB
  n3_dens2 <- MRC[nearest3,]$gene_density100kB
  n3_dens3 <- MRC[nearest3,]$gene_density10kB
  n3_dens4 <- MRC[nearest3,]$gene_density1kB
  densities1 <- c(exp_dens1,n1_dens1,n2_dens1,n3_dens1)
  densities2 <- c(exp_dens2,n1_dens2,n2_dens2,n3_dens2)
  densities3 <- c(exp_dens3,n1_dens3,n2_dens3,n3_dens3)
  densities4 <- c(exp_dens4,n1_dens4,n2_dens4,n3_dens4)
  names(densities1) <- c("exp","MRC1","MRC2","MRC3")
  names(densities2) <- c("exp","MRC1","MRC2","MRC3")
  names(densities3) <- c("exp","MRC1","MRC2","MRC3")
  names(densities4) <- c("exp","MRC1","MRC2","MRC3")
  densities1 <- sort(densities1,decreasing = T)
  densities2 <- sort(densities2,decreasing = T)
  densities3 <- sort(densities3,decreasing = T)
  densities4 <- sort(densities4,decreasing = T)
  rank1 = which(names(densities1)=="exp")
  rank2 = which(names(densities2)=="exp")
  rank3 = which(names(densities3)=="exp")
  rank4 = which(names(densities4)=="exp")
  annotated_df[i,]$ROC_rank_dens_1MB <- rank1
  annotated_df[i,]$ROC_rank_dens_100kB <- rank2
  annotated_df[i,]$ROC_rank_dens_10kB <- rank3
  annotated_df[i,]$ROC_rank_dens_1kB <- rank4
}

write.csv(annotated_df, "savepoint_roccurve.csv")

# calculate the number of matched random sites ranked below each one of the experimetal ones
# note this moves the default value for the MRCs into -4

annotated_df[,MRC_below_1MB:=4-ROC_rank_dens_1MB]
annotated_df[,MRC_below_100kB:=4-ROC_rank_dens_100kB]
annotated_df[,MRC_below_10kB:=4-ROC_rank_dens_10kB]
annotated_df[,MRC_below_1kB:=4-ROC_rank_dens_1kB]

test <- annotated_df[group!="MRC",] %>%  # have to exclude MRC here to not get fucky wucky with the data
  mutate(MRC_below_div_1MB = MRC_below_1MB /3,
         MRC_below_div_100kB = MRC_below_100kB /3,
         MRC_below_div_10kB = MRC_below_10kB /3,
         MRC_below_div_1kB = MRC_below_1kB /3) %>% 
  group_by(group) %>%
  mutate(counts = n()) %>%
  ungroup() %>%
  group_by(group) %>%
  mutate(ROC_value_1MB = sum(MRC_below_div_1MB)/counts,
         ROC_value_100kB = sum(MRC_below_div_100kB)/counts,
         ROC_value_10kB = sum(MRC_below_div_10kB)/counts,
         ROC_value_1kB = sum(MRC_below_div_1kB)/counts)







pdf("Histograms_ROC_test2.pdf", height=6, width=8)

p2
p3
p4
proc
dev.off()


pdf("STP_04_2025_v2.pdf", height=6, width=8)

plot34
plot35
dev.off()






# continue with comarisons in relation to R loops, H3K36me3, H3K4me1, H3K20 , Pol II H3K27Ac H3K27me3 

# 1. load the rloops and other data, lets try to do distance to rloop

rloops <- fread("C:/Users/38598/Desktop/hivint/peakCalling/ConsensusPeaksFromTriplicates.txt", header = F)
rloops <- makeGRangesFromDataFrame(rloops, keep.extra.columns = TRUE, start.field = "V2", end.field = "V3", seqnames.field = "V1")
seqinfo(rloops) <- seqinfo(txdb)[seqlevels(rloops)]
h3k4me1 <- fread("C:/Users/38598/Desktop/hivint/Miseq_part2/H3K4me1_SE.bed")
h3k4me1 <- makeGRangesFromDataFrame(h3k4me1, keep.extra.columns = TRUE, start.field = "V2", end.field = "V3", seqnames.field = "V1")
seqinfo(h3k4me1) <- seqinfo(txdb)[seqlevels(h3k4me1)]
H3K27ac <- fread("C:/Users/38598/Desktop/hivint/Miseq_part2/H3K27ac_SE.bed")
H3K27ac <- makeGRangesFromDataFrame(H3K27ac, keep.extra.columns = TRUE, start.field = "V2", end.field = "V3", seqnames.field = "V1")
seqinfo(H3K27ac) <- seqinfo(txdb)[seqlevels(H3K27ac)]
# step 1 fid the distance to nearest for each integration site to these freatures
# it is here that we will have to limit ourselves to cannonical chromosomes because this data is sparser then the annotated_df!!!

annotated_df2 <- annotated_df[annotated_df$seqnames %in% levels(seqnames(H3K27ac))[-length(levels(seqnames(H3K27ac)))]]
coords_gr <- makeGRangesFromDataFrame(annotated_df2, keep.extra.columns = TRUE, start.field = "start", end.field = "end", seqnames.field = "seqnames")
seqinfo(coords_gr) <- seqinfo(txdb)[seqlevels(coords_gr)]

windows_gr <- resize(coords_gr, width = 1, fix = "center") # were doing distance so to avoid problems with the way i made MRCs

rloopdist <- distanceToNearest(windows_gr, rloops)
h3k4me1dist <- distanceToNearest(windows_gr, h3k4me1)
h3k27acdist <- distanceToNearest(windows_gr, H3K27ac)

annotated_df2$nearest_rloop <- rloopdist@elementMetadata$distance
annotated_df2$nearest_H3K4me1 <- h3k4me1dist@elementMetadata$distance
annotated_df2$nearest_H3K27ac <- h3k27acdist@elementMetadata$distance

MRC <- annotated_df2[group == "MRC",]
MRCr <- IRanges(start = MRC$start, end = MRC$end, names = MRC$seqnames)
annotated_df2$ROC_rank_rloop <- 0
annotated_df2$ROC_rank_H3K4me1 <- 0
annotated_df2$ROC_rank_H3K27ac <- 0
for (i in 1:nrow(annotated_df2[group!="MRC",])){
  rdist <- annotated_df2[i,]$nearest_rloop
  h3k4dist <- annotated_df2[i,]$nearest_H3K4me1
  h3k27dist <- annotated_df2[i,]$nearest_H3K27ac
  
  exp <- IRanges(start = annotated_df2[i,]$start, end = annotated_df2[i,]$end, names = annotated_df2[i,]$seqnames)
  nearest1 <- nearest(exp,MRCr)
  
  n1_dens1 <- MRC[nearest1,]$nearest_rloop
  n1_dens2 <- MRC[nearest1,]$nearest_H3K4me1
  n1_dens3 <- MRC[nearest1,]$nearest_H3K27ac
  
  nearest2 <- nearest(exp,MRCr[-nearest1])
  n2_dens1 <- MRC[nearest2,]$nearest_rloop
  n2_dens2 <- MRC[nearest2,]$nearest_H3K4me1
  n2_dens3 <- MRC[nearest2,]$nearest_H3K27ac

  nearest3 <- nearest(exp,MRCr[-c(nearest1,nearest2)])
  n3_dens1 <- MRC[nearest3,]$nearest_rloop
  n3_dens2 <- MRC[nearest3,]$nearest_H3K4me1
  n3_dens3 <- MRC[nearest3,]$nearest_H3K27ac
  
  densities1 <- c(rdist,n1_dens1,n2_dens1,n3_dens1)
  densities2 <- c(h3k4dist,n1_dens2,n2_dens2,n3_dens2)
  densities3 <- c(h3k27dist,n1_dens3,n2_dens3,n3_dens3)
  
  names(densities1) <- c("exp","MRC1","MRC2","MRC3")
  names(densities2) <- c("exp","MRC1","MRC2","MRC3")
  names(densities3) <- c("exp","MRC1","MRC2","MRC3")

  densities1 <- sort(densities1,decreasing = F)
  densities2 <- sort(densities2,decreasing = F)
  densities3 <- sort(densities3,decreasing = F)

  
  rank1 = which(names(densities1)=="exp")
  rank2 = which(names(densities2)=="exp")
  rank3 = which(names(densities3)=="exp")

  
  annotated_df2[i,]$ROC_rank_rloop <- rank1
  annotated_df2[i,]$ROC_rank_H3K4me1 <- rank2
  annotated_df2[i,]$ROC_rank_H3K27ac <- rank3
  
}

# managed to fuck up the previous columns so this will take some finagling / we should do this anyway since the 
# first graph with the densities used all chromosomal integrations; not just the cannonical ones

nrow(annotated_df2)
annotated_df <- annotated_df[annotated_df$seqnames %in% levels(seqnames(H3K27ac))[-length(levels(seqnames(H3K27ac)))]]

annotated_df <- cbind(annotated_df,annotated_df2[,c(26:31)] )
annotated_df[,MRC_below_rloop:=4-ROC_rank_rloop]
annotated_df[,MRC_below_h3k4me1:=4-ROC_rank_H3K4me1]
annotated_df[,MRC_below_h3k27ac:=4-ROC_rank_H3K27ac]
# save this

saveRDS(annotated_df, file = "df_density_rl_h3k4_h3k27.rds")

annotated_df <- readRDS("C:/Users/38598/Documents/df_density_rl_h3k4_h3k27.rds")

# adding the rest of the comaprison values
H3K36me3 <- fread("C:/Users/38598/Desktop/hivint/H3K36me3.bed", header = F)
H3K36me3 <- makeGRangesFromDataFrame(H3K36me3, keep.extra.columns = TRUE, start.field = "V2", end.field = "V3", seqnames.field = "V1")
seqinfo(H3K36me3) <- seqinfo(txdb)[seqlevels(H3K36me3)]

H3K27me3.bed <- fread("C:/Users/38598/Desktop/hivint/H3K27me3.bed", header = F)
H3K27me3.bed <- makeGRangesFromDataFrame(H3K27me3.bed, keep.extra.columns = TRUE, start.field = "V2", end.field = "V3", seqnames.field = "V1")
seqinfo(H3K27me3.bed) <- seqinfo(txdb)[seqlevels(H3K27me3.bed)]

pol2.bed <- fread("C:/Users/38598/Desktop/hivint/pol2.bed", header = F)
pol2.bed <- makeGRangesFromDataFrame(pol2.bed, keep.extra.columns = TRUE, start.field = "V2", end.field = "V3", seqnames.field = "V1")
seqinfo(pol2.bed) <- seqinfo(txdb)[seqlevels(pol2.bed)]

coords_gr <- makeGRangesFromDataFrame(annotated_df, keep.extra.columns = TRUE, start.field = "start", end.field = "end", seqnames.field = "seqnames")
seqinfo(coords_gr) <- seqinfo(txdb)[seqlevels(coords_gr)]

windows_gr <- resize(coords_gr, width = 1, fix = "center") # were doing distance so to avoid problems with the way i made MRCs

H3K36me3 <- distanceToNearest(windows_gr, H3K36me3)
H3K27me3 <- distanceToNearest(windows_gr, H3K27me3.bed)
pol2 <- distanceToNearest(windows_gr, pol2.bed)

annotated_df$nearest_H3K36me3 <- H3K36me3@elementMetadata$distance
annotated_df$nearest_H3K27me3 <- H3K27me3@elementMetadata$distance
annotated_df$nearest_pol2<- pol2@elementMetadata$distance

MRC <- annotated_df[group == "MRC",]
MRCr <- IRanges(start = MRC$start, end = MRC$end, names = MRC$seqnames)
annotated_df$ROC_rank_H3K36me3 <- 0
annotated_df$ROC_rank_H3K27me3 <- 0
annotated_df$ROC_rank_pol2 <- 0
for (i in 1:nrow(annotated_df[group!="MRC",])){
  rdist <- annotated_df[i,]$nearest_H3K36me3
  h3k4dist <- annotated_df[i,]$nearest_H3K27me3
  h3k27dist <- annotated_df[i,]$nearest_pol2
  
  exp <- IRanges(start = annotated_df[i,]$start, end = annotated_df[i,]$end, names = annotated_df[i,]$seqnames)
  nearest1 <- nearest(exp,MRCr)
  
  n1_dens1 <- MRC[nearest1,]$nearest_H3K36me3
  n1_dens2 <- MRC[nearest1,]$nearest_H3K27me3
  n1_dens3 <- MRC[nearest1,]$nearest_pol2
  
  nearest2 <- nearest(exp,MRCr[-nearest1])
  n2_dens1 <- MRC[nearest2,]$nearest_H3K36me3
  n2_dens2 <- MRC[nearest2,]$nearest_H3K27me3
  n2_dens3 <- MRC[nearest2,]$nearest_pol2
  
  nearest3 <- nearest(exp,MRCr[-c(nearest1,nearest2)])
  n3_dens1 <- MRC[nearest3,]$nearest_H3K36me3
  n3_dens2 <- MRC[nearest3,]$nearest_H3K27me3
  n3_dens3 <- MRC[nearest3,]$nearest_pol2
  
  densities1 <- c(rdist,n1_dens1,n2_dens1,n3_dens1)
  densities2 <- c(h3k4dist,n1_dens2,n2_dens2,n3_dens2)
  densities3 <- c(h3k27dist,n1_dens3,n2_dens3,n3_dens3)
  
  names(densities1) <- c("exp","MRC1","MRC2","MRC3")
  names(densities2) <- c("exp","MRC1","MRC2","MRC3")
  names(densities3) <- c("exp","MRC1","MRC2","MRC3")
  
  densities1 <- sort(densities1,decreasing = F)
  densities2 <- sort(densities2,decreasing = F)
  densities3 <- sort(densities3,decreasing = F)
  
  
  rank1 = which(names(densities1)=="exp")
  rank2 = which(names(densities2)=="exp")
  rank3 = which(names(densities3)=="exp")
  
  
  annotated_df[i,]$ROC_rank_H3K36me3 <- rank1
  annotated_df[i,]$ROC_rank_H3K27me3 <- rank2
  annotated_df[i,]$ROC_rank_pol2 <- rank3
  
}

saveRDS(annotated_df, file = "df_density_all_comps.rds")


annotated_df[,MRC_below_H3K36me3:=4-ROC_rank_H3K36me3]
annotated_df[,MRC_below_H3K27me3:=4-ROC_rank_H3K27me3]
annotated_df[,MRC_below_pol2:=4-ROC_rank_pol2]

annotated_df <- readRDS("df_density_all_comps.rds")

# also add GC content to all this
# also do ALU repeats


coords_gr <- makeGRangesFromDataFrame(annotated_df, keep.extra.columns = TRUE, start.field = "start", end.field = "end", seqnames.field = "seqnames")
seqinfo(coords_gr) <- seqinfo(txdb)[seqlevels(coords_gr)]

windows_gr1 <- resize(coords_gr, width = 1e6, fix = "center")
windows_gr2 <- resize(coords_gr, width = 1e5, fix = "center")
windows_gr3 <- resize(coords_gr, width = 1e4, fix = "center")
windows_gr4 <- resize(coords_gr, width = 1e3, fix = "center")

windows_gr1 <- GenomicRanges::trim(windows_gr1)
windows_gr2 <- GenomicRanges::trim(windows_gr2)
windows_gr3 <- GenomicRanges::trim(windows_gr3)
windows_gr4 <- GenomicRanges::trim(windows_gr4)

#have to make it chunky

chunk_size <- 10000  # Tune this depending on your available memory
chunks <- split(windows_gr4, ceiling(seq_along(windows_gr4) / chunk_size))

get_gc_content <- function(gr_chunk, genome) {
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, gr_chunk)
  gc <- letterFrequency(seqs, letters = c("G", "C"), as.prob = TRUE)
  rowSums(gc)
}
gc_content <- unlist(lapply(chunks, get_gc_content, genome = genome))

annotated_df$gc_content_1kB <- gc_content

saveRDS(annotated_df, file = "df_density_all_comps_gc.rds")


#now lets deal with alu repeats
# get an alu bedfile from ucsc
alu_gr <- import("rmsk_alu_hg19.bed")
chunk_size <- 10000

count_alus_in_chunk <- function(gr_chunk, alu_gr) {
  hits <- findOverlaps(gr_chunk, alu_gr)
  tab <- table(queryHits(hits))
  counts <- integer(length(gr_chunk))
  counts[as.integer(names(tab))] <- as.integer(tab)
  counts
}

# the paper did alus in 100bp to 10kb
windows_gr1 <- resize(coords_gr, width = 1e4, fix = "center")
windows_gr2 <- resize(coords_gr, width = 1e3, fix = "center")
windows_gr3 <- resize(coords_gr, width = 100, fix = "center")


windows_gr1 <- GenomicRanges::trim(windows_gr1)
windows_gr2 <- GenomicRanges::trim(windows_gr2)
windows_gr3 <- GenomicRanges::trim(windows_gr3)

#have to make it chunky

chunk_size <- 10000  # Tune this depending on your available memory
chunks1 <- split(windows_gr1, ceiling(seq_along(windows_gr1) / chunk_size))
chunks2 <- split(windows_gr2, ceiling(seq_along(windows_gr2) / chunk_size))
chunks3 <- split(windows_gr3, ceiling(seq_along(windows_gr3) / chunk_size))



alu_counts1 <- unlist(lapply(chunks1, count_alus_in_chunk, alu_gr = alu_gr))
alu_counts2 <- unlist(lapply(chunks2, count_alus_in_chunk, alu_gr = alu_gr))
alu_counts3 <- unlist(lapply(chunks3, count_alus_in_chunk, alu_gr = alu_gr))


annotated_df$alu_content_10kB <- alu_counts1
annotated_df$alu_content_1kB <- alu_counts2
annotated_df$alu_content_100B <- alu_counts3

saveRDS(annotated_df, file = "df_density_all_comps_gc_alu.rds")

MRC <- annotated_df[group == "MRC",]
MRCr <- IRanges(start = MRC$start, end = MRC$end, names = MRC$seqnames)
annotated_df$ROC_rank_gc_content_1MB <- 0
annotated_df$ROC_rank_gc_content_100kB <- 0
annotated_df$ROC_rank_gc_content_10kB <- 0
annotated_df$ROC_rank_gc_content_1kB <- 0
annotated_df$ROC_rank_alu_content_10kB <- 0
annotated_df$ROC_rank_alu_content_1kB <- 0
annotated_df$ROC_rank_alu_content_100B <- 0


#ranking Gc and alu content
for (i in 1:nrow(annotated_df[group!="MRC",])){
  gc_content_1MB <- annotated_df[i,]$gc_content_1MB
  gc_content_100kB <- annotated_df[i,]$gc_content_100kB
  gc_content_10kB <- annotated_df[i,]$gc_content_10kB
  gc_content_1kB <- annotated_df[i,]$gc_content_1kB
  alu_content_10kB <- annotated_df[i,]$alu_content_10kB
  alu_content_1kB <- annotated_df[i,]$alu_content_1kB
  alu_content_100B <- annotated_df[i,]$alu_content_100B
  
  exp <- IRanges(start = annotated_df[i,]$start, end = annotated_df[i,]$end, names = annotated_df[i,]$seqnames)
  nearest1 <- nearest(exp,MRCr)
  
  n1_dens1 <- MRC[nearest1,]$gc_content_1MB
  n1_dens2 <- MRC[nearest1,]$gc_content_100kB
  n1_dens3 <- MRC[nearest1,]$gc_content_10kB
  n1_dens4 <- MRC[nearest1,]$gc_content_1kB
  n1_dens5 <- MRC[nearest1,]$alu_content_10kB
  n1_dens6 <- MRC[nearest1,]$alu_content_1kB
  n1_dens7 <- MRC[nearest1,]$alu_content_100B
  
  nearest2 <- nearest(exp,MRCr[-nearest1])
  n2_dens1 <- MRC[nearest2,]$gc_content_1MB
  n2_dens2 <- MRC[nearest2,]$gc_content_100kB
  n2_dens3 <- MRC[nearest2,]$gc_content_10kB
  n2_dens4 <- MRC[nearest2,]$gc_content_1kB
  n2_dens5 <- MRC[nearest2,]$alu_content_10kB
  n2_dens6 <- MRC[nearest2,]$alu_content_1kB
  n2_dens7 <- MRC[nearest2,]$alu_content_100B
  
  nearest3 <- nearest(exp,MRCr[-c(nearest1,nearest2)])
  n3_dens1 <- MRC[nearest3,]$gc_content_1MB
  n3_dens2 <- MRC[nearest3,]$gc_content_100kB
  n3_dens3 <- MRC[nearest3,]$gc_content_10kB
  n3_dens4 <- MRC[nearest3,]$gc_content_1kB
  n3_dens5 <- MRC[nearest3,]$alu_content_10kB
  n3_dens6 <- MRC[nearest3,]$alu_content_1kB
  n3_dens7 <- MRC[nearest3,]$alu_content_100B
  
  densities1 <- c(gc_content_1MB,n1_dens1,n2_dens1,n3_dens1)
  densities2 <- c(gc_content_100kB,n1_dens2,n2_dens2,n3_dens2)
  densities3 <- c(gc_content_10kB,n1_dens3,n2_dens3,n3_dens3)
  densities4 <- c(gc_content_1kB,n1_dens3,n2_dens3,n3_dens3)
  densities5 <- c(alu_content_10kB,n1_dens3,n2_dens3,n3_dens3)
  densities6 <- c(alu_content_1kB,n1_dens3,n2_dens3,n3_dens3)
  densities7 <- c(alu_content_100B,n1_dens3,n2_dens3,n3_dens3)
  
  names(densities1) <- c("exp","MRC1","MRC2","MRC3")
  names(densities2) <- c("exp","MRC1","MRC2","MRC3")
  names(densities3) <- c("exp","MRC1","MRC2","MRC3")
  names(densities4) <- c("exp","MRC1","MRC2","MRC3")
  names(densities5) <- c("exp","MRC1","MRC2","MRC3")
  names(densities6) <- c("exp","MRC1","MRC2","MRC3")
  names(densities7) <- c("exp","MRC1","MRC2","MRC3")
  
  densities1 <- sort(densities1,decreasing = T)
  densities2 <- sort(densities2,decreasing = T)
  densities3 <- sort(densities3,decreasing = T)
  densities4 <- sort(densities4,decreasing = T)
  densities5 <- sort(densities5,decreasing = T)
  densities6 <- sort(densities6,decreasing = T)
  densities7 <- sort(densities7,decreasing = T)
  
  
  rank1 = which(names(densities1)=="exp")
  rank2 = which(names(densities2)=="exp")
  rank3 = which(names(densities3)=="exp")
  rank4 = which(names(densities4)=="exp")
  rank5 = which(names(densities5)=="exp")
  rank6 = which(names(densities6)=="exp")
  rank7 = which(names(densities7)=="exp")
  
  annotated_df[i,]$ROC_rank_gc_content_1MB <- rank1
  annotated_df[i,]$ROC_rank_gc_content_100kB <- rank2
  annotated_df[i,]$ROC_rank_gc_content_10kB <- rank3
  annotated_df[i,]$ROC_rank_gc_content_1kB <- rank4
  annotated_df[i,]$ROC_rank_alu_content_10kB <- rank5
  annotated_df[i,]$ROC_rank_alu_content_1kB <- rank6
  annotated_df[i,]$ROC_rank_alu_content_100B <- rank7
  
}


annotated_df[,MRC_below_gc_1MB:=4-ROC_rank_gc_content_1MB]
annotated_df[,MRC_below_gc_100kB:=4-ROC_rank_gc_content_100kB]
annotated_df[,MRC_below_gc_10kB:=4-ROC_rank_gc_content_10kB]
annotated_df[,MRC_below_gc_1kB:=4-ROC_rank_gc_content_1kB]
annotated_df[,MRC_below_alu_10kB:=4-ROC_rank_alu_content_10kB]
annotated_df[,MRC_below_alu_1kB:=4-ROC_rank_alu_content_1kB]
annotated_df[,MRC_below_alu_100B:=4-ROC_rank_alu_content_100B]




annotated_df <- readRDS("df_density_all_comps_gc_alu2.rds")


# making it into a heatmap
test <- annotated_df[group!="MRC",] %>%  # have to exclude MRC here to not get fucky wucky with the data
  mutate(MRC_below_div_1MB = MRC_below_1MB /3,
         MRC_below_div_100kB = MRC_below_100kB /3,
         MRC_below_div_10kB = MRC_below_10kB /3,
         MRC_below_div_1kB = MRC_below_1kB /3,
         MRC_below_div_rloop = MRC_below_rloop/3,
         MRC_below_div_h3k4me1 = MRC_below_h3k4me1 /3,
         MRC_below_div_h3k27ac = MRC_below_h3k27ac /3,
         MRC_below_div_H3K36me3 = MRC_below_H3K36me3/3,
         MRC_below_div_H3K27me3 = MRC_below_H3K27me3 /3,
         MRC_below_div_pol2 = MRC_below_pol2 /3,
         MRC_below_div_gc_1MB = MRC_below_gc_1MB /3,
         MRC_below_div_gc_100kB = MRC_below_gc_100kB /3,
         MRC_below_div_gc_10kB = MRC_below_gc_10kB /3,
         MRC_below_div_gc_1kB = MRC_below_gc_1kB /3,
         MRC_below_div_alu_10kB = MRC_below_alu_10kB /3,
         MRC_below_div_alu_1kB = MRC_below_alu_1kB /3,
         MRC_below_div_alu_100B = MRC_below_alu_100B /3) %>% 
  group_by(group) %>%
  mutate(counts = n()) %>%
  ungroup() %>%
  group_by(group) %>%
  mutate(ROC_value_1MB = sum(MRC_below_div_1MB)/counts,
         ROC_value_100kB = sum(MRC_below_div_100kB)/counts,
         ROC_value_10kB = sum(MRC_below_div_10kB)/counts,
         ROC_value_1kB = sum(MRC_below_div_1kB)/counts,
         ROC_value_rloop = sum(MRC_below_div_rloop)/counts,
         ROC_value_H3K4me1 = sum(MRC_below_div_h3k4me1)/counts,
         ROC_value_H3K27ac = sum(MRC_below_div_h3k27ac)/counts,
         ROC_value_H3K36me3 = sum(MRC_below_div_H3K36me3)/counts,
         ROC_value_H3K27me3 = sum(MRC_below_div_H3K27me3)/counts,
         ROC_value_pol2 = sum(MRC_below_div_pol2)/counts,
         ROC_value_gc_1MB = sum(MRC_below_div_gc_1MB)/counts,
         ROC_value_gc_100kB = sum(MRC_below_div_gc_100kB)/counts,
         ROC_value_gc_10kB = sum(MRC_below_div_gc_10kB)/counts,
         ROC_value_gc_1kB = sum(MRC_below_div_gc_1kB)/counts,
         ROC_value_alu_10kB = sum(MRC_below_div_alu_10kB)/counts,
         ROC_value_alu_1kB = sum(MRC_below_div_alu_1kB)/counts,
         ROC_value_alu_100B = sum(MRC_below_div_alu_100B)/counts
         )



# draw an actual heatmap


df <- test[!test$group %in% c("Active mock","Resting mock"),] %>% 
  group_by(group) %>% 
  summarise(Gene_Density_1MB = mean(ROC_value_1MB),
            Gene_density_100kB = mean(ROC_value_100kB),
            Gene_density_10kB = mean(ROC_value_10kB),
            Gene_density_1kB = mean(ROC_value_1kB),
            Distance_to_rloop = mean(ROC_value_rloop),
            Distance_to_H3K4me1 = mean(ROC_value_H3K4me1),
            Distance_to_H3K27ac = mean(ROC_value_H3K27ac),
            Distance_to_H3K36me3 = mean(ROC_value_H3K36me3),
            Distance_to_H3K27me3 = mean(ROC_value_H3K27me3),
            Distance_to_pol2 = mean(ROC_value_pol2),
            GC_content_1MB = mean(ROC_value_gc_1MB),
            GC_content_100kB = mean(ROC_value_gc_100kB),
            GC_content_10kB = mean(ROC_value_gc_10kB),
            GC_content_1kB = mean(ROC_value_gc_1kB),
            Alu_repeats_10kB = mean(ROC_value_alu_10kB),
            Alu_repeats_1kB = mean(ROC_value_alu_1kB),
            Alu_repeats_100B = mean(ROC_value_alu_100B)
            )

# fuck it now w have to melt this

df2 <- melt(df)

# do the stats

annotated_df <- annotated_df[group!="MRC",]
stats_df <- data.table()


dt <- annotated_df[!group %in% c("Active mock","Resting mock"),c("group","ROC_rank_H3K27ac")]
current <- colnames(dt)[2]
colnames(dt) <- c("group","rank")
dt[, rank := factor(rank, levels=1:4)]
alphas <- c(0.05, 0.01, 0.001)
df     <- length(levels(dt$rank)) - 1  # here 4 ranks ⇒ df = 3
results <- dt[, {
  # 1) observed counts and χ² test
  obs  <- as.numeric(table(rank))
  test <- chisq.test(obs, p = rep(1/4, 4))
  chi2 <- unname(test$statistic)
  pval <- test$p.value
  pval <- p.adjust(pval, method="holm")
  # 2) observed Cramér’s V
  n    <- .N
  Vobs <- sqrt(chi2 / (n * df))
  # 3) critical V thresholds for each alpha
  #    χ²_{crit} = qchisq(1 - α, df) ⇒ Vcrit = sqrt(χ²_{crit} / (n * df))
  chi2_crit <- qchisq(1 - alphas, df)
  Vcrit     <- sqrt(chi2_crit / (n * df))
  # return a one-row data.table
  .(
    chi2    = chi2,
    p.value = pval,
    V.obs   = Vobs,
    V.0.05  = Vcrit[1],
    V.0.01  = Vcrit[2],
    V.0.001 = Vcrit[3],
    sig.0.05  = Vobs > Vcrit[1],
    sig.0.01  = Vobs > Vcrit[2],
    sig.0.001 = Vobs > Vcrit[3]
  )
}, by = group]
results$variable <- current
stats_df <- rbind(stats_df,results)
colnames(annotated_df)
table(stats_df$variable)

# okay now get a graphing column in here
#stats_df$padjust <- p.adjust(stats_df$p.value, method = "holm")
stats_df <- as.data.table(stats_df)
stats_df[,label := ifelse(sig.0.001, "***",
                          ifelse(sig.0.01,"**",
                                 ifelse(sig.0.05,"*","")))]

stats_df[,merge_column:=ifelse(variable == "ROC_rank_dens_100kB","Gene_density_100kB",
                               ifelse(variable == "ROC_rank_dens_10kB","Gene_density_10kB",
                                      ifelse(variable == "ROC_rank_dens_1kB","Gene_density_1kB",
                                             ifelse(variable == "ROC_rank_dens_1MB","Gene_Density_1MB",
                                                    ifelse(variable== "ROC_rank_alu_content_100B","Alu_repeats_100B",
                                                           ifelse(variable == "ROC_rank_alu_content_1kB","Alu_repeats_1kB",
                                                                  ifelse(variable == "ROC_rank_alu_content_10kB","Alu_repeats_10kB",
                                                                         ifelse(variable == "ROC_rank_gc_content_1MB","GC_content_1MB",
                                                                                ifelse(variable == "ROC_rank_gc_content_100kB","GC_content_100kB",
                                                                                       ifelse(variable == "ROC_rank_gc_content_10kB","GC_content_10kB",
                                                                                              ifelse(variable == "ROC_rank_gc_content_1kB","GC_content_1kB",
                                                                                                     ifelse(variable =="ROC_rank_H3K27ac","Distance_to_H3K27ac",
                                                                                                            ifelse(variable == "ROC_rank_H3K27me3","Distance_to_H3K27me3",
                                                                                                                   ifelse(variable == "ROC_rank_H3K36me3","Distance_to_H3K36me3",
                                                                                                                          ifelse(variable == "ROC_rank_H3K4me1","Distance_to_H3K4me1",
                                                                                                                                 ifelse(variable == "ROC_rank_pol2","Distance_to_pol2",
                                                                                                                                        ifelse(variable == "ROC_rank_rloop","Distance_to_rloop","mistake"))))))))))))))))
)]




# the one made using hci square test goodnes of fit and holm correction
saveRDS(stats_df, "Stats_of_Miseq_groups_Xsq.rds")
stats_df <- readRDS("Stats_of_Miseq_groups_Xsq_kramer.rds")

# calculate kramer's V

tbl <- table(annotated_df[!group %in% c("Active mock","Resting mock","MRC"),]$group)
gsizedf <- data.table(group = names(tbl), gsize = tbl)

stats_df <- merge(stats_df,gsizedf, by.x= c("group"), by.y = c("group"
))
stats_df[,effect_size:= sqrt(chi2/(gsize.N * 3))]  # kramaers V : 4 columns so its 4-1 = 3

#  0.06  small, 0.17 medium 0.29 large effect
stats_df[,label2 := ifelse(effect_size >= 0.29, "***",
                           ifelse(effect_size <= 0.29 & effect_size >= 0.17 ,"**",
                                  ifelse(effect_size <= 0.17 & effect_size >= 0.06,"**","")))]

# Here we do the statistics concerning activated vs resting comparisons
annotated_df <- annotated_df[group!="MRC",]
stats_df2 <- data.table()

# 
dt <- annotated_df[!group %in% c("Active mock","Resting mock"),c("group","alu_content_100B")]
current <- colnames(dt)[2]
colnames(dt) <- c("group","var")
#  normality test
shapiro.test(dt[group == "Active Nup153KO"]$var) # def cant assume normality

# run tests 
t1 <- wilcox.test(dt[group == "Active NT"]$var,dt[group == "Resting NT"]$var)
t2 <- wilcox.test(dt[group == "Active Nup153KO"]$var,dt[group == "Resting Nup153KO"]$var)
t3 <- wilcox.test(dt[group == "Active Nup62KO"]$var,dt[group == "Resting Nup62KO"]$var)
p_vals <- p.adjust(p=c(t1$p.value,t2$p.value,t3$p.value))
results <- data.table(pval = p_vals,
                      group = c("Resting NT","Resting Nup153KO","Resting Nup62KO"),
                      # these are added together and used only for effect size calculation
                      gsize =c(nrow(dt[group == "Resting NT"]) + nrow(dt[group == "Active NT"]),
                               nrow(dt[group == "Resting Nup153KO"])+nrow(dt[group == "Active Nup153KO"]),
                               nrow(dt[group == "Resting Nup62KO"]))+nrow(dt[group == "Active Nup62KO"]),
                      statistic=c(t1$statistic,t2$statistic,t3$statistic),
                      variable = current)
# calc effect size?
results[,effsize:=statistic / (sqrt(gsize))]

stats_df2 <- rbind(stats_df2,results)

# visualize
ggplot(dt, aes(x=group, y=var)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4)
colnames(annotated_df)

# forming the labels for plotting
stats_df2[,label := ifelse(pval <= 0.001, "***",
                          ifelse(pval <= 0.01,"**",
                                 ifelse(pval <= 0.05,"*","")))]
#thats all extremely significant;
# perhaps I need to consider effec size labels insted- but what does that eve mean for signed rank tests?

# just saving for now until requested!
saveRDS(stats_df2, "stats_miseq_rest_vs_act.rds")



df4 <- merge(df2, stats_df, by.x= c("group","variable"), by.y = c("group","merge_column"))

# fix stupid df4 gene density misprint
df4 <- as.data.table(df4)
df4[df4$variable == "Gene_Density_1MB",]$variable <- "Gene_density_1MB"

df4$variable <- factor(df4$variable)
df4$variable <- forcats::fct_relevel(df4$variable,"Gene_density_1MB")
levels(df4$variable)


proc2 <- ggplot(df4, aes(x = group, y = variable, fill = value)) +
  geom_tile(color = "white", width = 1, height = 1) +
  geom_text(aes(label = label2), size = 3) +
  ylab("")+
  xlab("")+
  scale_fill_gradient2(
    low = "blue",      # Color at 0
    mid = "white",     # Color at 0.5
    high = "red",      # Color at 1
    midpoint = 0.5,    # Set center point
    limits = c(0, 1),  # Fix scale from 0 to 1
    #oob = squish,                   # prevent color bleed
    breaks = c(0, 0.5, 1),
    labels = c("Depleted compared to random","", "Enriched compared to random")
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=0.8)
    
  ) +
  labs(
    title = "Comparison to random",
    fill = ""
  )

proc2



#split by sample
# new grouping

annotated_df[,group2:= ifelse(V7%in% c("D183_mock","D184_mock"),"Resting mock",
                             ifelse(V7%in% c("D183A_mock","D185A_mock"),"Active mock",
                                    ifelse(V7%in% c("D183_NT","D184_NT"),"Resting NT",
                                           ifelse(V7%in% c("D183A_NT","D185A_NT"),"Active NT",
                                                  ifelse(V7 ==  "D183_Nup153KO","D183_Nup153KO",
                                                         ifelse(V7 ==  "D184_Nup153KO","D184_Nup153KO",
                                                                ifelse(V7 ==  "D185A_Nup153KO","D185A_Nup153KO",
                                                                       ifelse(V7 ==  "D183A_Nup153KO","D183A_Nup153KO",
                                                                              ifelse(V7 ==  "D185A_Nup62KO","D185A_Nup62KO",
                                                                                     ifelse(V7 ==  "D183A_Nup62KO","D183A_Nup62KO",
                                                                                            ifelse(V7 ==  "D184_Nup62KO","D184_Nup62KO",
                                                                                                   ifelse(V7 ==  "D183_Nup62KO","D183_Nup62KO","MRC"
                                                                       ))))))))))))]



# making it into a heatmap
test <- annotated_df[group!="MRC",] %>%  # have to exclude MRC here to not get fucky wucky with the data
  mutate(MRC_below_div_1MB = MRC_below_1MB /3,
         MRC_below_div_100kB = MRC_below_100kB /3,
         MRC_below_div_10kB = MRC_below_10kB /3,
         MRC_below_div_1kB = MRC_below_1kB /3,
         MRC_below_div_rloop = MRC_below_rloop/3,
         MRC_below_div_h3k4me1 = MRC_below_h3k4me1 /3,
         MRC_below_div_h3k27ac = MRC_below_h3k27ac /3,
         MRC_below_div_H3K36me3 = MRC_below_H3K36me3/3,
         MRC_below_div_H3K27me3 = MRC_below_H3K27me3 /3,
         MRC_below_div_pol2 = MRC_below_pol2 /3,
         MRC_below_div_gc_1MB = MRC_below_gc_1MB /3,
         MRC_below_div_gc_100kB = MRC_below_gc_100kB /3,
         MRC_below_div_gc_10kB = MRC_below_gc_10kB /3,
         MRC_below_div_gc_1kB = MRC_below_gc_1kB /3,
         MRC_below_div_alu_10kB = MRC_below_alu_10kB /3,
         MRC_below_div_alu_1kB = MRC_below_alu_1kB /3,
         MRC_below_div_alu_100B = MRC_below_alu_100B /3) %>% 
  group_by(group) %>%
  mutate(counts = n()) %>%
  ungroup() %>%
  group_by(group) %>%
  mutate(ROC_value_1MB = sum(MRC_below_div_1MB)/counts,
         ROC_value_100kB = sum(MRC_below_div_100kB)/counts,
         ROC_value_10kB = sum(MRC_below_div_10kB)/counts,
         ROC_value_1kB = sum(MRC_below_div_1kB)/counts,
         ROC_value_rloop = sum(MRC_below_div_rloop)/counts,
         ROC_value_H3K4me1 = sum(MRC_below_div_h3k4me1)/counts,
         ROC_value_H3K27ac = sum(MRC_below_div_h3k27ac)/counts,
         ROC_value_H3K36me3 = sum(MRC_below_div_H3K36me3)/counts,
         ROC_value_H3K27me3 = sum(MRC_below_div_H3K27me3)/counts,
         ROC_value_pol2 = sum(MRC_below_div_pol2)/counts,
         ROC_value_gc_1MB = sum(MRC_below_div_gc_1MB)/counts,
         ROC_value_gc_100kB = sum(MRC_below_div_gc_100kB)/counts,
         ROC_value_gc_10kB = sum(MRC_below_div_gc_10kB)/counts,
         ROC_value_gc_1kB = sum(MRC_below_div_gc_1kB)/counts,
         ROC_value_alu_10kB = sum(MRC_below_div_alu_10kB)/counts,
         ROC_value_alu_1kB = sum(MRC_below_div_alu_1kB)/counts,
         ROC_value_alu_100B = sum(MRC_below_div_alu_100B)/counts
  )
df <- test[!test$group2 %in% c("Active mock","Resting mock"),] %>% 
  group_by(group2) %>% 
  summarise(Gene_Density_1MB = mean(ROC_value_1MB),
            Gene_density_100kB = mean(ROC_value_100kB),
            Gene_density_10kB = mean(ROC_value_10kB),
            Gene_density_1kB = mean(ROC_value_1kB),
            Distance_to_rloop = mean(ROC_value_rloop),
            Distance_to_H3K4me1 = mean(ROC_value_H3K4me1),
            Distance_to_H3K27ac = mean(ROC_value_H3K27ac),
            Distance_to_H3K36me3 = mean(ROC_value_H3K36me3),
            Distance_to_H3K27me3 = mean(ROC_value_H3K27me3),
            Distance_to_pol2 = mean(ROC_value_pol2),
            GC_content_1MB = mean(ROC_value_gc_1MB),
            GC_content_100kB = mean(ROC_value_gc_100kB),
            GC_content_10kB = mean(ROC_value_gc_10kB),
            GC_content_1kB = mean(ROC_value_gc_1kB),
            Alu_repeats_10kB = mean(ROC_value_alu_10kB),
            Alu_repeats_1kB = mean(ROC_value_alu_1kB),
            Alu_repeats_100B = mean(ROC_value_alu_100B)
  )

# fuck it now w have to melt this

df2 <- melt(df)

df2$group2 <- factor(df2$group2,
                                     levels = c("Resting NT","Active NT",
                                                "D183_Nup153KO","D184_Nup153KO","D183_Nup62KO","D184_Nup62KO",
                                                "D183A_Nup153KO","D185A_Nup153KO","D183A_Nup62KO","D185A_Nup62KO"
))

# now the stats of these as well
annotated_df <- annotated_df[group!="MRC",]
stats_df <- data.table()


dt <- annotated_df[!group2 %in% c("Active mock","Resting mock"),c("group2","ROC_rank_alu_content_100B")]
current <- colnames(dt)[2]
colnames(dt) <- c("group","rank")
dt[, rank := factor(rank, levels=1:4)]
alphas <- c(0.05, 0.01, 0.001)
df     <- length(levels(dt$rank)) - 1  # here 4 ranks ⇒ df = 3
results <- dt[, {
  # 1) observed counts and χ² test
  obs  <- as.numeric(table(rank))
  test <- chisq.test(obs, p = rep(1/4, 4))
  chi2 <- unname(test$statistic)
  pval <- test$p.value
  pval <- p.adjust(pval, method="holm")
  # 2) observed Cramér’s V
  n    <- .N
  Vobs <- sqrt(chi2 / (n * df))
  # 3) critical V thresholds for each alpha
  #    χ²_{crit} = qchisq(1 - α, df) ⇒ Vcrit = sqrt(χ²_{crit} / (n * df))
  chi2_crit <- qchisq(1 - alphas, df)
  Vcrit     <- sqrt(chi2_crit / (n * df))
  # return a one-row data.table
  .(
    chi2    = chi2,
    p.value = pval,
    V.obs   = Vobs,
    V.0.05  = Vcrit[1],
    V.0.01  = Vcrit[2],
    V.0.001 = Vcrit[3],
    sig.0.05  = Vobs > Vcrit[1],
    sig.0.01  = Vobs > Vcrit[2],
    sig.0.001 = Vobs > Vcrit[3]
  )
}, by = group]
results$variable <- current
stats_df <- rbind(stats_df,results)
colnames(annotated_df)
table(stats_df$variable)

# okay now get a graphing column in here
#stats_df$padjust <- p.adjust(stats_df$p.value, method = "holm")
stats_df <- as.data.table(stats_df)
stats_df[,label := ifelse(sig.0.001, "***",
                          ifelse(sig.0.01,"**",
                                 ifelse(sig.0.05,"*","")))]

stats_df[,merge_column:=ifelse(variable == "ROC_rank_dens_100kB","Gene_density_100kB",
                               ifelse(variable == "ROC_rank_dens_10kB","Gene_density_10kB",
                                      ifelse(variable == "ROC_rank_dens_1kB","Gene_density_1kB",
                                             ifelse(variable == "ROC_rank_dens_1MB","Gene_Density_1MB",
                                                    ifelse(variable== "ROC_rank_alu_content_100B","Alu_repeats_100B",
                                                           ifelse(variable == "ROC_rank_alu_content_1kB","Alu_repeats_1kB",
                                                                  ifelse(variable == "ROC_rank_alu_content_10kB","Alu_repeats_10kB",
                                                                         ifelse(variable == "ROC_rank_gc_content_1MB","GC_content_1MB",
                                                                                ifelse(variable == "ROC_rank_gc_content_100kB","GC_content_100kB",
                                                                                       ifelse(variable == "ROC_rank_gc_content_10kB","GC_content_10kB",
                                                                                              ifelse(variable == "ROC_rank_gc_content_1kB","GC_content_1kB",
                                                                                                     ifelse(variable =="ROC_rank_H3K27ac","Distance_to_H3K27ac",
                                                                                                            ifelse(variable == "ROC_rank_H3K27me3","Distance_to_H3K27me3",
                                                                                                                   ifelse(variable == "ROC_rank_H3K36me3","Distance_to_H3K36me3",
                                                                                                                          ifelse(variable == "ROC_rank_H3K4me1","Distance_to_H3K4me1",
                                                                                                                                 ifelse(variable == "ROC_rank_pol2","Distance_to_pol2",
                                                                                                                                        ifelse(variable == "ROC_rank_rloop","Distance_to_rloop","mistake"))))))))))))))))
)]


# the one made using hci square test goodnes of fit and holm correction
saveRDS(stats_df, "Stats_of_Miseq_groups_Xsq_perg.rds")
stats_df <- readRDS("Stats_of_Miseq_groups_Xsq_perg.rds")



tbl <- table(annotated_df[!group2 %in% c("Active mock","Resting mock","MRC"),]$group2)
gsizedf <- data.table(group = names(tbl), gsize = tbl)

stats_df <- merge(stats_df,gsizedf, by.x= c("group"), by.y = c("group"
))
stats_df[,effect_size:= sqrt(chi2/(gsize.N * 3))]  # kramaers V : 4 columns so its 4-1 = 3

#  0.1 is small, 0.3 medium 0.5 large efect
stats_df[,label2 := ifelse(effect_size >= 0.29, "***",
                           ifelse(effect_size <= 0.29 & effect_size >= 0.17 ,"**",
                                  ifelse(effect_size <= 0.17 & effect_size >= 0.06,"**","")))]


df4 <- merge(df2, stats_df, by.x= c("group2","variable"), by.y = c("group","merge_column"))

# fix stupid df4 gene density misprint
df4 <- as.data.table(df4)
df4[df4$variable == "Gene_Density_1MB",]$variable <- "Gene_density_1MB"

df4$variable <- factor(df4$variable)
df4$variable <- forcats::fct_relevel(df4$variable,"Gene_density_1MB")
levels(df4$variable)

proc3 <- ggplot(df4, aes(x = group2, y = variable, fill = value)) +
  geom_tile(color = "white", width = 1, height = 1) +
  geom_text(aes(label = label2), size = 3) +
  ylab("")+
  xlab("")+
  scale_fill_gradient2(
    low = "blue",      # Color at 0
    mid = "white",     # Color at 0.5
    high = "red",      # Color at 1
    midpoint = 0.5,    # Set center point
    limits = c(0, 1),  # Fix scale from 0 to 1
    #oob = squish,                   # prevent color bleed
    breaks = c(0, 0.5, 1),
    labels = c("Depleted compared to random","", "Enriched compared to random")
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=0.8)
    
  ) +
  labs(
    title = "Comparison to random",
    fill = ""
  )

proc3


# adding stats




pdf("STP&Miseq_ver4.pdf", height=6, width=8)
proc2
proc3
proc
dev.off()

