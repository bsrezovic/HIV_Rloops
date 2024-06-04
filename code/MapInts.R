# file for calling the itegration sites
# based on dunjas previous work

library(data.table)
library(stringr)
library(Biostrings)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)  #fails when downloading???

# list all .bed files obtained at the end of procedure in 1_Read_mapping_and_filtering.md and extract sample names:
setwd("C:/Users/38598/Desktop/hivint/IS_mapping/new_data")
bedfiles <- list.files(pattern="*.bed$")
samplenames <- substr(bedfiles, 1, 11)

# construct the sample table
treatment <- c("AQR_inf_S7","NTC_inf_S6","NTC_mock_S5","AQR_M_S3","AQR_MY_S4","AQR_S2","siCTRL_S1")
#UNFINISHED

#Load in the data from the bed files and distribute it somehow
all_samples_islands <- list()
all_samples_peaks <- data.table()
mean_cov_in_covered_areas <- c()

for(k in 1:length(bedfiles)) {   # process each sample separately
  
  # import .bed file with positions of mapped R1 reads after filtering:
  s <- import(bedfiles[k], format="bed")
  
  # for samples in which filtering by linker sequence match was possible, keep only mappings by R1 reads whose R2 pair contained a match with linker sequence:
  #if(samplenames[k] %in% c("SRR12322275", "SRR12322277", "SRR12322278", "SRR12330757", "SRR12330758", "SRR12330759", "SRR12330760", "SRR12330761", "SRR12330762")) {
  #  linkermatch <- fread(file=paste0(samplenames[k], "_pair_matches_linker.txt"), header=FALSE)
  #  s <- s[s$name %in% linkermatch$V1]
  #}
  
  # since integration is staggered, starts and ends have to be moved for 2 or 3 bp to match the central bp of integration site:
  s_dt <- as.data.table(s)
  s_dt <- s_dt[, start := ifelse(strand == "+", start+2, start-2)
  ][, end := ifelse(strand == "+", end+3, end-3)]
  s <- GRanges(s_dt[, score := NULL])
  seqlevels(s, pruning.mode="coarse") <- seqlevels(BSgenome.Hsapiens.UCSC.hg19)
  seqlengths(s) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
  
  # get coverage, extract areas where there is something mapped (islands), and calculate average coverage on that areas:
  cov <- coverage(s)
  isl <- reduce(s)
  isl$sample <- samplenames[k]
  #isl$group <- treatment[sample == samplenames[k]]$samplegroup
  #isl$treatment <- treatment[sample == samplenames[k]]$conditions
  mean_cov_in_covered_areas[k] <- sum(sum(cov)) / sum(width(isl))
  
  # find most likely IS position for each island:
  isl_cov <- RleViewsList(rleList = cov, rangesList = isl)
  max_position_left <- unlist(start(viewRangeMaxs(isl_cov)))
  max_position_right <- unlist(end(viewRangeMaxs(isl_cov)))
  isl$IS_position <- ifelse(strand(isl) == "+", max_position_left, max_position_right)
  
  # calculate some coverage stats for each island:
  isl$isl_width <- unlist(width(isl_cov))
  isl$isl_mean_cov <- unlist(viewMeans(isl_cov))
  isl$cov_at_max <- unlist(viewMaxs(isl_cov))
  isl$enrichment <- unlist(isl$isl_mean_cov / mean_cov_in_covered_areas[k])
  
  
  # islands wider than average read length are potentially problematic in downstream analysis (during duplicate removal), IF they are saturated (meaning high coverage throughout) - mark them:
  isl$note <- ifelse(width(isl) > quantile(isl$isl_width)[4], "wide", "narrow")   # mark wide peaks
  highcov_cutoff <- quantile(isl$isl_mean_cov)[4]            # decide on cutoff for "high coverage"
  highcoverage <- GRanges(IRanges::slice(cov, highcov_cutoff))        # get high coverage areas
  wideandhigh <- subsetByOverlaps(isl[isl$note == "wide"], highcoverage)   # overlap the two
  
  saturated <- lapply(cov[wideandhigh], function(x) {
    high_cov_area <- sum(as.numeric(x) > highcov_cutoff)   # if high-coverage area...
    length_cutoff <- 0.85*length(as.numeric(x))            # ...spans more than 85% of bases in peak...
    high_cov_area >= length_cutoff                         # ...return true
  })
  saturated_ranges <- wideandhigh[unlist(saturated)]         # mark saturated islands
  isl$note[subjectHits(findOverlaps(saturated_ranges, isl))] <- "saturated"
  
  # since IS positions are more or less approximate, we can define them as peaks including +/- 5 around IS:
  peaks <- isl
  ranges(peaks) <- IRanges(start=isl$IS_position, end=isl$IS_position)
  peaks <- trim(resize(peaks, width=11, fix = "center"))
  all_samples_peaks <- rbind(all_samples_peaks, as.data.table(peaks))
  
  # save intermediary data for later use:
  saveRDS(all_samples_peaks, file="HIV1_all_peaks_new_data.rds") 
}



#analysis

# read in all peaks: 
all_samples_peaks <- readRDS(file="HIV1_all_peaks_new_data.rds")
peaks <- GRanges(all_samples_peaks)
seqlevels(peaks, pruning.mode="coarse") <- seqlevels(BSgenome.Hsapiens.UCSC.hg19)
seqlengths(peaks) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)

# bin the genome and mark in which bins peaks/islands are located:
binsize <- 100
set.seed(1234)
# strange error fixing
tile <- tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19), 
                   tilewidth = binsize, cut.last.tile.in.chrom = TRUE)
#####
bins <- trim(IRanges::shift(tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19), 
                              tilewidth = binsize, cut.last.tile.in.chrom = TRUE), sample(1:(binsize/2), 1)  ))
binoverlap <- findOverlaps(peaks, bins)
group <- as.data.table(binoverlap)[, .SD[1], by=queryHits]$subjectHits   # takes care of peaks which are in multiple bins

peakindices <- 1:length(peaks)
peakindicesfound <- unique(queryHits(binoverlap))
peaks <- peaks[peakindices %in% peakindicesfound]                        # takes care of peaks which are not in any bins (e.g. at last <100bp of chromosome)

# mark bins which have multiple peaks assigned to them:
peaks$bin <- group
peaks <- as.data.table(peaks)
binindices <- peaks[, unique(bin), by=sample]$V1   # unique(bin) stops it from counting in-sample duplicates
duplicatedbins <- unique(binindices[duplicated(binindices)])

# remove in-sample duplicates:
peaks <- peaks[order(sample, bin, enrichment)][, .SD[.N], by=c("sample", "bin")]

# remove all duplicated peaks if one of them is saturated:
uniquepeaks <- peaks[!(bin %in% duplicatedbins)]   
duplicatedpeaks <- peaks[bin %in% duplicatedbins]
removebins <- duplicatedpeaks[note == "saturated"]$bin
duplicatedpeaks_nosat <- duplicatedpeaks[!(bin %in% removebins)]

# out of all remaining duplicates, keep the one with biggest enrichment:
duplicatedpeaks_clean <- duplicatedpeaks_nosat[order(bin, enrichment)][, .SD[.N], by=bin]

# format the end result into a nice-looking table:
peaks_final <- rbind(uniquepeaks, duplicatedpeaks_clean)[order(sample, bin)]

# this adds celltypes and treatment groups; which we dont have we have samples and shit
#peaks_final <- 
#  peaks_final[, celltype := gsub("_.*", "", gsub("Jurkat_", "", treatment))
#              ][, virus := ifelse(celltype == "A77V" | celltype == "N74D", celltype, "WT")
#                ][, celltype := ifelse(celltype == "A77V" | celltype == "N74D", "WT", celltype)
#                  ][, group := gsub(" ", "_", group)
#                    ][, group := gsub("/", "", group)
#                      ][, IS_position := start + 5   # center position is considered to be the IS
#                        ][, .(sample, treatment, group, celltype, virus, seqnames, IS_position, strand, isl_width, isl_mean_cov, cov_at_max, enrichment)]
peaks_final <- peaks_final[,IS_position := start+5]# center position is considered to be the IS
fwrite(peaks_final, file="HIV1_IS_clean_new_data.txt", sep="\t")   # save the table
peaks_final <- fread("HIV1_IS_clean_new_data.txt")
# optionally - if you need the insertion sites in .bed format:
peaks_bed <- fread("HIV1_IS_clean_new_data.txt")
peaks_bed <- peaks_bed[, .(seqnames, IS_position, strand, isl_mean_cov,sample)
                       ][strand == "+", ':=' (start = IS_position, end = IS_position+1)
                         ][strand == "-", ':=' (start = IS_position-1, end = IS_position)
                           ][, .(seqnames, start, end, strand, isl_mean_cov,sample)]
fwrite(peaks_bed, "HIV1_IS_new_data.bed", sep="\t", col.names = FALSE)

# group treatments from different experiment groups together (e.g. WT from capsid mutant experiments and LEDGFp75 knockout experiments are counted together):
peaks_group <- peaks_bed[, experiment := ifelse(virus == "A77V" | virus == "N74D", "capsid mutant", celltype)
                         ][, experiment := ifelse(celltype == "LKO" | celltype == "LEDGFKO", "LKO", experiment)
                           ][, experiment := ifelse(celltype == "NC" | celltype == "NT", "WT", experiment)
                             ][, experiment := ifelse(celltype == "IBD-/-", "IBD--", experiment)]

# optionally - if you want to split data by experimant conditions in separate .bed files:
peaks_group_list <- split(peaks_group, by="experiment")
for(i in 1:length(peaks_group_list)) {
  peaks_out <- peaks_group_list[[i]][, .(seqnames, start, end, strand, isl_mean_cov)]
  fwrite(peaks_out, str_c("HIV1_", names(h1_list)[i], "_IS.bed"), sep="\t", col.names = FALSE)
}




# Own work, without dicta dunja

library(biomaRt)
# define biomart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host    = 'grch37.ensembl.org')

atts <- c("ensembl_gene_id", "ensembl_transcript_id", "description", "chromosome_name", "start_position", 
          "end_position", "strand", 
          "gene_biotype",
          "hgnc_symbol")
genes <- as.data.table(getBM(attributes = atts, mart = mart))

# find overlaps between thse two things
overs <- findOverlaps(IRanges(start = peaks_bed$start, end = peaks_bed$end),
             IRanges(start=genes$start_position, end =genes$end_position))
# 1706 oit of 1720 are in genes? ist this to be expexted?

genes[subjectHits(overs)]
peaks_final$IS2 <- peaks_final$IS_position+1

#try with foverlaps
genes <- genes[genes$chromosome_name%in% c(1:22, "X","Y")]

genes$seqnames <- paste("chr",genes$chromosome_name, sep = "")
setkey(genes,seqnames, start_position, end_position)
setkey(peaks_final,seqnames, IS_position,IS2)

ols <- foverlaps(peaks_final, genes, mult="all")
ols  # there is a lot of them 

results <- ols[,c("sample","seqnames","bin","IS_position","isl_width","isl_mean_cov","cov_at_max","enrichment",
                  "strand","start","end","note","ensembl_gene_id","hgnc_symbol","gene_biotype","start_position","end_position")]
results[is.na(results$ensembl_gene_id),]$hgnc_symbol <- "Not in gene"
results[is.na(results$ensembl_gene_id),]$ensembl_gene_id <- "Not in gene"
results <- unique(results)
# this is just if it overlaps gene itself;


#check 
results[results$IS_position==42713271]
unique(results[results$IS_position==189041699])

# add the rloops information
genidf <- readRDS("C:/Users/38598/Desktop/hivint/MajasCode/Rloops/R-loops-code/peaks-analysis/genidf.RDS")

genidf_ex <- genidf[,c("names","rloop","RLOOPS_dist")] #just the name and rloop information
results4 <- merge(results, genidf_ex, by.x = "ensembl_gene_id", by.y = "names", all.x = T)

# distance to nearest rloop
# we have to find nearest rloop 
# above code gets us the RLOOPs dist but just for integration sites that are in genes, how bout the other ones

fwrite(results4, file="HIV_rloops_new_data", sep="\t")

#attempt to use nearest from Iranges

rloops <- genidf[rloop=="rloop",c("names","seqnames","start","end","RloopOverlaps","gene_biotype","Names")]

# attempt to set up a statistical idea of what gets integrated where?

results2 <- fread("HIV_rloops_new_data")

# First lets add the nearest r loop for those not in genes
# then additionally add nubmer of rloops in various ranges 10k 20k 30k 50k
# maybe add rloopoverlaps information and gene biotype somehow?


# 1 distance to nearest for all, but especially the non gene ones ( use all tho)
# so we have is_position from results2 
# but we dont have specific rloop peak position, only gene with rloop position, maybe use the raw peaks?

peaks <-readRDS("C:/Users/38598/Desktop/hivint/MajasCode/Rloops/R-loops-code/Consensus-peaks/consensusPeaks_withMeanSignal.RDS")
grrloops<-GRanges(peaks)

# using distance to nearest
# lets first filter out the non standard chromosomes, basically without this we just have two gnes thah make trouble
# both of them in NTC inf S6 so no reason to e alarmed yet 
results2 <- results_sample
chrnames <- paste("chr",c(1:22, "X","Y"), sep = "")
results2 <- results2[results2$seqnames%in% chrnames]
rgrange <- GRanges(seqnames = results2$seqnames, ranges = IRanges(start = results2$IS_position, 
                                                                  end = results2$IS_position))

dist_to_nearest <- distanceToNearest(rgrange, grrloops)
results2$dist_to_first <- mcols(dist_to_nearest)$distance

# add gene width to compare distnac ecalc numbers with 
results2$gwidth <- results2$end_position - results2$start_position
results2$RLOOPS_dist - results2$dist_to_first   + results2$gwidth

# add information of how many are in the nearest section
# first lets graps the distances
median(results2$RLOOPS_dist, na.rm = T)
median(results2$dist_to_first)


mean(results2$RLOOPS_dist, na.rm = T)
mean(results2$dist_to_first)

# add "no rloop" to the ones not in gene, as your rloop info is from dripcseq which should be genes only
      # so basically werre saying that if they aint in  a gene they aint gor a loop
      # this coincide with my current interpretation of what no rloop means which is "this gene
          # doesent overlapp with a peak"
results2[results2$rloop==""]$rloop <- "no rloop"

# based on this i guess 20 , 50 and 100+ seem like solid breakpoints for this, lets plot this

plot1 <- ggplot(results2, aes(x=rloop, y=RLOOPS_dist)) + 
  geom_boxplot()+ theme_bw() + ylab("")
plot1
plot2 <- ggplot(results2, aes(x=rloop, y=dist_to_first)) + 
  geom_boxplot()+ theme_bw() + ylab("")
plot2
plot3 <- ggplot(results2, aes(x=sample, y=dist_to_first)) + 
  geom_boxplot()+ theme_bw() + ylab("")
plot3
plot4 <- ggplot(results2, aes(x=sample, y=RLOOPS_dist)) + 
  geom_boxplot()+ theme_bw() + ylab("")
plot4

# tbe plots are too big to plot
# select randomly like 80

rownums <- sample(1:nrow(results2), 80)
# make 10 samples of 80 random points between them
results_x<- results2[rownums]
results_x$sampling_x <- 10
results_sample <- rbind(results_sample, results_x)
#saveRDS(results_sample, file="simulated_sampled_10x80.RDS")
# the plots make my distance calculation seem sensible enough
# now we try to add iforatio nof how many R loops are in the regions near this 
# for the new data
setwd("C:/Users/38598/Desktop/hivint/IS_mapping")
results_real <- readRDS("distances_intToRloop.RDS")

results_sample <- simulated1

results_sample$rl10k <- 0
# we might have to do for loop


maxdist <- 10000

for (x in 1:nrow(results_sample)){
  # find distance to nearest, remove that dist from search, continue until none remain
  #print("Doing row:")
 #print(x)
  results3 <- results_sample[x,]
  grr <- grrloops
  rg <- GRanges(seqnames = results3$seqnames, ranges = IRanges(start = results3$IS_position, 
                                                               end = results3$IS_position))
  dist_to_nearest <- distanceToNearest(rg, grr)
  distance <- mcols(dist_to_nearest)$distance
  #print(distance)
  #print(subjectHits(dist_to_nearest))
  # update a counter for how mamy there are
  i = 0
  if(length(distance) == 0){
    distance = 2*maxdist
   }
  while(distance<maxdist){
    i = i + 1
    print(i)
    grr <- grr[-(subjectHits(dist_to_nearest))]
    dist_to_nearest <- distanceToNearest(rg, grr)
    distance <- mcols(dist_to_nearest)$distance
    if (distance>maxdist){
      break
    }
  }
  results_sample[x,]$rl10k <- i
}

# add them all

saveRDS(results_sample, file="distances_sim_n274.RDS")
# this behaves as expected
setwd("C:/Users/38598/Desktop/hivint/IS_mapping/new_data")
results_real <- readRDS("distances_new_data.RDS")

results_real <- results_sample
results_real$sample <- factor(results_real$sample, levels = c("NTC","AQR", "siCTRL","siAQR","siAQR+wtAQR-R","siAQR+Y1196A-R"))
separate= c("NTC","AQR", "siCTRL","siAQR","siAQR+wtAQR-R","siAQR+Y1196A-R")
results_real[sample=="AQR_S2_inf_",]$sample <- "D403 AQR"
results_real[sample=="AQR_S4_inf_",]$sample <- "D410 AQR"
results_real[sample=="NTC_inf_hg1",]$sample <- "D403 NTC"
results_real[sample=="NTC_S3_inf_",]$sample <- "D410 NTC"
results_real[sample=="siAQR_S6_hg",]$sample <- "siAQR+M"
results_real[sample=="siAQRplusM_",]$sample <- "siAQR"
results_real[sample=="siCTRL_hg19",]$sample <- "siCTRL"

# kicking out mock NTC_S% from the reulst thtats usually not represented
results2 <- results2[sample!="NTC mock S5",]
separate <- c("NTC inf S6","AQR inf S7")


# analyisis of the data as above
give.n <- function(x){
  return(c(y = mean(x)-0.1, label = length(x)))
}
# Distance to nearest Rloop from Integration site
plot4 <- ggplot(results_s, aes(x=sample, y=dist_to_first)) + 
  geom_boxplot()+ theme_bw() + ylab("Distance to nearest Rloop")+
  theme(axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1))+
  stat_summary(fun.data = give.n, geom = "text",aes(y=-104000.5))
plot4


# percentage of sites that are within a certain distance of rloop

          # lets start with 
results_real <- results2
results_real[,within5:=ifelse(dist_to_first<=5000, "Closer","Farther")]
results_real[,within10:=ifelse(dist_to_first<=10000, "Closer","Farther")]
results_real[,within20:=ifelse(dist_to_first<=20000, "Closer","Farther")]
# get rekative frequencies of these
results5 <- results_real %>%
  group_by(sample, within5) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
results5$within5 <- factor(results5$within5, levels=c("Closer","Farther"))

results10 <- results_real %>%
  group_by(sample, within10) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
results10$within10 <- factor(results10$within10, levels=c("Closer","Farther"))

results20 <- results_real %>%
  group_by(sample, within20) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
results20$within20 <- factor(results20$within20, levels=c("Closer","Farther"))

# plot these
separate = unique(results_real$sample)
ggplot(results5[results5$sample%in% separate,], aes(fill=within5, y=freq, x=sample)) + 
  geom_bar(position="fill", stat="identity")+ theme_bw() + ylab("Percent Int within 5k bp of rloop")+
  theme(axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1))+ scale_fill_brewer(palette="Paired", direction = -9)

results10$sample <- factor(results10$sample, levels =c("sample 1","sample 2","sample 3","sample 4", "sample 5",
                                                       "sample 6","sample 7","sample 8","sample 9","sample 10"))
plot1 <- ggplot(results10[,], aes(fill=within10, y=freq, x=sample)) + 
  geom_bar(position="fill", stat="identity")+ theme_bw() + ylab("Percent Int within 10k bp of rloop")+
  theme(axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1))+ scale_fill_brewer(palette="Paired", direction = -1)
plot1


ggplot(results20[results20$sample%in% separate,], aes(fill=within20, y=freq, x=sample)) + 
  geom_bar(position="fill", stat="identity")+ theme_bw() + ylab("Percent Int within 20k bp of rloop")+
  theme(axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1))+ scale_fill_brewer(palette="Paired", direction = -1)



# Number of loops within some distance of Integration site


pdf("simulated_within10_n274.pdf", height=6, width=8)
plot1
dev.off()
plot4 <- ggplot(results_real[,], aes(x=sample, y=rl10k)) + 
  geom_boxplot()+ theme_bw() + ylab("Num of Rloops within 10k bp of Integration site")+
  theme(axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1))+
  stat_summary(fun.data = give.n, geom = "text")
plot4
plot4 <- ggplot(results_real[sample%in% separate,], aes(x=sample, y=rl10k)) + 
  geom_boxplot()+ theme_bw() + ylab("Num ofRloops within 10k bp of Integration site")+
  theme(axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1))+
  stat_summary(fun.data = give.n, geom = "text")
plot4

plot4 <- ggplot(results_real[sample%in% separate,], aes(x=sample, y=rl20k)) + 
  geom_boxplot()+ theme_bw() + ylab("Num ofRloops within 20k bp of Integration site")+
  theme(axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1))+
  stat_summary(fun.data = give.n, geom = "text")
plot4
ggplot(results5[results5$sample%in% separate,], aes(fill=within5, y=freq, x=sample)) + 
  geom_bar(position="fill", stat="identity")+ theme_bw() + ylab("Percent Int within 5k bp of rloop")+
  theme(axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1))+ scale_fill_brewer(palette="Paired", direction = -1)

ggplot(results10[results10$sample%in% separate,], aes(fill=within10, y=freq, x=sample)) + 
  geom_bar(position="fill", stat="identity")+ theme_bw() + ylab("Percent Int within 10k bp of rloop")+
  theme(axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1))+ scale_fill_brewer(palette="Paired", direction = -1)

ggplot(results20[results20$sample%in% separate,], aes(fill=within20, y=freq, x=sample)) + 
  geom_bar(position="fill", stat="identity")+ theme_bw() + ylab("Percent Int within 20k bp of rloop")+
  theme(axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1))+ scale_fill_brewer(palette="Paired", direction = -1)


dev.off()

# some stats of these
# are distances to the loops normally distributed

shapiro.test(results2$dist_to_first)  # hell no theyre not
# yes they bbe differing (there shuld be no problem for the group sizes)
kruskal.test(dist_to_first~sample, data = results2)

# Any post hocers?
pairwise.wilcox.test(results2$dist_to_first,results2$sample,
                     p.adjust.method = "BH")



# analyizing the simulated data
results_sample <- simulated1
results_sample$sampling_x <- as.factor(results_sample$sampling_x)

results_sample[,within20:=ifelse(dist_to_first<=20000, "Closer","Farther")]
results_sample[,within50:=ifelse(dist_to_first<=50000, "Closer","Farther")]
results_sample[,within100:=ifelse(dist_to_first<=100000, "Closer","Farther")]

# plotting simulated ones
plot1 <- ggplot(results_sample, aes(x=rloop, y=RLOOPS_dist)) + 
  geom_boxplot()+ theme_bw() + ylab("")
plot1
plot2 <- ggplot(results_sample, aes(x=rloop, y=dist_to_first)) + 
  geom_boxplot()+ theme_bw() + ylab("")
plot2
plot3 <- ggplot(results_sample, aes(x=sampling_x, y=dist_to_first)) + 
  geom_boxplot()+ theme_bw() + ylab("")
plot3
plot4 <- ggplot(results_sample, aes(x=sample, y=RLOOPS_dist)) + 
  geom_boxplot()+ theme_bw() + ylab("")
plot4

# proportion of Rloop vs no rloop 

results10 <- results_sample %>%
  group_by(sampling_x, within10) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
results50 <- results_sample %>%
  group_by(sampling_x, within50) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
results100 <- results_sample %>%
  group_by(sampling_x, within100) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
# plot these
#pdf("Simulated_rloop_distance_x10.pdf", height=6, width=8)
ggplot(results20, aes(fill=within20, y=freq, x=sampling_x)) + 
  geom_bar(position="fill", stat="identity")+ theme_bw() + ylab("Percent Int within 20k bp of rloop")+
  theme(axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1))

ggplot(results50, aes(fill=within50, y=freq, x=sampling_x)) + 
  geom_bar(position="fill", stat="identity")+ theme_bw() + ylab("Percent Int within 50k bp of rloop")+
  theme(axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1))

ggplot(results100, aes(fill=within100, y=freq, x=sampling_x)) + 
  geom_bar(position="fill", stat="identity")+ theme_bw() + ylab("Percent Int within 100k bp of rloop")+
  theme(axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1))
dev.off()


