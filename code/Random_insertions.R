

library(BSgenome.Hsapiens.UCSC.hg19) 
library(data.table)
library(stringr)
library(Biostrings)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19) 
library(parallel)
save_fasta <- function( reads, directory.path, prefix, compress = TRUE ){
  
  if( isTRUE( compress ) ){
    end <- ".fa.gz"
  } else {
    end <- ".fa"
  }
  
  R1 <- writeXStringSet( x = reads[[1]],
                         filepath = paste0( directory.path, "/", prefix, "_R1", end ),
                         compress = compress,
                         format = "fasta",
                         width = 1000)
  
  R2 <- Biostrings::writeXStringSet( x = reads[[2]],
                                     filepath = paste0( directory.path, "/", prefix, "_R2", end ),
                                     compress = compress,
                                     format = "fasta",
                                     width = 1000)
}

trim_seqs <- function( fragments, min.width = 14, max.distance = 1000, max.bp = 150 ){
  max.width <- max.distance + max.bp*2
  filtered <- fragments[ width( fragments ) >= min.width & width(fragments) <= max.width ]
  left <- subseq( x = filtered, start=1, width = pmin( width( filtered ), max.bp ) )
  right <- subseq(x = filtered,
                  start = pmax( 1, width( filtered ) - max.bp + 1 ), width = pmin( width( filtered ), max.bp ) )
  right <- reverseComplement( right )
  names( left ) <- paste0( "sequence_", seq_along( left ) )
  names( right ) <- paste0( "sequence_", seq_along( right ) )
  return( list( left, right ) )
}
bound_check <- function( fragments, genome.obj, include.lower = FALSE ){
  frag.seqnames <- as.character( GenomicRanges::seqnames( fragments ) )
  genome.seqlengths <- GenomeInfoDb::seqlengths( genome.obj )
  if( include.lower ){
    which( GenomicRanges::end( fragments ) > genome.seqlengths[ c( frag.seqnames ) ] | start( fragments ) < 1 )
  } else{
    which( GenomicRanges::end( fragments ) > genome.seqlengths[ c( frag.seqnames ) ] )
  }
}
# bedwell function that calls in numbe of sites and a genom object?
make_fragments <- function( insert.sites,
                            frag.sites = NULL,
                            random = TRUE,
                            mean = 500,
                            sd = 250,
                            genome.obj ){
  
  genome.seqlengths <- seqlengths( genome.obj )
  
  if( !isTRUE( random ) ){
    if( is.null( frag.sites ) ){
      stop( "frag.sites cannot be NULL for non-random fragmentation.", call.=FALSE )
    }
    
    matches <- lapply( X = frag.sites,
                       FUN = function(x){
                         precede( x = insert.sites,
                                  subject = x,
                                  select = "all",
                                  ignore.strand = FALSE ) } )
    
    coordinates <- lapply( X=seq_along (matches ),
                           FUN=function(x){
                             matrix(c( start( insert.sites[ queryHits( matches[[ x ]] ) ] ),
                                       start( frag.sites[[ x ]][ subjectHits( matches[[ x ]] ) ] ) ),
                                    ncol = 2) } )
    
    frag.ranges <- lapply( X = seq_along( coordinates ),
                           FUN=function(x) {
                             names <- seqnames( insert.sites[ queryHits( matches[[ x ]] ) ] )
                             strand <- strand( insert.sites[ queryHits( matches[[ x ]] ) ] )
                             GRanges( seqnames = names,
                                      ranges = IRanges( start = apply( X = coordinates[[ x ]], MARGIN = 1, FUN = min ),
                                                        end = apply( X = coordinates[[ x ]], MARGIN = 1, FUN = max ) ),
                                      strand = strand ) } )
    
    frag.ranges <- do.call(c, frag.ranges)
    
    unmatched <- insert.sites[ !insert.sites %over% frag.ranges, ]
    
    unmatched.pos <- unmatched[ strand( unmatched ) == "+" ]
    pos.seqnames <- as.character( seqnames( unmatched.pos ) )
    end( unmatched.pos ) <- genome.seqlengths[ c( pos.seqnames ) ]
    
    unmatched.neg <- unmatched[ strand( unmatched ) == "-" ]
    neg.seqnames <- as.character( seqnames( unmatched.neg ) )
    start( unmatched.neg ) <- 1
    
    unmatched.frags <- do.call( c, list( unmatched.pos, unmatched.neg ) )
    
    frag.ranges <- do.call( c, list( frag.ranges, unmatched.frags ) )
  }
  
  else {
    lnorm.loc <- log( mean^2 / sqrt( sd^2 + mean^2 ) )
    lnorm.shape <- sqrt( log( 1 + ( sd^2 / mean^2 ) ) )
    
    plus.sites <- insert.sites[ strand( insert.sites ) == "+" ]
    plus.widths <- rlnorm( n = length(plus.sites), meanlog = lnorm.loc, sdlog = lnorm.shape)
    plus.sites <- GRanges( seqnames = seqnames(plus.sites),
                           ranges = IRanges(start = pmax( start( plus.sites ), 1 ),
                                            end = start( plus.sites ) + plus.widths ),
                           strand = "+" )
    
    minus.sites <- insert.sites[ strand( insert.sites ) == "-" ]
    minus.widths <- rlnorm( n = length( minus.sites ), meanlog = lnorm.loc, sdlog = lnorm.shape )
    minus.sites <- GRanges( seqnames = seqnames( minus.sites ),
                            ranges = IRanges(start = pmax( start( minus.sites ) - minus.widths, 1),
                                             end = end( minus.sites ) ),
                            strand = "-" )
    
    frag.ranges <- do.call( c, list( plus.sites, minus.sites ) )
    
    outliers <- bound_check( fragments = frag.ranges,
                             genome.obj = genome.obj,
                             include.lower = FALSE )
    
    if( length( outliers ) != 0 ){
      extreme.frags <- frag.ranges[ outliers ]
      extreme.seqnames <- as.character( seqnames( extreme.frags ) )
      end( extreme.frags ) <- genome.seqlengths[ c( extreme.seqnames ) ]
      frag.ranges <- frag.ranges[ -outliers ]
      frag.ranges <- do.call( c, list( frag.ranges, extreme.frags ) )
    }
  }
  
  frag.ranges <- sortSeqlevels( frag.ranges )
  frag.ranges <- sort( frag.ranges, ignore.strand = TRUE )
  
  return( frag.ranges )
  
}
random_sites <- function( n.sites, genome.obj ){
  
  chr.names <- seqlevels( genome.obj )
  chr.lengths <- seqlengths( genome.obj )
  chr.weights <- chr.lengths / sum( chr.lengths )
  chr.vec <- sample( chr.names, n.sites, prob = chr.weights, replace = TRUE ) # sample chromosome names weighted by size
  site.pos <- sapply( X = chr.vec,
                      FUN = function(x) {
                        sample(x = 1:chr.lengths[[x]], size = 1 )
                      }
  )
  site.strand <- sample( x = c( "+", "-" ), size = n.sites, replace = TRUE )
  sitesG <- GRanges( seqnames = chr.vec,
                     ranges = IRanges( start=site.pos, width=1 ),
                     strand = site.strand )
  return( sitesG )
  
}

simulate_random_data <- function( genome.obj,
                                  re.sites = NULL,
                                  cut.after = 1,
                                  n.sites,
                                  mean,
                                  sd,
                                  min.width = 14,
                                  max.distance = 1000,
                                  max.bp = 150,
                                  iterations = 1,
                                  n.cores = 1,
                                  write.ranges = FALSE,
                                  prefix = NULL,
                                  directory.path = NULL,
                                  compress = TRUE,
                                  collapse = TRUE ){

  cat( "\n" )
  cat( "Getting chromosome sequences...", "\n\n" )

  chr.seqs <- get_chromosome_seqs( genome.obj = genome.obj )

  if ( is.null( re.sites ) ){
    random <- TRUE
    re.cuts <- NULL
  } else{
    cat( "Finding restriction enzyme target sequences...", "\n\n" )
    re.cuts <- digest( string.list = chr.seqs,
                       re.sites = re.sites,
                       cut.after = cut.after )
    random <- FALSE
  }

  if ( isTRUE( random ) ){
    if ( missing( mean ) | missing( sd ) )
      stop( "The mean and sd of the sampling distribution must be defined to generate fragments via random fragmentation.",
            call. = FALSE )
  }

  frags <- mclapply( X = 1:iterations,
                     FUN = function(x){

                       cat( "Set", x, "\n" )
                       cat( "-----", "\n\n" )

                       cat( "Generating random site positions...", "\n\n" )

                       rand.sites <- random_sites( n.sites = n.sites,
                                                   genome.obj = genome.obj )

                       cat( "Getting fragment ranges...", "\n\n" )

                       rand.fragments <- make_fragments( insert.sites = rand.sites,
                                                         frag.sites = re.cuts,
                                                         random = random,
                                                         genome.obj = genome.obj,
                                                         mean = mean,
                                                         sd = sd )

                       cat( "Extracting fragment sequences...", "\n\n" )

                       frag.seqs <- getSeq( x = genome.obj,
                                            names = rand.fragments,
                                            as.character = FALSE )

                       cat( "Trimming fragment sequences...", "\n\n" )

                       frag.trim <- trim_seqs( fragments = frag.seqs,
                                               min.width = min.width,
                                               max.distance = max.distance,
                                               max.bp = max.bp )

                       cat( "Saving FASTA files...", "\n\n" )

                       if ( is.null( directory.path ) ){
                         stop( "directory.path and prefix must be defined to save FASTA files.",
                               call. = FALSE)
                       }

                       if( is.null( prefix ) ){
                         prefix <- paste0( "set_", x )
                       } else{
                         prefix <- paste0( prefix, "_", x )
                       }

                       save_fasta( reads = frag.trim,
                                   directory.path = directory.path,
                                   prefix = prefix,
                                   compress = compress )

                       return( rand.fragments )

                       },

                     mc.cores = n.cores

                     )

  if( isTRUE( collapse ) ){

    generated.fragments <- do.call( c, frags )

  } else{

    generated.fragments <- frags

  }

  if( isTRUE( write.ranges ) ){

    cat( "Writing generated fragment ranges...", "\n\n" )

    save( generated.fragments,
          file = paste0( directory.path, "/", "generated_fragments.RData.gz" ),
          compression_level = 6 )

    cat( "Done!", "\n\n" )

    } else {

      cat( "Done!", "\n\n" )

      return( generated.fragments )

    }

}
get_chromosome_seqs <- function( genome.obj ){
  seqs <- lapply( X = seqnames( genome.obj ),
                  FUN = function(x){
                    string <- DNAString( genome.obj[[ x ]] )
                    metadata( string )$name <- x
                    return( string )
                  } )
  names( seqs ) <- seqnames( genome.obj )
  return( seqs )
}

t2 <- simulate_random_data(BSgenome.Hsapiens.UCSC.hg19, n.sites = 9000000, mean = 234, sd = 50000,
                     directory.path = "C:/Users/38598/Desktop/hivint/IS_mapping/simulated",
                     prefix = "test2")


saveRDS(t2, file="9mil_random_sites_1.RDS")


# comapring various integration distances etc.

setwd("C:/Users/38598/Desktop/hivint/IS_mapping/simulated")
#simulated1 <- readRDS("distances_simulated_Rloopdists.RDS") 
simulated_full <- fread("HIV_rloops_1_simulated")
# make a full table out of this?
# making another 10 samples of full data but this time sample size is 2
uzorak <- sample(1:nrow(simulated_full), 274)

simulated_x <- simulated_full[uzorak,]
simulated_x$sample <- "sample 10"
simulated1 <- rbind(simulated1, simulated_x)

saveRDS(simulated1, file="10randoms_samples_n274.RDS")

setwd("C:/Users/38598/Desktop/hivint/IS_mapping")
results_real <- readRDS("distances_intToRloop.RDS")
genidf <- readRDS("C:/Users/38598/Desktop/hivint/MajasCode/Rloops/R-loops-code/peaks-analysis/genidf.RDS")

# all the genes from simulated1 are also in simulated full here, so no errors yet appear
#results_real <-combined
# merge some of the info from genidf (rigs rags) into all three tables
results_real <- merge(x = results_real, y = genidf[,c("names","rloop","RIG","rig","rag")],
                      by.x = "ensembl_gene_id", by.y = "names", all.x = T)
simulated1 <- merge(x = simulated1, y = genidf[,c("names","rloop","RIG","rig","rag")],
                      by.x = "ensembl_gene_id", by.y = "names", all.x = T)
simulated_full <- merge(x = simulated_full, y = genidf[,c("names","rloop","RIG","rig","rag")],
                      by.x = "ensembl_gene_id", by.y = "names", all.x = T)



# draft some comparisons of these things; 

#start with rloop vs no rloop graph
library(ggplot2)
library(dplyr)
# make frequency tables for any conceivable selection of differences; then calculate
# we have a problem cuz you set results real to no rloop if not in gene!!! fixing this now
results_real[results_real$ensembl_gene_id == "Not in gene",]$rloop.x <- "Not in gene"

results_real[sample=="AQR_inf_S7_",]$sample <- "AQR"
results_real[sample=="NTC_inf_S6_",]$sample <- "NTC"
results_real[sample=="NTC_mock_S5",]$sample <- "mock"
results_real[sample=="siAQR_M_S3_",]$sample <- "siAQR+wtAQR-R"
results_real[sample=="siAQR_MY_S4",]$sample <- "siAQR+Y1196A-R"
results_real[sample=="siAQR_S2_hg",]$sample <- "siAQR"
results_real[sample=="siCTRL_S1_h",]$sample <- "siCTRL"


simulated1[simulated1$ensembl_gene_id == "Not in gene",]$rloop.x <- "Not in gene"
simulated_full[simulated_full$ensembl_gene_id=="Not in gene",]$rloop.x <- "Not in gene"
# here again we decide that genes that have no info are genes with n RLoop
simulated_full[simulated_full$rloop.x=="",]$rloop.x <- "No info"
simulated1[is.na(simulated1$rloop.x),]$rloop.x <- "No info"
separate <- c("NTC","AQR")

results_real[is.na(rloop.x),]$rloop.x <- "No info"

rr <- results_real%>%
  group_by(sample) %>%
  count(rloop.x) %>%
  mutate(prop = prop.table(n))
#rr$rloop.x <- rr$rloop
#rr$rloop <- NULL
sf <- simulated_full %>%
  count(rloop.x) %>%
  mutate(prop = prop.table(n))

s1 <- simulated1 %>%
  count(rloop.x) %>%
  mutate(prop = prop.table(n))

sf$sample <- "All simulated integrations"
s1$sample <- "10 simulated samples n=80"

all_plot <- rbind(rr,sf,s1)
#reorder the naming factors
all_plot$sample <- factor(all_plot$sample, levels = c("10 simulated samples n=80","All simulated integrations",
                                                         "NTC","AQR", "siCTRL","siAQR","siAQR+wtAQR-R","siAQR+Y1196A-R"))
#organizing big and small letters
all_plot[all_plot$rloop.x=="rloop",]$rloop.x <- "Rloop"
all_plot[all_plot$rloop.x=="no rloop",]$rloop.x <- "No rloop"
all_plot$rloop.x <- factor(all_plot$rloop.x, levels = c( "Rloop","No rloop","Not in gene","No info"))

separate2 <- c(separate, "All simulated integrations","10 simulated samples n=80")

p1 <- ggplot(data=all_plot, aes(y=prop ,x=rloop.x)) +
  geom_bar(stat="identity", aes(fill = rloop.x))+facet_wrap(~pool) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="Rloop presence")+ scale_fill_brewer(palette="Blues", direction = -1)
p1
p15 <- ggplot(data=all_plot[!all_plot$sample %in% separate,], aes(y=prop ,x=rloop.x)) +
  geom_bar(stat="identity", aes(fill = rloop.x))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="Rloop presence")+ scale_fill_brewer(palette="Blues", direction = -1)
p15

p16 <- ggplot(data=all_plot[all_plot$rloop.x!="No info" &!all_plot$sample %in% separate,], aes(y=prop ,x=rloop.x)) +
  geom_bar(stat="identity", aes(fill = rloop.x))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="Rloop presence")+ scale_fill_brewer(palette="Blues", direction = -1)
p16

p17 <- ggplot(data=all_plot[all_plot$rloop.x!="No info" &all_plot$sample %in% separate2,], aes(y=prop ,x=rloop.x)) +
  geom_bar(stat="identity", aes(fill = rloop.x))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="Rloop presence")+ scale_fill_brewer(palette="Blues", direction = -1)
p17

p18 <- ggplot(data=all_plot[all_plot$sample %in% separate2 &all_plot$sample!="10 simulated samples n=80",], aes(y=prop ,x=rloop.x)) +
  geom_bar(stat="identity", aes(fill = rloop.x))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="Rloop presence")+ scale_fill_brewer(palette="Blues", direction = -1)
p18
p19 <- ggplot(data=all_plot[!all_plot$sample %in% separate&all_plot$sample!="10 simulated samples n=80",,], aes(y=prop ,x=rloop.x)) +
  geom_bar(stat="identity", aes(fill = rloop.x))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="Rloop presence")+ scale_fill_brewer(palette="Blues", direction = -1)
p19

p20 <- ggplot(data=all_plot[all_plot$rloop.x!="No info" &!all_plot$sample %in% separate&all_plot$sample!="10 simulated samples n=80",,], aes(y=prop ,x=rloop.x)) +
  geom_bar(stat="identity", aes(fill = rloop.x))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="Rloop presence")+ scale_fill_brewer(palette="Blues", direction = -1)
p20

p21 <- ggplot(data=all_plot[all_plot$rloop.x!="No info" &all_plot$sample %in% separate2&all_plot$sample!="10 simulated samples n=80",,], aes(y=prop ,x=rloop.x)) +
  geom_bar(stat="identity", aes(fill = rloop.x))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="Rloop presence")+ scale_fill_brewer(palette="Blues", direction = -1)
p21



# next leths do rig/RAG info 
results_real[results_real$ensembl_gene_id == "Not in gene",]$RIG <- "Not in gene"

simulated1[simulated1$ensembl_gene_id == "Not in gene",]$RIG <- "Not in gene"
simulated_full[simulated_full$ensembl_gene_id=="Not in gene",]$RIG <- "Not in gene"
# here again we decide that genes that have no info are genes with n RLoop
simulated_full[is.na(simulated_full$RIG),]$RIG <- "No info"
simulated1[is.na(simulated1$RIG),]$RIG <- "No info"
results_real[is.na(results_real$RIG),]$RIG <- "No info"



rr <- results_real %>%
  group_by(pool) %>%
  count(RIG) %>%
  mutate(prop = prop.table(n))

sf <- simulated_full %>%
  count(RIG) %>%
  mutate(prop = prop.table(n))

s1 <- simulated1 %>%
  count(RIG) %>%
  mutate(prop = prop.table(n))
sf$pool <- "All simulated integrations"
s1$pool <- "10 simulated samples n=80"
all_plot <- rbind(rr,sf,s1)

#all_plot$RIG <- as.character(all_plot$RIG)
#all_plot[is.na(all_plot$RIG),]$RIG <- "no info"
all_plot$pool <- factor(all_plot$pool, levels = c("10 simulated samples n=80","All simulated integrations",
                                                      "NTC","AQR"))
all_plot$RIG <- factor(all_plot$RIG, levels = c("RIG","1","not targeted","Not in gene", "No info"))





p2 <- ggplot(data=all_plot, aes(y=prop ,x=RIG)) +
  geom_bar(stat="identity", aes(fill = RIG))+facet_wrap(~pool) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="RIG info")+ scale_fill_brewer(palette="Blues", direction = -1)
p2

p25 <- ggplot(data=all_plot[all_plot$pool!="10 simulated samples n=80" & all_plot$RIG!="No info",], aes(y=prop ,x=RIG)) +
  geom_bar(stat="identity", aes(fill = RIG))+facet_wrap(~pool) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="RIG info")+ scale_fill_brewer(palette="Blues", direction = -1)
p25

p26 <- ggplot(data=all_plot[all_plot$sample %in% separate2& all_plot$RIG!="No info",], aes(y=prop ,x=RIG)) +
  geom_bar(stat="identity", aes(fill = RIG))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="RIG info")+ scale_fill_brewer(palette="Blues", direction = -1)
p26

p27 <- ggplot(data=all_plot[!all_plot$sample %in% separate & all_plot$RIG!="No info",], aes(y=prop ,x=RIG)) +
  geom_bar(stat="identity", aes(fill = RIG))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="RIG info")+ scale_fill_brewer(palette="Blues", direction = -1)
p27


p28 <- ggplot(data=all_plot[all_plot$sample %in% separate2&all_plot$sample!="10 simulated samples n=80",], aes(y=prop ,x=RIG)) +
  geom_bar(stat="identity", aes(fill = RIG))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="RIG info")+ scale_fill_brewer(palette="Blues", direction = -1)
p28

p29 <- ggplot(data=all_plot[!all_plot$sample %in% separate&all_plot$sample!="10 simulated samples n=80",], aes(y=prop ,x=RIG)) +
  geom_bar(stat="identity", aes(fill = RIG))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="RIG info")+ scale_fill_brewer(palette="Blues", direction = -1)
p29

p30 <- ggplot(data=all_plot[all_plot$sample %in% separate2& all_plot$RIG!="No info"&all_plot$sample!="10 simulated samples n=80",], aes(y=prop ,x=RIG)) +
  geom_bar(stat="identity", aes(fill = RIG))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="RIG info")+ scale_fill_brewer(palette="Blues", direction = -1)
p30

p31 <- ggplot(data=all_plot[!all_plot$sample %in% separate & all_plot$RIG!="No info"&all_plot$sample!="10 simulated samples n=80",], aes(y=prop ,x=RIG)) +
  geom_bar(stat="identity", aes(fill = RIG))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="RIG info")+ scale_fill_brewer(palette="Blues", direction = -1)
p31



# do the spads? meaning the pads followed by the exon info etc.

spads <- fread("spad_rigs_rloops.csv")



results_real <- merge(x = results_real, y = spads,
                      by.x = "ensembl_gene_id", by.y = "gene_id", all.x = T)
simulated1 <- merge(x = simulated1, y = spads,
                    by.x = "ensembl_gene_id", by.y = "gene_id", all.x = T)
simulated_full <- merge(x = simulated_full, y = spads,
                        by.x = "ensembl_gene_id", by.y = "gene_id", all.x = T)

results_real[results_real$ensembl_gene_id == "Not in gene",]$cat <- "Not in gene"

simulated1[simulated1$ensembl_gene_id == "Not in gene",]$cat <- "Not in gene"
simulated_full[simulated_full$ensembl_gene_id=="Not in gene",]$cat<- "Not in gene"
# here again we decide that genes that have no info are genes with n RLoop
simulated_full[is.na(simulated_full$cat),]$cat <- "No info"
simulated1[is.na(simulated1$cat),]$cat <- "No info"
results_real[is.na(results_real$cat),]$cat<- "No info"



# first the spads
rr <- results_real[] %>%
  group_by(pool) %>%
  count(cat) %>%
  mutate(prop = prop.table(n))

sf <- simulated_full %>%
  count(cat) %>%
  mutate(prop = prop.table(n))

s1 <- simulated1 %>%
  count(cat) %>%
  mutate(prop = prop.table(n))

sf$pool <- "All simulated integrations"
s1$pool <- "10 simulated samples n=80"

all_plot <- rbind(rr,sf,s1)
#all_plot[is.na(all_plot$cat),]$cat <- "no info"
all_plot$pool<- factor(all_plot$pool, levels =c("10 simulated samples n=80","All simulated integrations",
                                                     "NTC","AQR"))
all_plot$cat <- factor(all_plot$cat, levels = c("SPAD","Not SPAD","Not in gene", "No info"))


p3 <- ggplot(data=all_plot[], aes(y=prop ,x=cat)) +
  geom_bar(stat="identity", aes(fill = cat))+facet_wrap(~pool) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="SPADs")+ scale_fill_brewer(palette="Blues", direction = -1)
p3
p4 <- ggplot(data=all_plot[all_plot$pool!="10 simulated samples n=80"& all_plot$cat!="No info",], aes(y=prop ,x=cat)) +
  geom_bar(stat="identity", aes(fill = cat))+facet_wrap(~pool) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="SPADs")+ scale_fill_brewer(palette="Blues", direction = -1)
p4

p5 <- ggplot(data=all_plot[all_plot$sample %in% separate2& all_plot$cat!="No info",], aes(y=prop ,x=cat)) +
  geom_bar(stat="identity", aes(fill = cat))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="SPADs")+ scale_fill_brewer(palette="Blues", direction = -1)
p5
p6 <- ggplot(data=all_plot[!all_plot$sample %in% separate& all_plot$cat!="No info",], aes(y=prop ,x=cat)) +
  geom_bar(stat="identity", aes(fill = cat))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="SPADs")+ scale_fill_brewer(palette="Blues", direction = -1)
p6

p7 <- ggplot(data=all_plot[all_plot$sample %in% separate2&all_plot$sample!="10 simulated samples n=80",], aes(y=prop ,x=cat)) +
  geom_bar(stat="identity", aes(fill = cat))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="SPADs")+ scale_fill_brewer(palette="Blues", direction = -1)
p7
p8 <- ggplot(data=all_plot[!all_plot$sample %in% separate&all_plot$sample!="10 simulated samples n=80",], aes(y=prop ,x=cat)) +
  geom_bar(stat="identity", aes(fill = cat))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="SPADs")+ scale_fill_brewer(palette="Blues", direction = -1)
p8

p9 <- ggplot(data=all_plot[all_plot$sample %in% separate2& all_plot$cat!="No info"&all_plot$sample!="10 simulated samples n=80",], aes(y=prop ,x=cat)) +
  geom_bar(stat="identity", aes(fill = cat))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="SPADs")+ scale_fill_brewer(palette="Blues", direction = -1)
p9
p9 <- ggplot(data=all_plot[!all_plot$sample %in% separate& all_plot$cat!="No info"&all_plot$sample!="10 simulated samples n=80",], aes(y=prop ,x=cat)) +
  geom_bar(stat="identity", aes(fill = cat))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="SPADs")+ scale_fill_brewer(palette="Blues", direction = -1)
p9



# then the gene type, exon type and 


rr <- results_real[sample!="NTC mock S5"] %>%
  group_by(sample) %>%
  count(gene_type) %>%
  mutate(prop = prop.table(n))

sf <- simulated_full %>%
  count(gene_type) %>%
  mutate(prop = prop.table(n))

s1 <- simulated1 %>%
  count(gene_type) %>%
  mutate(prop = prop.table(n))

sf$sample <- "All simulated integrations"
s1$sample <- "10 simulated samples n=80"

all_plot <- rbind(rr,sf,s1)
all_plot[is.na(all_plot$gene_type),]$gene_type <- "no info"
#check whether those not in genes overlap totally with those not in ovrlap with spad
#  they dont, this isnt comprehensive data; some of these genes are not in the SPAD 

# here drop all rows that are smaller than 0.01
all_plot <- all_plot[all_plot$prop>0.01,]

p6 <- ggplot(data=all_plot[all_plot$sample %in% separate2,], aes(y=prop ,x=gene_type)) +
  geom_bar(stat="identity", aes(fill = gene_type))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="gene_type")+ggtitle("Gene types with > 1% frequency")
p6
p7 <- ggplot(data=all_plot[!all_plot$sample %in% separate,], aes(y=prop ,x=gene_type)) +
  geom_bar(stat="identity", aes(fill = gene_type))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="gene_type")+ggtitle("Gene types with > 1% frequency")
p7

p8 <- ggplot(data=all_plot, aes(y=prop ,x=gene_type)) +
  geom_bar(stat="identity", aes(fill = gene_type))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="gene_type")+ggtitle("Gene types with > 1% frequency")
p8


rr <- results_real[sample!="NTC mock S5"] %>%
  group_by(sample) %>%
  count(exon_type) %>%
  mutate(prop = prop.table(n))

sf <- simulated_full %>%
  count(exon_type) %>%
  mutate(prop = prop.table(n))

s1 <- simulated1 %>%
  count(exon_type) %>%
  mutate(prop = prop.table(n))

sf$sample <- "All simulated integrations"
s1$sample <- "10 simulated samples n=80"

all_plot <- rbind(rr,sf,s1)
all_plot[is.na(all_plot$exon_type),]$exon_type <- "no info"
#check whether those not in genes overlap totally with those not in ovrlap with spad
#  they dont, this isnt comprehensive data; some of these genes are not in the SPAD 


p9 <- ggplot(data=all_plot[all_plot$sample %in% separate2,], aes(y=prop ,x=exon_type)) +
  geom_bar(stat="identity", aes(fill = exon_type))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="exon_type")
p9
pa <- ggplot(data=all_plot[!all_plot$sample %in% separate,], aes(y=prop ,x=exon_type)) +
  geom_bar(stat="identity", aes(fill = exon_type))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="exon_type")
pa

pb <- ggplot(data=all_plot, aes(y=prop ,x=exon_type)) +
  geom_bar(stat="identity", aes(fill = exon_type))+facet_wrap(~sample) +theme_bw() +
  ggtitle("")+ ylab("Percentage of integrations")+
  xlab("")  + labs(fill="exon_type")
pb


pdf("Simulated_vs_real_plusSPADS_v2.pdf", height=6, width=8)
p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p15 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p16 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p17 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p18+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p19+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p20+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p21+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p25 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p26 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p27 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p28 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p29 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p30 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p31 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p4 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p5 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p6 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p7 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p8 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p9 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#p10 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


pdf("Simulated_vs_real_plusSPADS_new_data.pdf", height=6, width=8)
p1+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("Simulated_vs_real_plusSPADS_new_data.pdf", height=6, width=8)
p1+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# figuring out the statistics of all of this


# trying the seqpare method of similar

v1 <- results_real[results_real$sample == "siAQR_MY_S4"]
v2 <- results_real[results_real$sample == "siAQR_S2_hg"]
simulated1
setkey(v1, seqnames ,start,end)
setkey(v2, seqnames ,start,end)

v1 <- v1[v1$ensembl_gene_id !="Not in gene",]
v2 <- v2[v2$ensembl_gene_id !="Not in gene",]
# so this has no intersect with the insertions themselves but has some regarding the genes 

intersekcija <- intersect(IRanges(start = v1$start_position, end = v1$end_position, names = v1$seqnames ),
  IRanges(start = v2$start_position, end = v2$end_position, names = v2$seqnames ))


intersekcija <- intersect(IRanges(start = simulated1$start, end = simulated1$end, names = simulated1$seqnames ),
                          IRanges(start = v2$start, end = v2$end, names = v2$seqnames ))


