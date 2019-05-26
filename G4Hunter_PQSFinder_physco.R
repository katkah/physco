###############################################################################
#upload physco sequence
###############################################################################
library("Biostrings")
my_path <- "physcomitrella_patens"
setwd(my_path)
s <- readDNAStringSet("NC_037272.1[124466..135981]18S_5.8S_25S_5S_IGS.fa")
my_seq <- s[[1]]
s_string <- toString(my_seq)
putative_TSS <- regexpr("TATGTGGGGG",s_string)

###############################################################################
#pqsfinder
###############################################################################
source("https://bioconductor.org/biocLite.R")
biocLite("pqsfinder")
biocLite("rtracklayer")
biocLite("Gviz")
library(pqsfinder)
library(rtracklayer)
library(Gviz)

pqs <- pqsfinder(my_seq) 

##export the result to .gff and .fa
setwd(paste(my_path,"/G4_results", sep="")) 
gr <- as(pqs, "GRanges")
export(gr, "pqs_result.gff", version = "3")
dss <- as(pqs, "DNAStringSet")
writeXStringSet(dss, file = "pqs_result.fa", format = "fasta")

###############################################################################
#G4hunter
###############################################################################
source(paste(my_path,"/function_G4Hunter.r", sep = "")) # Functions created by Bedrat et al.
Refinment <- function(res, y) {
  ################################################################################
  #### function to refine the ends of extracted regions
  #### y is a long string from which we want the G4 regions to be extracted
  #### res is the not-yet-refined result of extracted regions from function G4extract in form "start end sequence G4Hscore"

  if (res[1] == -1){
    return(res)
  }
  RefinedStarts<- c()
  RefinedEnds<- c()
  RefinedSequences <- c()
  RefinedScores <- c()
  for (s in 1:nrow(res)) {
    #initializing step for each unrefined region
    toBeRefined <- as.character(res[s,3])#has to be converted to character due to the header of res
    letter <- "G" 
    #plus strand
    if (as.numeric(res[s,4]) < 0 ){ 
      #minus strand
      letter <- "C"
    }
    pos <- as.vector(gregexpr(letter,toBeRefined)[[1]])
    #positions of all G/C in toBeRefined region
    fstG <- pos[1] + as.numeric(res[s,1])-1
    #position of first G/C plus the absolute position of toBeRefined region in the original string y 
    lastG <- pos[length(pos)] + as.numeric(res[s,1])-1
    #position of last G/C plus the absolute position of toBeRefined region in the original string y 
    letter1 <- letter
    if (substr(toBeRefined,nchar(toBeRefined), nchar(toBeRefined) ) == letter )  {
      #G/C at the end => checks if the toBeRefined region can be elongated
      while (letter1 == letter && lastG < nchar(y)){ 
        #enlongates the END of extracted regions if there are neighboring G/C
        #does not allow for overflow of the original sequence y
        letter1 <- substr(y, lastG+1, lastG+1)
        if(letter1 == letter) {
          lastG <- lastG+1            
        }
      }
    }
    letter1 <- letter    
    if (substr(toBeRefined,1,1) == letter){
      #G/C at the start => checks if the toBeRefined region can be elongated 
      while (letter1 == letter && fstG > 1){ 
        #enlongates the start of extracted regions if there are neighboring G/C
        #does not allow for overflow of the original sequence y
        letter1 <- substr(y, fstG-1, fstG-1)
        if(letter1 == letter) {
          fstG <- fstG-1            
        }
      }      
    }  
    #saving the refined item
    RefinedStarts <- c(RefinedStarts,fstG)
    RefinedEnds <- c(RefinedEnds,lastG)
    RefinedSequences <- c(RefinedSequences,substr(y,fstG,lastG) )
    RefinedScores <- c(RefinedScores,signif(G4Hscore(substr(y,fstG,lastG)),2))
  }
  #returns table of refined regions
  newRes<-cbind(RefinedStarts,RefinedEnds,RefinedSequences,RefinedScores)
  newRes
}
OverlapMerge <- function(Borders){
  ##################################################
  ####merging of overlapping k-long regions
  checkingBorder <- Borders[1,]
  newBorders <- c()
  if (nrow(Borders) != 1){#in case nrow(Borders)==1 there are no regions to merge together
    for ( i in 2:nrow(Borders)) { 
      if (checkingBorder[2] >= Borders[i,1]){#merging regions
        checkingBorder[2] <- Borders[i,2] 
      } else {
        newBorders <- c(newBorders,checkingBorder)#saving fused regions which can no longer be elongated
        checkingBorder <- c(Borders[i,1],Borders[i,2])
      }
    }
  }
  newBorders = matrix(newBorders, ncol = 2, byrow = TRUE)
  newBorders
}
G4extract <- function(header, y, k = 25, hl = 1)
{
  ##################################################
  #identifies k-long regions in y which has higher G4Hscore than threshold hl
  #@param y input sequence
  #@param k length of region
  #@param hl threshold

  ##################################################
  #identifies regions with the score higher than hl
  #separate computation for minus strand and plus strand G4 (C-rich and G-rich quadruplexes)
  runmeanScore <- as.vector(G4runmean(y, k))

  myBoolRle_minus <- Rle(runmeanScore <= -hl) 
  startTrue_minus <- start(myBoolRle_minus)[ runValue(myBoolRle_minus) == 'TRUE']
  endTrue_minus <- end(myBoolRle_minus)[ runValue(myBoolRle_minus) == 'TRUE'] + k -1
  Borders_minus <- cbind(startTrue_minus, endTrue_minus)
  
  myBoolRle_plus <- Rle(runmeanScore >= hl) 
  startTrue_plus <- start(myBoolRle_plus)[ runValue(myBoolRle_plus) == 'TRUE']
  endTrue_plus <- end(myBoolRle_plus)[ runValue(myBoolRle_plus) == 'TRUE'] + k -1
  Borders_plus <- cbind(startTrue_plus, endTrue_plus)
  
  #check if there are any G4s otherwise it falls down
  #myBoolRle <- cbind(myBoolRle_minus, myBoolRle_plus)
  #if (nrun(myBoolRle) == 1 & runValue(myBoolRle)[1] == FALSE){
  #  res <- c(-1,-1,"no G4 found",-1)
  #  return(res)
  #}
  
  ##################################################
  ####merging of overlapping k-long regions
  newBorders <- rbind(OverlapMerge(Borders_minus),OverlapMerge(Borders_plus))
    
  if (length(newBorders) != 0 ){
   # 
    extract = c()
    for (i in 1:nrow(newBorders) ){
      extract <- c(extract, c(substr(y,newBorders[i,1],newBorders[i,2]) ) )   
    }
    G4score = c()
    for ( sequence in extract) {
      G4score <- c(G4score, signif(G4Hscore(sequence),2))  
    }
    res <- cbind(newBorders, extract, G4score) 
    colnames(res) <- c("start", "end","sequence", "G4score")
  }else {
    res <- c(-1,-1,"no G4 found",-1) 
    return(res)
  }
  
  ##################################################
  ####extremity repair and final G4Hscore computation   
  ref <- Refinment(res,y)  
  ref_coor <- cbind(as.numeric(ref[,1]), as.numeric(ref[,2]))
  ref

}

res <- G4extract("rDNA", s_string)
res <- as.data.frame(res)
hunter_scores = as.numeric(levels(res$RefinedScores)[as.integer(res$RefinedScores)])
score_to_strand <- function(hunter_score) {
  if (sign(hunter_score) < 0){
    return("-")
  }else{
    return("+")
  }
}

#####################################################
#export result
#prepare hunter_gr suitable for exporting the result
hunter_df1 <- data.frame(
  chrom = "chr1", 
  start = as.numeric(levels(res$RefinedStarts)[as.integer(res$RefinedStarts)]),
  end = as.numeric(levels(res$RefinedEnds)[as.integer(res$RefinedEnds)]),
  strand = sapply(hunter_scores, score_to_strand),
  scores = hunter_scores,
  sequence = res$RefinedSequences
)

hunter_gr1 <- as(hunter_df1, "GRanges")
hunter_gr1 <- sortSeqlevels(hunter_gr1)
hunter_gr1 <- sort(hunter_gr1)
export(hunter_gr1, "G4hunter_result.gff", version = "3")


#prepare hunter_gr in a form suitable for Gviz DataTracks
hunter_df <- data.frame(
  chrom = "chr1", 
  start = as.numeric(levels(res$RefinedStarts)[as.integer(res$RefinedStarts)]),
  end = as.numeric(levels(res$RefinedEnds)[as.integer(res$RefinedEnds)]),
  strand = sapply(hunter_scores, score_to_strand),
  plus_strand = hunter_scores,
  minus_strand = hunter_scores
)
hunter_df[which(hunter_df$strand == "+"),"minus_strand"]  <- NA
hunter_df$minus_strand <- abs(hunter_df$minus_strand)
hunter_df[which(hunter_df$strand == "-"),"plus_strand"]  <- NA
hunter_df[,"strand"] <- "*"

hunter_gr <- as(hunter_df, "GRanges")

###############################################################################
#rDNA sequence annotation
###############################################################################
ir <- IRanges(start = c(416,2176,2655,6779,8747,11427), end = c(1864,2328,6100,6897,8757,11516), names = c("18S","5.8S","28S", "5S", "putative TSS", "18S") )
gr_annotation <- GRanges(seqnames = "chr1", strand = c("+", "+", "+", "+","+","+"),
              ranges = ir)

###############################################################################
#plot the result
###############################################################################
gtrack <- AnnotationTrack(gr_annotation, 
                          name = "rDNA", 
                          #id = names(gr_annotation), 
                          cex = 0.7, 
                          fontcolor.item="black",
                          #showFeatureId=TRUE,
                          shape = "box",
                          group = c("18S","5.8S","28S","5S", "TSS", "18S next unit"))
atrack <- GenomeAxisTrack()

pqs_df <- data.frame(chrom = "chr1", start = start(pqs), end = end(pqs), strand = strand(pqs), plus_strand = score(pqs), minus_strand = score(pqs))
pqs_df[which(pqs_df$strand == "+"),"minus_strand"]  <- NA
pqs_df[which(pqs_df$strand == "-"),"plus_strand"]  <- NA
pqs_df[,"strand"] <- "*"
gr <- as(pqs_df, "GRanges")
pqs_track <- DataTrack(
  gr,
  name = "pqsfinder G4 score",
  col=c("red","cornflowerblue"),
  groups=c("plus_strand","minus_strand")
)

G4hunter_track <- DataTrack(
  hunter_gr,
  name = "G4hunter score",
  col = c("red","cornflowerblue"),
  groups=c("plus_strand","minus_strand")
)

tiff("Plot3.tiff", width = 6, height = 4, units = 'in', res = 300)
suppressWarnings(plotTracks(c(gtrack, pqs_track, G4hunter_track, atrack), type = "h", labelPos="below", groupAnnotation = "group", just.group = "below"))
dev.off()

suppressWarnings(plotTracks(c(gtrack, pqs_track, G4hunter_track, genomeAxis), type = "h"))




###############################################################################
#plant/human telomere sequence
###############################################################################
plant <- paste(replicate(30,"TTTAGGG"), collapse = "")
plant <- DNAString(plant)
pqsfinder(plant)
G4extract("rDNA", plant)


human <- paste(replicate(30,"TTAGGG"), collapse = "")
human <- DNAString(human)
pqsfinder(human)
G4extract("rDNA", human)

###############################################################################
#zoom 5S-18S
###############################################################################
my_path <- c("physcomitrella_patens")
setwd(my_path)
s = readDNAStringSet("physco_reference_IGS.fasta")
my_seq <- s[[1]]
setwd(paste(my_path,"/G4_results/zoom_5S-18S", sep=""))
pqs <- pqsfinder(my_seq)
gr <- as(pqs, "GRanges")
export(gr, "5S-18S_pqs_result.gff", version = "3")
dss <- as(pqs, "DNAStringSet")
writeXStringSet(dss, file = "5S-18S_pqs_result.fa", format = "fasta")
ir <- IRanges(start = c(1,4649), end = c(119,4738), names = c("5S", "18S") )
gr_annotation <- GRanges(seqnames = "chr1", strand = c("+", "+"),
                         ranges = ir)
gtrack <- AnnotationTrack(gr_annotation, name = "rDNA", id = names(gr_annotation), cex = 0.7, 
                          fontcolor.item="black",
                          showFeatureId=TRUE )
atrack <- GenomeAxisTrack()
strack <- DataTrack(
  start = start(pqs), end = end(pqs),
  data = score(pqs), genome = "rDNA_unit", chr = "chr1", name = "pqsfinder G4 score")

suppressWarnings(plotTracks(c(gtrack, strack, atrack), type = "h"))




