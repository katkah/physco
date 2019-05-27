
################################################################################
################# Function used for the G4Hunter paper #########################
################################################################################
#### L. Lacroix, laurent.lacroix@inserm.fr , 20150928

################################################################################
###### G4translate change the DNA code into G4Hunter code.
###### Only G or C are taken into account. non G/C bases are translated in 0
###### It is OK if N or U are in the sequence
###### but might be a problem if other letters or numbers are present
###### lowercase ARE not welcome 
G4translate <- function(x)		
# x is a Rle of a sequence
	{
	xres=x
	runValue(xres)[runValue(x)=='C' & runLength(x)>3] <- -4
	runValue(xres)[runValue(x)=='C' & runLength(x)==3] <- -3
	runValue(xres)[runValue(x)=='C' & runLength(x)==2] <- -2
	runValue(xres)[runValue(x)=='C' & runLength(x)==1] <- -1
	runValue(xres)[runValue(x)=='G' & runLength(x)>3] <- 4
	runValue(xres)[runValue(x)=='G' & runLength(x)==3] <- 3
	runValue(xres)[runValue(x)=='G' & runLength(x)==2] <- 2
	runValue(xres)[runValue(x)=='G' & runLength(x)==1] <- 1
	runValue(xres)[runValue(x)!='C' & runValue(x)!='G'] <- 0
	Rle(as.numeric(xres))
	}
################################################################################


################################################################################
###### G4Hscore return the G4Hscore of a sequence
###### The output is a number
G4Hscore <- function(y)		
# y can be DNAString or a DNAStringSet or a simple char string.
	{
	require(S4Vectors)
	y2 <- Rle(strsplit(as.character(y),NULL)[[1]])
	y3 <- G4translate(y2)
	mean(y3)
	}
################################################################################


################################################################################
###### G4runmean return the G4Hscore in a window of 25 using runmean.
###### The output is a Rle of the runmeans
G4runmean <- function(y,k=25)		
#y is a DNAString or a DNAStringSet, k is the window for the runmean
	{
	require(S4Vectors)
	y2 <- Rle(strsplit(as.character(y),NULL)[[1]])
	y3 <- G4translate(y2)
	runmean(y3,k)
	}
################################################################################





################################################################################
#### function to extract sequences with abs(G4Hscore)>=hl (threshold)
#### from a chromosome (i) of a genome (gen) 
#### return a GRanges
#### use masked=5 for unmasked genome
#### need a genome to be in the memory in the genome variable
#### k is the window size for the runmean
#### i is the chromosome number
#### hl is the threshold

G4hunt <- function(i,k=25,hl=1.5,gen=genome,masked=5)
	{
	require(GenomicRanges)
	chr <- gen[[i]]
	if (masked==2) {chr=injectHardMask(chr)}
	if (masked==3) {active(masks(chr))['RM']=T;chr=injectHardMask(chr)}
	if (masked==0) {active(masks(chr))=F;chr=injectHardMask(chr)}
	if (masked==4) {active(masks(chr))=T;chr=injectHardMask(chr)}
	
	chr_G4hk <- G4runmean(chr,k)
	tgenome <- G4translate(Rle(strsplit(as.character(chr),NULL)[[1]]))
	if (class(gen)=="DNAStringSet")
		{seqname <- names(gen)[i]
		}else{
		seqname <- seqnames(gen)[i]
		}	
	
		j <- hl
		chrCh <- Views(chr_G4hk, chr_G4hk<=(-j))
		chrGh <- Views(chr_G4hk, chr_G4hk>=j)
	
		IRC <- IRanges(start=start(chrCh),end=(end(chrCh)+k-1))
		IRG <- IRanges(start=start(chrGh),end=(end(chrGh)+k-1))
		nIRC <- reduce(IRC)
		nIRG <- reduce(IRG)
		# overlapping results on the same strand are fused
		nchrCh <- Views(chr_G4hk,start=start(nIRC),end=(end(nIRC)-k+1))
		nchrGh <- Views(chr_G4hk,start=start(nIRG),end=(end(nIRG)-k+1))
		nscoreC <- signif(mean(Views(tgenome,nIRC)),3)
		nscoreG <- signif(mean(Views(tgenome,nIRG)),3)
		nstraC <- Rle(rep('-',length(nIRC)))
		nstraG <- Rle(rep('+',length(nIRG)))
		nhlC <- Rle(rep(j,length(nIRC)))
		nhlG <- Rle(rep(j,length(nIRG)))
		nkC <- Rle(rep(k,length(nIRC)))
		nkG <- Rle(rep(k,length(nIRG)))
		maskC <- Rle(rep(masked,length(nIRC)))
		maskG <- Rle(rep(masked,length(nIRG)))
	
		if (length(nIRC)==0)
			{
			nxC <- GRanges()                                
			}else{      
			nxC <- GRanges(seqnames=Rle(seqname),
		        ranges=nIRC,
		        strand=nstraC,
		        score=nscoreC,hl=nhlC,k=nkC,mask=maskC)
			}
	
		if (length(nIRG)==0)
			{        
			nxG <- GRanges()
			}else{          
			nxG <- GRanges(seqnames=Rle(seqname),
		        ranges=nIRG,
		        strand=nstraG,
		        score=nscoreG,hl=nhlG,k=nkG,mask=maskG)
			}
		nx <- c(nxC,nxG)
		return(nx)
	}
################################################################################


################################################################################
##### like G4hunt but for several thresholds
##### warning, if the list of threshold is created with the seq function
##### rounding error can occur
##### return a GRangesList
G4huntlist <- function(i,k=25,hl=c(1,1.2,1.5,1.75,2),gen=genome,masked=5)
	{
	require(GenomicRanges)
	nx <- list()
	chr <- gen[[i]]
	if (masked==2) {chr=injectHardMask(chr)}
	if (masked==3) {active(masks(chr))['RM']=T;chr=injectHardMask(chr)}
	if (masked==0) {active(masks(chr))=F;chr=injectHardMask(chr)}
	if (masked==4) {active(masks(chr))=T;chr=injectHardMask(chr)}

	chr_G4hk <- G4runmean(chr,k)
	tgenome <- G4translate(Rle(strsplit(as.character(chr),NULL)[[1]]))
	if (class(gen)=="DNAStringSet")
		{seqname <- names(gen)[i]
		}else{
		seqname <- seqnames(gen)[i]
		}

	hl <- round(hl,2)		

	for (j in hl)
		{
		chrCh <- Views(chr_G4hk, chr_G4hk<=(-j))
		chrGh <- Views(chr_G4hk, chr_G4hk>=j)
	
		IRC <- IRanges(start=start(chrCh),end=(end(chrCh)+k-1))
		IRG <- IRanges(start=start(chrGh),end=(end(chrGh)+k-1))
		nIRC <- reduce(IRC)
		nIRG <- reduce(IRG)
		# overlapping results on the same strand are fused
		nchrCh <- Views(chr_G4hk,start=start(nIRC),end=(end(nIRC)-k+1))
		nchrGh <- Views(chr_G4hk,start=start(nIRG),end=(end(nIRG)-k+1))
		nscoreC <- signif(mean(Views(tgenome,nIRC)),3)
		nscoreG <- signif(mean(Views(tgenome,nIRG)),3)
		nstraC <- Rle(rep('-',length(nIRC)))
		nstraG <- Rle(rep('+',length(nIRG)))
		nhlC <- Rle(rep(j,length(nIRC)))
		nhlG <- Rle(rep(j,length(nIRG)))
		nkC <- Rle(rep(k,length(nIRC)))
		nkG <- Rle(rep(k,length(nIRG)))
		maskC <- Rle(rep(masked,length(nIRC)))
		maskG <- Rle(rep(masked,length(nIRG)))
	
		if (length(nIRC)==0)
			{
			nxC <- GRanges()                                
			}else{      
			nxC <- GRanges(seqnames=Rle(seqname),
		        ranges=nIRC,
		        strand=nstraC,
		        score=nscoreC,hl=nhlC,k=nkC,mask=maskC)
			}
	
		if (length(nIRG)==0)
			{        
			nxG <- GRanges()
			}else{          
			nxG <- GRanges(seqnames=Rle(seqname),
		        ranges=nIRG,
		        strand=nstraG,
		        score=nscoreG,hl=nhlG,k=nkG,mask=maskG)
			}
		nx[[which(hl==j)]]=c(nxC,nxG)
		}
		names(nx) <- as.character(hl)
		nx <- GRangesList(nx)
		return(nx)
	}
################################################################################


################################################################################

##### function to refine G4hunt results

G4startrun <- function(y,chrom=chr,letter='C')	#y is a START
	{
		if (y!=1)
			{if (letter(chrom,y)==letter)
				{
				while (letter(chrom,y-1)==letter) {y <- y-1}
				}
				else
				{
				y <- y+1
				while (letter(chrom,y)!=letter) {y <- y+1}
				}
			}	
		return(y)
	}

G4endrun <- function(y,chrom=chr,letter='C')	#y is a END
	{
		if (y!=length(chrom))
			{if (letter(chrom,y)==letter)
				{
				while (letter(chrom,y+1)==letter) {y <- y+1}
				}
				else
				{
				y <- y-1
				while (letter(chrom,y)!=letter) {y <- y-1}
				}
			}
		return(y)
	}	
################################################################################



################################################################################



G4huntrefined <- function(Res_hl,gen=genome,i=25)		
# take in input the result from the function G4hunt, a GRanges
# it also require a genome and chromosome number (i) to extract the seq 
# and  to build the GRanges
{
chr <- gen[[i]]
if (class(gen)=="DNAStringSet")
	{seqname <- names(gen)[i]
	}else{
	seqname <- seqnames(gen)[i]
	}
Reschr <- Res_hl[seqnames(Res_hl)==seqname]		
tgenome <- G4translate(Rle(strsplit(as.character(chr),NULL)[[1]]))

### C treatment 
ReschrC <- Reschr[which(score(Reschr)<0)]
IRC <- IRanges(start=start(ReschrC),end=end(ReschrC))
if (length(IRC)==0)
{
totC <- NULL                                
}else{
chr_transC <- Views(tgenome,IRC)
	
nnIRC <- IRC
start(nnIRC) <- sapply(start(IRC),G4startrun,letter='C',chrom=chr)
end(nnIRC) <- sapply(end(IRC),G4endrun,letter='C',chrom=chr)
nG4scoreC <- signif(mean(Views(tgenome,nnIRC)),3)
nnseqC <- as.character(Views(chr,nnIRC))

totC <- cbind.data.frame(start(nnIRC),end(nnIRC),nnseqC,nG4scoreC)
names(totC) <- c('newstart','newend','newseq','newG4h')
}

### G treatment
ReschrG <- Reschr[which(score(Reschr)>0)]
IRG=IRanges(start=start(ReschrG),end=end(ReschrG))
if (length(IRG)==0)
{
  totG <- NULL                            
}else{
chr_transG <- Views(tgenome,IRG)
nnIRG <- IRG
start(nnIRG) <- sapply(start(IRG),G4startrun,letter='G',chrom=chr)
end(nnIRG) <- sapply(end(IRG),G4endrun,letter='G',chrom=chr)
nG4scoreG <- signif(mean(Views(tgenome,nnIRG)),3)
nnseqG <- as.character(Views(chr,nnIRG))

totG <- cbind.data.frame(start(nnIRG),end(nnIRG),nnseqG,nG4scoreG)
names(totG) <- c('newstart','newend','newseq','newG4h')
}
if (is.null(totC) & is.null(totG))
	{
	tot2 <- NULL
	tot2GR <- NULL
	}else{
	tot <- rbind.data.frame(totC,totG)
	tot2 <- tot[order(tot[,'newstart']),]
	stra <- sign(tot2$newG4h)
	stra[stra==1] <- '+'
	stra[stra=='-1'] <- '-'
	tot2GR <- GRanges(seqnames=Rle(seqname,length(tot2$newstart)), ranges=IRanges(start=tot2$newstart,end=tot2$newend),strand=Rle(stra),score=tot2$newG4h,sequence=tot2$newseq,seqinfo=seqinfo(gen))
	}
return(tot2GR)
}
################################################################################


################################################################################
#### Home made function to test if a sequence matchs with a Quadparser type G4
#### this function return 1 if there is a seqence matching the pattern 
#### "G{gmin,gmax}.{lmin,lmax}G{gmin,gmax}.{lmin,lmax}Ggmin,gmax}.{lmin,lmax}G{gmin,gmax}"
#### by default, the pattern is "G{3,}.{1,7}G{3,}.{1,7}G{3,}.{1,7}G{3,}"
#### the function return -1 for a match on the complementary strand
QPtest <- function(sequence,gmin=3,gmax='',lmin=1,lmax=7)

	{ 
		Gpat <- paste('G{',gmin,',',gmax,'}.{',lmin,',',lmax,'}G{',gmin,',',gmax,'}.{',lmin,',',lmax,'}G{',gmin,',',gmax,'}.{',lmin,',',lmax,'}G{',gmin,',',gmax,'}',sep='')
		Cpat <- paste('C{',gmin,',',gmax,'}.{',lmin,',',lmax,'}C{',gmin,',',gmax,'}.{',lmin,',',lmax,'}C{',gmin,',',gmax,'}.{',lmin,',',lmax,'}C{',gmin,',',gmax,'}',sep='')
		Gres <- gregexpr(Gpat, sequence, perl=T)
		Cres <- gregexpr(Cpat, sequence, perl=T)
		output <- 0
		if (Gres[[1]][1]>0)
			{output <- 1}
		if (Cres[[1]][1]>0)
			{output <- -1}
		return(output)	
	}
################################################################################


################################################################################
#### Homemade function for Quadparser
quadparser <- function(i,gen=genome,gmin=3,gmax='',lmin=1,lmax=7)
# this function extract form a chromosome 'i' in a genome from a BSGenome all the sequence matching the Gpattern
# "G{gmin,gmax}.{lmin,lmax}G{gmin,gmax}.{lmin,lmax}Ggmin,gmax}.{lmin,lmax}G{gmin,gmax}"
# or its complement pattern
# by default, the pattern is "G{3,}.{1,7}G{3,}.{1,7}G{3,}.{1,7}G{3,}"
# the Gpat2/Cpat2 pattern are here to compensate the issue with overlapping sequences excluded from gregexpr
	{
	require(GenomicRanges)	
	Gpat <- paste('G{',gmin,',',gmax,'}.{',lmin,',',lmax,'}G{',gmin,',',gmax,'}.{',lmin,',',lmax,'}G{',gmin,',',gmax,'}.{',lmin,',',lmax,'}G{',gmin,',',gmax,'}',sep='')
	Cpat <- paste('C{',gmin,',',gmax,'}.{',lmin,',',lmax,'}C{',gmin,',',gmax,'}.{',lmin,',',lmax,'}C{',gmin,',',gmax,'}.{',lmin,',',lmax,'}C{',gmin,',',gmax,'}',sep='')
	Gpat2 <- paste('G{',gmin,',',gmax,'}.{',lmin,',',lmax,'}G{',gmin,',',gmax,'}.{',lmin,',',lmax,'}G{',gmin,',',gmax,'}.{',lmin,',',lmax,'}G{',gmin,',',gmax,'}','(.{',lmin,',',lmax,'}G{',gmin,',',gmax,'})?',sep='')
	Cpat2 <- paste('C{',gmin,',',gmax,'}.{',lmin,',',lmax,'}C{',gmin,',',gmax,'}.{',lmin,',',lmax,'}C{',gmin,',',gmax,'}.{',lmin,',',lmax,'}C{',gmin,',',gmax,'}','(.{',lmin,',',lmax,'}C{',gmin,',',gmax,'})?',sep='')
		
	if (class(gen)=="DNAStringSet")
		{
		seqname <- names(gen)[i]
		}else{
		seqname <- seqnames(gen)[i]
		}	
	seqinf=seqinfo(gen)	
	Gres <- gregexpr(Gpat, gen[[i]], perl=T)
	Cres <- gregexpr(Cpat, gen[[i]], perl=T)
	Gres2 <- gregexpr(Gpat2, gen[[i]], perl=T)
	Cres2 <- gregexpr(Cpat2, gen[[i]], perl=T)
	if (Gres[[1]][1]>0) 
		{
		Gwidth <- attr(Gres[[1]],"match.length")
		Gstart <- unlist(Gres)
		G4Ranges <- GRanges(seqnames=seqname,IRanges(start=Gstart,width=Gwidth),strand='+',seqinfo=seqinf)
		}else{
		G4Ranges <- GRanges()
		}
	if (Cres[[1]][1]>0)		
		{
		Cwidth <- attr(Cres[[1]],"match.length")
		Cstart <- unlist(Cres)
		C4Ranges <- GRanges(seqnames=seqname,IRanges(start=Cstart,width=Cwidth),strand='-',seqinfo=seqinf)
		}else{
		C4Ranges <- GRanges()
		}	
	if (Gres2[[1]][1]>0) 
		{
		Gwidth2 <- attr(Gres2[[1]],"match.length")
		Gstart2 <- unlist(Gres2)
		G4Ranges2 <- GRanges(seqnames=seqname,IRanges(start=Gstart2,width=Gwidth2),strand='+',seqinfo=seqinf)
		}else{
		G4Ranges2 <- GRanges()
		}
	if (Cres2[[1]][1]>0)		
		{
		Cwidth2 <- attr(Cres2[[1]],"match.length")
		Cstart2 <- unlist(Cres2)
		C4Ranges2 <- GRanges(seqnames=seqname,IRanges(start=Cstart2,width=Cwidth2),strand='-',seqinfo=seqinf)
		}else{
		C4Ranges2 <- GRanges()
		}	
	res <- sort(reduce(c(G4Ranges,C4Ranges,G4Ranges2,C4Ranges2)),ignore.strand=T)
	res$seq <- as.character(Views(gen[[i]],ranges(res)))
	return(res)
	}	
################################################################################


################################################################################
#### a function to transform a GRanges of G4FS to G4FS per kb, used to generate bigwig
G4GR2bigwig <- function(gr,gen=genome,mcore=mcore)
	{
	chrlist <- 1:length(seqlevels(gr))
	tt=mclapply(chrlist,function(y) 
						{
						a=GRanges(seqnames=seqnames(gen)[y],ranges=breakInChunks(totalsize=seqlengths(gen)[y],chunksize=1000),strand='*',seqinfo=seqinfo(gr))
						a$score=countOverlaps(a,resize(gr[seqnames(gr)==seqnames(gen)[y]],fix='center',width=1),ignore.strand=T)
						return(a)
						},mc.cores=mcore)
	zc=do.call(c,tt)
	return(zc)
	}	
################################################################################
