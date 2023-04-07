#Copywrite Valentina Boeva, 2017

# >>> SOURCE LICENSE >>>
#   This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation (www.fsf.org); either version 2 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# A copy of the GNU General Public License is available at
# http://www.fsf.org/licensing/licenses
# 
# >>> END OF LICENSE >>>

#################################################################################
# to run:
# cat runLILY.R | R --slave --args SAMPLE_NAME OUTPUT_DIR Distance_toStitch distFromTSS transcriptome_GFF_file faiFileToCreateBigWig
#################################################################################

cl_args <- commandArgs()

if(
  sum(grep("RStudio",cl_args,ignore.case = TRUE))>0 || #rstudio
  sum(grep("RGui",cl_args,ignore.case = TRUE))>0 #rgui
){
  
  if (sum(c("sampleName",
            "OUTPUT_DIR",
            "maxDistanceToStitch",
            "distFromTSS",
            "transcriptomeFile",
            "faiFileToCreateBigWig") %in% ls())<6){
    stop("We have started from a GUI and we cannot use args.\n",
         "Please, define sampleName, OUTPUT_DIR, maxDistanceToStitch, distFromTSS,\n",
         "transcriptomeFile, and faiFileToCreateBigWig before we start.\n")
  }
} else { #Rscript cli run. Let's parse
  sampleName=cl_args[4]
  OUTPUT_DIR=cl_args[5]
  maxDistanceToStitch=as.numeric(cl_args[6])
  distFromTSS=as.numeric(cl_args[7])
  transcriptomeFile=cl_args[8]
  faiFileToCreateBigWig=cl_args[9]
}
cat (paste("..Working with sample",sampleName,"\n"))
cat (paste("..Will stitch enhancers at maximal distance of",maxDistanceToStitch,"\n"))
cat (paste("..Will consider promoter regions as regions around +-",distFromTSS,"from gene TSS\n"))


#################check that all files and directories exist ####################

dir.create(OUTPUT_DIR, showWarnings = FALSE)
setwd(OUTPUT_DIR)
cat (paste("..Will write the output into ",OUTPUT_DIR,"\n"))

if(!file.exists(transcriptomeFile)) {
  cat ("Error:Please prodive a valid path to a file with transcriptome information\nFor example from here: https://github.com/linlabbcm/rose2/tree/master/rose2/annotation\nEXIT!\n")
  quit()
}

if(!file.exists(paste0(sampleName,"_peaks.narrowPeak"))) {
  cat ("Error:Please prodive a valid path to output files of HMCan\n")
  quit()
}


#####################################################

suppressMessages(library(rtracklayer))

################# FUNCTIONS #########################


import.bed3<- function(filename){
  peaks=read.table(filename)
  return(GRanges(peaks[,1],IRanges(peaks[,2],peaks[,3])))
}

import.bed6<- function(filename){
  peaks=read.table(filename)
  return(GRanges(peaks[,1],IRanges(peaks[,2],peaks[,3]),score=peaks[,5]))
}

import.ucsc<- function(filename){
  peaks=read.table(filename,header = T,comment.char = "!")
  peaks=GRanges(peaks$chrom,IRanges(peaks$txStart,peaks$txEnd),strand = peaks$strand)
  return(peaks)
}


getThresholdOnPeakScore <- function (filename) {
  dataTable <-read.table(filename, header=F);  
  lengths = dataTable$V3-dataTable$V2  
  lengths = (lengths - 1)/50 + 1  
  scores=dataTable$V5
  
  tranch=200
  ascendingScores = rev(scores)
  lengthForAscendingScores = rev(lengths)
  
  x = NULL
  y = NULL
  z=NULL
  for (i in c(0:(floor(length(lengths)/tranch)-1))) {
    tt=c((tranch*i+1):(tranch*i+tranch))
    m=mean(ascendingScores[tt]/lengthForAscendingScores[tt])
    s=sd(ascendingScores[tt]/lengthForAscendingScores[tt])
    x=c(x,max(ascendingScores[tt]))
    y=c(y,m)
    z=c(z,s)
  }
  
  thresholdToReturn=x[min(which(y>=min(y[which(x>5)])))]
  if (thresholdToReturn<2) {thresholdToReturn=2} #warning!!!
  return(thresholdToReturn)
}

output_peaks_density <- function(enhancersStitched,enhancers,promoters,bwFile,outputFile){
  
  bedRegions = enhancersStitched
  strand(bedRegions)="*"; 
  bedRegions$score=0; 
  densities = import.bw(bwFile)  
  gc()
  overlaps=findOverlaps(bedRegions,densities)  
  ensids <- densities$score[subjectHits(overlaps)]
  x <- tapply(ensids, queryHits(overlaps), sum)
  bedRegions$score[unique(queryHits(overlaps))]=x
  rm(overlaps)
  rm(densities)
  gc()  
  strand(bedRegions)="+"; 
  bedRegions=bedRegions[order(bedRegions$score,decreasing=T)]
  cutoff=calculate_cutoff(bedRegions$score,F)
  bedRegions$name="enhancer"
  bedRegions$name[which(bedRegions$score>cutoff)]="SE"
  
  #unstitch enhancers:
  strand(enhancers)="+";enhancers$score=0;enhancers$name="enhancer"; 
  tt1=which(countOverlaps(enhancers,bedRegions[which(bedRegions$name=="enhancer")])>0)
  tt2=which(countOverlaps(bedRegions,enhancers)>0 & bedRegions$name=="enhancer")
  
  bedRegions=bedRegions[-tt2]
  bedRegions=c(bedRegions,enhancers[tt1])
  
  #find promoters intersecting with enhancers (not SEs) and call these enhancers: promoters
  promoters$score=0
  promoters$name="promoter"
  strand(promoters)="+"
  
  enh_promoters=c(promoters,bedRegions[which(bedRegions$name=="enhancer")])
  enh_promoters=resize(enh_promoters,20+width(enh_promoters))
  enh_promoters=reduce(enh_promoters)
  enh_promoters=resize(enh_promoters,-20+width(enh_promoters))
  
  tt1=which(countOverlaps(enh_promoters,promoters)>0)
  tt2=which(countOverlaps(enh_promoters,promoters)==0)
  largePromoters=enh_promoters[tt1]
  largeEnh=enh_promoters[tt2]
  
  largePromoters$score=0;largePromoters$name="promoter"
  largeEnh$score=0;largeEnh$name="enhancer"
  tt1=which(countOverlaps(largeEnh,bedRegions[which(bedRegions$name=="enhancer")])>0)
  tt2=which(countOverlaps(bedRegions,largeEnh)>0 & bedRegions$name=="enhancer")
  bedRegions=bedRegions[-tt2]
  bedRegions=c(bedRegions,largeEnh[tt1])
  
  finalSet=c(bedRegions,largePromoters)
  strand(finalSet)="*"; 
  finalSet$score=0; 
  densities = import.bw(bwFile)  
  gc()
  overlaps=findOverlaps(finalSet,densities)  
  ensids <- densities$score[subjectHits(overlaps)]
  x <- tapply(ensids, queryHits(overlaps), sum)
  finalSet$score[unique(queryHits(overlaps))]=x
  rm(overlaps)
  rm(densities)
  gc()  
  strand(finalSet)="+"; 
  
  export.bed(finalSet,outputFile)
}

#============================================================================
#==============SUPER-ENHANCER CALLING AND PLOTTING FUNCTIONS=================
#======   from the original ROSE package developed by Charles Lin  (c) ======
#======   ROSE licence:    http://younglab.wi.mit.edu/ROSE/LICENSE.txt ======
#============================================================================

#this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
numPts_below_line <- function(myVector,slope,x){
  yPt <- myVector[x]
  b <- yPt-(slope*x)
  xPts <- 1:length(myVector)
  return(sum(myVector<=(xPts*slope+b)))
}

#This function calculates the cutoff by sliding a diagonal line and finding where it is tangential (or as close as possible)
calculate_cutoff <- function(inputVector, drawPlot=FALSE,...){
  print("this version will try to get more than 600 SEs")
  t1=600
  
  inputVector <- sort(inputVector)
  inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
  
  numberOfSE=0
  while (numberOfSE<t1) {
    slope <- (max(inputVector)-min(inputVector))/(length(inputVector)) #This is the slope of the line we want to slide. This is the diagonal.
   
    xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
    y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
    
    numberOfSE = length(which(inputVector>y_cutoff))
    inputVector=inputVector[-length(inputVector)]
  }
  
  if(drawPlot){  #if TRUE, draw the plot
    plot(1:length(inputVector), inputVector,type="l",...)
    b <- y_cutoff-(slope* xPt)
    abline(v= xPt,h= y_cutoff,lty=2,col=8)
    points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
    abline(coef=c(b,slope),col=2)
    title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
    axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
  }
  return(y_cutoff)
}

######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## select threshold based on Narrow peaks:

threshold_score = getThresholdOnPeakScore(paste0(sampleName,"_peaks.narrowPeak"))
cat (paste("..Minimal peak score set to ",threshold_score,"\n"))

####### read regions and select regions with high score:

regionsToStitch=import.bed6(paste0(sampleName,"_regions.bed"))
tt=which(regionsToStitch$score>=threshold_score)
cat (paste("will select ",length(tt), " regions out of ",length(regionsToStitch)," initial peaks\n"))

regionsToStitch=regionsToStitch[tt]

######## read gene info and extract TSS:
if (substr(transcriptomeFile,start = nchar(transcriptomeFile)-2,nchar(transcriptomeFile))=="gff") {
  transcripts=import(transcriptomeFile)
}
if  (substr(transcriptomeFile,start = nchar(transcriptomeFile)-3,nchar(transcriptomeFile))=="ucsc"){
  transcripts=import.ucsc(transcriptomeFile)
}

cat(paste("..Read file",transcriptomeFile,"\n"))

TSSRegions=GRanges(chrom(transcripts),IRanges(start(transcripts)-distFromTSS,start(transcripts)+distFromTSS-1))
tt=which(strand(transcripts)=='-')
TSSRegions[tt]=GRanges(chrom(transcripts[tt]),IRanges(end(transcripts[tt])-distFromTSS,end(transcripts[tt])+distFromTSS-1))

cat(paste("..Created",length(TSSRegions),"TSS regions\n"))

########## remove TSS peaks from enhancers:
enhancers=setdiff(regionsToStitch,TSSRegions,ignore.strand=T)
promoters=intersect(regionsToStitch,TSSRegions,ignore.strand=T)

########## stitch enhancers:
enhancersStitched=resize(enhancers,maxDistanceToStitch+width(enhancers))
enhancersStitched=reduce(enhancersStitched)
enhancersStitched=resize(enhancersStitched,-maxDistanceToStitch+width(enhancersStitched))

cat(paste("..Created",length(enhancersStitched),"stitched regions\n"))
#export(enhancersStitched,paste0(sampleName,".stitched_regions.bed"))

######### create big wig file out of wig if .bw does not exists ######
hasBW=F
bwFile=paste0(sampleName,".wig.bw")
if(file.exists(bwFile)) {
  hasBW = TRUE
}

if (!hasBW) {
  bwFile=paste0(sampleName,".bw")
  if(file.exists(bwFile)) {
    hasBW = TRUE
  }
}

if (!hasBW) {
  cat ("Will create bw file using 'wigToBigWig'\nCheck what wigToBigWig is added to your PATH!\nEx: PATH=$PATH:/usr/local/bin/ucsc_tools/wigToBigWig\n")
  outSys=system(paste("wigToBigWig -clip", paste0(sampleName,".wig"), faiFileToCreateBigWig, paste0(sampleName,".bw")))
  bwFile=paste0(sampleName,".bw")
  if (outSys==0) {
    cat (paste(bwFile,"was created\n"))
  }
  if(file.exists(bwFile)) {
    hasBW = TRUE
  }else {
    cat(paste0("Could not create bigwig file out of ",paste0(sampleName,".wig"),"\nPlease create this file yourself and rerun\n"))
    cat (paste0("Ex: wigToBigWig -clip ", paste0(sampleName,".wig"), " YOUR_PATH_TO/Human/hg19/hg19.fa.fai ", paste0(sampleName,".bw\n")))
    cat ("wigToBigWig can be downloaded here: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/\n")
    quit()
  }
}

outputFile=paste0(sampleName,".scores.bed")
outputFile=basename(outputFile)
cat (paste("..Printing SEs, enhancers and promoters with their scores into",outputFile,"\n"))
output_peaks_density(enhancersStitched,enhancers,promoters,bwFile,outputFile)



