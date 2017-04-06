#Copywrite Valentina Boeva, 2016
#This code cannot be made public


########################################################################################################
#Will create a BED file with scores; will sort regions from top score to the lowest score
#To run this script: Rscript normalizeHMCanDensity.R <gff or bed - input type> <BED of GFF1 peak file> <bw file> <output file> 

########################################################################################################



library(rtracklayer)


#============================================================================
#==============SUPER-ENHANCER CALLING AND PLOTTING FUNCTIONS=================
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
 	inputVector <- sort(inputVector)
	inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
	slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
	xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
	y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
	
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


output_peaks_density_bed <- function(info){
  bedRegions = import.bed(info$peak)
  strand(bedRegions)="*"
  bedRegions$score=0
  densities = import.bw(info$density)  
  gc()
  overlaps=findOverlaps(bedRegions,densities)  
  
  ensids <- densities$score[subjectHits(overlaps)]
  x <- tapply(ensids, queryHits(overlaps), sum)
  bedRegions$score[unique(queryHits(overlaps))]=x
  rm(densities)
  rm(overlaps)
  gc()  
  strand(bedRegions)="+"
  bedRegions=bedRegions[order(bedRegions$score,decreasing=T)]  
  cutoff=calculate_cutoff(bedRegions$score,F)
  print(cutoff)
  bedRegions$name="enhancer"
  bedRegions$name[which(bedRegions$score>cutoff)]="SE"
  export.bed(bedRegions,info$output)
}

output_peaks_density_gff <- function(info){
  bedRegions = import.gff1(info$peak)
  strand(bedRegions)="*"
  bedRegions$score=0
  densities = import.bw(info$density)  
  gc()
  overlaps=findOverlaps(bedRegions,densities)  
  
  ensids <- densities$score[subjectHits(overlaps)]
  x <- tapply(ensids, queryHits(overlaps), sum)
  bedRegions$score[unique(queryHits(overlaps))]=x
  rm(densities)
  rm(overlaps)
  gc()  
  strand(bedRegions)="+"
  bedRegions=bedRegions[order(bedRegions$score,decreasing=T)]
  cutoff=calculate_cutoff(bedRegions$score,F)
  bedRegions$name="enhancer"
  bedRegions$name[which(bedRegions$score>cutoff)]="SE"
  export.bed(bedRegions,info$output)
}

args<- commandArgs()
files_info = NULL
files_info$inputType=args[6]
files_info$peak=args[7]
files_info$density=args[8]
files_info$output=args[9]

print(files_info) 

if (files_info$inputType=="bed") {
  output_peaks_density_bed(files_info)
}else {
  output_peaks_density_gff(files_info)
}
