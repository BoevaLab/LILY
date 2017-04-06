#Copywrite Valentina Boeva, 2016

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
# cat makeGraph.R | R --slave --args < directory with bw and narrowPeak files > < mark >
#################################################################################
args <- commandArgs()

MYDIRWITHWIGS=args[4]

cat (paste("working with .bw and .peaks.narrowPeak files in",MYDIRWITHWIGS))

setwd(MYDIRWITHWIGS)

type = args[5] # e.g "K27ac"

files=dir(MYDIRWITHWIGS)
library(rtracklayer)

import.bed3<- function(filename){
  peaks=read.table(filename)
  return(GRanges(peaks[,1],IRanges(peaks[,2],peaks[,3])))
}

get_peaks_density <- function(info){
  peaks = import.bed3(info["peak"])
  densities = import.bw(info["density"])
  rangeToConsider=c(100:5000)
  peak_density = subsetByOverlaps(densities,peaks[rangeToConsider])  
  total_density = median(peak_density$score)
  gc()
  return(total_density)
}


scale_wig_file <- function(info){
  gc()
  densities = import.bw(info["density"])
  print("densities read")
  gc()
  densities$score = as.numeric(info["scalingFactor"])*densities$score
  print("densities recalculated")  
  export.bw(densities,BigWigFile(info["output"]))
  print("densities printed")  
  gc()
}

#here it tries to get filenames with narrow peaks transformed to the BED format:
files.regions=files[which(  (sapply(as.character(files),FUN=function(x) {strsplit(x,split="_",fixed = F)[[1]][2]})=="peaks.narrowPeak" | sapply(as.character(files),FUN=function(x) {strsplit(x,split="_",fixed = F)[[1]][3]})=="peaks.narrowPeak")
		& (sapply(as.character(files),FUN=function(x) {strsplit(x,split=".",fixed = T)[[1]][2]})==type)  )]

#here it tries to get corresponding sample names:
samplesToProcess=sapply(as.character(files.regions),FUN=function(x) {substr(x,1,nchar(x)-nchar("peaks.narrowPeak")-1)})

#here it tries to get bigWig files to process:
files.Wig=paste0(samplesToProcess,".wig.bw")

files_info=cbind(files.regions,files.Wig,paste0(samplesToProcess,".renorm.bw"))
colnames(files_info) = c("peak","density","output")

total_densities = apply(files_info,1,get_peaks_density)

reference = median(total_densities)
scalingFactor = reference/total_densities
cbind(scalingFactor,samplesToProcess)

files_info = cbind(files_info,scalingFactor)
colnames(files_info) = c("peak","density","output","scalingFactor")


for (i in c(1:nrow(files_info))) {
  print (files_info[i,3])
  scale_wig_file(files_info[i,])
}

