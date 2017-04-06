#Copywrite Valentina Boeva, 2016
#This code cannot be made public

args <- commandArgs()

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
  return(thresholdToReturn)
}

cat(getThresholdOnPeakScore(args[4]))

