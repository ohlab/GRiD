library(dplyr)
library(getopt)
library(ggplot2)
library(gsubfn)

spec = matrix(c('cfile','i',1,"character",'pfix','p',1,"character",'cov','c',1,"double"),byrow=TRUE,ncol=4)
opt=getopt(spec)
inputfile <- opt$cfile
prefix <- opt$pfix
coverage <- opt$cov
title<-strapplyc(prefix, "(.*).........", simplify = TRUE)
slidingwindowplot <- function(windowsize, inputseq)
{
  starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
  n <- length(starts)
  chunkbps <- numeric(n)
  chunkstats<- numeric(n)
  for (i in 1:n) {
    chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)] 
    chunkmean <- mean(chunk)
    chunkstdv<-sd(chunk)
    chunkbps[i] <- chunkmean
    chunkstats[i]<-chunkstdv
  }
  return (list(starts,chunkbps,chunkstats))
}
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
binSize<-10000
column.types <- c("character", "numeric", "numeric")
###########################################################################################################
#  Calculate GRiD 
###########################################################################################################
all.data <- read.csv(inputfile, header=FALSE, sep="\t",colClasses=column.types)
myvector_all<-as.vector(as.matrix(all.data[3]))
windowAll<-slidingwindowplot(binSize,myvector_all)
df<-data.frame(windowAll[[1]],windowAll[[2]],windowAll[[3]])
colname<-c("x","mean","sd")
colnames(df)<-colname

if(median(df$mean) > 0.15){
df[,-1] <- remove_outliers(df$mean)
df[, -1] <- lapply( df[, -1], function(x){ (x/sum(x, na.rm=TRUE))*100} )
df[, -1] <- log2(df$mean)
Plot <- ggplot(data = df, aes(x = x, y = mean))
P2 <- Plot + stat_smooth(aes(outfit=fit<<-..y..),alpha=0,method ="loess",span=1,method.args = list(family="symmetric")) + geom_jitter()+theme_bw()+xlab("Genome Location (bp)")+ylab("Coverage")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+ggtitle(title)
V <-ggplot_build(P2)$data[[1]]

V[V==""] <- NA
V <- na.omit(V)

ter<-filter(V, y == min(y)) %>% select(ymax)
ori<-filter(V, y == max(y)) %>% select(ymin)

fit[fit==""] <- NA
fit <- na.omit(fit)

   GRiD<-(2^ori)/(2^ter)
   maxF<-max(fit)
   minF<-min(fit)
   GRiD_unrefined<- (2^maxF)/(2^minF)
   GRiD_unrefined <- round(GRiD_unrefined, 2)
GRiD_rounded<-round(GRiD, 2)
if(GRiD_rounded < 1){
GRiD_rounded<-1.00
} else {
GRiD_rounded<-GRiD_rounded
}
}else {
GRiD_rounded<-1.00
GRiD_unrefined<-1.00
}

#########################################
species_heterogeneity <- (1-(GRiD_rounded/GRiD_unrefined))
###################################################

merge_data <- paste(title, GRiD_rounded, GRiD_unrefined, species_heterogeneity, coverage, sep = '\t')

write(merge_data, file = prefix)

