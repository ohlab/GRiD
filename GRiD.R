library(dplyr)
library(getopt)
library(ggplot2)
library(gsubfn)

spec = matrix(c('cfile','i',1,"character",'cfile2','x',1,"character",'cfile3','y',1,"character",'npdf','o',1,"character",'pfix','p',1,"character"),byrow=TRUE,ncol=4)
opt=getopt(spec)
inputfile <- opt$cfile
inputfile2 <- opt$cfile2
inputfile3 <- opt$cfile3
image <- opt$npdf
prefix <- opt$pfix
title<-strapplyc(prefix, "(.*)....", simplify = TRUE)

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
############################################################################################################################
# GRiD calculation using first round of subsampled reads in order to estimate 95% confidence interval 
###########################################################################################################################
all.data <- read.csv(inputfile2, header=FALSE, sep="\t",colClasses=column.types)
myvector_all<-as.vector(as.matrix(all.data[3]))
windowAll<-slidingwindowplot(binSize,myvector_all)
df<-data.frame(windowAll[[1]],windowAll[[2]],windowAll[[3]])
colname<-c("x","mean","sd")
colnames(df)<-colname

df[,-1] <- remove_outliers(df$mean)
df[, -1] <- lapply( df[, -1], function(x){ (x/sum(x, na.rm=TRUE))*100} )
df[, -1] <- log2(df$mean)
P <- ggplot(data = df, aes(x = x, y = mean))
P2 <- P + stat_smooth(aes(outfit=fit<<-..y..),alpha=0,method ="loess",span=1,method.args = list(family="symmetric")) + geom_jitter()+theme_bw()+xlab("Genome Location (bp)")+ylab("Coverage")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+ggtitle(title)
V <-ggplot_build(P2)$data[[1]]

V[V==""] <- NA
V <- na.omit(V)

ter<-filter(V, y == min(y)) %>% select(ymax)
ori<-filter(V, y == max(y)) %>% select(ymin)

fit[fit==""] <- NA
fit <- na.omit(fit)

   GRiD2<-(2^ori)/(2^ter)
   maxF<-max(fit)
   minF<-min(fit)
   #GRiD_unrefined2 <- (2^maxF)/(2^minF)
GRiD_rounded2 <-round(GRiD2, 2)
if(GRiD_rounded2 < 1){
GRiD_rounded2 <- 1.00
} else {
GRiD_rounded2 <- GRiD_rounded2
}
###########################################################################################################################
# GRiD calculation using second round of subsampled reads in order to estimate 95% confidence interval
###########################################################################################################################
all.data <- read.csv(inputfile3, header=FALSE, sep="\t",colClasses=column.types)
myvector_all<-as.vector(as.matrix(all.data[3]))
windowAll<-slidingwindowplot(binSize,myvector_all)
df<-data.frame(windowAll[[1]],windowAll[[2]],windowAll[[3]])
colname<-c("x","mean","sd")
colnames(df)<-colname

df[,-1] <- remove_outliers(df$mean)
df[, -1] <- lapply( df[, -1], function(x){ (x/sum(x, na.rm=TRUE))*100} )
df[, -1] <- log2(df$mean)
P <- ggplot(data = df, aes(x = x, y = mean))
P2 <- P + stat_smooth(aes(outfit=fit<<-..y..),alpha=0,method ="loess",span=1,method.args = list(family="symmetric")) + geom_jitter()+theme_bw()+xlab("Genome Location (bp)")+ylab("Coverage")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+ggtitle(title)
V <-ggplot_build(P2)$data[[1]]

V[V==""] <- NA
V <- na.omit(V)

ter<-filter(V, y == min(y)) %>% select(ymax)
ori<-filter(V, y == max(y)) %>% select(ymin)

fit[fit==""] <- NA
fit <- na.omit(fit)

   GRiD3<-(2^ori)/(2^ter)
   maxF<-max(fit)
   minF<-min(fit)
   #GRiD_unrefined3 <- (2^maxF)/(2^minF)
GRiD_rounded3 <-round(GRiD3, 2)
if(GRiD_rounded3 < 1){
GRiD_rounded3 <- 1.00
} else {
GRiD_rounded3 <- GRiD_rounded3
}

#########################################
GRiD_data <- rbind(GRiD_rounded,GRiD_rounded2,GRiD_rounded3)
GRiD_average <- mean(GRiD_data)
GRiD_average <- round(GRiD_average, 2)
GRiD_CI <- (1.96*(sd(GRiD_data)))/sqrt(nrow(GRiD_data))
GRiD_CI <- round(GRiD_CI, 2)

GRiD_CI_lower <- GRiD_average - GRiD_CI
GRiD_CI_upper <- GRiD_average + GRiD_CI

species_heterogeneity <- (1-(GRiD_rounded/GRiD_unrefined))
###################################################
var1 <- paste(GRiD_CI_lower," - ",GRiD_CI_upper, sep="")
image_title<-paste(title, " =", GRiD_rounded, ", 95% CI = ", var1, sep = ' ')

merge_data <- paste(title, GRiD_rounded, var1, GRiD_unrefined, species_heterogeneity, sep = '\t')
colname<-paste("Sample","GRiD","95% CI","GRiD unrefined","Species heterogeneity", sep = '\t')
output_results <- rbind(colname,merge_data)

write(output_results, file = prefix)

pdf(image)
Plot + stat_smooth(aes(outfit=fit<<-..y..),alpha=0,method ="loess",span=1,method.args = list(family="symmetric")) + geom_jitter()+theme_bw()+xlab("Genome Location (bp)")+ylab("Log2 (% Coverage)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+ggtitle(image_title)
dev.off()
