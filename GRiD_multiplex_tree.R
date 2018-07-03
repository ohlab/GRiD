library(getopt)
library(gsubfn)
library(gplots)

spec = matrix(c('ifile','i',1,"character",'ofile','o',1,"character"),byrow=TRUE,ncol=4)
opt=getopt(spec)
input <- opt$ifile
outputfile <- opt$ofile

data <- read.csv(input, sep="\t", header=TRUE)
data <- data[,-(3:5)]
df = setNames(data.frame(t(data[,-1])), data[,1]) #transform the data. Here, the rowname is reomved
df <- cbind(Sample = rownames(df), df) #duplicate the 1st row and add a row name called "Sample"
rownames(df) <- NULL #remove the duplicate
row.names(df) <- df$Sample
df <- df[, -1]
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)
df <- as.matrix(df)
df1 <- rbind(df,df)

lmat = rbind(c(0,3),c(2,1),c(0,4))
lwid = c(0.5,9)
lhei = c(4,9,2.7)
pdf(outputfile)
heatmap.2(df1,margin=c(23,0.4),denscol=NA,key.ytickfun = NA,key.ylab=NA,dendrogram="column",key.title=NA,key.xlab="GRiD",cexCol=1,lmat = lmat, lwid = lwid, lhei = lhei,labRow = "",trace="none",col = scaleyellowred)
dev.off()
