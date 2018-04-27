library(getopt)
library(gsubfn)
library(gplots)

spec = matrix(c('ifile','i',1,"character",'cfile','c',1,"character",'ofile','o',1,"character",'cut','x',1,"double"),byrow=TRUE,ncol=4)
opt=getopt(spec)
input1 <- opt$ifile
input2 <- opt$cfile
outputfile <- opt$ofile
cutoff <- opt$cut

depth <- read.csv(input1,sep = '\t', header = F)
length <- read.csv(input2,sep = '\t', header = F)
colname<-c("Genome","base_cov")
colnames(depth)<-colname
colname<-c("Genome","length")
colnames(length)<-colname

cov_temp <- aggregate(base_cov~Genome,depth,sum)
cov <- merge(cov_temp, length, by="Genome")
cov$coverage <- cov$base_cov/cov$length
cov<-cov[which(cov[,4]>cutoff),]
passed <- cov[,1, drop=FALSE]
write.table(passed, file = outputfile,quote = F, col.names = F,row.names = F)
