library(getopt)

spec = matrix(c('ifile','i',1,"character",'ofile','o',1,"character",'cut','x',1,"double"),byrow=TRUE,ncol=4)
opt=getopt(spec)
input1 <- opt$ifile
outputfile <- opt$ofile
cutoff <- opt$cut

depth <- read.csv(input1, sep = '\t', header = F)

colname<-c("Genome","base_cov","prop")
colnames(depth)<-colname
cov <- aggregate(prop~Genome,depth,sum)
cov <- cov[-(nrow(cov)),]
cov$prop <- cov$prop - 1
cov<-cov[which(cov[,2]>cutoff),]
passed <- cov[,1, drop=FALSE]
write.table(passed, file = outputfile,quote = F, col.names = F,row.names = F)
