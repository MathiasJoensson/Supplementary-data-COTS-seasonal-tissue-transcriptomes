library(clusterProfiler)

#########################################full genome
GOinfo <- read.delim(file.choose(),header=TRUE,sep="\t", stringsAsFactors=FALSE) #GO.tb

TERM2GENE<- read.table(file.choose(),header=T,sep="\t") #gene2go_full_blank.txt
GenesofInterest  <- read.table(file.choose(),header=F,sep=",")

GO <- enricher(GenesofInterest$V1,
               TERM2GENE=TERM2GENE,
               TERM2NAME=GOinfo[1:2],
               pAdjustMethod = "fdr",
               pvalueCutoff = 0.05, 
               qvalueCutoff = 0.05,
               minGSSize = 10,
               maxGSSize = 100000)

GO@result<-merge(GO@result,go2ont(GO@result$ID),by.x= "ID",by.y="go_id")
GO@result<- GO@result[order(GO@result$pvalue),]

sig <- subset(GO@result, p.adjust < 0.05)

write.table(as.data.frame(sig),file="Output.txt", quote=FALSE, sep="\t", row.names =F)





