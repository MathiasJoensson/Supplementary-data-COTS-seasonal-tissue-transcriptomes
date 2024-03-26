library(pheatmap)
library(RColorBrewer)
#install.packages("colorRamps")
library(colorRamps)
hmcol <- c(colorRampPalette(c("blue", "black", "red"))(n = 9))
hmcol = blue2red (400)
#hmcol = colorRampPalette("YlOrBr")

coul <- c(colorRampPalette(c("blue", "white", "red"))(n=25))

a <- read.table("All.txt", sep="\t", row.names = 1, header = TRUE)

#test <- a[apply(a, 1, function(x) sd(x)!=0),]
pheatmap(as.matrix(a), 
         #color = coul,
         color = colorRampPalette(brewer.pal(n=7, name = "PuBuGn"))(11),
         cluster_cols=F, 
         cluster_rows=T, 
         show_rownames=T,
         show_colnames=TRUE, 
         scale = "row", 
         fontsize_row=3, 
         cellheight = 4,
         cellwidth = 15,
         #gaps_col = c(7,13),
         clustering_method = "complete",
         #breaks = NULL,
         border_color = "NA",
         #gaps_row = c(3,17,59,66,88,104,108,110,117),
         filename="All_smr_heatmap.pdf")

#Quartile script
a <- read.table("orexin receptors quartile.txt", sep="\t", row.names = 1, header = TRUE)
library(colorRamps)
hmqcol <- colorRampPalette(c("blue", "white", "firebrick2"))(n=100)
hmqcol = blue2red(5)
qcoul <- c(colorRampPalette(c("blue", "white", "red"))(n=5))
#test <- a[apply(a, 1, function(x) sd(x)!=0),]
pheatmap(as.matrix(a), 
         color = hmqcol,
         #color = colorRampPalette(rev(brewer.pal(n=11, name = "RdBu")))(11),
         cluster_cols=F, 
         cluster_rows=F, 
         show_rownames=T,
         show_colnames=TRUE, 
         #scale = "row", 
         #border_color = "NA", 
         fontsize_row=3, 
         #cellheight = 4,
         cellwidth = 15,
         #gaps_col = c(3,6,9,12,15,18,21),
         #clustering_method = "complete",
         #breaks = NULL,
         border_color = "NA",
         #gaps_row = c(3,17,59,66,88,104,108,110,117),
         filename="quartile.pdf")