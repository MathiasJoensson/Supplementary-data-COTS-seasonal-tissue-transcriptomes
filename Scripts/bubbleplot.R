#setwd("C:/Users/marie/OneDrive/Bureau/Pfam_Enrichment/GO enrichment/Seasonal")
library(ggplot2)
library(plyr)
library(RColorBrewer)
#install.packages("colorRamps")


test <- read.delim("bubbleplot_KEGG_testes.txt")
pdf("KEGGenrichment_testes.pdf") 
bubble<-ggplot(test, aes(x=Enrich_ratio, y= Term, size= Input_number, color=qvalue)) +
  geom_point(alpha=0.5) +
  scale_size(range = c(.1, 6), name = "Number of genes") +
  labs(
    y="Pathway",
    fill=NA,
    color="qvalue",
    size="Number of genes"
  ) +
  theme(axis.text.y = element_text(size = 8), # x axis (pfam) domain
        axis.ticks = element_blank(), 
        #panel.grid.major.y = element_line(colour = "grey84"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        #panel.background = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.title.y = element_text(size = 6, angle = 90), #rel(0.7) ie x
        axis.title.x = element_text(size = 6),
        #title = element_text(size = 12),
        legend.text=element_text(size=6), # legend
        legend.justification = "top",
        legend.key.width=unit(0.5, "lines"), # legend width
        legend.key.height=unit(0.5, "lines"), # legend height
        panel.background = element_rect(fill = "white", colour = "grey50"),
        #plot.background=element_blank()
  )
print(bubble + scale_color_gradient(low = "red", high = "blue"))
dev.off()

test <- read.delim("test.txt")
pdf("GO_enrichment.pdf") 
bubble<-ggplot(test, aes(x=Enrichment_ratio, y= GO_Name, size= Nr_Test, color=qvalue)) +
  geom_point(alpha=0.5) +
  scale_size(range = c(.1, 6), name = "Number of genes") +
  labs(
    y="Pathway",
    fill=NA,
    color="qvalue",
    size="Number of genes"
  ) +
  theme(axis.text.y = element_text(size = 8), # x axis (pfam) domain
        axis.ticks = element_blank(), 
        #panel.grid.major.y = element_line(colour = "grey84"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        #panel.background = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.title.y = element_text(size = 6, angle = 90), #rel(0.7) ie x
        axis.title.x = element_text(size = 6),
        #title = element_text(size = 12),
        legend.text=element_text(size=6), # legend
        legend.justification = "top",
        legend.key.width=unit(0.5, "lines"), # legend width
        legend.key.height=unit(0.5, "lines"), # legend height
        panel.background = element_rect(fill = "white", colour = "grey50"),
        #plot.background=element_blank()
  )
print(bubble + scale_color_gradient(low = "red", high = "blue"))
dev.off()