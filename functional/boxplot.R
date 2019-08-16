library("ggpubr")

INPUT="C:/Users/David/Desktop/Sequencing/R_data/boxplotdata-closed.txt"
IMAGE="C:/Users/David/Desktop/Sequencing/R_data/vboxplot_T0-closed.eps"
YLABEL="log10 Distance to T0-closed\n or T0-closed-in-LPS Enhancer (bp)"
TOP=10

data<-read.delim(INPUT,header=TRUE,stringsAsFactors=TRUE)
data$Geneset <- factor(data$Geneset,levels=c("T0-induced","T0-repressed","Random"))

#eps(IMAGE,width=4.5,height=4.5)
#p <- ggboxplot(data,x="Geneset",y="Distance",color="Geneset",
#              palette=get_palette("Dark2",3),add="jitter",
#              shape="group",outlier.shape=NA)
p <- ggviolin(data,x="Geneset",y="Distance",fill="Geneset",
              palette=get_palette("Dark2",3),
              add="boxplot",add.params=list(fill="white"))

my_comparisons <- list( c("T0-induced","Random"),
                        c("T0-repressed","Random"))

p <- p + stat_compare_means(comparisons = my_comparisons,label="p.signif",
                     label.y = c(TOP-1,TOP)) + labs(y=YLABEL)
p <- p +  theme(text=element_text(size=12),
                axis.title.x=element_blank(),
                legend.position="none") + scale_y_continuous(breaks=0:TOP)
  
ggsave(IMAGE,p,width=4,height=4)