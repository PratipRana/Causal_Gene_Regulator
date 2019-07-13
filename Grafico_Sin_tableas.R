library(canvasXpress)
library(ggplot2)
y=read.table("http://www.canvasxpress.org/data/cX-dumbbell-dat.txt", header=TRUE, sep="\t", quote="", row.names=1, fill=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
z=read.table("http://www.canvasxpress.org/data/cX-dumbbell-var.txt", header=TRUE, sep="\t", quote="", row.names=1, fill=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
canvasXpress(
  data=Result_mirPath_Merge_Paper_T,
  
  axisAlgorithm="wilkinson",
  connectBy="Connect",
  dotplotType="stacked",
  graphType="Dotplot",
  showTransition=TRUE,
  smpTitle="MicroRNA",
  sortDir="descending",
  title="Microrna pathway analysis",
    xAxisMinorTicks=FALSE,
  xAxisShow=FALSE,
  xAxisTickFormat="",
  afterRender=list(list("sortSamplesByVariable", list("Men")))
)


x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
y <- as.data.frame(x)
head(y)


Result_mirPath_Merge_Paper=read.csv('Result_mirPath_Merge_Paper_T2.csv', header = TRUE)
ggplot(Result_mirPath_Merge_Paper, # you can replace the numbers to the row number of pathway of your interest
       aes(x =miRNAs , y = KEGG.pathway)) + 
  geom_point(aes(size = genes, color = P_value)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab(NULL) + scale_x_continuous(breaks = c(2,3,4,5), limits = c(2,5))
  ggtitle("MicroRNA pathway")




          