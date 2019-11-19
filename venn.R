library(VennDiagram)

sigOE <- dplyr::filter(res_ids, KO_DMSO_PL_3A._vs_WT_DMSO_PL_1A._p.adj< 0.05)
sigOE_genes1 <- as.character(sigOE$GENEID)

sigOE <- dplyr::filter(res_ids, WT_DMSO_PL_1A._vs_WT_PDZ1i_PL_2A._p.adj< 0.05)
sigOE_genes2 <- as.character(sigOE$GENEID)

sigOE <- dplyr::filter(res_ids, KO_DMSO_PL_3A._vs_KO_PDZ1i_PL_4A._p.adj< 0.01)
sigOE_genes3 <- as.character(sigOE$GENEID)

x<-list(sigOE_genes1, sigOE_genes2, sigOE_genes3)

v0<- venn.diagram(
  x,
  category.names = c("KO_DMSO_PL_3A._vs_WT_DMSO_PL_1A" , " WT_DMSO_PL_1A._vs_WT_PDZ1i_PL_2A" , "WT_DMSO_PL_1A._vs_WT_PDZ1i_PL_2A"),
  filename = NULL,
  fill = c("red", "blue", "green"),
  sub.cex=0.5,
  alpha = 0.50,
  col = "transparent")

grid.draw(v0)

overlaps <- calculate.overlap(x)

# extract indexes of overlaps from list names
indx <- as.numeric(substr(names(overlaps),2,2))


# labels start at position 7 in the list for Venn's with 3 circles
for (i in 1:length(overlaps)){
  v0[[6 + indx[i] ]]$label <- paste(overlaps[[i]], collapse = "\n") 
}


grid.newpage()
grid.draw(v0)