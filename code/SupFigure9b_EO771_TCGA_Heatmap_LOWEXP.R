## Amy Olex
## 3/13/2020
## Script to create Supplementary Figure 9b in mansucript by Clark et al 
## "Regulatory T cells support breast cancer progression by opposing IFN-γ -dependent 
##     functional reprogramming of myeloid cells"


library(NMF)

center.palette <- function(data, palette_length = 100, color1 = "blue", color2 = "red"){
  
  my_range <- range(data)
  if((my_range[1] >= 0) | (my_range[2] <= 0)){
    stop("Range does not cross zero. Cannot center color palette.")
  }
  
  diff <- my_range[2] - my_range[1]
  length1 <- floor((abs(my_range[1])/diff)*palette_length)
  length2 <- floor((abs(my_range[2])/diff)*palette_length)
  
  
  color1=colorRampPalette(colors=c(color1, "white"))(length1)
  color2=colorRampPalette(colors=c("white", color2))(length2)
  #my_colors <- c(color1, color2)
  
  #my_breaks <- seq(from = my_range[1], to = my_range[2], length.out = length(my_colors)+1)
  
  #return(list(colors = my_colors, breaks = my_breaks))
  return(c(color1, color2))
}



##############
## Import Data and metadata
##############
## TCGA data from the BRCA project was obtained from a TCGABioLinks Query.
## The scaled_estimate data was extracted, multiplied by 1 million, then log2 transformed to obtain TPM values.
## TPM values were then upper quantile normalized and filtered for only those genes in the biomarker signature.
## This processed and filtered TCGA data is provided in the TCGA_data file.

datadir = "../data/"
TCGA_data_LOW = "SupFigure9b_DATA_TCGA-BRCA_Expression_92-Gene_Signature_LOWEXP_log2TPMUpperQuantNorm.txt"

## Print out LOW gene sug for Sup Fig 9b
biomarker_exp <- read.delim(file=paste0(datadir,TCGA_data_LOW))

## Row Median Center Genes 
rowMedians_biomarker_exp <- rowMedians(as.matrix(biomarker_exp))
biomarker_exp_scaled <- biomarker_exp - rowMedians_biomarker_exp


## Get a centered color map
centered_colors <- center.palette(biomarker_exp_scaled, palette_length = 100)

jpeg(filename=paste("tmp.jpeg", sep=""))
my_tree <- aheatmap(biomarker_exp_scaled, distfun="euclidean", hclustfun="ward", scale="none", color=centered_colors, treeheight=100)#annCol=annot_data, annColors=my_colors, fontsize=10, cexRow=4, treeheight=150)
dev.off()


##########
## Extract Clusters
##########

clusters <- cutree(as.hclust(my_tree$Colv),k=3)

category <- clusters
category[names(clusters)[clusters == 3]] <- "DT-like"
category[names(clusters)[clusters == 2]] <- "Control-like"
category[names(clusters)[clusters == 1]] <- "Intermediate"

annot_data <- data.frame(Cluster=category)
my_colors = list(Cluster=c("#000000", "#008000", "#7F7F7F"))


png(filename=paste("SuppFigure9b_Heatmap_TCGA-BRCA_Expression_92-Gene_Signature_16LOWEXP_log2TPMUpperQuantNorm.png", sep=""), width = 2000, height=1000, res=200)
aheatmap(biomarker_exp_scaled, distfun="euclidean", hclustfun="ward", scale="none", color=centered_colors, 
         Rowv = F,treeheight=20, annColors=my_colors, annCol=annot_data, fontsize=10, cexRow=.7, cexCol=0)
dev.off()

## Save Expression Data To file
## order to be same as heatmap:
ordered_clusters <- clusters[my_tree$colInd]

all(names(biomarker_exp) == names(biomarker_exp_scaled))
all(row.names(biomarker_exp) == row.names(biomarker_exp_scaled))

biomarker_exp_t <- as.data.frame(t(biomarker_exp))
biomarker_exp_t <- biomarker_exp_t[names(ordered_clusters),,drop=F]
all(row.names(biomarker_exp_t) == names(ordered_clusters))

## re-order genes
biomarker_exp_t <- biomarker_exp_t[,my_tree$rowInd,drop=F]

biomarker_exp_t$EO771.geneSig.clusters <- ordered_clusters

## Save all cluster data to a file
write.table(biomarker_exp_t, file="SuppFigure9b_DATA_TCGA-BRCA_Expression_92-Gene_Signature_16LOWEXP_log2TPMUpperQuantNorm.txt", sep="\t", quote=FALSE)


