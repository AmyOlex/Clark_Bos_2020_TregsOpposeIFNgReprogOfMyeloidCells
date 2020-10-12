## Amy Olex
## 3/13/2020
## Script to create Figure 7c in mansucript by Clark et al 
## "Regulatory T cells support breast cancer progression by opposing IFN-γ -dependent 
##     functional reprogramming of myeloid cells"
## This derivation of the script adds in the Machrophage signatures from Thorsson, Vésteinn, David L. Gibbs, Scott D. Brown, Denise Wolf, Dante S. Bortone, Tai-Hsien Ou Yang, Eduard Porta-Pardo, et al. “The Immune Landscape of Cancer.” Immunity, April 2018. https://doi.org/10.1016/j.immuni.2018.03.023.



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
immune_file <- "Figure7b_mmc2.txt"
TCGA_data = "Figure7b_DATA_TCGA-BRCA_Expression_92-Gene_Signature_HIGHEXP_log2TPMUpperQuantNorm.txt"

biomarker_exp <- read.delim(file=paste0(datadir,TCGA_data))

## Row Median Center Genes 
rowMedians_biomarker_exp <- rowMedians(as.matrix(biomarker_exp))
biomarker_exp_scaled <- biomarker_exp - rowMedians_biomarker_exp

## import the immune cell signatures.
immune_data <- read.delim(paste0(datadir,immune_file), header=TRUE, row.names=1, stringsAsFactors = FALSE)
row.names(immune_data) <- make.names(row.names(immune_data))
## Extract TCGA patients that match our list.
immune_data_filtered <- immune_data[names(biomarker_exp),,drop=FALSE]
all(names(biomarker_exp) == row.names(immune_data_filtered))
## Extract M1 high/low/NA
M1 <- immune_data_filtered$Macrophages.M1
quantile(na.exclude(M1))
lowM1 = which(M1<.032)
highM1 = which(M1>=.088)
naM1 = which(is.na(M1))
midM1 = setdiff(which(M1==M1), union(union(lowM1,highM1),naM1))

M1[lowM1] <- "low"
M1[highM1] <- "high"
M1[midM1] <- "mid"
M1[naM1] <- "NA"



## Extract M2 high/low/NA
M2 <- immune_data_filtered$Macrophages.M2
quantile(na.exclude(M2))
lowM2 = which(M2<.189)
highM2 = which(M2>=.354)
naM2 = which(is.na(M2))
midM2 = setdiff(which(M2==M2), union(union(lowM2,highM2),naM2))

M2[lowM2] <- "low"
M2[highM2] <- "high"
M2[midM2] <- "mid"
M2[naM2] <- "NA"

## Get a centered color map
centered_colors <- center.palette(biomarker_exp_scaled, palette_length = 100)

jpeg(filename=paste("tmp.jpeg", sep=""))
my_tree <- aheatmap(biomarker_exp_scaled, distfun="euclidean", hclustfun="ward", scale="none", color=centered_colors, treeheight=100)
dev.off()


##########
## Extract Clusters
##########

clusters <- cutree(as.hclust(my_tree$Colv),k=3)

category <- clusters
category[names(clusters)[clusters == 3]] <- "DT-like"
category[names(clusters)[clusters == 2]] <- "Control-like"
category[names(clusters)[clusters == 1]] <- "Intermediate"

annot_data <- data.frame(Cluster=category, Macrophage1=M1, Macrophage2=M2)
my_colors = list(Cluster=c("#000000", "#008000", "#7F7F7F"), Macrophage1=c("red", "blue", "grey","black"),Macrophage2=c("red", "blue", "grey","black"))


png(filename=paste("Figure7b_Heatmap_TCGA-BRCA_Expression_92-Gene_Signature_55HIGHEXP_log2TPMUpperQuantNorm_M1sig.png", sep=""), width = 2000, height=1000, res=200)
aheatmap(biomarker_exp_scaled, distfun="euclidean", hclustfun="ward", scale="none", color=centered_colors, 
         Rowv = F,treeheight=20, annColors=my_colors, annCol=annot_data, fontsize=10, cexRow=.4, cexCol=0)
dev.off()


