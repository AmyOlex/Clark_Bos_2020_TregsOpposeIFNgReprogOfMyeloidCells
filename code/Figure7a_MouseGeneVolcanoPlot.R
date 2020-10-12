## Amy Olex
## 3/13/2020
## Script to create Figure 7a in mansucript by Clark et al 
## "Regulatory T cells support breast cancer progression by opposing IFN-γ -dependent 
##     functional reprogramming of myeloid cells"

library("ggplot2")
library("ggrepel")

## import DEG analysis file
degs <- read.delim("../data/Figure7a_Data_Mouse_DESeq2_EO771-CodingOnly_DEGs_CTL.vs.DT.txt")

## filter to only those genes with an adjusted p-value (FDR)
degs_filt <- degs[!is.na(degs$padj),,drop=F]

## Set colors for points
degs_filt$padj.neglog10 <- -log10(degs_filt$padj)
degs_filt$threshold <- rep("GREY", length(degs_filt$padj))
degs_filt$threshold[degs_filt$padj <= 0.05 & abs(degs_filt$log2FoldChange) >= 1.4] <- "RED"

## Create volcano plot

png(filename="Figure7a_MouseExpression_Volcano.png", width=1000, height=1000, res=200)

  ggplot(degs_filt, aes(x=log2FoldChange, y=padj.neglog10)) +
  geom_point(aes(colour = threshold), size=1.3) +
  scale_colour_manual(values = c("GREY"= "grey", "RED"="red")) +
  #scale_y_continuous(trans = log10_trans(), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = .3) +
  geom_vline(xintercept = 1.4, linetype = "dashed", color = "black", size = .3) +
  geom_vline(xintercept = -1.4, linetype = "dashed", color = "black", size = .3) +
  labs(x = "Log2 FC: DT-treated vs Control", y="-Log10(p-value)") +
  geom_text_repel(
  data = subset(degs_filt, (gene.symbol %in% c("Cxcl10", "Cxcl9", "Cd274", "Ifng", "Ido1", "Cxcl11"))),
  aes(label = gene.symbol), size = 4, nudge_x=2, force=5, segment.size=.2) +
  theme(legend.position="none")

dev.off()  


  