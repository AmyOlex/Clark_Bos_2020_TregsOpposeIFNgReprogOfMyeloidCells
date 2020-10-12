## Amy Olex
## 3/13/2020
## Script to create Figure 7d and 7e in mansucript by Clark et al 
## "Regulatory T cells support breast cancer progression by opposing IFN-γ -dependent 
##     functional reprogramming of myeloid cells"
library(survival)
library(survminer)

datadir <- "../data/"
survival_file <- "Figure7cde_BRCA_Survival_AllPatients.tsv"
cluster_file <- "Figure7cde_EO771_55_HIGHEXPGeneSig_TCGA-BRCA_Clusters.txt"
er_file <- "Figure7de_TCGA_BRCA_ERstatus.txt"

## read in survival data
survival_data <- read.delim(paste(datadir, survival_file,sep=""), sep="\t", stringsAsFactors = F)
row.names(survival_data) <- survival_data$submitter_id

## read in EO771 Clusters
cluster_data <- read.delim(paste(datadir, cluster_file,sep=""), sep="\t", stringsAsFactors = F)
row.names(cluster_data) <- cluster_data$Patient.ID

## import the TCGA ER status clinical data:
er_status <- read.delim(paste(datadir, er_file,sep=""), header=T, sep="\t")
row.names(er_status) <- er_status$patient.id

## get patient ids common to all three files
common_ids <- intersect(intersect(row.names(survival_data), row.names(cluster_data)), row.names(er_status))

## crop survival data by columns and filter by common ids
survival_data <- survival_data[common_ids,c("time", "censored"),drop=F]

## crop cluster data to the common_ids
cluster_data <- cluster_data[common_ids,,drop=F]

## crop er_status to the common IDs
er_status <- er_status[common_ids,,drop=F]


### Sanity check
all(row.names(cluster_data) == row.names(survival_data))
all(row.names(cluster_data) == row.names(er_status))
all(row.names(er_status) == row.names(survival_data))




## Augment the survival data
survival_data$event <- survival_data$censored
survival_data[row.names(survival_data)[survival_data$censored == "true"],"event"] <- 0
survival_data[row.names(survival_data)[survival_data$censored == "false"],"event"] <- 1
survival_data$event <- as.numeric(survival_data$event)

survival_data$cluster <- cluster_data$EO771.geneSig.clusters
survival_data$er.status <- er_status$er.status.ihc

## Filter to only those with ER+ status
survival_data_ERpos <- survival_data[survival_data$er.status=="Positive",,drop=F]
survival_data_ERneg <- survival_data[survival_data$er.status=="Negative",,drop=F]


############### Using all patients
########## ER Positive Survival
## Filter ER+ Survival Data to 10-year Survival (3655 days)
ten_yr_pos <- survival_data_ERpos
ten_yr_pos$event[ten_yr_pos$time >= 3655] <- 0
ten_yr_pos$time[ten_yr_pos$time >= 3655] <- 3655

fit_pos <- survfit(Surv(ten_yr_pos$time, ten_yr_pos$event) ~ cluster, data = ten_yr_pos)

png(filename=paste("Figure7d_alt_10yr_Survival_ER-Positive_EO771_92GeneSignature_HIGHEXP_TCGA-BRCA_normalized_All_1079_Patients_withIntermediates.png", sep=""), width=1800, height=1200, res=200)

ggsurvplot(fit_pos, pval = TRUE, pval.method = TRUE, surv.median.line = "hv", conf.int = FALSE, palette = c("#000000", "#008000", "#7F7F7F"), 
           legend = "bottom", fontsize = 10, xlab = "Number of Days", legend.title = "", legend.labs = c("Control-Like", "DT-Like", "Intermediate"))

dev.off()

########## ER Negative Survival
## Filter ER+ Survival Data to 10-year Survival (3655 days)
ten_yr_neg <- survival_data_ERneg
ten_yr_neg$event[ten_yr_neg$time >= 3655] <- 0
ten_yr_neg$time[ten_yr_neg$time >= 3655] <- 3655

fit_neg <- survfit(Surv(ten_yr_neg$time, ten_yr_neg$event) ~ cluster, data = ten_yr_neg)

png(filename=paste("Figure7e_alt_10yr_Survival_ER-Negative_EO771_92GeneSignature_HIGHEXP_TCGA-BRCA_normalized_All_1079_Patients_withIntermediates.png", sep=""), width=1800, height=1200, res=200)

ggsurvplot(fit_neg, pval = TRUE, pval.method = TRUE, surv.median.line = "hv", conf.int = FALSE, palette = c("#000000", "#008000", "#7F7F7F"), 
           legend = "bottom", fontsize = 10, xlab = "Number of Days", legend.title = "", legend.labs = c("Control-Like", "DT-Like", "Intermediate"))

dev.off()


########## no Intermediates
########## ER Positive Survival
## Filter ER+ Survival Data to 10-year Survival (3655 days)
ten_yr_pos <- survival_data_ERpos
ten_yr_pos$event[ten_yr_pos$time >= 3655] <- 0
ten_yr_pos$time[ten_yr_pos$time >= 3655] <- 3655

ten_yr_pos <- ten_yr_pos[ten_yr_pos$cluster %in% c("DT-like","Control-like"),,drop=F]
ten_yr_pos_surv <- Surv(ten_yr_pos$time, ten_yr_pos$event)
ten_yr_pos_hr <- exp(coxph(ten_yr_pos_surv ~ cluster, ten_yr_pos)$coefficients)

fit_pos <- survfit(Surv(ten_yr_pos$time, ten_yr_pos$event) ~ cluster, data = ten_yr_pos)

png(filename=paste("Figure7d_10yr_Survival_ER-Positive_EO771_92GeneSignature_HIGHEXP_TCGA-BRCA_normalized_All_1079_Patients_noIntermediate.png", sep=""), width=1800, height=1200, res=200)

ggsurvplot(fit_pos, pval = TRUE, pval.method = TRUE, surv.median.line = "hv", conf.int = FALSE, palette = c("#000000", "#008000"), 
           legend = "bottom", fontsize = 10, xlab = "Number of Days", legend.title = paste("HR = ", round(ten_yr_pos_hr, digits=2), sep=""), legend.labs = c("Control-Like", "DT-Like"))

dev.off()

########## ER Negative Survival
## Filter ER+ Survival Data to 10-year Survival (3655 days)
ten_yr_neg <- survival_data_ERneg
ten_yr_neg$event[ten_yr_neg$time >= 3655] <- 0
ten_yr_neg$time[ten_yr_neg$time >= 3655] <- 3655

ten_yr_neg <- ten_yr_neg[ten_yr_neg$cluster %in% c("DT-like","Control-like"),,drop=F]
ten_yr_neg_surv <- Surv(ten_yr_neg$time, ten_yr_neg$event)
ten_yr_neg_hr <- exp(coxph(ten_yr_neg_surv ~ cluster, ten_yr_neg)$coefficients)

fit_neg <- survfit(Surv(ten_yr_neg$time, ten_yr_neg$event) ~ cluster, data = ten_yr_neg)

png(filename=paste("Figure7e_10yr_Survival_ER-Negative_EO771_92GeneSignature_HIGHEXP_TCGA-BRCA_normalized_All_1079_Patients_noIntermediate.png", sep=""), width=1800, height=1200, res=200)

ggsurvplot(fit_neg, pval = TRUE, pval.method = TRUE, surv.median.line = "hv", conf.int = FALSE, palette = c("#000000", "#008000"), 
           legend = "bottom", fontsize = 10, xlab = "Number of Days", legend.title = paste("HR = ", round(ten_yr_neg_hr, digits=2), sep=""), legend.labs = c("Control-Like", "DT-Like"))

dev.off()
