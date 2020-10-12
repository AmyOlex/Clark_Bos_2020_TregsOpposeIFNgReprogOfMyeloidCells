## Amy Olex
## 3/13/2020
## Script to create Figure 7c in mansucript by Clark et al 
## "Regulatory T cells support breast cancer progression by opposing IFN-γ -dependent 
##     functional reprogramming of myeloid cells"
library(survival)
library(survminer)

datadir <- "../data/"
survival_file <- "Figure7cde_BRCA_Survival_AllPatients.tsv"
cluster_file <- "Figure7cde_EO771_55_HIGHEXPGeneSig_TCGA-BRCA_Clusters.txt"

## read in survival data
survival_data <- read.delim(paste(datadir, survival_file,sep=""), sep="\t", stringsAsFactors = F)
row.names(survival_data) <- survival_data$submitter_id

## read in EO771 Clusters
cluster_data <- read.delim(paste(datadir, cluster_file,sep=""), sep="\t", stringsAsFactors = F)
row.names(cluster_data) <- cluster_data$Patient.ID


## crop survival data by columns and filter to only those in our dataset
survival_data <- survival_data[row.names(survival_data) %in% row.names(cluster_data),c("time", "censored"),drop=F]

## crop cluster data to those data point we have survivial data for
cluster_data <- cluster_data[row.names(cluster_data) %in% row.names(survival_data),,drop=F]

## re-order survival data to same row order as cluster data
survival_data <- survival_data[row.names(cluster_data),,drop=F]

all(row.names(cluster_data) == row.names(survival_data))


survival_data$event <- survival_data$censored
survival_data[row.names(survival_data)[survival_data$censored == "true"],"event"] <- 0
survival_data[row.names(survival_data)[survival_data$censored == "false"],"event"] <- 1
survival_data$event <- as.numeric(survival_data$event)


survival_data$cluster <- cluster_data$EO771.geneSig.clusters

######### Using all patients

ten_yr <- survival_data
ten_yr$event[ten_yr$time >= 3655] <- 0
ten_yr$time[ten_yr$time >= 3655] <- 3655

ten_yr_surv <- Surv(ten_yr$time, ten_yr$event)
fit <- survfit(Surv(ten_yr$time, ten_yr$event) ~ cluster, data = ten_yr)

png(filename=paste("2018.07.12_10yr_Survival_EO771_92GeneSignature_HIGHEXP_TCGA-BRCA_normalized_All_1079_Patients.png", sep=""), width=1800, height=1200, res=200)

ggsurvplot(fit, pval = TRUE, pval.method = TRUE, surv.median.line = "hv", conf.int = FALSE, palette = c("#000000", "#008000", "#7F7F7F"), 
           legend = "bottom", fontsize = 10, xlab = "Number of Days", legend.title = "", legend.labs = c("Control-Like", "DT-Like", "Intermediate"))

dev.off()

### Remove Intermediates
ten_yr <- ten_yr[ten_yr$cluster %in% c("DT-like","Control-like"),,drop=F]
ten_yr_surv <- Surv(ten_yr$time, ten_yr$event)
ten_yr_hr <- exp(coxph(ten_yr_surv ~ cluster, ten_yr)$coefficients)

fit <- survfit(Surv(ten_yr$time, ten_yr$event) ~ cluster, data = ten_yr)

png(filename=paste("2018.07.12_10yr_Survival_EO771_92GeneSignature_HIGHEXP_TCGA-BRCA_normalized_All_1079_Patients_noIntermediate.png", sep=""), width=1800, height=1200, res=200)

ggsurvplot(fit, pval = TRUE, pval.method = TRUE, surv.median.line = "hv", conf.int = FALSE, palette = c("#000000", "#008000"), 
           legend = "bottom", fontsize = 10, xlab = "Number of Days", legend.title = paste("HR = ", round(ten_yr_hr, digits=2), sep=""), legend.labs = c("Control-Like", "DT-Like"))

dev.off()

