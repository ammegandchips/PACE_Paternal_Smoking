######################################################
# PACE PATERNAL SMOKING PHASE ONE ANALYSIS PLAN CODE # 
#                    GEMMA SHARP                     #
#                     28/12/2017                     #
######################################################

##################################################################################################################################################################################
# The following R code will allow you to complete all the EWAS requested in the PACE Paternal Smoking analysis plan.
# If you have insufficient data to complete one or more of the EWAS, you can just skip those models.
# The code also produces .csv files summarising the variables included in the EWASs.
# You shouldn't have to rewrite or add to the following code, unless otherwise stated.
# There are just two inputs required for this analysis:
# 1) pheno: a dataframe containing all the "phenotype" data needed for this project. 
#    Each row is a sample(individual) and each column is a different variable. 
#    Necessary variable names are: "pat.active.smoking", "mat.active.smoking","mat.passive.smoking","sex","ses","pat.age","mat.age","pat.bmi","mat.bmi","parity"
#    If these columns are named differently in your dataset, please rename the columns accordingly
#    Details on how to code these variables are provided in the analysis plan.
# 2) meth: a matrix of methylation illumina beta values. Each column is a sample and each row is a probe on the array (450k or EPIC). 
#    Column names must correspond to the sample.id column in pheno.
##################################################################################################################################################################################

# Load required packages (if these are not already installed, you will have to install them as the first step)
library(sva)
library(tableone)
library(matrixStats)
library(limma)

# Setup the necessary functions
## Function to remove outliers using the IQR*3 (Tukey) method
IQR.removal <- function(meth.matrix){
  rowIQR <- rowIQRs(meth.matrix, na.rm = T)
  row2575 <- rowQuantiles(meth.matrix, probs = c(0.25, 0.75), na.rm = T)
  maskL <- meth.matrix < row2575[,1] - 3 * rowIQR 
  maskU <- meth.matrix > row2575[,2] + 3 * rowIQR 
  meth.matrix[maskL] <- NA
  meth.matrix[maskU] <- NA
  meth.matrix
}
## Function to generate surrogate variables and merge them with the phenotype data (used to adjust for batch)
SVA.generate <- function(meth.matrix, pheno.data, variable.of.interest, model.covariates,n.sv){
  intersecting.samples <- intersect(pheno.data$sample.id,colnames(meth.matrix))
  pheno.data <- na.omit(pheno.data[which(pheno.data$sample.id %in% intersecting.samples),unique(c("sample.id",variable.of.interest,model.covariates))])
  meth.matrix <- meth.matrix[,match(pheno.data$sample.id,colnames(meth.matrix))]
  k = which(is.na(meth.matrix), arr.ind=TRUE)
  meth.matrix[k] = rowMedians(meth.matrix, na.rm=TRUE)[k[,1]]
  mod = model.matrix(reformulate(paste0("pheno.data$",colnames(pheno.data[-1]))))
  mod0 = mod[,-grep(paste0(variable.of.interest,collapse="|"),colnames(mod))]
  sva.ret = sva(meth.matrix, mod=mod, mod0=mod0, n.sv=n.sv)
  SVs = as.data.frame(sva.ret$sv)
  colnames(SVs) <-paste0("sv",1:ncol(SVs))
  cbind(pheno.data,SVs)
}
## Function to run EWAS
ewas.function <-  function(meth.matrix, pheno.data, variable.of.interest){   
  meth.matrix <- meth.matrix[,match(pheno.data$sample.id,colnames(meth.matrix))]
  model.covariates <- colnames(pheno.data)[-which(colnames(pheno.data) %in% c(variable.of.interest,"sample.id"))]
  des = model.matrix(reformulate(paste0("pheno.data$",c(variable.of.interest,model.covariates))))
  fit = lmFit(meth.matrix, des)
  fit.ebayes = eBayes(fit)
  n = rowSums(!is.na(meth.matrix))
  se = (sqrt(fit.ebayes$s2.post) * fit.ebayes$stdev.unscaled[,grep(paste0(variable.of.interest,collapse="|"),colnames(fit.ebayes$stdev.unscaled))])
  res = data.frame(n=n,
                 coef=fit.ebayes$coefficient[,grep(paste0(variable.of.interest,collapse="|"),colnames(fit.ebayes$coefficient))],
                 se=se,
                 p=fit.ebayes$p.value[,grep(paste0(variable.of.interest,collapse="|"),colnames(fit.ebayes$p.value))])
  res
}

# Set initial parameters
study <- "ALSPAC" #change to your study identifier
timepoint <- "birth" #change depending on the age of the children with methylation samples. Can be "birth", "early_childhood", "late_childhood", "adolescence" or "adult"
cell.names <- if(timepoint=="birth"){
  c("nk","gran","bcell","cd8t","cd4t","mono","nrbc")
    }else{
     c("nk","gran","bcell","cd8t","cd4t","mono")
   }
traits.and.covariates <- c("pat.active.smoking", "mat.active.smoking","mat.passive.smoking","sex","ses","pat.age","mat.age","pat.bmi","mat.bmi","parity")
covariates <- c("pat.age", "pat.bmi", "ses", "mat.age", "mat.bmi", "parity", cell.names)

# Load and check phenotype data
pheno <- read.csv("EWAS/pat_smoke/phenofile.alspac.csv",header=TRUE,stringsAsFactors=FALSE) #change filename/location to point to your phenotype file

for(i in 1:length(c("sample.id",traits.and.covariates,cell.names))) {
print(ifelse(c("sample.id",traits.and.covariates,cell.names)[i] %in% colnames(pheno)==FALSE,
  paste("CAUTION: the variable called",c("sample.id",traits.and.covariates,cell.names)[i],"is missing from pheno"),
  paste("variable called",c("sample.id",traits.and.covariates,cell.names)[i],"is present in pheno")))
}

table(pheno$mat.passive.smoking,(pheno$mat.active.smoking|pheno$pat.active.smoking)) #checking number of maternal passive smokers with no active smoking for either mother or father
pheno$mat.passive.smoking[which(pheno$pat.active.smoking==1 |pheno$mat.active.smoking==1)] <- NA #setting maternal passive smoking to missing if mother or father is an active smoker

# Load methylation data and perform QC (e.g. filter probes with high detection P-values)
load("/panfs/panasas01/dedicated-mrcieu/studies/latest/alspac/epigenetic/methylation/450k/aries/released/2016-05-03/data/betas/data.Robj") #change filename/location to point to your methylation data (meth)
# Perform any cohort-specific QC at this point (e.g. you might want to remove probes with high detection p-values)

# IQR*3 method to remove outliers (if this has not already been applied to your data)
log.iqr <- data.frame(cpgs = row.names(meth),NAs.before.IQR3 = rowSums(is.na(meth)))
meth <- IQR.removal(meth)
log.iqr$NAs.after.IQR3 <- rowSums(is.na(meth))
save(log.iqr, file=paste0(study,".patsmoking.logIQR.",timepoint,".Rdata"))

# Generate surrogate variables for technical batch and merge with pheno data to prepare pheno data frames for each EWAS
pheno.min.mutual <- SVA.generate(meth, pheno, variable.of.interest = c("pat.active.smoking","mat.active.smoking"), model.covariates = NULL,n.sv=20)
pheno.min.pat <-pheno.min.mutual[,-which(colnames(pheno.min.mutual) %in% "mat.active.smoking")]
pheno.min.mat <-pheno.min.mutual[,-which(colnames(pheno.min.mutual) %in% "pat.active.smoking")]
pheno.covs.mutual <- SVA.generate(meth, pheno, variable.of.interest = c("pat.active.smoking","mat.active.smoking"), model.covariates = c(covariates,"sex"),n.sv=20)
pheno.covs.pat <-pheno.covs.mutual[,-which(colnames(pheno.covs.mutual) %in% "mat.active.smoking")]
pheno.covs.mat <-pheno.covs.mutual[,-which(colnames(pheno.covs.mutual) %in% "pat.active.smoking")]
pheno.covs.passive <- SVA.generate(meth, pheno, variable.of.interest = "mat.passive.smoking", model.covariates = c(covariates,"sex"),n.sv=20)
pheno.covs.pat.only <-  pheno.covs.mutual[which(pheno.covs.mutual$mat.active.smoking == 0),]
pheno.covs.pat.only <-  pheno.covs.pat.only[, -which(colnames(pheno.covs.pat.only) %in% c("mat.active.smoking"))]
pheno.covs.mutual.boys.only <- pheno.covs.mutual[which(pheno.covs.mutual$sex == 0),]
pheno.covs.mutual.girls.only <- pheno.covs.mutual[which(pheno.covs.mutual$sex == 1),]

# Summarise pheno data and save summaries as .csv files
min.pat.tableone <- print(CreateTableOne(data=pheno.min.mutual[,-1],strata="pat.active.smoking",factorVars=c("pat.active.smoking","mat.active.smoking","ses","parity","sex")))
min.mat.tableone <- print(CreateTableOne(data=pheno.min.mutual[,-1],strata="mat.active.smoking",factorVars=c("pat.active.smoking","mat.active.smoking","ses","parity","sex")))
covs.pat.tableone <- print(CreateTableOne(data=pheno.covs.mutual[,-1],strata="pat.active.smoking",factorVars=c("pat.active.smoking","mat.active.smoking","ses","parity","sex")))
covs.mat.tableone <- print(CreateTableOne(data=pheno.covs.mutual[,-1],strata="mat.active.smoking",factorVars=c("pat.active.smoking","mat.active.smoking","ses","parity","sex")))
covs.pat.only.tableone <- print(CreateTableOne(data=pheno.covs.pat.only[,-1],strata="pat.active.smoking",factorVars=c("pat.active.smoking","ses","parity","sex")))
covs.passive.tableone <- print(CreateTableOne(data=pheno.covs.passive[,-1],strata="mat.passive.smoking",factorVars=c("ses","parity","sex")))

write.csv(min.pat.tableone,file=paste0(study,".patsmoking.min.pat.summary.",timepoint,".csv"))
write.csv(min.mat.tableone,file=paste0(study,".patsmoking.min.mat.summary.",timepoint,".csv"))
write.csv(covs.pat.tableone,file=paste0(study,".patsmoking.covs.pat.summary.",timepoint,".csv"))
write.csv(covs.mat.tableone,file=paste0(study,".patsmoking.covs.mat.summary.",timepoint,".csv"))
write.csv(pat.only.tableone,file=paste0(study,".patsmoking.pat.only.summary.",timepoint,".csv"))
write.csv(passive.tableone,file=paste0(study,".patsmoking.passive.summary.",timepoint,".csv"))

# Run each EWAS
ewas.res.min.pat <- ewas.function(meth, pheno.min.pat[,!colnames(pheno.min.pat) =="sex"], variable.of.interest = "pat.active.smoking")
ewas.res.min.mat <- ewas.function(meth, pheno.min.mat[,!colnames(pheno.min.mat) =="sex"], variable.of.interest = "mat.active.smoking")
ewas.res.min.mutual <- ewas.function(meth, pheno.min.mutual[,!colnames(pheno.min.mutual) =="sex"], variable.of.interest = c("pat.active.smoking","mat.active.smoking"))
ewas.res.covs.pat <- ewas.function(meth, pheno.covs.pat[,!colnames(pheno.covs.pat) =="sex"], variable.of.interest = "pat.active.smoking")
ewas.res.covs.mat <- ewas.function(meth, pheno.covs.mat[,!colnames(pheno.covs.mat) =="sex"], variable.of.interest = "mat.active.smoking")
ewas.res.covs.mutual <- ewas.function(meth, pheno.covs.mutual[,!colnames(pheno.covs.mutual) =="sex"], variable.of.interest = c("pat.active.smoking","mat.active.smoking"))
ewas.res.covs.mutual.boys.only <- ewas.function(meth, pheno.covs.mutual.boys.only[,!colnames(pheno.covs.mutual.boys.only) =="sex"], variable.of.interest = c("pat.active.smoking","mat.active.smoking"))
ewas.res.covs.mutual.girls.only <- ewas.function(meth, pheno.covs.mutual.girls.only[,!colnames(pheno.covs.mutual.girls.only) =="sex"], variable.of.interest = c("pat.active.smoking","mat.active.smoking"))
ewas.res.covs.pat.only <- ewas.function(meth, pheno.covs.pat.only[,!colnames(pheno.covs.pat.only) =="sex"], variable.of.interest = "pat.active.smoking")
ewas.res.covs.passive <- ewas.function(meth, pheno.covs.passive[,!colnames(pheno.covs.passive) =="sex"], variable.of.interest = "mat.passive.smoking")

# Save EWAS results as an Rdata file
save(list=intersect(ls(),
            c("ewas.res.min.pat",
            "ewas.res.min.mat",
            "ewas.res.min.mutual",
            "ewas.res.covs.pat",
            "ewas.res.covs.mat",
            "ewas.res.covs.mutual",
            "ewas.res.mutual.boys.only",
            "ewas.res.mutual.girls.only",
            "ewas.res.pat.only",
            "ewas.res.passive")),
     file=paste0(study,".patsmoking.ewasresults.",timepoint,".Rdata"))
