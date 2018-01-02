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
#    Necessary variable names are: "pat.active.smoking", "mat.active.smoking","mat.passive.smoking","sex","pat.ses","mat.ses","pat.age","mat.age","pat.bmi","mat.bmi","parity"
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
traits.and.covariates <- c("pat.active.smoking", "mat.active.smoking","mat.passive.smoking","sex","pat.ses","mat.ses","pat.age","mat.age","pat.bmi","mat.bmi","parity")
covariates <- c("pat.age", "pat.bmi", "pat.ses", "mat.age", "mat.bmi", "mat.ses" , "parity", cell.names)

# Load and check phenotype data
pheno <- read.csv("EWAS/pat_smoke/phenofile.alspac.csv",header=TRUE,stringsAsFactors=FALSE) #change filename/location to point to your phenotype file

for(i in 1:length(c("sample.id",traits.and.covariates,cell.names))) {
print(ifelse(c("sample.id",traits.and.covariates,cell.names)[i] %in% colnames(pheno)==FALSE,
  paste("CAUTION: the variable called",c("sample.id",traits.and.covariates,cell.names)[i],"is missing from pheno"),
  paste("variable called",c("sample.id",traits.and.covariates,cell.names)[i],"is present in pheno")))
}

table(pheno$mat.passive.smoking,(pheno$mat.active.smoking|pheno$pat.active.smoking)) #checking number of maternal passive smokers with no active smoking for either mother or father
pheno$mat.passive.smoking[which(pheno$pat.active.smoking==1 |pheno$mat.active.smoking==1)] <- NA #setting maternal passive smoking to missing if mother or father is an active smoker
table(pheno$bio.dad) #checking number of partners that are not biological fathers
pheno <- pheno[which(pheno$bio.dad==1),] #removing partners that are not biological fathers

# Load methylation data and perform QC (e.g. filter probes with high detection P-values)
load("/panfs/panasas01/dedicated-mrcieu/studies/latest/alspac/epigenetic/methylation/450k/aries/released/2016-05-03/data/betas/data.Robj") #change filename/location to point to your methylation data (meth)
# Perform any cohort-specific QC at this point (e.g. you might want to remove probes with high detection p-values)

# IQR*3 method to remove outliers (if this has not already been applied to your data)
meth <- IQR.removal(meth)

# Generate surrogate variables for technical batch and merge with pheno data to prepare pheno data frames for each EWAS
pheno.mutual <- SVA.generate(meth, pheno, variable.of.interest = c("pat.active.smoking","mat.active.smoking"), model.covariates = c(covariates,"sex"),n.sv=20)
pheno.minimal.pat <-pheno.mutual[,-which(colnames(pheno.mutual) %in% "mat.active.smoking")]
pheno.minimal.mat <-pheno.mutual[,-which(colnames(pheno.mutual) %in% "pat.active.smoking")]
pheno.passive <- SVA.generate(meth, pheno, variable.of.interest = "mat.passive.smoking", model.covariates = c(covariates,"sex"),n.sv=20)
pheno.pat.only <-  pheno.mutual[which(pheno.mutual$mat.active.smoking == 0),]
pheno.pat.only <-  pheno.pat.only[, -which(colnames(pheno.pat.only) %in% c("mat.active.smoking","pat.active.smoking"))]
pheno.mutual.boys.only <- pheno.mutual[which(pheno.mutual$sex == 0),]
pheno.mutual.girls.only <- pheno.mutual[which(pheno.mutual$sex == 1),]

# Summarise pheno data and save summaries as .csv files
mutual.pat.tableone <- print(CreateTableOne(data=pheno.mutual[,-1],strata="pat.active.smoking",factorVars=c("pat.active.smoking","mat.active.smoking","mat.ses","pat.ses","parity","sex")))
mutual.mat.tableone <- print(CreateTableOne(data=pheno.mutual[,-1],strata="mat.active.smoking",factorVars=c("pat.active.smoking","mat.active.smoking","mat.ses","pat.ses","parity","sex")))
pat.only.tableone <- print(CreateTableOne(data=pheno.pat.only[,-1],strata="pat.active.smoking",factorVars=c("pat.active.smoking","mat.ses","pat.ses","parity","sex")))
passive.tableone <- print(CreateTableOne(data=pheno.passive[,-1],strata="mat.passive.smoking",factorVars=c("mat.ses","pat.ses","parity","sex")))

write.csv(mutual.pat.tableone,file=paste0(study,".patsmoking.mutual.pat.summary.",timepoint,".csv"))
write.csv(mutual.mat.tableone,file=paste0(study,".patsmoking.mutual.mat.summary.",timepoint,".csv"))
write.csv(pat.only.tableone,file=paste0(study,".patsmoking.pat.only.summary.",timepoint,".csv"))
write.csv(passive.tableone,file=paste0(study,".patsmoking.passive.summary.",timepoint,".csv"))

# Run each EWAS
ewas.res.minimal.pat <- ewas.function(meth, pheno.minimal.pat[,!colnames(pheno.minimal.pat) =="sex"], variable.of.interest = "pat.active.smoking")
ewas.res.minimal.mat <- ewas.function(meth, pheno.minimal.mat[,!colnames(pheno.minimal.mat) =="sex"], variable.of.interest = "mat.active.smoking")
ewas.res.mutual <- ewas.function(meth, pheno.mutual[,!colnames(pheno.mutual) =="sex"], variable.of.interest = c("pat.active.smoking","mat.active.smoking"))
ewas.res.mutual.boys.only <- ewas.function(meth, pheno.mutual.boys.only[,!colnames(pheno.mutual.boys.only) =="sex"], variable.of.interest = c("pat.active.smoking","mat.active.smoking"))
ewas.res.mutual.girls.only <- ewas.function(meth, pheno.mutual.girls.only[,!colnames(pheno.mutual.girls.only) =="sex"], variable.of.interest = c("pat.active.smoking","mat.active.smoking"))
ewas.res.pat.only <- ewas.function(meth, pheno.pat.only[,!colnames(pheno.pat.only) =="sex"], variable.of.interest = "pat.active.smoking")
ewas.res.passive <- ewas.function(meth, pheno.passive[,!colnames(pheno.passive) =="sex"], variable.of.interest = "mat.passive.smoking")

# Save EWAS results as an Rdata file
save(list=intersect(ls(),
            c("ewas.res.minimal.pat",
            "ewas.res.minimal.mat",
            "ewas.res.mutual",
            "ewas.res.mutual.boys.only",
            "ewas.res.mutual.girls.only",
            "ewas.res.pat.only",
            "ewas.res.passive")),
     file=paste0(study,".patsmoking.ewasresults.",timepoint,".Rdata"))
