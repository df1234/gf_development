################################################################################

# This script reproduces key analyses for the paper

# Fuhrmann, D., Simpson-Kent, I.L., Bathelt, J., the CALM team & Kievit, R. A.
#    (2019). A hierarchical watershed model of fluid intelligence in childhood and adolescence.
#    Cerebral Cortex, bhz091, https://doi.org/10.1093/cercor/bhz091

# The script uses the covariance matrix fromt he CALM sample. 
# Estimates and test statistics will differ from values reported in the paper due to missingness. 
# To reproduce results exactly, please request raw data from http://fcon_1000.projects.nitrc.org/indi/enhanced/data.html 
# for NKI-RS. Raw data for CALM has not yet been released due to ethics restrictions. Please check
# http://calm.mrc-cbu.cam.ac.uk/ for updates.

# Script author: Delia Fuhrmann, August 2018, using RStudio 1.0.143 and R 
# version 3.4.3 (Kite-Eating Tree)

################################################################################

rm(list=ls()) # clear all

### Load lavaan
# install.packages("lavaan") # this only needs to be done once
library(lavaan)

### Load covariance matrix
cov.matrix.CALM = read.csv("CovMat_CALM.csv", row.names = 1) # CALM data

################################################################################
### Measurement model of the watershed model for both samples

model_ThreeFactor_CALM<-
  '
# latent variables: gf (fluid intelligence), WM (working memory) and PS (processing speed)
gf =~ gf_MatrixReasoning
WM =~ WM_DigitRecall + WM_DotMatrix + WM_BackwardDigit + WM_MrX
PS =~ PS_PhAB_il + PS_Teach_il + PS_DKEFS_il     
'
fit_ThreeFactor_CALM <- cfa(model_ThreeFactor_CALM, sample.cov=as.matrix(cov.matrix.CALM), 
                            sample.nobs=551, std.lv=T)
summary(fit_ThreeFactor_CALM, standardized = T)

################################################################################
### MIMIC models (modelling cognitive levels of the watershed model)

model_MIMIC_CALM <-
  '
# latent variables
gf =~ gf_MatrixReasoning
WM =~ WM_DigitRecall + WM_DotMatrix + WM_BackwardDigit + WM_MrX
PS =~ PS_PhAB_il + PS_Teach_il + PS_DKEFS_il

# regressions
gf ~ PS + WM
'
fit_MIMIC_CALM <- sem(model_MIMIC_CALM, sample.cov=as.matrix(cov.matrix.CALM), 
                          sample.nobs=551, std.lv=T)
summary(fit_MIMIC_CALM, standardized = T)

################################################################################
### Full watershed models including neuroimaging data (FA in 10 different white matter tracts)

model_watershed_CALM <-
  '
# latent variables
gf =~ gf_MatrixReasoning
WM =~ WM_DigitRecall + WM_DotMatrix + WM_BackwardDigit + WM_MrX
PS =~ PS_PhAB_il + PS_Teach_il + PS_DKEFS_il

# regressions
gf ~ PS + WM
PS ~ fa_UF + fa_SLF + fa_IFOF + fa_ATR + fa_CST + fa_FMaj + fa_FMin + fa_CG + fa_CH + fa_ILF
WM ~ fa_UF + fa_SLF + fa_IFOF + fa_ATR + fa_CST + fa_FMaj + fa_FMin + fa_CG + fa_CH + fa_ILF

# covariances
PS ~~ WM
'
fit_watershed_CALM <- sem(model_watershed_CALM, sample.cov=as.matrix(cov.matrix.CALM), 
                          sample.nobs=551, std.lv=T, fixed.x=F)
summary(fit_watershed_CALM, standardized = T)
