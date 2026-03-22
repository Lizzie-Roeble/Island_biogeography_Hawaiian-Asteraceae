###########################
# DAISIE Model Comparison #
# Constant                #
###########################

# 1. Load  pacakges and data

## load libraries
library(DAISIE)  
library(readxl)
library(writexl)

## update for specific run:
file_path <- "/scratch/p300503/DAISIE_Hawaii/outputs/Models_Hawaii_constant/DAISIE_ML_results_cleaned.xlsx"
output_path <- "/scratch/p300503/DAISIE_Hawaii/outputs/Models_Hawaii_constant/"

## load and format data
models_all <- read_excel(file_path)
# subset to columns needed
models_all <- models_all[, c("model", "clado", "mu", "K", "gamma", "ana", "loglik", "n_pars")]
# Ensure it's a dataframe
models_all <- as.data.frame(models_all)
# find unique models
models <- sort(unique(models_all[,'model']))
# count number of models
length(models)


# 2. Define IC functions

#LogLik- vector of logliks e.g. LogLik = c(-100.1,-102,-105.2)
# k - number of parameters e.g. k =c(4,5,3)
# n - sample size, should be the same as M e.g. n = 7000
# IC_type - for the weights function only - 'AIC', 'BIC' or 'AICc'

# AIC Akaike Information Criterion
myAIC <- function(LogLik,k){
  aic <- (2*k)-(2*LogLik)
  return(aic)
}

## AIC for small sample sizes
myAICc <- function(LogLik,k,n){
  aicc <- myAIC(LogLik,k)+((2*k*(k+1))/(n-k-1))
  return(aicc)
}

# BIC Bayesian Information Criterion
myBIC <- function(LogLik,k,n){
  bic <- -2*LogLik+k*(log(n)+log(2*pi))
  return(bic)
}

# IC weights
ICweights <- function(LogLik,k,n,IC_type){
  if(IC_type == 'AIC'){
    IC <- myAIC(LogLik,k)
  }
  if(IC_type == 'AICc'){
    IC <- myAICc(LogLik,k,n)
  }
  if(IC_type == 'BIC'){
    IC <- myBIC(LogLik,k,n)
  }
  bestmodelIC <- min(IC)
  weights <- exp(-0.5*(IC-bestmodelIC))
  weights <- weights/sum(weights)
  return(weights)
}

# 3. Selecting the Best Constant Model

## Initializes empty vectors to store best log-likelihoods and parameter counts
the_logliks <- c()
the_npars <- c()
the_clado <- c()
the_mu <- c()
the_K <- c()
the_gamma <- c()
the_ana <- c()

for (i in 1:length(models)) {
  the_model <- models[i]
  the_sub <- subset(models_all, model == the_model)
  the_row <- which(the_sub[, 'loglik'] == max(the_sub[, 'loglik']))[1]
  
  the_logliks <- append(the_logliks, the_sub[the_row, 'loglik'])
  the_npars <- append(the_npars, the_sub[the_row, 'n_pars'])
  the_clado <- append(the_clado, the_sub[the_row, 'clado'])
  the_mu <- append(the_mu, the_sub[the_row, 'mu'])
  the_K <- append(the_K, the_sub[the_row, 'K'])
  the_gamma <- append(the_gamma, the_sub[the_row, 'gamma'])
  the_ana <- append(the_ana, the_sub[the_row, 'ana'])
}

## Compute BIC and AIC for the selected 4 models
bics <- myBIC(the_logliks, the_npars, 7000)
aics <- myAIC(the_logliks, the_npars)
## Compute weights for BIC, AIC, and AICc
bic_weights <- ICweights(the_logliks, the_npars, 7000, 'BIC')
aic_weights <- ICweights(the_logliks, the_npars, 7000, 'AIC')
aicc_weights <- ICweights(the_logliks, the_npars, 7000, 'AICc')

## Identify the best constant model for each criterion
best_constant_model_BIC <- paste('The best model is: M', models[which(bic_weights == max(bic_weights))], sep = '')
best_constant_model_AIC <- paste('The best model is: M', models[which(aic_weights == max(aic_weights))], sep = '')
best_constant_model_AICc <- paste('The best model is: M', models[which(aicc_weights == max(aicc_weights))], sep = '')

## Create a comparison table
constant_model_comparison <- data.frame(
  Model = models,
  Cladogenesis = the_clado,
  Extinction = the_mu,
  Carrying_Capacity = the_K,
  Colonization = the_gamma,
  Anagenesis = the_ana,
  Num_pars = the_npars,
  Loglik = the_logliks,
  BIC = bics,
  BIC_weight = bic_weights,
  AIC = aics,
  AIC_weight = aic_weights,
  AICc_weight = aicc_weights
)

## Sort by BIC weight (highest to lowest)
constant_model_comparison[order(-constant_model_comparison$BIC_weight), ]


# 5. Review and save model comparison results

best_constant_model_BIC
best_constant_model_AIC
best_constant_model_AICc

write_xlsx(constant_model_comparison, paste0(output_path, "Model_comparison_constant.xlsx"))