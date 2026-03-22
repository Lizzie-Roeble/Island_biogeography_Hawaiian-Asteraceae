######################
# DAISIE Model Array #
# Constant           #
######################

# Load DAISIE
library(DAISIE)  

# Update paths at the start
param_file <- "model_parameters_constant.txt"
type1_data_path <- "../data/Hawaii_DAISIE_datalist.Rdata"
output_dir <- "../outputs/Models_Hawaii_constant/Individual_models/"

# 1. Loading Model Parameters from File

## Capture command-line arguments (replicate number)
args = commandArgs(trailingOnly=TRUE)
replicate_number = as.numeric(args[[1]])  # Convert argument to numeric

## Read input parameter file containing model settings
a_vals <- read.table(param_file, header=F)

the_dataset <- a_vals[replicate_number,1]  # Dataset index
model <- a_vals[replicate_number,2]  # Model number
the_init <- a_vals[replicate_number,3]  # Seed for reproducibility
cond <- a_vals[replicate_number,4]  # Condition parameter
res <- a_vals[replicate_number,5]  # Resolution setting
methode <- a_vals[replicate_number,6]  # Optimization method

set.seed(the_init)  # Set random seed for reproducibility

print(sessionInfo())

# 2. Load DAISIE dataset

## load DAISIE type1 (only need type1 for constant models)
type1_data <- load(type1_data_path)
type1_datalist <- get(type1_data)

## Extract island age and mainland pool
island_age <- type1_datalist[[1]]$island_age
M <- type1_datalist[[1]]$not_present + length(type1_datalist) - 1

# Use type1 datalist for all models since we're only running constant models
datalist <- type1_datalist
P_type2 <- NA  # Keep NA for type1 models

# 3. Generate starting parameter values 

## Generate random initial parameter values for optimization
r_lamc <- runif(1, min=0, max=2)   # Cladogenesis speciation rate
r_mu <- runif(1, min=0, max=2)     # Extinction rate
r_K <- runif(1, min=35, max=200)  # Carrying capacity - set higher than the largest clade
r_gam <- runif(1, min=0, max=0.0005)  # Immigration rate
r_ana <- runif(1, min=0.1, max=5)  # Anagenesis speciation rate

# 4. Define DAISIE Models

## M1: CS - DD
if(model == "1") {
  ddmodel = 11  # Density-dependent model
  idparsopt <- 1:5  # Parameters to optimize
  parsfix <- NULL  # No fixed parameters
  idparsfix <- NULL
  idparsnoshift = 6:10  # Parameters unaffected by shifts
  initparsopt <- c(r_lamc, r_mu, r_K, r_gam, r_ana)
  type <- 'type1'  # Uses type1 dataset
}
## M2: CS - DI
if(model=="2") {
  ddmodel= 0
  idparsopt<-c(1,2,4,5)
  parsfix<-Inf
  idparsfix<-3
  idparsnoshift = 6:10
  initparsopt<-c(r_lamc,r_mu,r_gam,r_ana)
  type <- 'type1'
}
## M3: CS - DD no anagenesis
if(model=="3") {
  ddmodel= 11
  idparsopt<-c(1,2,3,4)
  parsfix<-0
  idparsfix<-5
  idparsnoshift = 6:10
  initparsopt<-c(r_lamc,r_mu,r_K,r_gam)
  type <- 'type1'
}
## M4: CS - DI no anagenesis
if(model=="4") {
  ddmodel= 0
  idparsopt<-c(1,2,4)
  parsfix<-c(Inf,0)
  idparsfix<-c(3,5)
  idparsnoshift = 6:10
  initparsopt<-c(r_lamc,r_mu,r_gam)
  type <- 'type1'
}

## Round initial parameter values for consistency
initparsopt <- round(initparsopt, 5)

# 5. Run DAISIE Maximum Likelihood (ML) optimization

## Run1
lik_res <- DAISIE_ML(datalist= datalist, initparsopt=initparsopt, 
                     idparsnoshift=idparsnoshift, idparsopt=idparsopt, 
                     parsfix=parsfix, idparsfix=idparsfix, 
                     ddmodel=ddmodel, methode=methode, cond=cond, res=res)
## Save first run results
lik_res<-as.matrix(lik_res)
file_name <- file.path(output_dir, paste0("results_R1_", replicate_number, ".txt"))
to.write <- c(replicate_number, the_dataset, model, the_init, island_age, M, type, P_type2, res, methode, "R1", as.matrix(lik_res), initparsopt, "\n")
cat(to.write, file=file_name, append=TRUE)

## Run2
## Run second optimization with updated parameters if Run1 succeeds
if(!is.na(lik_res[idparsopt][1])) { 
  initparsopt <- lik_res[idparsopt]  # Update parameters based on first run
}

lik_res <- DAISIE_ML(datalist= datalist, initparsopt=initparsopt, 
                     idparsnoshift=idparsnoshift, idparsopt=idparsopt, 
                     parsfix=parsfix, idparsfix=idparsfix, 
                     ddmodel=ddmodel, methode=methode, cond=cond, res=res)
## Save second run results
lik_res<-as.matrix(lik_res)
file_name <- file.path(output_dir, paste0("results_R2_", replicate_number, ".txt"))
to.write <- c(replicate_number, the_dataset, model, the_init, island_age, M, type, P_type2, res, methode, "R2", as.matrix(lik_res), initparsopt, "\n")
cat(to.write, file=file_name, append=TRUE)