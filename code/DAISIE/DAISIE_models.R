######################
# DAISIE Model Array #
######################

## Capture command-line arguments (replicate number)
args <- commandArgs(trailingOnly=TRUE)
replicate_number <- as.numeric(args[1])  # Convert argument to numeric

.libPaths('/home1/p223208/R/library')

# Load DAISIE
library(DAISIE)  

# Update paths at the start
output_dir <- "/home1/p223208/R/Lizzie"
type1_data_path_vec <- c("Hawaii_DAISIE_datalist_precise_6.15myr.Rdata",
                         "Hawaii_DAISIE_datalist_precise_5.14myr.Rdata",
                         "Hawaii_DAISIE_datalist_maxage_6.15myr.Rdata",
                         "Hawaii_DAISIE_datalist_maxage_5.14myr.Rdata")

# 1. Loading Model Parameters from File

## Read input parameter file containing model settings

model <- 1  # Model number
the_init <- replicate_number  # Seed for reproducibility
cond <- 0  # Condition parameter
res <- 200  # Resolution setting
methode <- 'odeint::runge_kutta_cash_karp54'

set.seed(the_init)  # Set random seed for reproducibility

print(sessionInfo())

# 2. Load DAISIE dataset
for(the_dataset in 1:4) {

type1_data_path <- type1_data_path_vec[the_dataset]

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

## Round initial parameter values for consistency
initparsopt <- round(initparsopt, 5)

# 5. Run DAISIE Maximum Likelihood (ML) optimization

## Run1
lik_res <- DAISIE_ML(datalist= datalist, initparsopt=initparsopt, 
                     idparsnoshift=idparsnoshift, idparsopt=idparsopt, 
                     parsfix=parsfix, idparsfix=idparsfix, 
                     ddmodel=ddmodel, methode=methode, cond=cond, res=res, num_cycles = 5)
## Save first run results
lik_res<-as.matrix(lik_res)
file_name <- file.path(output_dir, paste0(the_dataset,"_results_R1_", replicate_number, ".txt"))
to.write <- c(replicate_number, model, the_init, island_age, M, type, P_type2, res, methode, as.matrix(lik_res), initparsopt, "\n")
cat(to.write, file=file_name, append=TRUE)
}