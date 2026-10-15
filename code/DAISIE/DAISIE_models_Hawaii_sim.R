DAISIE_models_Hawaii_sim <- function(datasetno, row_index) {
  # 1. Loading Data and Model Parameters from File and setting initial parameters

  if(datasetno == 1) {
     load('DAISIE_Hawaii_sim_precise_6.15myr.Rdata')
     results <- read.table('results_precise_6.15myr.txt', header = FALSE)
     outputfile <- "DAISIE_Hawaii_sim_precise_6.15myr.txt"
  }
  if(datasetno == 2) {
     load('DAISIE_Hawaii_sim_maxage_6.15myr.Rdata')
     results <- read.table('results_maxage_6.15myr.txt', header = FALSE)
     outputfile <- "DAISIE_Hawaii_sim_maxage_6.15myr.txt"
  }
  if(datasetno == 3) {
     load('DAISIE_Hawaii_sim_precise_5.14myr.Rdata')
     results <- read.table('results_precise_5.14myr.txt', header = FALSE)
     outputfile <- "DAISIE_Hawaii_sim_precise_5.14myr.txt"
  }
  if(datasetno == 4) {
     load('DAISIE_Hawaii_sim_maxage_5.14myr.Rdata')
   results <- read.table('results_maxage_5.14myr.txt', header = FALSE)
     outputfile <- "DAISIE_Hawaii_sim_maxage_5.14myr.txt"
  }
  results <- results[order(-results[, 15]), ]
  initparsopt <- as.numeric(results[1,10:14])

  print(sessionInfo())
  cat(paste0('Dataset: ',datasetno, '\n'))
  cat(paste0('Simulation: ',row_index, '\n'))

  # 2. Load DAISIE simulation
  datalist <- sim_results[[row_index]]

  ## Extract island age and mainland pool
  island_age <- datalist[[1]]$island_age
  M <- datalist[[1]]$not_present + length(datalist) - 1
  res <- 200
  cond <- 1
  methode <- 'odeint::runge_kutta_cash_karp54'
  ddmodel = 11  # Density-dependent model

  stac4 <- FALSE
  for(j in 2:length(datalist)) {
     if(datalist[[j]]$stac == 4) stac4 <- TRUE
  }

  if(stac4) {
  idparsopt <- 1:5  # Parameters to optimize
  parsfix <- NULL  # No fixed parameters
  idparsfix <- NULL
  initpars <- initparsopt

  #onlystac2 <- TRUE
  #onlyclades <- TRUE
  #for(j in 2:length(datalist)) {
  #  if(datalist[[j]]$stac != 2) onlystac2 <- FALSE
  #  if(datalist[[j]]$stac == 2 & length(datalist[[j]]$branching_times) == 2) onlyclades <- FALSE 
  #}
  #if(onlystac2 == FALSE | (onlystac2 == TRUE & onlyclades == TRUE)) {
  #   idparsopt <- 1:5  # Parameters to optimize
  #   parsfix <- NULL  # No fixed parameters
  #   idparsfix <- NULL
  #   initpars <- initparsopt
  #} else {
  #   initparsopt <- initparsopt[1:4]
  #   idparsopt <- 1:4
  #   parsfix <- 500 #anagenesis to 1000
  #   idparsfix <- 5
  #   initpars <- c(initparsopt,parsfix)
  #}
  idparsnoshift = 6:10  # Parameters unaffected by shifts
  type <- 'type1'  # Uses type1 dataset

  # 5. Run DAISIE Maximum Likelihood (ML) optimization
  lik_res <- DAISIE_ML(datalist = datalist,
                       initparsopt = initparsopt,
                       idparsnoshift = idparsnoshift,
                       idparsopt = idparsopt,
                       parsfix = parsfix,
                       idparsfix = idparsfix,
                       verbose = 1,
                       ddmodel = ddmodel,
                       methode = methode,
                       cond = cond,
                       res = res,
                       CS_version = list(model = 1, function_to_optimize = "DAISIE"),
                       equal_extinction = TRUE,
                       num_cycles = 5)
  ## Save results
  lik_res<-as.matrix(lik_res)
  file_name <- file.path(outputfile)
  to.write <- c(row_index, island_age, M, type, res, cond, methode, lik_res, initpars, "\n")
  cat(to.write, file = file_name, append = TRUE)
}
}

