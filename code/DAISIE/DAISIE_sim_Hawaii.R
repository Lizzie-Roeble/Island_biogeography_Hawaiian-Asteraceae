args <- commandArgs(trailingOnly = TRUE)
home_dir <- substring(getwd(),1,14)
lib_dir <- paste(home_dir,'/R/library',sep = '')
.libPaths(lib_dir)
library(DDD)
library(DAISIE)
set.seed(42)
if(args == 1) {
   load('Hawaii_DAISIE_datalist_precise_6.15myr.Rdata');
   datalist <- Hawaii_DAISIE_datalist_precise;
   results <- read.table('results_precise_6.15myr.txt', header = FALSE)
   results <- results[order(-results[, 15]), ]
   pars <- as.numeric(results[1,10:14])
   outputfile <- 'DAISIE_Hawaii_sim_precise_6.15myr.Rdata'
}
if(args == 2) {
   load('Hawaii_DAISIE_datalist_maxage_6.15myr.Rdata');
   datalist <- Hawaii_DAISIE_datalist;
   results <- read.table('results_maxage_6.15myr.txt', header = FALSE)
   results <- results[order(-results[, 15]), ]
   pars <- as.numeric(results[1,10:14])
   outputfile <- 'DAISIE_Hawaii_sim_maxage_6.15myr.Rdata'
}
if(args == 3) {
   load('Hawaii_DAISIE_datalist_precise_5.14myr.Rdata');
   datalist <- Hawaii_DAISIE_datalist_precise;
   results <- read.table('results_precise_5.14myr.txt', header = FALSE)
   results <- results[order(-results[, 15]), ]
   pars <- as.numeric(results[1,10:14])
    outputfile <- 'DAISIE_Hawaii_sim_precise_5.14myr.Rdata'
}
if(args == 4) {
   load('Hawaii_DAISIE_datalist_maxage_5.14myr.Rdata');
   datalist <- Hawaii_DAISIE_datalist;
   results <- read.table('results_maxage_5.14myr.txt', header = FALSE)
   results <- results[order(-results[, 15]), ]
   pars <- as.numeric(results[1,10:14])
   outputfile <- 'DAISIE_Hawaii_sim_maxage_5.14myr.Rdata'
}

print(args)
print(pars)

sim_results <- DAISIE_sim(time = datalist[[1]]$island_age,
                          M = 7000, #datalist[[1]]$not_present + length(datalist) - 1,
                          pars = pars,
                          replicates = 10000,
                          plot_sims = FALSE,
                          files_to_write = 0)
save(sim_results, file = outputfile)

#old
#load('Hawaii_DAISIE_datalist.Rdata'); datalist <- Hawaii_DAISIE_datalist; pars <- c(1.55906359,0.818905529,87.70998181,0.000554282,10.69543004); outputfile <- 'DAISIE_Hawaii_sim_maxage2.Rdata'
#load('Hawaii_DAISIE_datalist_precise.Rdata'); datalist <- Hawaii_DAISIE_datalist_precise; pars <- c(1.197907, 0.5638870, 158.369260903136, 0.0004939022, 2.058786); outputfile <- 'DAISIE_Hawaii_sim_preciseage2.Rdata'
