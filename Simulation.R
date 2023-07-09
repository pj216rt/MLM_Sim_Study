source("functions.R")
#CMDSTANR
library(cmdstanr)
library(posterior)
library(bayesplot)
library(parallel) # to run parallel

color_scheme_set("brightblue")
check_cmdstan_toolchain()
cmdstan_path()

#install_cmdstan() IF NEEDED

#Set up grid
num_con <- 3   #Number of data conditions
priors <- c("uninform", "ridge", "studentst", "lasso") #List the priors to check
cond <- 1:num_con
#Create grid to search
sim_conditions <- expand.grid(prior = priors, condition=cond)
sim_conditions$prior<- as.character(sim_conditions$prior) #To help with some of the output, not treat it as a factor anymore

#SIM
#Because we are on a Windows machine for now, we need to use the socket method
#of paralleization.  This launches a new version of R on each core being used.  
nworkers <- detectCores() # number of cores to use
cl <- makePSOCKcluster(nworkers) # create cluster
# clusterCall(cl, function() library(cmdstanr))
# clusterCall(cl, function() library(bayesplot))
# clusterCall(cl, function() { source("Dat_gen.R") })
clusterEvalQ(cl, {
  library(cmdstanr)
  library(bayesplot)
  library(posterior)
  })
clusterExport(cl, list("out", "out1", "out2", "out3", "out4", "out5"), envir = environment())
simulation <- clusterApplyLB(cl, 1:nrow(sim_conditions), simulate_study, cond=sim_conditions) # run simulation
stopCluster(cl) # shut down the nodes
