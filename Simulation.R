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
sim_conditions$prior<- as.character(sim_conditions$prior) #TO help with some of the output, not treat it as a factor anymore

#test <- simulate_study(pos=1, cond = sim_conditions)
#test <- apply(1:nrow(sim_conditions), function(x) simulate_study(pos = x, cond = sim_conditions))

if(.Platform$OS.type == "unix") {
  cl = makeCluster(numcores, type = "FORK")
} else {
  cl <- makeCluster(numcores, type = "PSOCK")
  # If you're using libraries in each node, you can call this.
}

#SIM
nworkers <- detectCores() # number of cores to use
cl <- makePSOCKcluster(nworkers) # create cluster
clusterCall(cl, function() library(cmdstanr))
clusterCall(cl, function() library(bayesplot))
clusterCall(cl, function() { source("Dat_gen.R") })
simulation <- clusterApplyLB(cl, 1:nrow(sim_conditions), simulate_study, cond=sim_conditions) # run simulation
stopCluster(cl) # shut down the nodes
