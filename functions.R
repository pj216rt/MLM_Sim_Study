#Functions for simulation
require(MBESS)
require(mnormt)
require(lattice)
require(nlme)
require(ggplot2)
require(reshape2)
require(rstan)
require(dplyr)
require(tidyverse)
require(lme4)
#require(brms)

genData <- function(nSubjs = 100, sdErr = 1, 
                    # intercept and slope fixed effects
                    coef1 = c(4, 3),
                    # types of level 2 covariates
                    level2Binary = c(.2, 0.7),
                    level2Continuous = list(c(mu = 0, sd = 1),
                                            c(mu = 5, sd = 1)),
                    # corr between errors on subject-specific int and slope
                    corrRE = 0.20,
                    # sd of errors on subject-specific int and slope
                    sdRE = c(1, 1),
                    # for each predictor in level 2, (int2, slope2) 
                    # specify effect on level 2 intercept and level 2 slope
                    coef2Binary = list(c(int2 = 2.0, slope2 = 1.0),
                                       c(int2 = 1.0, slope2 = 3.0)),
                    coef2Continuous = list(c(int2 = 1.0, slope2 = 1.5),
                                           c(int2 = 0.2, slope2 = 3.0))
){
  # subject IDs
  ids <- 1:nSubjs
  
  # number of observation times for each subject
  ntimes <- rpois(n = nSubjs, lambda = 3) + 1
  # subject observation times (stacked)
  subjTimes <- sequence(ntimes, from = 0, by = 1)
  
  # reshape ids for long format
  ids <- rep(ids, ntimes)
  
  # generate level 2 predictor variables
  # binary predictors
  nBinaryVars <- length(level2Binary)
  B <- rbinom(nSubjs*nBinaryVars, size = 1, 
              prob = rep(level2Binary, each = nSubjs))
  Bmat <- matrix(B, ncol = nBinaryVars)
  # continuous predictors
  nContVars <- length(level2Continuous)
  means <- sapply(level2Continuous, function(x) x["mu"])
  sds <- sapply(level2Continuous, function(x) x["sd"])
  Z <- rnorm(nSubjs*nContVars, 
             mean = rep(means, each = nSubjs),
             sd = rep(sds, each = nSubjs))
  Zmat <- matrix(Z, ncol = nContVars)
  # all level 2 predictors
  Xmat <- cbind(Bmat, Zmat)
  colnames(Xmat) <- paste0("X", 1:(nBinaryVars + nContVars))
  
  # var-cov of the random effects at level 2
  covRE <- cor2cov(matrix(c(1, corrRE, corrRE, 1), nrow = 2), sdRE)
  cat("Population corr matrix of level 2 errors is:\n")
  print(cov2cor(covRE))
  
  # random errors at level 2: columns are Intercept, Slope
  REs <- rmnorm(nSubjs, mean = rep(0, 2), varcov = covRE)
  cat("Sample corr matrix of level 2 errors is:\n")
  print(cor(REs), digits = 3)
  
  # subject-specific (level 2) intercept and slope
  allcoefs <- matrix(c(unlist(coef2Binary), unlist(coef2Continuous)), nrow = 2)
  #print(allcoef)
  rownames(allcoefs) <- c("int2", "slope2")
  coef2int <- allcoefs[1,]
  coef2slope <- allcoefs[2,]
  REintslope <- Xmat %*% cbind(coef2int, coef2slope)
  REcoefs <- REintslope + REs
  # reshape to allow for longitudinal data
  REcoefs <- matrix(rep(REcoefs, rep(ntimes, 2)), ncol = 2)
  
  # iid error for the response
  respError <- rnorm(sum(ntimes), 0, sd = sdErr)
  
  Y <- (coef1[1] + REcoefs[,1]) +
    (coef1[2] + REcoefs[,2]) * subjTimes + respError
  
  # construct long format data frame
  covars <- matrix(rep(Xmat, rep(ntimes, nBinaryVars+nContVars)), ncol = ncol(Xmat))
  colnames(covars) <- colnames(Xmat)
  ans <- cbind(id = ids, time = subjTimes, covars, Y)
  ans <- data.frame(ans)
  
  attributes(ans) <- c(attributes(ans),
                       list(coef1 = coef1, coef2 = allcoefs,
                            REcoefs = REcoefs,
                            sdErr = sdErr, corrRE = corrRE, sdRE = sdRE
                       )
  )
  ans
}

###Want to generate data where every subject is observed a constant number of times
genData_balanced <- function(nSubjs = 100, num_obs = 5, sdErr = 10, 
                             # intercept and slope fixed effects
                             coef1 = c(4, 3),
                             # types of level 2 covariates
                             level2Binary = c(.2, 0.7),
                             level2Continuous = list(c(mu = 0, sd = 1),
                                                     c(mu = 5, sd = 1)),
                             # corr between errors on subject-specific int and slope
                             corrRE = 0.20,
                             # sd of errors on subject-specific int and slope
                             sdRE = c(1, 1),
                             # for each predictor in level 2, (int2, slope2) 
                             # specify effect on level 2 intercept and level 2 slope
                             coef2Binary = list(c(int2 = 2.0, slope2 = 1.0),
                                                c(int2 = 1.0, slope2 = 3.0)),
                             coef2Continuous = list(c(int2 = 1.0, slope2 = 1.5),
                                                    c(int2 = 0.2, slope2 = 3.0))
){
  # subject IDs
  ids <- 1:nSubjs
  
  # number of observation times for each subject
  ntimes <- replicate(n=nSubjs, num_obs)
  # subject observation times (stacked)
  subjTimes <- sequence(ntimes, from = 0, by = 1)
  
  # reshape ids for long format
  ids <- rep(ids, ntimes)
  
  # generate level 2 predictor variables
  # binary predictors
  nBinaryVars <- length(level2Binary)
  B <- rbinom(nSubjs*nBinaryVars, size = 1, 
              prob = rep(level2Binary, each = nSubjs))
  Bmat <- matrix(B, ncol = nBinaryVars)
  # continuous predictors
  nContVars <- length(level2Continuous)
  means <- sapply(level2Continuous, function(x) x["mu"])
  sds <- sapply(level2Continuous, function(x) x["sd"])
  Z <- rnorm(nSubjs*nContVars, 
             mean = rep(means, each = nSubjs),
             sd = rep(sds, each = nSubjs))
  Zmat <- matrix(Z, ncol = nContVars)
  # all level 2 predictors
  Xmat <- cbind(Bmat, Zmat)
  colnames(Xmat) <- paste0("X", 1:(nBinaryVars + nContVars))
  
  # var-cov of the random effects at level 2
  covRE <- cor2cov(matrix(c(1, corrRE, corrRE, 1), nrow = 2), sdRE)
  cat("Population corr matrix of level 2 errors is:\n")
  print(cov2cor(covRE))
  
  # random errors at level 2: columns are Intercept, Slope
  REs <- rmnorm(nSubjs, mean = rep(0, 2), varcov = covRE)
  cat("Sample corr matrix of level 2 errors is:\n")
  print(cor(REs), digits = 3)
  
  # subject-specific (level 2) intercept and slope
  allcoefs <- matrix(c(unlist(coef2Binary), unlist(coef2Continuous)), nrow = 2)
  rownames(allcoefs) <- c("int2", "slope2")
  coef2int <- allcoefs[1,]
  coef2slope <- allcoefs[2,]
  REintslope <- Xmat %*% cbind(coef2int, coef2slope)
  REcoefs <- REintslope + REs
  # reshape to allow for longitudinal data
  REcoefs <- matrix(rep(REcoefs, rep(ntimes, 2)), ncol = 2)
  
  # iid error for the response
  respError <- rnorm(sum(ntimes), 0, sd = sdErr)
  
  Y <- (coef1[1] + REcoefs[,1]) +
    (coef1[2] + REcoefs[,2]) * subjTimes + respError
  
  # construct long format data frame
  covars <- matrix(rep(Xmat, rep(ntimes, nBinaryVars+nContVars)), ncol = ncol(Xmat))
  colnames(covars) <- colnames(Xmat)
  ans <- cbind(id = ids, time = subjTimes, covars, Y)
  ans <- data.frame(ans)
  
  attributes(ans) <- c(attributes(ans),
                       list(coef1 = coef1, coef2 = allcoefs,
                            REcoefs = REcoefs,
                            sdErr = sdErr, corrRE = corrRE, sdRE = sdRE
                       )
  )
  ans
}

#Balanced datasets generation, lots
gen_balanced_datasets <- function(nreps = 10, nSubjs = 100, num_obs = 5, sdErr = 10, 
                                  # intercept and slope fixed effects
                                  coef1 = c(4, 3),
                                  # types of level 2 covariates
                                  level2Binary = c(.2, 0.7),
                                  level2Continuous = list(c(mu = 0, sd = 1),
                                                          c(mu = 5, sd = 1)),
                                  # corr between errors on subject-specific int and slope
                                  corrRE = 0.20,
                                  # sd of errors on subject-specific int and slope
                                  sdRE = c(1, 1),
                                  # for each predictor in level 2, (int2, slope2) 
                                  # specify effect on level 2 intercept and level 2 slope
                                  coef2Binary = list(c(int2 = 2.0, slope2 = 1.0),
                                                     c(int2 = 1.0, slope2 = 3.0)),
                                  coef2Continuous = list(c(int2 = 1.0, slope2 = 1.5),
                                                         c(int2 = 0.2, slope2 = 3.0))){
  
  #generate nreps datasets
  simdata <- list()
  for(i in 1:nreps){
    simdata[[i]] <- genData_balanced(nSubjs, num_obs, sdErr, coef1, level2Binary, 
                                     level2Continuous, corrRE, sdRE,coef2Binary,coef2Continuous)
  }
  return(simdata)
}


#Function to automatically select the continuous variable columns, and scale them
scale_cont_data <- function(datasets){
  for(i in seq_along(datasets)){
    dat <- datasets[[i]]
    for(j in seq_along(dat)){
      #Basically, don't select id, time, or the response variable
      temp <- subset(dat[[j]], select = -c(id, time, Y))
      temp1 <- ncol(temp)
      temp2 <- ncol(temp)/2
      temp3 <- temp[, (temp2+1):temp1] #Taking advantage of the way that we have structured the data
      
      #get names of this
      temp4 <- colnames(temp3)
      #Scale ONLY the continuous level 2 variables
      dat[[j]] <- dat[[j]] %>% mutate_at(c(temp4), ~(scale(.) %>% as.vector))
    }
    datasets[[i]] <- dat
  }
  return(datasets)
}

#Split data into test and train sets.
#Split data into test and train
tt_split <- function(datasets, var_to_select = id ,percent_train = 0.80){
  train_dat <- list()
  test_dat <- list()
  for(i in seq_along(datasets)){
    data <- datasets[[i]]
    
    #Essentially create a new column that specifies what group each participant is in
    groups <- data %>% dplyr::select(id) %>% distinct(id) %>% rowwise() %>%
      mutate(group = sample(
        c("train", "test"), 1,
        replace = TRUE,
        prob = c(percent_train, (1-percent_train)) #weights for each group
      ))
    
    #left joined groups to the data
    data <- data %>% left_join(groups)
    
    #Split data into test and train
    #Train data
    train_dat[[i]] <- filter(data, group == "train")
    train_dat[[i]]$id.new <- match(train_dat[[i]]$id, unique(train_dat[[i]]$id))
    
    #Test data
    test_dat[[i]] <- filter(data, group == "test")
    test_dat[[i]]$id.new <- match(test_dat[[i]]$id, unique(test_dat[[i]]$id))
  }
  out <- list(Training = train_dat, Testing = test_dat)
  return(out)
}

###If for some reason, this code doesn't work, we can try the following:

# #Now need to split into test and train sets
# #Split into test and train
# split.sim1 <- list()
# for(i in seq_along(dataset)){
#   temp <- dataset[[i]]
#   temp <- tt_split(datasets = temp, percent_train = 0.80)
#   split.sim1[[i]] <- temp
# }

#Perhaps there is some nicer way to do this, but this seems to work
extract_lev2 <- function(dat, cols_to_drop = c("time", "subject", "group")){
  f <- dat[!duplicated(dat$id), ] #Select the first row of each id "group"
  df <- f[, ! names(f) %in% cols_to_drop, drop = F]
  return(df)
}

#Extract level 2 variables across multiple datasets
multiple_extract_lev2_var <- function(datasets, 
                                      cols_drop = c("id", "time", "Y", 
                                                    "group", "id.new")){
  output <- list()
  for(i in seq_along(datasets)){
    print("hello")
    lev2_vars <- extract_lev2(datasets[[i]], cols_to_drop=cols_drop)
    output[[i]] <- lev2_vars
  }
  return(output)
}

#Simple stan data loop construction
stan_data_loop <- function(training_datasets, testing_datasets){
  level2_vars <- multiple_extract_lev2_var(training_datasets)
  #Creating formula to plug into the test data model matrix later
  names_of_lev2_vars <- colnames(level2_vars[[1]])
  test_dat_lev2_vars <- paste0(names_of_lev2_vars, collapse = "+") #get level 2 variables and combine them into a formula for later
  test_dat_lev2_vars <- paste0("Y~(",test_dat_lev2_vars, ")", "*time")
  test_dat_form <- as.formula(paste0(test_dat_lev2_vars))
  print(test_dat_lev2_vars)
  stan_dat <- list()
  
  #Creating lists of data to feed into STAN sampler
  for(i in seq_along(training_datasets)){
    holder <- as.data.frame(level2_vars[[i]])
    temp <- list(
      N_obs_train = nrow(training_datasets[[i]]),
      N_pts_train = n_distinct(training_datasets[[i]]$id.new),
      L = 2, K = ncol(level2_vars[[i]])+1,
      pid_train = training_datasets[[i]]$id.new,
      x_train = cbind(1, training_datasets[[i]]$time),
      x2_train = cbind(1, holder),
      y_train = training_datasets[[i]]$Y,
      N_obs_test = nrow(testing_datasets[[i]]),
      test_data = model.matrix(test_dat_form, data = testing_datasets[[i]])
    )
    stan_dat[[i]] <- temp
  }
  #print(stan_dat[[1]])
  return(stan_dat)
}

#big function for scaling, splitting, and construction
# ready_dat <- function(sim_dat){
#   #Scaling
#   temp <- scale_cont_data(sim_dat)
#   
#   #Splitting
#   split.sim1 <- list()
#   for(i in seq_along(temp)){
#     temp1 <- temp[[i]]
#     temp2 <- tt_split(datasets = temp1, percent_train = 0.80)
#     split.sim1[[i]] <- temp2
#   }
#   #Constructing
#   split.dat <- list()
#   for(i in seq_along(split.sim1)){
#     temp <- split.sim1[[i]]
#     temp <- stan_data_loop(training_datasets = temp$Training, testing_datasets = temp$Testing)
#     split.dat[[i]] <- temp
#   }
#   return(split.dat)
# }

#Trying to see if we can simply return a list with both the split and the stan data

ready_dat1 <- function(sim_dat){
  #Scaling
  temp <- scale_cont_data(sim_dat)
  
  #Splitting
  split.sim1 <- list()
  for(i in seq_along(temp)){
    temp1 <- temp[[i]]
    temp2 <- tt_split(datasets = temp1, percent_train = 0.80)
    split.sim1[[i]] <- temp2
  }
  #Constructing
  split.dat <- list()
  for(i in seq_along(split.sim1)){
    temp <- split.sim1[[i]]
    temp <- stan_data_loop(training_datasets = temp$Training, testing_datasets = temp$Testing)
    split.dat[[i]] <- temp
  }
  results <- setNames(list(split.sim1, split.dat), c('splitdat', 'standat'))
  return(results)
}


#function to extract all of the STAN output
stan_out <- function(stan_data_collection, stan_file = "pred_error_uninform.stan", iterations = 3000,
                     chains_to_run = 4){
  output <- list()
  for(i in seq_along(stan_data_collection)){
    print("Hello")
    #print(stan_data_collection[[i]])
    stan_fit <- stan(file = stan_file, data = stan_data_collection[[i]], iter = iterations, 
                     chains = chains_to_run)
    df_of_draws <- stan_fit
    output[[i]] <- df_of_draws
  }
  return(output)
}

#Function for Simulation
simulate_study <- function(pos, cond){
  ###select condition
  prior <- cond$prior[pos]
  conditione <- cond$condition[pos]
  print(prior)
  print(conditione)
  
  #If Else stuff to select the correct dataset from the generated sets
  if(conditione==1){ #If we are looking at condition 1, we want to use the the out dataset
    #and use the portion of the list for the standata
    dataset <- out$standat
  }
  else if(conditione==2){
    dataset <- out1$standat
  }
  else if(conditione==3){
    dataset <- out2$standat
  }
  else if(conditione==4){
    dataset <- out3$standat
  }
  else if(conditione==5){
    dataset <- out4$standat
  }
  else if(conditione==6){
    dataset <- out5$standat
  }
  
  print(dataset)
  #Get STAN results from the stan data
  file_name <- paste0("new_", prior, ".stan")
  print(file_name)
  
  #CMDSTAN to SAMPLE
  mod <- cmdstan_model(file_name)
  
  #Storing the STAN data
  stan.data <- list()
  #Loop across the ~500 datasets for each condition that we generate
  for(i in seq_along(dataset)){
    temp <- dataset[[i]]
    stan.data_1 <- list()
    for(j in seq_along(temp)){
      #Perform STAN sampling using CmdstanR
      temp1 <- mod$sample(data = temp[[j]], seed = 123, chains = 4, parallel_chains = 4, refresh = 500)
      #Convert the list to a df of the MC draws
      draws_temp <- as_draws_df(temp1)
      #Store the complete list of draws in a list in case we need it later
      stan.data_1[[j]] <- draws_temp
    }
    stan.data[[i]] <- stan.data_1 
  }
  return(stan.data)
}

#Specify function to extract the gamma draws
gamma_study <- function(list_of_complete_draws){
  gamma_draws <- list()
  for(i in seq_along(list_of_complete_draws)){
    temp <- list_of_complete_draws[[i]]
    gam_d <- list()
    for(j in seq_along(temp)){
      temp1 <- subset_draws(temp[[j]], variable = "gamma")
      gam_d[[j]] <- temp1
    }
    gamma_draws[[i]] <- gam_d
  }
  return(gamma_draws)
}

#Function to get point estimate of gamma draws
gamma_point_estimates <- function(list_of_gamma_draws){
  point_list <- list()
  for(i in seq_along(list_of_gamma_draws)){
    temp <- list_of_gamma_draws[[i]]
    pe_d <- list()
    for(j in seq_along(temp)){
      temp1 <- summarise_draws(temp[[j]], "mean", "median") #get both the mean and median
      pe_d[[j]] <- temp1
    }
    point_list[[i]] <- pe_d
  }
  return(point_list)
}

#Do the same with the y generated 
y_gen_study <- function(list_of_complete_draws){
  y_gen_draws <- list()
  for(i in seq_along(list_of_complete_draws)){
    temp <- list_of_complete_draws[[i]]
    y_d <- list()
    for(j in seq_along(temp)){
      temp1 <- subset_draws(temp[[j]], variable = "y_new") #Select only y_new variables
      y_d[[j]] <- temp1
    }
    y_gen_draws[[i]] <- y_d
  }
  return(y_gen_draws)
}

#Need point estimates of the y generated draws
y_gen_point_estimates <- function(list_generated_draws){
  point_list <- list()
  for(i in seq_along(list_generated_draws)){
    temp <- list_generated_draws[[i]]
    pe_d <- list()
    for(j in seq_along(temp)){
      temp1 <- summarise_draws(temp[[j]], "mean", "median") #get both the mean and median
      pe_d[[j]] <- temp1
    }
    point_list[[i]] <- pe_d
  }
  return(point_list)
}

#Get upper and lower bounds on CI based on an 80% CI
gamm_interval_construct <- function(list_of_gamma_draws){
  confidence_int <- list()
  for(i in seq_along(list_of_gamma_draws)){
    temp <- list_of_gamma_draws[[i]]
    ci_d <- list()
    for(j in seq_along(temp)){
      temp1 <- summarise_draws(temp[[j]], function(x) quantile(x, probs = c(0.1, 0.9)))
      ci_d[[j]] <- temp1
    }
    confidence_int[[i]] <- ci_d
  }
  return(confidence_int)
}

#Count the number of times that 0 is inside the interval
#This should get us a Boolean list of TRUE/FALSE
inlc_excl_gam <- function(gamma_ci_list){
  value <- 0 #we care if this value is inside the interval
  in_out <- list()
  for(i in seq_along(gamma_ci_list)){
    temp <- gamma_ci_list[[i]]
    temp1 <- list()
    for(j in seq_along(temp)){
      #Basically see if the 10% bound and the 90% bound contain 0.
      #Would ordinarily return 1/TRUE if 0 is included and 0/FALSE if 0 is not in the range
      #We flip that with the ! operator b/c we want 0 if the variable is not significant
      #and 1 if the variables CI excludes 0
      temp2 <- !data.table::between(value, temp[[j]]$'10%', temp[[j]]$'90%')
      temp1[[j]] <- temp2
    }
    in_out[[i]] <- temp1
  }
  return(in_out)
}

#Now need a function to count the number of times that TRUE/FALSE
#check <- do.call(rbind, test4[[1]])
sum_inlc_excl <- function(boolean_list){
  out <- list()
  for(i in seq_along(boolean_list)){
    temp <- data.frame(do.call(rbind, boolean_list[[i]]))
    out[[i]] <- temp
  }
  return(out)
}


###WORKING HERE
#Sum each column in the combined df
percent_in_ex <- function(list_of_dfs){
  out <- list()
  for(i in seq_along(list_of_dfs)){
    temp <- list_of_dfs[[i]]
    temp1 <- colSums(temp)
    #print(temp1)
    out[[i]] <- temp1/nrow(temp)
  }
  #get number of rows in df
  return(out)
}

# 
# names(ah4) <- names
# for (i in seq_along(ah4)){
# 
#   print(names(temp))
#   for (j in seq_along(temp)){
#     
#     print(length(temp1))
#     names(temp1) <- names1
#     print(names(temp1))
#   }
# }

#Need to build a function to compare prediction error
#Include option to use the mean or the median for calculations
prediction_error_func <- function(gen_data, test_data, moc_to_use = "mean"){
  rmse_vals <- list()
  for(i in seq_along(gen_data)){
    #access the same generated/test data for each condition
    temp <- gen_data[[i]]
    temp1 <- test_data$splitdat[[i]]
    #print(i)
    testing <- list()
    for(j in seq_along(temp)){
      print(j)
      actdat <- temp1$Testing[[j]]$Y
      if(moc_to_use == "mean"){
        generated_dat <- temp[[j]]$mean
      }
      else{
        generated_dat <- temp[[j]]$median
      }
      print(generated_dat)
    #   if(moc_to_use == "mean"){
    #     actdat <- temp1[[j]]$mean
    #   }
    #   else{
    #     actdat <- temp1[[j]]$median
    #   }
    #   print(generated_dat)
      val <- sqrt(mean((actdat - generated_dat)^2)) #Compute RMSE
      testing[[length(testing)+1]] <- val #Store RMSE value in list
    }
    rmse_vals[[length(rmse_vals)+1]] <- testing
  }
  return(rmse_vals)
}

#RMSE values are returned as a list.  We need to unlist these
rmse_unlist <- function(rmse_list){
  out <- list()
  for(i in seq_along(rmse_list)){
    temp1 <- rmse_list[[i]]
    temp2 <- unlist(temp1)
    out[[i]] <- temp2
  }
  return(out)
}
