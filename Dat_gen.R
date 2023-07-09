source("functions.R")

#Set seed
set.seed(12345)

#Generating data for the study
#Conditions to simulate

num_datasets <- 4
sample_size <- c(50, 100)
corr_sizes <- c(0.0)
stand_error <- c(5)

#Using expand.grid to come up with all combinations fo sample size and correlation
conditions <- expand.grid(sample.size = sample_size, correlations = corr_sizes, noise = stand_error)

#Datasets
#First case
test <- apply(conditions, MARGIN = 1, 
              function(x) gen_balanced_datasets(nreps = num_datasets, 
                                                nSubjs = x[1], corrRE = x[2], sdErr = x[3],
                                                coef1 = c(4, 3),
                                                # types of level 2 covariates
                                                level2Binary = c(.2, 0.7),
                                                level2Continuous = list(c(mu = 0, sd = 1),
                                                                        c(mu = 5, sd = 1)),
                                                # sd of errors on subject-specific int and slope
                                                sdRE = c(1, 1),
                                                # for each predictor in level 2, (int2, slope2) 
                                                # specify effect on level 2 intercept and level 2 slope
                                                coef2Binary = list(c(int2 = 4.0, slope2 = 3.0),
                                                                   c(int2 = 2.0, slope2 = 0.0)),
                                                coef2Continuous = list(c(int2 = 0.0, slope2 = 0.0),
                                                                       c(int2 = 0.0, slope2 = 0.0))))

#Getting Data Ready for Stan Compilation
out <- ready_dat(test)
#Plotting Data
xyplot(Y ~ time, data = test[[1]][[1]], type = "b", groups = id, 
       xlab="Time",
       ylab="Response Variable",
       main="Simulated Data: Response Variable vs. Time")

###Second case###
test1 <- apply(conditions, MARGIN = 1, 
              function(x) gen_balanced_datasets(nreps = num_datasets, 
                                                nSubjs = x[1], corrRE = x[2], sdErr = x[3],
                                                coef1 = c(4, 3),
                                                # types of level 2 covariates
                                                level2Binary = c(.2, 0.7, 0.5, 0.3),
                                                level2Continuous = list(c(mu = 0, sd = 1),
                                                                        c(mu = 5, sd = 1),
                                                                        c(mu = 0, sd = 1),
                                                                        c(mu = 5, sd = 1)),
                                                # sd of errors on subject-specific int and slope
                                                sdRE = c(1, 1),
                                                coef2Binary = list(c(int2 = 4.0, slope2 = 3.0),
                                                                   c(int2 = 2.0, slope2 = 0.0),
                                                                   c(int2 = 0.0, slope2 = 0.0),
                                                                   c(int2 = 0.0, slope2 = 0.0)),
                                                coef2Continuous = list(c(int2 = 0.0, slope2 = 0.0),
                                                                       c(int2 = 0.0, slope2 = 0.0),
                                                                       c(int2 = 0.0, slope2 = 5.0),
                                                                       c(int2 = 2.0, slope2 = 2.0))))
#Getting Data Ready for Stan Compilation
out1 <- ready_dat(test1)
#Plotting Data
xyplot(Y ~ time, data = test1[[1]][[1]], type = "b", groups = id, 
       xlab="Time",
       ylab="Response Variable",
       main="Simulated Data: Response Variable vs. Time")


###Third Case###
test2 <- apply(conditions, MARGIN = 1, 
              function(x) gen_balanced_datasets(nreps = num_datasets, 
                                                nSubjs = x[1], corrRE = x[2], sdErr = x[3],
                                                coef1 = c(2, 0),
                                                # types of level 2 covariates
                                                level2Binary = c(.2, 0.7),
                                                level2Continuous = list(c(mu = 0, sd = 1),
                                                                        c(mu = 5, sd = 1)),
                                                # sd of errors on subject-specific int and slope
                                                sdRE = c(1, 1),
                                                # for each predictor in level 2, (int2, slope2) 
                                                # specify effect on level 2 intercept and level 2 slope
                                                coef2Binary = list(c(int2 = 2.0, slope2 = 0.0),
                                                                   c(int2 = 0.0, slope2 = 0.0)),
                                                coef2Continuous = list(c(int2 = 0.0, slope2 = 0.0),
                                                                       c(int2 = 0.0, slope2 = 0.0))))

#Getting Data Ready for Stan Compilation
out2 <- ready_dat(test2)
#Plotting Data
xyplot(Y ~ time, data = test2[[1]][[1]], type = "b", groups = id, 
       xlab="Time",
       ylab="Response Variable",
       main="Simulated Data: Response Variable vs. Time")


##Fourth Case###
test3 <- apply(conditions, MARGIN = 1, 
              function(x) gen_balanced_datasets(nreps = num_datasets, 
                                                nSubjs = x[1], corrRE = x[2], sdErr = x[3],
                                                coef1 = c(4, 3),
                                                # types of level 2 covariates
                                                level2Binary = rep(c(.2, 0.7, 0.5, 0.3),  3),
                                                level2Continuous = rep(list(c(mu = 0, sd = 1),
                                                                            c(mu = 5, sd = 1),
                                                                            c(mu = 0, sd = 1),
                                                                            c(mu = 5, sd = 1)), 3),
                                                # sd of errors on subject-specific int and slope
                                                sdRE = c(1, 1),
                                                # for each predictor in level 2, (int2, slope2) 
                                                # specify effect on level 2 intercept and level 2 slope
                                                coef2Binary = rep(list(c(int2 = 4.0, slope2 = 3.0),
                                                                       c(int2 = 2.0, slope2 = 0.0),
                                                                       c(int2 = 0.0, slope2 = 0.0),
                                                                       c(int2 = 0.0, slope2 = 0.0)),3),
                                                coef2Continuous = rep(list(c(int2 = 0.0, slope2 = 0.0),
                                                                           c(int2 = 0.0, slope2 = 0.0),
                                                                           c(int2 = 0.0, slope2 = 5.0),
                                                                           c(int2 = 2.0, slope2 = 2.0)),3)))

#Getting Data Ready for Stan Compilation
out3 <- ready_dat(test3)
#Plotting Data
xyplot(Y ~ time, data = test3[[1]][[1]], type = "b", groups = id, 
       xlab="Time",
       ylab="Response Variable",
       main="Simulated Data: Response Variable vs. Time")


#Fifrth Case#
test4 <- apply(conditions, MARGIN = 1, 
               function(x) gen_balanced_datasets(nreps = num_datasets, 
                                                 nSubjs = x[1], corrRE = x[2], sdErr = x[3],
                                                 coef1 = c(2, 0),
                                                 # types of level 2 covariates
                                                 level2Binary = rep(c(.2, 0.7),  3),
                                                 level2Continuous = rep(list(c(mu = 0, sd = 1),
                                                                             c(mu = 5, sd = 1)), 3),
                                                 # sd of errors on subject-specific int and slope
                                                 sdRE = c(1, 1),
                                                 # for each predictor in level 2, (int2, slope2) 
                                                 # specify effect on level 2 intercept and level 2 slope
                                                 coef2Binary = rep(list(c(int2 = 2.0, slope2 = 0.0),
                                                                        c(int2 = 0.0, slope2 = 0.0)),3),
                                                 coef2Continuous = rep(list(c(int2 = 0.0, slope2 = 0.0),
                                                                            c(int2 = 0.0, slope2 = 0.0)),3)))

#Getting Data Ready for Stan Compilation
out4 <- ready_dat(test4)
#Plotting Data
xyplot(Y ~ time, data = test4[[1]][[1]], type = "b", groups = id, 
       xlab="Time",
       ylab="Response Variable",
       main="Simulated Data: Response Variable vs. Time")


#Sixth Case#
test5 <- apply(conditions, MARGIN = 1, 
               function(x) gen_balanced_datasets(nreps = num_datasets, 
                                                 nSubjs = x[1], corrRE = x[2], sdErr = x[3],
                                                 coef1 = c(2, 0),
                                                 # types of level 2 covariates
                                                 level2Binary = rep(c(.2, 0.7, 0.5, 0.3, 0.4),  20),
                                                 level2Continuous = rep(list(c(mu = 0, sd = 1),
                                                                             c(mu = 5, sd = 1),
                                                                             c(mu = 0, sd = 1),
                                                                             c(mu = 5, sd = 1),
                                                                             c(mu = 0, sd = 1)), 20),
                                                 # sd of errors on subject-specific int and slope
                                                 sdRE = c(1, 1),
                                                 # for each predictor in level 2, (int2, slope2) 
                                                 # specify effect on level 2 intercept and level 2 slope
                                                 coef2Binary = rep(list(c(int2 = 4.0, slope2 = 3.0),
                                                                        c(int2 = 2.0, slope2 = 0.0),
                                                                        c(int2 = 0.0, slope2 = 0.0),
                                                                        c(int2 = 0.0, slope2 = 0.0),
                                                                        c(int2 = 0.0, slope2 = 0.0)), 20),
                                                 coef2Continuous = rep(list(c(int2 = 0.0, slope2 = 0.0),
                                                                            c(int2 = 0.0, slope2 = 0.0),
                                                                            c(int2 = 0.0, slope2 = 0.0),
                                                                            c(int2 = 0.0, slope2 = 5.0),
                                                                            c(int2 = 2.0, slope2 = 2.0)),20)))

#Getting Data Ready for Stan Compilation
out5 <- ready_dat(test5)
#Plotting Data
xyplot(Y ~ time, data = test5[[1]][[1]], type = "b", groups = id, 
       xlab="Time",
       ylab="Response Variable",
       main="Simulated Data: Response Variable vs. Time")
