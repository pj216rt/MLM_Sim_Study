library(cmdstanr)
library(posterior)
library(bayesplot)

check_cmdstan_toolchain()
cmdstan_path()

#Testing new stan files
mod <- cmdstan_model("new_uniform_mixture.stan")
fit <- mod$sample(data = out$standat[[1]][[1]], seed = 123, chains = 4, parallel_chains = 4, refresh = 500)
fit <- as_draws_df(fit)

test <- simulate_study(pos = 2, cond = sim_conditions)

#Applies to the test case
test1 <- gamma_study(test)
test2 <- y_gen_study(test)
check <- gamma_point_estimates(test1)
check1 <- y_gen_point_estimates(test2)

test3 <- gamm_interval_construct(test1)
test4 <- inlc_excl_gam(test3)
test5 <- sum_inlc_excl(test4)
test6 <- percent_in_ex(test5)

#We have PRMSE values here
test7 <- prediction_error_func(gen_data = check1, test_data = out)
test8 <- rmse_unlist((test7))

#Boxplot of test case
boxplot(test8, col = "orange",border = "brown",
        horizontal = TRUE,
        notch = TRUE)

#Can we turn this list of RMSE values into a df?
test9 <- as.data.frame(do.call(cbind, test8))


###IMPORTANT###
library(stats) #For setnames function
#Testing naming different levels of a list in a list
varA = paste0("varA", 1:10)
varB = paste0("varB", 1:3)

library(foreach)
tabs = foreach(j = 1:length(varA)) %do% {
  main = varA[j]
  mytabs = lapply(1:length(varB), class)
}
#setnames, then lapply(nested list, call setNames, call the inner names), then finally outer names
tabs <- setNames(lapply(tabs, setNames, varB), varA)
