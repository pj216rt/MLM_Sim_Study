source("functions.R")

#Applies to the test case
# test1 <- gamma_study(test)
# test2 <- y_gen_study(test)
# test3 <- gamm_interval_construct(test1)
# test4 <- inlc_excl_gam(test3)
# test5 <- sum_inlc_excl(test4)
# test6 <- percent_in_ex(test5)
# 
# check <- bind_rows(test6)
# row.names(check) <- sample_size
# df <- melt(as.matrix(check))

#Need to use lapply
ah <- lapply(simulation, gamma_study)
ah1 <- lapply(ah, gamm_interval_construct)
ah2 <- lapply(ah1, inlc_excl_gam)
ah3 <- lapply(ah2, sum_inlc_excl)
ah4 <- lapply(ah3, percent_in_ex)

#name the output list to match up with the conditions of the simulation
names <- apply(sim_conditions, MARGIN = 1, function(x) paste0(x[1], " condition:", x[2]))
names(ah4) <- names 

#name the inside lists to match the data generation conditions
names1 <- apply(conditions, MARGIN = 1, function(x) paste0("n: ", x[1], " r: ", x[2], " SE: ", x[3]))
names1

for (i in seq_along(ah4)){
  temp <- ah4[i]
  print(names(temp))
  for (j in seq_along(temp)){
    temp1 <- temp[[j]]
    print(length(temp1))
    names(temp1) <- names1
    print(names(temp1))
  }
}

  
#Plotting probability of inclusion.
p <- ggplot(df, aes(x=X1, y=value, color = X2)) + geom_point() +
  guides(color = guide_legend(title = "Title")) + 
  scale_color_discrete(labels = paste("Gamma", 1:10)) + 
  labs(x = "Sample Size", y = "Probability of inclusion", title = "Probability of Variable Inclusion vs. Sample Size")
p


#Tables?
library(gt)
check %>% gt(rownames_to_stub = T)
