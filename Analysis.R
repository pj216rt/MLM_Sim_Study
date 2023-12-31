source("functions.R")
library(stats) #For setnames function
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

#Need to use lapply across these rows
ah <- lapply(simulation, gamma_study)
ah1 <- lapply(ah, gamm_interval_construct)
ah2 <- lapply(ah1, inlc_excl_gam)
ah3 <- lapply(ah2, sum_inlc_excl)
ah4 <- lapply(ah3, percent_in_ex)


#Naming lists.  Outer list
outer_names <- apply(sim_conditions, MARGIN = 1, function(x) paste0(x[1], " condition:", x[2]))
#Naming inner list
inner_names <- apply(conditions, MARGIN = 1, function(x) paste0(" n=", x[1], " p=", x[2], " SE=", x[3]))

#setnames, then lapply(nested list, call setNames, call the inner names), then finally outer names
named_list <- setNames(lapply(ah4, setNames, inner_names), outer_names)

#Saving this output
#Saving output of simulation to enviroment
saveRDS(named_list, file = "small simulation results")


library(purrr)
library(data.table)
library(stringr)
library(reshape2)
#Extract condition 1:
cond1_subset <- named_list[grepl("condition:1", names(named_list))]
cond1_subset

do1 <- unlist(cond1_subset, recursive = FALSE)
names(do1)
do2 <- map_dfr(do1, ~as.data.frame(t(.x)))
do2$situation <- names(do1)
do2 <- do2 %>% relocate(situation)
do3 <- str_split_fixed(do2$situation, " ", 5)
do3
do2$prior <- do3[,1]
do2$condition <- do3[,2]
do2$n <- do3[,3]
do2$p <- do3[,4]
do2$SE <- do3[,5]
do2 <- subset(do2, select = -(situation))
do2 <- do2 %>% relocate(prior, condition, n, p, SE)




#Plotting probability of inclusion.
p <- ggplot(df, aes(x=X1, y=value, color = X2)) + geom_point() +
  guides(color = guide_legend(title = "Title")) + 
  scale_color_discrete(labels = paste("Gamma", 1:10)) + 
  labs(x = "Sample Size", y = "Probability of inclusion", title = "Probability of Variable Inclusion vs. Sample Size")
p

#Tables?
library(gt)
check %>% gt(rownames_to_stub = T)
