source("functions.R")

library(reshape)

test1 <- gamma_study(test)
test2 <- y_gen_study(test)
test3 <- gamm_interval_construct(test1)
test4 <- inlc_excl_gam(test3)
test5 <- sum_inlc_excl(test4)
test6 <- percent_in_ex(test5)

check <- bind_rows(test6)
row.names(check) <- sample_size
df <- melt(as.matrix(check))


#Plotting probability of inclusion.
p <- ggplot(df, aes(x=X1, y=value, color = X2)) + geom_point() +
  guides(color = guide_legend(title = "Title")) + 
  scale_color_discrete(labels = paste("Gamma", 1:10)) + 
  labs(x = "Sample Size", y = "Probability of inclusion", title = "Probability of Variable Inclusion vs. Sample Size")
p


#Tables?
library(gt)
check %>% gt(rownames_to_stub = T)
