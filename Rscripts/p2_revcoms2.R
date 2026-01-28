# model based predictions and effect sizes

source(file = "scripts/00_set_up.R")
source("scripts/h0_fxn.R")

set.seed(seed=11)

################
######## data ## political complexity ~ alcohol + phylogeny + space + env
################

# wide
read.csv(file=paste0(outfol,"workingdata/thedata_wide.csv"), stringsAsFactors = F)->thedata
colnames(thedata)[colnames(thedata)%in%"samplewide"]<-"Alcohol"
table(thedata$Alcohol)

thedata$Alcohol<-as.factor(thedata$Alcohol) # binary predictor
thedata$PolC<-factor(thedata$PolC, levels=levels(as.factor(thedata$PolC)), ordered=T) # ordered response
thedata$PC1<-scale(thedata$PC1) # continuous predictor
thedata$PC2<-scale(thedata$PC2)
thedata$PC3<-scale(thedata$PC3)
thedata$Agriculture<-factor(thedata$Agriculture, levels=levels(as.factor(thedata$Agriculture)), ordered=T) # ordered predictor
table(thedata$PolC)

thedata$Language_ID <- thedata$SCCS_ID
thedata$Language_ID2 <- thedata$SCCS_ID
#checking for duplicated locations and jittering them
duplicate_coords = thedata[duplicated(thedata[,c("lon", "lat")]) | duplicated(thedata[,c("lon", "lat")], fromLast = TRUE),"Language_ID"]
duplicate_rowid = thedata$Language_ID %in% duplicate_coords
thedata$lat[duplicate_rowid] = jitter(thedata$lat[duplicate_rowid], factor = 1)
thedata$lon[duplicate_rowid] = jitter(thedata$lon[duplicate_rowid], factor = 1)


##################
######## model 5 #
##################

library(brms)
library(marginaleffects)
library(tidyverse)

load(file=paste0(outfol,"models/model5_pc1_R.rds")) ; model5_pc1 -> model ; rm(model5_pc1)

# predictions using mean model coefficients + values for PC1, Agr at each datapoint, A=0 and A=1
eff <- marginaleffects::predictions(model, newdata = thedata, variables = "Alcohol", type = "response",re_formula = NULL, resp="PolC") # re_formula=NULL is default anyway
grtest <- 5 # here modify group number to look at each group
filtered_preds <- eff %>%
  filter(group == grtest, Alcohol %in% c(0,1))
ggplot(filtered_preds, aes(x = estimate, fill = factor(Alcohol))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("red", "blue"), labels = c("A = 0", "A = 1")) +
  labs(
    title = "Posterior Distributions of Predicted Probability",
    subtitle = "Comparing A = 0 vs. A = 1 for group = 1",
    x = "Predicted Probability",
    y = "Density",
    fill = "A"
  ) +
  theme_minimal()

eff_diffs <- eff %>%
  select(SCCS_ID, group, Alcohol, estimate) %>%  # Keep only relevant columns
  pivot_wider(names_from = Alcohol, values_from = estimate, names_prefix = "A_") %>%  # Reshape
  mutate(diff = A_1 - A_0)  # Compute the difference

ggplot(eff_diffs, aes(x = diff, fill = factor(group))) +
  stat_halfeye() +
  scale_fill_manual(values = colorspace::lighten(c("red", "blue", "green", "purple", "orange"), 0.4)) +
  labs(x = "Change in Probability (A=1 - A=0)", y = "Density", fill = "Y category",
       title = "Marginal Effect of A on Probability of Y Categories",
       subtitle = "Random effects averaged / marginalized / integrated") +
  theme_minimal()


# predictions using posterior model coefficients + all values in df for PC1 & Agr from the data*
  # ie: each draw = one sample set of model coefficients, run across all rows in thedata (ie for values for PC1 and Agr), then average
marginal_preds <- predictions(
  model, 
  newdata = thedata, 
  type = "response",  # Ensures probabilities are returned
  by = c("Alcohol", "PolC"),  # Keep track of which Y category each row corresponds to
  resp="PolC",
  re_formula = NULL  # This averages over random effects, as opposed to NA
) %>% posterior_draws()

# average over the PolC column (undertainty in Agriculture) - this needs to be done only for wmodel5_pc1, not for model56_pc1
marginal_preds <- marginal_preds %>%
  group_by(drawid, group, Alcohol) %>%  # Aggregate over uncertainty in Ag
  summarize(
    draw = mean(draw),  # Compute the mean probability across all Ag samples
    estimate = mean(estimate),
    conf.low = mean(conf.low),
    conf.high = mean(conf.high),
    .groups = "drop")

grtest <- 5
filtered_preds <- marginal_preds %>%
  filter(group == grtest, Alcohol %in% c(0,1))
ggplot(filtered_preds, aes(x = draw, fill = factor(Alcohol))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("red", "blue"), labels = c("A = 0", "A = 1")) +
  labs(
    title = "Posterior Distributions of Predicted Probability",
    subtitle = "Comparing A = 0 vs. A = 1 for group = 1",
    x = "Predicted Probability",
    y = "Density",
    fill = "A"
  ) +
  theme_minimal()

# differences
marginal_diffs <- marginal_preds %>%
  select(drawid, group, Alcohol, draw) %>%  # Keep only relevant columns
  pivot_wider(names_from = Alcohol, values_from = draw, names_prefix = "A_") %>%  # Reshape
  mutate(diff = A_1 - A_0)  # Compute the difference

p1 <- ggplot(marginal_diffs, aes(x = diff, fill = factor(group))) +
  stat_halfeye() +
  scale_fill_manual(values = colorspace::lighten(c("red", "blue", "green", "purple", "orange"), 0.4)) +
  labs(x = "Change in Probability (A=1 - A=0)", y = "Density", fill = "Y category",
       #title = "Marginal Effect of A on Probability of Y Categories",
       title = "Strictly empirical grid",
       #subtitle = "Random effects averaged / marginalized / integrated") +
       subtitle = "") +
  theme_minimal()
plot(p1)

# what we don't have now is the values for societies for A=0 and A=1 (because it only did the A that the society actually has in thedata)
  # duplicate thedata so that each soc gets an Alcohol =1 and Alcohol=0 ie Counterfactual grid + Empirical grid
newdata2 <- thedata %>%
  mutate(Alcohol = 0) %>%  # Set A=0
  bind_rows(thedata %>% mutate(Alcohol = 1))  # Duplicate rows with A=1

marginal_preds2 <- predictions(
  model, 
  newdata = newdata2, 
  type = "response",  # Ensures probabilities are returned
  by = c("Alcohol", "PolC"),  # Keep track of which Y category each row corresponds to
  resp="PolC",
  re_formula = NULL  # This averages over random effects, as opposed to NA
) %>% posterior_draws()

marginal_preds2 <- marginal_preds2 %>%
  group_by(drawid, group, Alcohol) %>%  # Aggregate over uncertainty in Ag
  summarize(
    draw = mean(draw),  # Compute the mean probability across all Ag samples
    estimate = mean(estimate),
    conf.low = mean(conf.low),
    conf.high = mean(conf.high),
    .groups = "drop")

grtest <- 5
filtered_preds <- marginal_preds2 %>%
  filter(group == grtest, Alcohol %in% c(0,1))
ggplot(filtered_preds, aes(x = draw, fill = factor(Alcohol))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("red", "blue"), labels = c("A = 0", "A = 1")) +
  labs(
    title = "Posterior Distributions of Predicted Probability",
    subtitle = "Comparing A = 0 vs. A = 1 for group = 1",
    x = "Predicted Probability",
    y = "Density",
    fill = "A"
  ) +
  theme_minimal()

# differences
marginal_diffs2 <- marginal_preds2 %>%
  select(drawid, group, Alcohol, draw) %>%  # Keep only relevant columns
  pivot_wider(names_from = Alcohol, values_from = draw, names_prefix = "A_") %>%  # Reshape
  mutate(diff = A_1 - A_0)  # Compute the difference

p2 <- ggplot(marginal_diffs2, aes(x = diff, fill = factor(group))) +
  stat_halfeye() +
  scale_fill_manual(values = colorspace::lighten(c("red", "blue", "green", "purple", "orange"), 0.4)) +
  labs(x = "Change in Probability (A=1 - A=0)", y = "Density", fill = "Y category",
       #title = "Marginal Effect of A on Probability of Y Categories",
       title = "Empirical counterfactual grid",
       #subtitle = "Random effects averaged / marginalized / integrated") +
       subtitle = "") +
  theme_minimal()
plot(p2)

# or a more detailed plot (in ms)

library(ggplot2)
library(ggridges)
library(patchwork)
library(dplyr)

# Compute summary statistics for left panel (Predicted Probabilities)
summary_data_left <- marginal_diffs2 %>%
  group_by(group) %>%
  summarize(
    mean_A0 = mean(A_0, na.rm = TRUE),
    lower_A0_95 = quantile(A_0, 0.025, na.rm = TRUE),
    upper_A0_95 = quantile(A_0, 0.975, na.rm = TRUE),
    lower_A0_66 = quantile(A_0, 0.17, na.rm = TRUE),  # 66% CI
    upper_A0_66 = quantile(A_0, 0.83, na.rm = TRUE),  # 66% CI
    mean_A1 = mean(A_1, na.rm = TRUE),
    lower_A1_95 = quantile(A_1, 0.025, na.rm = TRUE),
    upper_A1_95 = quantile(A_1, 0.975, na.rm = TRUE),
    lower_A1_66 = quantile(A_1, 0.17, na.rm = TRUE),  # 66% CI
    upper_A1_66 = quantile(A_1, 0.83, na.rm = TRUE),  # 66% CI
    .groups = "drop"
  )
colnames(summary_data_left)

summary_data_right <- marginal_diffs2 %>%
  group_by(group) %>%
  summarize(
    mean_diff = mean(diff, na.rm = TRUE),
    lower_95 = quantile(diff, 0.025, na.rm = TRUE),
    upper_95 = quantile(diff, 0.975, na.rm = TRUE),
    lower_66 = quantile(diff, 0.17, na.rm = TRUE),  # 66% CI
    upper_66 = quantile(diff, 0.83, na.rm = TRUE),  # 66% CI
    .groups = "drop"
  )
colnames(summary_data_right)

# for plotting
summary_data_left_A0 <- summary_data_left %>%
  mutate(A = "A=0")  # Label for A=0

summary_data_left_A1 <- summary_data_left %>%
  mutate(A = "A=1")  # Label for A=1

# Left Panel: Ridge Plot of Predicted Probabilities for A=0 vs. A=1
p_left <- ggplot(marginal_diffs2, aes(x = A_0, y = factor(group), fill = factor("A=0"))) +
  geom_density_ridges(color = NA, alpha = 0.8, scale = 1, rel_min_height = 0.01) +
  geom_density_ridges(data = marginal_diffs2, aes(x = A_1, y = factor(group), fill = factor("A=1")), 
                      color = NA, alpha = 0.8, scale = 1, rel_min_height = 0.01) +
  
  geom_point(data = summary_data_left_A0, aes(x = mean_A0, y = factor(group), fill = factor(A)), 
             color = "black", size = 2, shape = 19) +  # Ensures black dot is inside the fill
  geom_point(data = summary_data_left_A1, aes(x = mean_A1, y = factor(group), fill = factor(A)), 
             color = "black", size = 2, shape = 19) +  # Black dot inside the filled shape
  
  # 95% CI (thinner)
  geom_segment(data = summary_data_left_A0, 
               aes(x = lower_A0_95, xend = upper_A0_95, y = factor(group), yend = factor(group)), 
               color = "black", linewidth = 0.8) +
  geom_segment(data = summary_data_left_A1, 
               aes(x = lower_A1_95, xend = upper_A1_95, y = factor(group), yend = factor(group)), 
               color = "black", linewidth = 0.8) +
  
  # 66% CI (thicker)
  geom_segment(data = summary_data_left_A0, 
               aes(x = lower_A0_66, xend = upper_A0_66, y = factor(group), yend = factor(group)), 
               color = "black", linewidth = 1.5) +
  geom_segment(data = summary_data_left_A1, 
               aes(x = lower_A1_66, xend = upper_A1_66, y = factor(group), yend = factor(group)), 
               color = "black", linewidth = 1.5) +
  
  scale_fill_manual(values = c("A=0" = "orange", "A=1" = "royalblue")) +
  scale_y_discrete(limits = rev(unique(marginal_diffs2$group))) +
  labs(x = "Predicted probability", y = "Political complexity categories", fill = "A",
       title = "Posterior Predictions") +
  theme_minimal() + theme(plot.title=element_text(face="bold",size=14), 
                          axis.title = element_text(size = 14),  # Increase axis labels
                          axis.text = element_text(size = 12))
p_left

# Right Panel: Ridge Plot of Differences (A_1 - A_0)
p_right <- ggplot(marginal_diffs2, aes(x = diff, y = factor(group), fill = factor(group))) +
  geom_density_ridges(color = NA, alpha = 0.8, rel_min_height = 0.01, scale = 1) +
  geom_point(data = summary_data_right, aes(x = mean_diff, y = factor(group)), color = "black", size = 2) +
  
  # 95% CI (thinner)
  geom_segment(data = summary_data_right, 
               aes(x = lower_95, xend = upper_95, y = factor(group), yend = factor(group)), 
               color = "black", linewidth = 0.8) +
  # 66% CI (thicker)
  geom_segment(data = summary_data_right, 
               aes(x = lower_66, xend = upper_66, y = factor(group), yend = factor(group)), 
               color = "black", linewidth = 1.5) +
  
  scale_fill_manual(values = c("deeppink4", "deeppink4", "deeppink4", "deeppink4", "deeppink4")) +
  scale_y_discrete(limits = rev(unique(marginal_diffs2$group))) +
  labs(x = "Probability difference", y = "", fill = "Y category",
       title = "Marginal effect of alcohol") +
  theme_minimal() + theme(plot.title=element_text(face="bold",size=14), legend.position="none",
                          axis.title = element_text(size = 14),  # Increase axis labels
                          axis.text = element_text(size = 12))
p_right

# Combine the two plots into a single figure
final_plot <- p_left + p_right + plot_layout(ncol = 2)

# Print the final figure
pdf(paste0(outfol,"plots_review/margeff.pdf"), width=10, height=7)
print(final_plot)
dev.off()


# set the other predictors to median/mode = balanced grid
newdata3 <- newdata2 # take the data with alcohol = 0 and alcohol = 1
newdata3$PC1 <- mean(thedata$PC1)
table(thedata$Agriculture) # most common one is 4
newdata3$Agriculture <- 4 
marginal_preds3 <- predictions(
  model, 
  newdata = newdata3, 
  type = "response",  # Ensures probabilities are returned
  by = c("Alcohol", "PolC"),  # Keep track of which Y category each row corresponds to
  resp="PolC",
  re_formula = NULL  # This averages over random effects, as opposed to NA
) %>% posterior_draws()

marginal_preds3 <- marginal_preds3 %>%
  group_by(drawid, group, Alcohol) %>%  # Aggregate over uncertainty in Ag
  summarize(
    draw = mean(draw),  # Compute the mean probability across all Ag samples
    estimate = mean(estimate),
    conf.low = mean(conf.low),
    conf.high = mean(conf.high),
    .groups = "drop")

grtest <- 5
filtered_preds <- marginal_preds3 %>%
  filter(group == grtest, Alcohol %in% c(0,1))
ggplot(filtered_preds, aes(x = draw, fill = factor(Alcohol))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("red", "blue"), labels = c("A = 0", "A = 1")) +
  labs(
    title = "Posterior Distributions of Predicted Probability",
    subtitle = "Comparing A = 0 vs. A = 1 for group = 1",
    x = "Predicted Probability",
    y = "Density",
    fill = "A"
  ) +
  theme_minimal()

# differences
marginal_diffs3 <- marginal_preds3 %>%
  select(drawid, group, Alcohol, draw) %>%  # Keep only relevant columns
  pivot_wider(names_from = Alcohol, values_from = draw, names_prefix = "A_") %>%  # Reshape
  mutate(diff = A_1 - A_0)  # Compute the difference

p3 <- ggplot(marginal_diffs3, aes(x = diff, fill = factor(group))) +
  stat_halfeye() +
  scale_fill_manual(values = colorspace::lighten(c("red", "blue", "green", "purple", "orange"), 0.4)) +
  labs(x = "Change in Probability (A=1 - A=0)", y = "Density", fill = "Y category",
       #title = "Marginal Effect of A on Probability of Y Categories",
       title = "Empirical counterfactual grid",
       #subtitle = "Random effects averaged / marginalized / integrated") +
       subtitle = "other predictors set at averages") +
  theme_minimal()
plot(p3)
# lastly if we use a sampling method over a distribution of PC1 and Ag values? ie not an empirical grid, but a more smoothed out one based still on our data
  # ie not sample over purely hypothetical PC1 and balanced Ag rather hypothetical PC1 and Ag inspired from the data

hist(thedata$PC1)
E_mean <- mean(thedata$PC1, na.rm = TRUE)
E_sd <- sd(thedata$PC1, na.rm = TRUE)
Ag_probs <- thedata %>%
  count(Agriculture) %>%
  mutate(prob = n / sum(n))  # Convert counts to probabilities
n_samples <- 160  # Same as your dataset size
newdata_grid <- tibble(
  Alcohol = rep(c(0, 1), each = n_samples),  # Each row appears with A=0 and A=1
  PC1 = rnorm(n_samples * 2, mean = E_mean, sd = E_sd),  # Sample E from normal distribution
  Agriculture= sample(Ag_probs$Ag, size = n_samples * 2, replace = TRUE, prob = Ag_probs$prob),  # Sample Ag based on observed distribution
  Language_ID = sample(thedata$Language_ID, size = n_samples * 2, replace = TRUE),  # Assign a Language_ID for random effects
  Language_ID2 = sample(thedata$Language_ID2, size = n_samples * 2, replace = TRUE)
)

expanded_newdata_grid <- newdata_grid %>%
  crossing(PolC = 1:5)  # Ensure each row has all Y categories
marginal_preds_grid <- predictions(
  model, 
  newdata = expanded_newdata_grid, 
  type = "response",  
  by = c("Alcohol", "PolC"),  
  resp="PolC",
  re_formula = NULL  
) %>% posterior_draws()

marginal_preds_grid <- marginal_preds_grid %>%
  group_by(drawid, group, Alcohol) %>%  # Aggregate over uncertainty in Ag
  summarize(
    draw = mean(draw),  # Compute the mean probability across all Ag samples
    estimate = mean(estimate),
    conf.low = mean(conf.low),
    conf.high = mean(conf.high),
    .groups = "drop")

marginal_diffs4 <- marginal_preds_grid %>%
  select(drawid, group, Alcohol, draw) %>%  # Keep only relevant columns
  pivot_wider(names_from = Alcohol, values_from = draw, names_prefix = "A_") %>%  # Reshape
  mutate(diff = A_1 - A_0)  # Compute the difference

p4 <- ggplot(marginal_diffs4, aes(x = diff, fill = factor(group))) +
  stat_halfeye() +
  scale_fill_manual(values = colorspace::lighten(c("red", "blue", "green", "purple", "orange"), 0.4)) +
  labs(x = "Change in Probability (A=1 - A=0)", y = "Density", fill = "Y category",
       #title = "Marginal Effect of A on Probability of Y Categories",
       title = "Hypothetical grid inspired by data distribution",
       #subtitle = "Random effects averaged / marginalized / integrated") +
       subtitle = "") +
  theme_minimal()
plot(p4)


library(patchwork)
pdf("MargEffects_totalmodel5.pdf", width=10)
p1 + p2 + 
  p3 + p4 + 
  plot_layout(ncol = 2, nrow = 2)
dev.off()


######################
#### other models ####
######################

library(ggplot2)
library(ggridges)
library(patchwork)
library(dplyr)


newdata2 <- thedata %>%
  mutate(Alcohol = 0) %>%  # Set A=0
  bind_rows(thedata %>% mutate(Alcohol = 1))  # Duplicate rows with A=1

load(file=paste0(outfol,"models/model1R.rds")) ; model1 -> wmodel1 ; rm(model1)
load(file=paste0(outfol,"models/model2R.rds")) ; model2 -> wmodel2 ; rm(model2)
load(file=paste0(outfol,"models/model3_pc1R.rds")) ; model3_pc1 -> wmodel3_pc1 ; rm(model3_pc1) 
load(file=paste0(outfol,"models/model4R.rds")) ; model4 -> wmodel4 ; rm(model4)
load(file=paste0(outfol,"models/model5_pc1_R.rds")) ; model5_pc1 -> wmodel5_pc1 ; rm(model5_pc1)


function_me <- function(model,namem,newdata2){
  marginal_preds2 <- predictions(
    model, 
    newdata = newdata2, 
    type = "response",  # Ensures probabilities are returned
    by = c("Alcohol", "PolC"),  # Keep track of which Y category each row corresponds to
    resp="PolC",
    re_formula = NULL  # This averages over random effects, as opposed to NA
  ) %>% posterior_draws()
  
  marginal_preds2 <- marginal_preds2 %>%
    group_by(drawid, group, Alcohol) %>%  # Aggregate over uncertainty in Ag
    summarize(
      draw = mean(draw),  # Compute the mean probability across all Ag samples
      estimate = mean(estimate),
      conf.low = mean(conf.low),
      conf.high = mean(conf.high),
      .groups = "drop")
  
  
  # differences
  marginal_diffs2 <- marginal_preds2 %>%
    select(drawid, group, Alcohol, draw) %>%  # Keep only relevant columns
    pivot_wider(names_from = Alcohol, values_from = draw, names_prefix = "A_") %>%  # Reshape
    mutate(diff = A_1 - A_0)  # Compute the difference
  
  p2 <- ggplot(marginal_diffs2, aes(x = diff, fill = factor(group))) +
    stat_halfeye() +
    scale_fill_manual(values = colorspace::lighten(c("red", "blue", "green", "purple", "orange"), 0.4)) +
    labs(x = "Change in Probability (A=1 - A=0)", y = "Density", fill = "Y category",
         #title = "Marginal Effect of A on Probability of Y Categories",
         title = "Empirical counterfactual grid",
         #subtitle = "Random effects averaged / marginalized / integrated") +
         subtitle = "") +
    theme_minimal()
  plot(p2)
  
  # or a more detailed plot (in ms)
  
  # Compute summary statistics for left panel (Predicted Probabilities)
  summary_data_left <- marginal_diffs2 %>%
    group_by(group) %>%
    summarize(
      mean_A0 = mean(A_0, na.rm = TRUE),
      lower_A0_95 = quantile(A_0, 0.025, na.rm = TRUE),
      upper_A0_95 = quantile(A_0, 0.975, na.rm = TRUE),
      lower_A0_66 = quantile(A_0, 0.17, na.rm = TRUE),  # 66% CI
      upper_A0_66 = quantile(A_0, 0.83, na.rm = TRUE),  # 66% CI
      mean_A1 = mean(A_1, na.rm = TRUE),
      lower_A1_95 = quantile(A_1, 0.025, na.rm = TRUE),
      upper_A1_95 = quantile(A_1, 0.975, na.rm = TRUE),
      lower_A1_66 = quantile(A_1, 0.17, na.rm = TRUE),  # 66% CI
      upper_A1_66 = quantile(A_1, 0.83, na.rm = TRUE),  # 66% CI
      .groups = "drop"
    )
  colnames(summary_data_left)
  
  summary_data_right <- marginal_diffs2 %>%
    group_by(group) %>%
    summarize(
      mean_diff = mean(diff, na.rm = TRUE),
      lower_95 = quantile(diff, 0.025, na.rm = TRUE),
      upper_95 = quantile(diff, 0.975, na.rm = TRUE),
      lower_66 = quantile(diff, 0.17, na.rm = TRUE),  # 66% CI
      upper_66 = quantile(diff, 0.83, na.rm = TRUE),  # 66% CI
      .groups = "drop"
    )
  colnames(summary_data_right)
  
  # for plotting
  summary_data_left_A0 <- summary_data_left %>%
    mutate(A = "A=0")  # Label for A=0
  
  summary_data_left_A1 <- summary_data_left %>%
    mutate(A = "A=1")  # Label for A=1
  
  # Left Panel: Ridge Plot of Predicted Probabilities for A=0 vs. A=1
  p_left <- ggplot(marginal_diffs2, aes(x = A_0, y = factor(group), fill = factor("A=0"))) +
    geom_density_ridges(color = NA, alpha = 0.8, scale = 1, rel_min_height = 0.01) +
    geom_density_ridges(data = marginal_diffs2, aes(x = A_1, y = factor(group), fill = factor("A=1")), 
                        color = NA, alpha = 0.8, scale = 1, rel_min_height = 0.01) +
    
    geom_point(data = summary_data_left_A0, aes(x = mean_A0, y = factor(group), fill = factor(A)), 
               color = "black", size = 2, shape = 19) +  # Ensures black dot is inside the fill
    geom_point(data = summary_data_left_A1, aes(x = mean_A1, y = factor(group), fill = factor(A)), 
               color = "black", size = 2, shape = 19) +  # Black dot inside the filled shape
    
    # 95% CI (thinner)
    geom_segment(data = summary_data_left_A0, 
                 aes(x = lower_A0_95, xend = upper_A0_95, y = factor(group), yend = factor(group)), 
                 color = "black", linewidth = 0.8) +
    geom_segment(data = summary_data_left_A1, 
                 aes(x = lower_A1_95, xend = upper_A1_95, y = factor(group), yend = factor(group)), 
                 color = "black", linewidth = 0.8) +
    
    # 66% CI (thicker)
    geom_segment(data = summary_data_left_A0, 
                 aes(x = lower_A0_66, xend = upper_A0_66, y = factor(group), yend = factor(group)), 
                 color = "black", linewidth = 1.5) +
    geom_segment(data = summary_data_left_A1, 
                 aes(x = lower_A1_66, xend = upper_A1_66, y = factor(group), yend = factor(group)), 
                 color = "black", linewidth = 1.5) +
    
    scale_fill_manual(values = c("A=0" = "orange", "A=1" = "royalblue")) +
    scale_y_discrete(limits = rev(unique(marginal_diffs2$group))) +
    labs(x = "Predicted probability", y = "Political complexity categories", fill = "A",
         title = "Posterior Predictions") +
    theme_minimal() + theme(plot.title=element_text(face="bold",size=14), 
                            axis.title = element_text(size = 14),  # Increase axis labels
                            axis.text = element_text(size = 12))
  
  # Right Panel: Ridge Plot of Differences (A_1 - A_0)
  p_right <- ggplot(marginal_diffs2, aes(x = diff, y = factor(group), fill = factor(group))) +
    geom_density_ridges(color = NA, alpha = 0.8, rel_min_height = 0.01, scale = 1) +
    geom_point(data = summary_data_right, aes(x = mean_diff, y = factor(group)), color = "black", size = 2) +
    
    # 95% CI (thinner)
    geom_segment(data = summary_data_right, 
                 aes(x = lower_95, xend = upper_95, y = factor(group), yend = factor(group)), 
                 color = "black", linewidth = 0.8) +
    # 66% CI (thicker)
    geom_segment(data = summary_data_right, 
                 aes(x = lower_66, xend = upper_66, y = factor(group), yend = factor(group)), 
                 color = "black", linewidth = 1.5) +
    
    scale_fill_manual(values = c("deeppink4", "deeppink4", "deeppink4", "deeppink4", "deeppink4")) +
    scale_y_discrete(limits = rev(unique(marginal_diffs2$group))) +
    labs(x = "Probability difference", y = "", fill = "Y category",
         title = "Marginal effect of alcohol") +
    theme_minimal() + theme(plot.title=element_text(face="bold",size=14), legend.position="none",
                            axis.title = element_text(size = 14),  # Increase axis labels
                            axis.text = element_text(size = 12))
  
  # Combine the two plots into a single figure
  final_plot <- p_left + p_right + plot_layout(ncol = 2)
  # Print the final figure
  pdf(paste0(outfol,"plots_review/margeff_model",namem,".pdf"), width=10, height=7)
  print(final_plot)
  dev.off()
  plot(final_plot)
  return(summary_data_right)
}

# run for all
function_me(wmodel1,namem="model1",newdata2)
function_me(wmodel2,namem="model2",newdata2)
function_me(wmodel3_pc1,namem="model3",newdata2)
function_me(wmodel4,namem="model4",newdata2)
function_me(wmodel5_pc1,namem="model5",newdata2)

# conservative
read.csv(file=paste0(outfol,"workingdata/thedata_cons.csv"), stringsAsFactors = F)->thedatac
colnames(thedatac)[colnames(thedatac)%in%"samplecons"]<-"Alcohol"
table(thedatac$Alcohol)
thedatac$Alcohol<-as.factor(thedatac$Alcohol) # binary predictor
thedatac$PolC<-factor(thedatac$PolC, levels=levels(as.factor(thedatac$PolC)), ordered=T) # ordered response
thedatac$Language_ID <- thedatac$SCCS_ID
thedatac$Language_ID2 <- thedatac$SCCS_ID
thedatac$Agriculture<-factor(thedatac$Agriculture, levels=levels(as.factor(thedatac$Agriculture)), ordered=T) # ordered predictor
duplicate_coords = thedatac[duplicated(thedatac[,c("lon", "lat")]) | duplicated(thedatac[,c("lon", "lat")], fromLast = TRUE),"Language_ID"]
duplicate_rowid = thedatac$Language_ID %in% duplicate_coords
thedatac$lat[duplicate_rowid] = jitter(thedatac$lat[duplicate_rowid], factor = 1)
thedatac$lon[duplicate_rowid] = jitter(thedatac$lon[duplicate_rowid], factor = 1)

newdata2c <- thedatac %>%
  mutate(Alcohol = 0) %>%  # Set A=0
  bind_rows(thedatac %>% mutate(Alcohol = 1))  # Duplicate rows with A=1

# conservative
load(file=paste0(outfol,"models_cons/cmodel1R.rds")) ; model1 -> cmodel1 ; rm(model1)
load(file=paste0(outfol,"models_cons/cmodel2R.rds")) ; model2 -> cmodel2 ; rm(model2)
load(file=paste0(outfol,"models_cons/cmodel3_pc1R.rds")) ; model3_pc1 -> cmodel3_pc1 ; rm(model3_pc1)
load(file=paste0(outfol,"models_cons/cmodel4R.rds")) ; model4 -> cmodel4 ; rm(model4)
load(file=paste0(outfol,"models_cons/cmodel5_pc1_R.rds")) ; model5_pc1 -> cmodel5_pc1 ; rm(model5_pc1)

# run for all
function_me(cmodel1,namem="cmodel1",newdata2c)
function_me(cmodel2,namem="cmodel2",newdata2c)
function_me(cmodel3_pc1,namem="cmodel3_pc1",newdata2c)
function_me(cmodel4,namem="cmodel4",newdata2c)
function_me(cmodel5_pc1,namem="cmodel5_pc1",newdata2c)
