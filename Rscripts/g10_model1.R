source(file = "scripts/00_set_up.R")

# output_files <- c(paste0(outfol,"models/model1R.rds"), 
#                   paste0(outfol,"models/model1_summary.txt"))
# 
# if(all(file.exists(output_files))){
#   cat(paste0("The output files for this model run already exists.\n", 
#              output_files), ".\n")
#   
# }else{
#   script_name <- "g10_model1.R"
#    start_time <- Sys.time()



################
######## data ## political complexity ~ alcohol
################


# wide
read.csv(file=paste0(outfol,"workingdata/thedata_wide.csv"), stringsAsFactors = F)->thedata
colnames(thedata)[colnames(thedata)%in%"samplewide"]<-"Alcohol"
table(thedata$Alcohol)

# thedata$Alcohol<-thedata$Alcohol+1
thedata$Alcohol<-as.factor(thedata$Alcohol) # binary predictor
thedata$PolC<-factor(thedata$PolC, levels=levels(as.factor(thedata$PolC)), ordered=T) # ordered response


#############
#### model ## ordered response ~ binary predictor
#############

library(phylopath)
library(tidyr)
library(tidyverse)

m <- define_model_set( model1 = c(PolC ~ Alcohol))
positions <- data.frame(name = c('Alcohol', 'PolC'), x = c(1, 2), y = c(1, 1))
plot_model_set(m,manual_layout = positions, edge_width = 0.5)

library(brms)
m1 <- brms:: bf(PolC ~ Alcohol, family = cumulative()) 
  # PC = ordered response so fam = cumulative - we can add probit or no
  # Alcohol = binary

model1<- brm (m1, data=thedata, cores=4, chains=2) # no prior
summary(model1)
  # intercepts = related to cutpoints
  # Alcohol1 = effect of alcohol
 save(model1, file=paste0(outfol,"models/model1R.rds"))

load(file=paste0(outfol,"models/model1R.rds"))
conditional_effects(model1, "Alcohol", categorical = TRUE)
  # probability of values in each PolC category for Alcohol absent (0) or present (1) 
  # is this equivalent to posterior predictions plot(?)

res <- as_draws_df(model1)
hist(res$b_Alcohol1, xlab="", main="Alcohol --> PolC", breaks=40, col="#eda31d", border="#e9b4a3")
median(res$b_Alcohol1)
sum(res$b_Alcohol1>0)/length(res$b_Alcohol1)

# EXPORT
capture.output(summary(model1), file=paste0(outfol,"models/model1_summary.txt"))

# inspect priors
get_prior(model1)

# end_time <- Sys.time()
# time <- end_time - start_time
# cat(
# paste0("I'm done with script ", script_name, ".\n",
#   "It took ", round(time[[1]], 2), " ",  units(time)),
# ".\n"
# )
# }

#########################
######## alternatives ### category specific effects 
#########################

m_seq <- brm(formula = PolC ~ Alcohol, data = thedata,family = sratio())
m_seq_cs <- brm(formula = PolC ~ cs(Alcohol), data = thedata,family = sratio())

#save(m_seq, file=paste0(outfol, "models_otherCAT/m_seq.rds"))
#save(m_seq_cs, file=paste0(outfol, "models_otherCAT/m_seq_cs.rds"))
modelc<- loo(m_seq, m_seq_cs) # For elpd (expected log predictive density), the larger the better; negative diff - second model = worse; see diff compared to SE difference
# save(modelc, file=paste0(outfol, "models_otherCAT/model1_comp.rds"))
# capture.output(print(modelc), file=paste0(outfol,"models_otherCAT/model1_comp_summary.txt"))

