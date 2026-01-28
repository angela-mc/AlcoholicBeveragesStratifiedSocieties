source(file = "scripts/00_set_up.R")

################
######## data ## political complexity ~ alcohol
################

# conservative
read.csv(file=paste0(outfol,"workingdata/thedata_cons.csv"), stringsAsFactors = F)->thedata
colnames(thedata)[colnames(thedata)%in%"samplecons"]<-"Alcohol"
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

model1<- brm (m1, data=thedata, cores=4, chains=2) # no prior
summary(model1)
 save(model1, file=paste0(outfol,"models_cons/cmodel1R.rds"))

load(file=paste0(outfol,"models_cons/cmodel1R.rds"))
conditional_effects(model1, "Alcohol", categorical = TRUE)
# probability of values in each PolC category for Alcohol absent (0) or present (1) 

res <- as_draws_df(model1)
hist(res$b_Alcohol1, xlab="", main="Alcohol --> PolC", breaks=40, col="#eda31d", border="#e9b4a3")
median(res$b_Alcohol1)

# EXPORT
capture.output(summary(model1), file=paste0(outfol,"models_cons/cmodel1_summary.txt"))

