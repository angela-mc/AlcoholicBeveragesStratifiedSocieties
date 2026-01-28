
rm(list=ls())
source(file = "scripts/00_set_up.R")
dir <-  paste0(outfol, "models_para")
if(!dir.exists(dir)) dir.create(dir)

library(brms)
library(readxl)


##############
## read data #
##############

para<-read_excel(paste0(projfol,"Dataset/PARAGRAPHS_dataset.xlsx")) 
read.csv(file=paste0(outfol,"workingdata/thedata_wide.csv"), stringsAsFactors = F)->thedata
colnames(thedata)[colnames(thedata)%in%"samplewide"]<-"Alcohol"
thedata$Alcohol<-as.factor(thedata$Alcohol) # binary response

para <- para[para$SCCS_ID%in%thedata$SCCS_ID,]
para<-para[match(thedata$SCCS_ID, para$SCCS_ID),]
thedata$Paragraphs_270 <- para$Paragraphs_270
thedata$Paragraphs_advanced <- para$Paragraphs_advanced

library(brms)
# m270 <- brms:: bf(Alcohol ~ Paragraphs_270, family = bernoulli()) 
# model270<- brm (m270, data=thedata, cores=4, chains=2) 
# summary(model270)
# 
# mAdv<- brms:: bf(Alcohol ~ Paragraphs_advanced, family = bernoulli()) 
# modelAdv<- brm (mAdv, data=thedata, cores=4, chains=2) 
# summary(modelAdv)


# Paragraphs_270 ~ Paragraphs_advanced
plot(sqrt(thedata$Paragraphs_270) ~ sqrt(thedata$Paragraphs_advanced), xlab="Advanced", ylab="270", pch=21)
abline(a=0,b=1)
correlation <- cor(thedata$Paragraphs_270, thedata$Paragraphs_advanced, method = "spearman") #0.6188093
correlation

# SQRT 
thedata$Paragraphs_270_s <- sqrt(thedata$Paragraphs_270)
thedata$Paragraphs_advanced_s <- sqrt(thedata$Paragraphs_advanced)
m270s <- brms:: bf(Alcohol ~ Paragraphs_270_s, family = bernoulli()) 
model270s<- brm (m270s, data=thedata, cores=4, chains=2) 
summary(model270s)

res <- as_draws_df(model270s)
hist(res$b_Paragraphs_270_s, xlab="", main="Para 270 sq--> Alcohol", breaks=40, col="#eda31d", border="#e9b4a3")
median(res$b_Paragraphs_270_s)
sum(res$b_Paragraphs_270_s>=0)/length(res$b_Paragraphs_270_s) # 0.814 ie most are higher than 0 but the effect is very small

outl <- list(summary(model270s),sum(res$b_Paragraphs_270_s>=0)/length(res$b_Paragraphs_270_s))
names(outl) <- c("model summary", "proportion of posterior coefficient >= 0") 
capture.output(outl, file=paste0(outfol,"models_para/model270_summary.txt"))


mAdvs<- brms:: bf(Alcohol ~ Paragraphs_advanced_s, family = bernoulli()) 
modelAdvs<- brm (mAdvs, data=thedata, cores=4, chains=2) 
summary(modelAdvs)

res <- as_draws_df(modelAdvs)
hist(res$b_Paragraphs_advanced_s, xlab="", main="Para adv sq--> Alcohol", breaks=40, col="#eda31d", border="#e9b4a3")
median(res$b_Paragraphs_advanced_s)
sum(res$b_Paragraphs_advanced_s>=0)/length(res$b_Paragraphs_advanced_s) # 0.8015 ie most are higher than 0 but the effect is very small

outl <- list(summary(modelAdvs),sum(res$b_Paragraphs_advanced_s>=0)/length(res$b_Paragraphs_advanced_s))
names(outl) <- c("model summary", "proportion of posterior coefficient >= 0") 
capture.output(outl, file=paste0(outfol,"models_para/modelAdv_summary.txt"))

save(model270s, file=paste0(outfol,"models_para/model270s.rds"))
save(modelAdvs, file=paste0(outfol,"models_para/modelAdvs.rds"))

source("scripts/h0_fxn.R")
load(paste0(outfol,"models_para/model270s.rds"))
res <- as_draws_df(model270s)
mediqr (res$b_Paragraphs_270_s)
  # 0.0276048; -0.004682758 -> 0.05989235
load(paste0(outfol,"models_para/modelAdvs.rds"))
res <- as_draws_df(modelAdvs)
mediqr (res$b_Paragraphs_advanced_s)
  # 0.02773574 -0.006714897 -> 0.06218638