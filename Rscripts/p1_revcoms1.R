
# model based predictions and effect sizes

source(file = "scripts/00_set_up.R")
source("scripts/h0_fxn.R")


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


#################### 
# load all models ##
####################

library(brms)
load(file=paste0(outfol,"models/model1R.rds")) ; model1 -> wmodel1 ; rm(model1)
load(file=paste0(outfol,"models/model2R.rds")) ; model2 -> wmodel2 ; rm(model2)
load(file=paste0(outfol,"models/model3_pc1R.rds")) ; model3_pc1 -> wmodel3_pc1 ; rm(model3_pc1) 
load(file=paste0(outfol,"models/model4R.rds")) ; model4 -> wmodel4 ; rm(model4)
load(file=paste0(outfol,"models/model5_pc1_R.rds")) ; model5_pc1 -> wmodel5_pc1 ; rm(model5_pc1)

load(file=paste0(outfol,"models_cons/cmodel1R.rds")) ; model1 -> cmodel1 ; rm(model1)
load(file=paste0(outfol,"models_cons/cmodel2R.rds")) ; model2 -> cmodel2 ; rm(model2)
load(file=paste0(outfol,"models_cons/cmodel3_pc1R.rds")) ; model3_pc1 -> cmodel3_pc1 ; rm(model3_pc1)
load(file=paste0(outfol,"models_cons/cmodel4R.rds")) ; model4 -> cmodel4 ; rm(model4)
load(file=paste0(outfol,"models_cons/cmodel5_pc1_R.rds")) ; model5_pc1 -> cmodel5_pc1 ; rm(model5_pc1)


###################### 
# model predictions ##
######################


# # use predict (less data wrangling, though annoying layout)
# dataA1 <- thedata
# dataA1$Alcohol <- 1
# dataA0 <- thedata
# dataA0$Alcohol <- 0
# pred1 <- predict(wmodel5_pc1, newdata = dataA1, type = "response", resp = "PolC")
# pred0 <- predict(wmodel5_pc1, newdata = dataA0, type = "response", resp = "PolC")
# diff in probabilities
# diffp <- pred1[,1:5]-pred0[,1:5]
# diffp$glottocode <- pred1$glottocode
# colMeans(diffp[,1:5]) # probability differences 

# E(Y)
# v<- runif(5)  # Random values between 0 and 1
# v <- v / sum(v)  # Normalize to sum to 1
# names(v) <- c(1,2,3,4,5)
# print(round(v, digits=2))
# sum(v)
# ey <-1*v[1] + 2*v[2] + 3*v[3] + 4*v[4] + 5*v[5]
# print(ey)


library(marginaleffects)
fxn_predictions<- function(df){
  # make two dfs, one for alcohol = 1 and one for alcohol = 0
  df1 <- data.frame("glottocode"=unique(df$glottocode), "P1"=-1, "P2"=-1, "P3"=-1, "P4"=-1,"P5"=-1, "MaxY"=-1, "EY"=-1, "Alcohol"=1, "lat"=-999, "lon"=-999)
  df0 <- data.frame("glottocode"=unique(df$glottocode), "P1"=-1, "P2"=-1, "P3"=-1, "P4"=-1,"P5"=-1, "MaxY"=-1, "EY"=-1, "Alcohol"=0, "lat"=-999, "lon"=-999)

  for(i in 1: dim(df1)[1]){
    dfi <- df[df$glottocode%in%df1$glottocode[i],]
    df1$lat[i] <- df0$lat[i] <- dfi$lat[1] # add coordinates for the possibility of plotting this on the map
    df1$lon[i] <- df0$lon[i] <- dfi$lon[1]
    df1$P1[i] <- dfi[dfi$Alcohol==1 & dfi$group==1,]$estimate
    df1$P2[i] <- dfi[dfi$Alcohol==1 & dfi$group==2,]$estimate
    df1$P3[i] <- dfi[dfi$Alcohol==1 & dfi$group==3,]$estimate
    df1$P4[i] <- dfi[dfi$Alcohol==1 & dfi$group==4,]$estimate
    df1$P5[i] <- dfi[dfi$Alcohol==1 & dfi$group==5,]$estimate
    
    df0$P1[i] <- dfi[dfi$Alcohol==0 & dfi$group==1,]$estimate
    df0$P2[i] <- dfi[dfi$Alcohol==0 & dfi$group==2,]$estimate
    df0$P3[i] <- dfi[dfi$Alcohol==0 & dfi$group==3,]$estimate
    df0$P4[i] <- dfi[dfi$Alcohol==0 & dfi$group==4,]$estimate
    df0$P5[i] <- dfi[dfi$Alcohol==0 & dfi$group==5,]$estimate
    
    # max category
    df1$MaxY[i] <- c(1:5) [ which(df1[i,2:6] %in% max(df1[i,2:6] ) )]
    df0$MaxY[i] <- c(1:5) [ which(df0[i,2:6] %in% max(df0[i,2:6] ) )]
    
    # expected Y: E(Y)=1×P(Y=1)+2×P(Y=2)+3×P(Y=3)+4×P(Y=4)+5×P(Y=5)
    df1$EY[i] <-  1* df1$P1[i] + 2* df1$P2[i] + 3* df1$P3[i]  + 4* df1$P4[i]  + 5* df1$P5[i]
    df0$EY[i] <-  1* df0$P1[i] + 2* df0$P2[i] + 3* df0$P3[i]  + 4* df0$P4[i]  + 5* df0$P5[i]
  }  
  toret <- list(df1, df0)
  names(toret) <- c("A1","A0")
  return(toret)
}

# model 5
effects_all <- marginaleffects::predictions(wmodel5_pc1, newdata = thedata, variables = "Alcohol", type = "response", resp="PolC") 
alld <- fxn_predictions(effects_all)
pred1 <- alld$A1
pred0 <- alld$A0
diffp <- pred1[,2:8]-pred0[,2:8]
diffp$glottocode <- pred1$glottocode
colMeans(diffp[,1:7]) # differences 
table(diffp$MaxY)
hist(diffp$EY)
abline(v=quantile(diffp$EY, probs = c(0.055, 0.945))) # where 89% of data is
abline(v=quantile(diffp$EY, probs = c(0.025, 0.975))) # where 95% of the data is
quantile(diffp$EY, probs = c(0.025, 0.975))

# if you want averages and looking at probabilities only:
# effects <- marginaleffects::avg_predictions(model, newdata = thedata, variables = "Alcohol", type = "response",resp = "PolC")
  # give predictions() thedata, ie the empirical grid + setting alcohol at A=0 and then at A=1 (~ counterfactual grid)

# model 4
effects_all <- marginaleffects::predictions(wmodel4, newdata = thedata, variables = "Alcohol", type = "response", resp="PolC") 
alld <- fxn_predictions(effects_all)
pred1 <- alld$A1
pred0 <- alld$A0
diffp <- pred1[,2:8]-pred0[,2:8]
diffp$glottocode <- pred1$glottocode
colMeans(diffp[,1:7]) # differences 
hist(diffp$EY)
abline(v=quantile(diffp$EY, probs = c(0.025, 0.975))) # where 95% of the data is
quantile(diffp$EY, probs = c(0.025, 0.975))

# model 3
effects_all <- marginaleffects::predictions(wmodel3_pc1, newdata = thedata, variables = "Alcohol", type = "response", resp="PolC") 
alld <- fxn_predictions(effects_all)
pred1 <- alld$A1
pred0 <- alld$A0
diffp <- pred1[,2:8]-pred0[,2:8]
diffp$glottocode <- pred1$glottocode
colMeans(diffp[,1:7]) # differences 
hist(diffp$EY)
abline(v=quantile(diffp$EY, probs = c(0.025, 0.975))) # where 95% of the data is
quantile(diffp$EY, probs = c(0.025, 0.975))

# model 2
effects_all <- marginaleffects::predictions(wmodel2, newdata = thedata, variables = "Alcohol", type = "response", resp="PolC") 
alld <- fxn_predictions(effects_all)
pred1 <- alld$A1
pred0 <- alld$A0
diffp <- pred1[,2:8]-pred0[,2:8]
diffp$glottocode <- pred1$glottocode
colMeans(diffp[,1:7]) # differences 
hist(diffp$EY)
abline(v=quantile(diffp$EY, probs = c(0.025, 0.975))) # where 95% of the data is
quantile(diffp$EY, probs = c(0.025, 0.975))

# model 1
effects_all <- marginaleffects::predictions(wmodel1, newdata = thedata, variables = "Alcohol", type = "response", resp="PolC") 
alld <- fxn_predictions(effects_all)
pred1 <- alld$A1
pred0 <- alld$A0
diffp <- pred1[,2:8]-pred0[,2:8]
diffp$glottocode <- pred1$glottocode
colMeans(diffp[,1:7]) # differences # with marginal effects it predicts the same

# with predict
dataA1 <- thedata
dataA1$Alcohol <- 1
dataA0 <- thedata
dataA0$Alcohol <- 0
effects_predict1 <- predict(wmodel1, newdata = dataA1, type = "response", resp = "PolC") 
effects_predict0 <- predict(wmodel1, newdata = dataA0, type = "response", resp = "PolC")
effects_A0_df <- as.data.frame(effects_predict0) %>% mutate(Alcohol = 0)
effects_A1_df <- as.data.frame(effects_predict1) %>% mutate(Alcohol = 1)
effects_A0_df$glottocode <- effects_A1_df$glottocode <- thedata$glottocode
effects_A0_df$Ey <- effects_A1_df$Ey <- -1
effects_A0_df$MaxY <- effects_A1_df$MaxY <- -1
for(i in 1:160){
  effects_A0_df$Ey[i] <- 1* effects_A0_df$`P(Y = 1)`[i] + 2* effects_A0_df$`P(Y = 2)`[i]+ 3* effects_A0_df$`P(Y = 3)`[i] + 4* effects_A0_df$`P(Y = 4)`[i] + 5* effects_A0_df$`P(Y = 5)`[i]
  effects_A1_df$Ey[i] <- 1* effects_A1_df$`P(Y = 1)`[i] + 2* effects_A1_df$`P(Y = 2)`[i]+ 3* effects_A1_df$`P(Y = 3)`[i] + 4* effects_A1_df$`P(Y = 4)`[i] + 5* effects_A1_df$`P(Y = 5)`[i]
  effects_A0_df$MaxY[i] <- c(1:5) [which(effects_A0_df[i,1:5] %in% max(effects_A0_df[i,1:5]))]
  effects_A1_df$MaxY[i] <- c(1:5) [which(effects_A1_df[i,1:5] %in% max(effects_A1_df[i,1:5]))]
}
diffE <- effects_A1_df$Ey - effects_A0_df$Ey
hist(diffE)
mean(diffE)
median(diffE)
diffMaxY <-effects_A1_df$MaxY - effects_A0_df$MaxY
hist(diffMaxY)
mean(diffMaxY)


################
######## data ## conservative
################

# conservative
read.csv(file=paste0(outfol,"workingdata/thedata_cons.csv"), stringsAsFactors = F)->thedata
colnames(thedata)[colnames(thedata)%in%"samplecons"]<-"Alcohol"
table(thedata$Alcohol)

thedata$Alcohol<-as.factor(thedata$Alcohol) # binary predictor
thedata$PolC<-factor(thedata$PolC, levels=levels(as.factor(thedata$PolC)), ordered=T) # ordered response
thedata$PC1<-scale(thedata$PC1) # continuous predictor
thedata$PC2<-scale(thedata$PC2)
thedata$PC3<-scale(thedata$PC3)
thedata$Agriculture<-factor(thedata$Agriculture, levels=levels(as.factor(thedata$Agriculture)), ordered=T) # ordered predictor

thedata$Language_ID <- thedata$SCCS_ID
thedata$Language_ID2 <- thedata$SCCS_ID
#checking for duplicated locations and jittering them
duplicate_coords = thedata[duplicated(thedata[,c("lon", "lat")]) | duplicated(thedata[,c("lon", "lat")], fromLast = TRUE),"Language_ID"]
duplicate_rowid = thedata$Language_ID %in% duplicate_coords
thedata$lat[duplicate_rowid] = jitter(thedata$lat[duplicate_rowid], factor = 1)
thedata$lon[duplicate_rowid] = jitter(thedata$lon[duplicate_rowid], factor = 1)

# model 5
effects_all <- marginaleffects::predictions(cmodel5_pc1, newdata = thedata, variables = "Alcohol", type = "response", resp="PolC") 
alld <- fxn_predictions(effects_all)
pred1 <- alld$A1
pred0 <- alld$A0
diffp <- pred1[,2:8]-pred0[,2:8]
diffp$glottocode <- pred1$glottocode
colMeans(diffp[,1:7]) # differences 
hist(diffp$EY)
abline(v=quantile(diffp$EY, probs = c(0.025, 0.975))) # where 95% of the data is
quantile(diffp$EY, probs = c(0.025, 0.975))

# model 4
effects_all <- marginaleffects::predictions(cmodel4, newdata = thedata, variables = "Alcohol", type = "response", resp="PolC") 
alld <- fxn_predictions(effects_all)
pred1 <- alld$A1
pred0 <- alld$A0
diffp <- pred1[,2:8]-pred0[,2:8]
diffp$glottocode <- pred1$glottocode
colMeans(diffp[,1:7]) # differences 
hist(diffp$EY)
abline(v=quantile(diffp$EY, probs = c(0.025, 0.975))) # where 95% of the data is
quantile(diffp$EY, probs = c(0.025, 0.975))

# model 3
effects_all <- marginaleffects::predictions(cmodel3_pc1, newdata = thedata, variables = "Alcohol", type = "response", resp="PolC") 
alld <- fxn_predictions(effects_all)
pred1 <- alld$A1
pred0 <- alld$A0
diffp <- pred1[,2:8]-pred0[,2:8]
diffp$glottocode <- pred1$glottocode
colMeans(diffp[,1:7]) # differences 
hist(diffp$EY)
abline(v=quantile(diffp$EY, probs = c(0.025, 0.975))) # where 95% of the data is
quantile(diffp$EY, probs = c(0.025, 0.975))

# model 2
effects_all <- marginaleffects::predictions(cmodel2, newdata = thedata, variables = "Alcohol", type = "response", resp="PolC") 
alld <- fxn_predictions(effects_all)
pred1 <- alld$A1
pred0 <- alld$A0
diffp <- pred1[,2:8]-pred0[,2:8]
diffp$glottocode <- pred1$glottocode
colMeans(diffp[,1:7]) # differences 
hist(diffp$EY)
abline(v=quantile(diffp$EY, probs = c(0.025, 0.975))) # where 95% of the data is
quantile(diffp$EY, probs = c(0.025, 0.975))

# model 1
effects_all <- marginaleffects::predictions(cmodel1, newdata = thedata, variables = "Alcohol", type = "response", resp="PolC") 
alld <- fxn_predictions(effects_all)
pred1 <- alld$A1
pred0 <- alld$A0
diffp <- pred1[,2:8]-pred0[,2:8]
diffp$glottocode <- pred1$glottocode
colMeans(diffp[,1:7]) # differences 
