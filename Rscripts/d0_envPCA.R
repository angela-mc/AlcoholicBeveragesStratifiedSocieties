source(file = "scripts/00_set_up.R")

dir <-  paste0(outfol, "envpca/")
if(!dir.exists(dir)){
  dir.create(dir)
}


######################################
######## ecoclimate folder dplace #### # and Moderate Resolution Imaging Spectroradiometer data (MODIS) & GMTED2010 (Global Multi-resolution Terrain Elevation Data 2010)
######################################

ecoclimate_folder<-paste0(projfol,"dplace-data-master/datasets/ecoClimate/")
read.csv(paste0(ecoclimate_folder,"data.csv"), stringsAsFactors = F)-> ECdata
read.csv(paste0(ecoclimate_folder,"variables.csv"), stringsAsFactors = F)-> ECvar

modis_folder<-paste0(projfol,"dplace-data-master/datasets/MODIS/")
read.csv(paste0(modis_folder,"data.csv"), stringsAsFactors = F)-> MODISdata
read.csv(paste0(modis_folder,"variables.csv"), stringsAsFactors = F)-> MODISvar

elev_folder<-paste0(projfol,"dplace-data-master/datasets/GMTED2010/")
read.csv(paste0(elev_folder,"data.csv"), stringsAsFactors = F)-> ELEVdata
read.csv(paste0(elev_folder,"variables.csv"), stringsAsFactors = F)-> ELEVvar

# subset to SCCS & relevant variables
socs<-paste0("SCCS", 1:186)
sum(socs%in%MODISdata$soc_id) # SCCS108 missing
sum(socs%in%ECdata$soc_id)
sum(socs%in%ELEVdata$soc_id)

ECdata<-ECdata[grepl("SCCS",ECdata$soc_id),]
MODISdata<-MODISdata[grepl("SCCS",MODISdata$soc_id),]
ELEVdata<-ELEVdata[grepl("SCCS",ELEVdata$soc_id),]


#################################
############ build dataset ######
#################################

alcohol<-read.csv(paste0(projfol,"Dataset/ALCOHOL_dataset.csv"),stringsAsFactors = F) 
envdata<-data.frame(alcohol$SCCS_ID)

envdata$AnnualMeanTemperature<- -999
envdata$AnnualTemperatureVariance<- -999
envdata$TemperaturePredictability<- -999

envdata$MonthlyMeanPrecipitation<- -999
envdata$PrecipitationPredictability<- -999
envdata$AnnualPrecipitationVariance<- -999

envdata$MonthlyMeanNetPrimaryProduction<- -999
envdata$AnnualNetPrimaryProductionVariance<- -999
envdata$NetPrimaryProductionPredictability<- -999

envdata$Elev<- -999
envdata$Slope<- -999

for(i in 1:length(envdata[,1])){
  envdata$AnnualMeanTemperature[i]<-ECdata[ECdata$var_id%in%"AnnualMeanTemperature" & ECdata$soc_id%in%envdata$alcohol.SCCS_ID[i],]$code
  envdata$AnnualTemperatureVariance[i]<-ECdata[ECdata$var_id%in%"AnnualTemperatureVariance" & ECdata$soc_id%in%envdata$alcohol.SCCS_ID[i],]$code
  envdata$TemperaturePredictability[i]<-ECdata[ECdata$var_id%in%"TemperaturePredictability" & ECdata$soc_id%in%envdata$alcohol.SCCS_ID[i],]$code
  envdata$MonthlyMeanPrecipitation[i]<-ECdata[ECdata$var_id%in%"MonthlyMeanPrecipitation" & ECdata$soc_id%in%envdata$alcohol.SCCS_ID[i],]$code
  envdata$PrecipitationPredictability[i]<-ECdata[ECdata$var_id%in%"PrecipitationPredictability" & ECdata$soc_id%in%envdata$alcohol.SCCS_ID[i],]$code
  envdata$AnnualPrecipitationVariance[i]<-ECdata[ECdata$var_id%in%"AnnualPrecipitationVariance" & ECdata$soc_id%in%envdata$alcohol.SCCS_ID[i],]$code
  
  if(envdata$alcohol.SCCS_ID[i]%in%MODISdata$soc_id){
    envdata$MonthlyMeanNetPrimaryProduction[i]<-MODISdata[MODISdata$var_id%in%"MonthlyMeanNetPrimaryProduction" & MODISdata$soc_id%in%envdata$alcohol.SCCS_ID[i],]$code
    envdata$AnnualNetPrimaryProductionVariance[i]<-MODISdata[MODISdata$var_id%in%"AnnualNetPrimaryProductionVariance" & MODISdata$soc_id%in%envdata$alcohol.SCCS_ID[i],]$code
    envdata$NetPrimaryProductionPredictability[i]<-MODISdata[MODISdata$var_id%in%"NetPrimaryProductionPredictability" & MODISdata$soc_id%in%envdata$alcohol.SCCS_ID[i],]$code
  }
  
  envdata$Elev[i]<-ELEVdata[ELEVdata$var_id%in%"Elevation" & ELEVdata$soc_id%in%envdata$alcohol.SCCS_ID[i],]$code
  envdata$Slope[i]<-ELEVdata[ELEVdata$var_id%in%"Slope" & ELEVdata$soc_id%in%envdata$alcohol.SCCS_ID[i],]$code

  if(verbose == TRUE){
  print(i)}
}

envdata[envdata$MonthlyMeanNetPrimaryProduction%in%-999,]$MonthlyMeanNetPrimaryProduction<-NA
envdata[envdata$AnnualNetPrimaryProductionVariance%in%-999,]$AnnualNetPrimaryProductionVariance<-NA
envdata[envdata$NetPrimaryProductionPredictability%in%-999,]$NetPrimaryProductionPredictability<-NA
envdata<-envdata[complete.cases(envdata),]

# normalize and scale variables (following Haynie et al)
envdata$SqrtMeanP <- scale(sqrt(envdata$MonthlyMeanPrecipitation-min(envdata$MonthlyMeanPrecipitation,na.rm=T)+1))
envdata$LogVarP <- scale(log(envdata$AnnualPrecipitationVariance-min(envdata$AnnualPrecipitationVariance,na.rm=T)+1))
envdata$PredictP <- scale(envdata$PrecipitationPredictability)

envdata$MeanT <- scale(envdata$AnnualMeanTemperature)
envdata$LogVarT <- scale(log(envdata$AnnualTemperatureVariance - min(envdata$AnnualTemperatureVariance,na.rm=T) + 1))
envdata$PredictT <- scale(envdata$TemperaturePredictability)

envdata$SqrtMeanN <- scale(sqrt(envdata$MonthlyMeanNetPrimaryProduction-min(envdata$MonthlyMeanNetPrimaryProduction,na.rm=T)+1))
envdata$LogVarN <- scale(log(envdata$AnnualNetPrimaryProductionVariance - min(envdata$AnnualNetPrimaryProductionVariance,na.rm=T) + 1))
envdata$PredictN <- scale(envdata$NetPrimaryProductionPredictability)

envdata$LogElev <- scale(log(envdata$Elev - min(envdata$Elev,na.rm=T) + 1))
envdata$LogSlope <- scale(log(envdata$Slope - min(envdata$Slope) + 1))


####################
############ PCA ###
####################

library(psych)
library(nFactors)

PCAdata <- envdata[,c('SqrtMeanP', 'LogVarP','PredictP','MeanT','LogVarT','PredictT','SqrtMeanN', 'PredictN', 'LogVarN', 'LogSlope', 'LogElev')]
pairs.panels(PCAdata, 
             method = "pearson", # correlation method
             hist.col = "#a3afd1",# "#a9d1a3","",""),
             density = TRUE,  # show density plots
             ellipses = F, # show correlation ellipses
             cex.labels= 0.75,cor=T,lm=T, ci = T, cex.cor = 0.9,stars = T) # smoother= T,

# Determine number of factors
ev <- eigen(cor(PCAdata))
ap <- parallel(subject=nrow(na.omit(PCAdata)),var=ncol(na.omit(PCAdata)),rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS) 
ev # suggests 3 based on values > 1

myPCA <- principal(PCAdata, nfactors = 3, rotate = "varimax", scores = T)
capture.output(myPCA, file=paste0(outfol,"envpca/mypca.txt"))
myPCA
  # PC1 = +meanP +varP +meanT -varT +PredictT +meanN -predictN : PC1 describes a gradient including variables associated with environmental productivity. Higher values of PC1 are associated with predictable and invariably warmer temperatures, seasonal and high levels of precipitation and high levels of NPP
  # PC2 = +predictP +varN : increasingly predictable and seasonal environments, including NPP variance and precipitation predictability
  # PC3 = +Elev +Slope : C3 describes a gradient including increasing slope and elevation. Higher values of topographic complexity may increase the patchiness of resources and may thus also contribute to defensibility (see Discussion for details).

# Add PCA scores to data
colnames(myPCA$scores)[which(colnames(myPCA$scores)=='RC1')] <- 'PC1'
colnames(myPCA$scores)[which(colnames(myPCA$scores)=='RC2')] <- 'PC2'
colnames(myPCA$scores)[which(colnames(myPCA$scores)=='RC3')] <- 'PC3'
envdata <- cbind(envdata, myPCA$scores) 

write.csv(envdata, file=paste0(outfol,"envpca/envdata.csv"))



