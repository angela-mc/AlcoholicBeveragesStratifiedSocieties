source(file = "scripts/00_set_up.R")

################
######## data ## political complexity ~ alcohol + agriculture + phylogeny + space
################


# conservative
read.csv(file=paste0(outfol,"workingdata/thedata_cons.csv"), stringsAsFactors = F)->thedata
colnames(thedata)[colnames(thedata)%in%"samplecons"]<-"Alcohol"
table(thedata$Alcohol)

thedata$Alcohol<-as.factor(thedata$Alcohol) # binary predictor
thedata$PolC<-factor(thedata$PolC, levels=levels(as.factor(thedata$PolC)), ordered=T) # ordered response
thedata$Agriculture[thedata$Agriculture==1]<-0 # binarize agriculture
thedata$Agriculture[thedata$Agriculture>0]<-1
thedata$Agriculture<-factor(thedata$Agriculture) # binary predictor
table(thedata$Agriculture)


###########################
#### space and phylogeny ##
###########################

thedata$Language_ID <- thedata$SCCS_ID
thedata$Language_ID2 <- thedata$SCCS_ID

#checking for duplicated locations and jittering them
duplicate_coords = thedata[duplicated(thedata[,c("lon", "lat")]) | duplicated(thedata[,c("lon", "lat")], fromLast = TRUE),"Language_ID"]
duplicate_rowid = thedata$Language_ID %in% duplicate_coords
thedata$lat[duplicate_rowid] = jitter(thedata$lat[duplicate_rowid], factor = 1)
thedata$lon[duplicate_rowid] = jitter(thedata$lon[duplicate_rowid], factor = 1)

thetree <- read.nexus(file=paste0(projfol,"EDGEtree/edgetree_sccs.nex"))
dtips <- thetree$tip.label[!thetree$tip.label%in%thedata$SCCS_ID]
thetree <- drop.tip(thetree, tip=dtips)
phylo_covar_mat<- ape::vcv.phylo(thetree, corr=TRUE)

source("scripts/space_correl.R")
spatial_covar_mat <- space_mat(thedata)



#############
#### model ## political complexity ~ alcohol + agriculture + space + phylogeny
#############

library(phylopath)
library(tidyr)
library(tidyverse)
library(tidybayes)
library(brms)
prior <- c(#set_prior("student_t(3, 0, 2.5)", class = "b"),
  set_prior("exponential(1)", class = "sd"))

m40 <- brms:: bf(PolC ~ Alcohol + Agriculture + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)),
                 family = cumulative())
m41 <- brms:: bf(Alcohol ~ Agriculture + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                 family = bernoulli()) 

model4_bin<- brm (m40 + m41 + set_rescor(FALSE), data=thedata, cores=4, chains=2, data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                  iter=4000,control = list(adapt_delta = 0.99), prior=c(prior(exponential(1), class = sd, resp = Alcohol),
                                                                        prior(exponential(1), class = sd, resp = PolC))) 
summary(model4_bin)
 save(model4_bin, file=paste0(outfol,"models_cons/binaryAgriculture/cmodel4R_bin.rds"))

# separate just to be able to plot conditional effects
model40_bin<- brm (m40 ,data=thedata, cores=4, chains=2, data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                   iter=7000,control = list(adapt_delta = 0.99), prior=prior) 
 save(model40_bin, file=paste0(outfol,"models_cons/binaryAgriculture/cmodel40R_bin.rds"))
conditional_effects(model40_bin, "Alcohol", categorical = TRUE)

load(file=paste0(outfol,"models_cons/binaryAgriculture/cmodel4R_bin.rds"))
res <- as_draws_df(model4_bin)
hist(res$b_PolC_Alcohol1, xlab="", main="Alcohol --> PolC, with agr, space & phylogeny", breaks=40, col="#eda31d", border="#e9b4a3")
median(res$b_PolC_Alcohol1)
sum(res$b_PolC_Alcohol1>=0)/length(res$b_PolC_Alcohol1)
hist(res$b_Alcohol_Agriculture1, xlab="", main="", breaks=40, col="#eda31d", border="#e9b4a3")

# # EXPORT
# capture.output(summary(model4_bin), file=paste0(outfol,"models_cons/binaryAgriculture/cmodel4_bin_summary.txt"))


#############
#### model ## political complexity ~ agriculture alone
#############

prior <- c(#set_prior("student_t(3, 0, 2.5)", class = "b"),
  set_prior("exponential(1)", class = "sd"))

m42 <- brms:: bf(PolC ~ Agriculture + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), family = cumulative())
model42_bin<- brm (m42, data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                   iter=7000,control = list(adapt_delta = 0.99)) 
summary(model42_bin)

 save(model42_bin, file=paste0(outfol,"models_cons/binaryAgriculture/cmodel42R_bin.rds"))
 capture.output(summary(model42_bin), file=paste0(outfol,"models_cons/binaryAgriculture/cmodel42_bin_summary.txt"))

load(file=paste0(outfol,"models_cons/binaryAgriculture/cmodel42R_bin.rds"))
conditional_effects(model42_bin, "Agriculture", categorical = TRUE)
