source(file = "scripts/00_set_up.R")

################
######## data ## political complexity ~ alcohol + phylogeny + space
################

# conservative
read.csv(file=paste0(outfol,"workingdata/thedata_cons.csv"), stringsAsFactors = F)->thedata
colnames(thedata)[colnames(thedata)%in%"samplecons"]<-"Alcohol"
table(thedata$Alcohol)

# thedata$Alcohol<-thedata$Alcohol+1
thedata$Alcohol<-as.factor(thedata$Alcohol) # binary predictor
thedata$PolC<-factor(thedata$PolC, levels=levels(as.factor(thedata$PolC)), ordered=T) # ordered response


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

library(brms)
prior <- c(#set_prior("student_t(3, 0, 2.5)", class = "b"),
  set_prior("exponential(1)", class = "sd"))


#############
#### model ## ordered response ~ binary predictor
#############

library(phylopath)
library(tidyr)
library(tidyverse)

m2 <- brms:: bf(PolC ~ Alcohol + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                family = cumulative()) 
model2<- brm (m2, data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
              iter=4000,control = list(adapt_delta = 0.99)) 
summary(model2)
# intercepts = related to cutpoints
# Alcohol1 = effect of alcohol
 save(model2, file=paste0(outfol,"models_cons/cmodel2R.rds"))

load(file=paste0(outfol,"models_cons/cmodel2R.rds"))
conditional_effects(model2, "Alcohol", categorical = TRUE)

res <- as_draws_df(model2)
hist(res$b_Alcohol1, xlab="", main="Alcohol --> PolC, with space & phylogeny", breaks=40, col="#eda31d", border="#e9b4a3")
median(res$b_Alcohol1)

# EXPORT
capture.output(summary(model2), file=paste0(outfol,"models_cons/cmodel2_summary.txt"))

