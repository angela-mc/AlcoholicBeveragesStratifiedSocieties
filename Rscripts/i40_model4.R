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
thedata$Agriculture<-factor(thedata$Agriculture, levels=levels(as.factor(thedata$Agriculture)), ordered=T) # ordered predictor


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

m40 <- brms:: bf(PolC ~ Alcohol + mo(Agriculture) + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)),
                 family = cumulative())
m41 <- brms:: bf(Alcohol ~ mo(Agriculture) + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                 family = bernoulli()) 
  # PC = ordered response so fam = cumulative - we can add probit or no
  # Alcohol = binary
  # Agriculture = ordered predictor so model as monotonic effect, 5 categories so 4 steps


# separate models (just for plotting)
model40<- brm (m40, data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
               iter=4000,control = list(adapt_delta = 0.99)) 
summary(m40)

model41<- brm (m41, data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
               iter=4000,control = list(adapt_delta = 0.99)) 
summary(model41)


# FULL MODEL
model4<- brm (m40 + m41 + set_rescor(FALSE), data=thedata, cores=4, chains=2, data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
              iter=4000,control = list(adapt_delta = 0.99), prior=c(prior(exponential(1), class = sd, resp = Alcohol),
                                                                    prior(exponential(1), class = sd, resp = PolC))) 
summary(model4)
  # intercepts = related to cutpoints
  # SimplexParameters 1 - 4 = δ (StatRe) or ζ (brms)
  # expected difference between categories i and (i-1) = proportion of the overall difference between lowest and highest categories
  # moAgr = b/beta = average change between two adjacent categories of Agr = expected average difference between two adjacent categories of x

 save(model4, file=paste0(outfol,"models_cons/cmodel4R.rds"))
 save(model40, file=paste0(outfol,"models_cons/cmodel40R.rds"))
 save(model41, file=paste0(outfol,"models_cons/cmodel41R.rds"))

load(file=paste0(outfol,"models_cons/cmodel4R.rds"))
load(file=paste0(outfol,"models_cons/cmodel40R.rds"))
load(file=paste0(outfol,"models_cons/cmodel41R.rds"))

# total effect of Agr = moAgriculture * 4 - this is when alcohol is kept constant
res <- as_draws_df(model4)
res$bE <- res$bsp_PolC_moAgriculture*4
median_qi(res$bE, .width = .89)[1:3]

res$bE <- res$bsp_Alcohol_moAgriculture*4
median_qi(res$bE, .width = .89)[1:3]

conditional_effects(model40, "Alcohol", categorical = TRUE)
  # mean estimate of the probability of values in each PolC category for each Alcohol category

# # EXPORT
# capture.output(summary(model4), file=paste0(outfol,"models_cons/cmodel4_summary.txt"))


#############
#### model ## political complexity ~ agriculture alone (no alcohol)
#############

library(brms)
prior <- c(#set_prior("student_t(3, 0, 2.5)", class = "b"),
  set_prior("exponential(1)", class = "sd"))

m42 <- brms:: bf(PolC ~ mo(Agriculture) + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), family = cumulative())
model42<- brm (m42, data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
               iter=7000,control = list(adapt_delta = 0.99)) 
summary(model42)

 save(model42, file=paste0(outfol,"models_cons/cmodel42R.rds"))
 capture.output(summary(model42), file=paste0(outfol,"models_cons/cmodel42_summary.txt"))

load(file=paste0(outfol,"models_cons/cmodel42R.rds"))

# total effect of Agr = moAgriculture * 4 - this is when alcohol is kept constant
res <- as_draws_df(model42)
res$bE <- res$bsp_moAgriculture*4
median_qi(res$bE, .width = .89)[1:3]
