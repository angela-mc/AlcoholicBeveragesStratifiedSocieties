source(file = "scripts/00_set_up.R")


# output_files <- c(
# paste0(outfol,"models/model4R.rds"),
# paste0(outfol,"models/model40R.rds"),
# paste0(outfol,"models/model41R.rds"),
# paste0(outfol,"models/model42R.rds"),
# paste0(outfol,"models/model42_summary.txt")
# )
# if(all(file.exists(output_files))){
#   cat(paste0("The output files for this model run already exists.\n", 
#              output_files), ".\n")
# }else{
#   script_name <- "g40_model4.R"
#   start_time <- Sys.time()
  


################
######## data ## political complexity ~ alcohol + agriculture + phylogeny + space
################


# wide
read.csv(file=paste0(outfol,"workingdata/thedata_wide.csv"), stringsAsFactors = F)->thedata
colnames(thedata)[colnames(thedata)%in%"samplewide"]<-"Alcohol"
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

m <- define_model_set( model3 = c(PolC ~ Alcohol + Agriculture + space + phylogeny, Alcohol ~ space + phylogeny + Agriculture, Agriculture ~ space + phylogeny))
positions <- data.frame(name = c('Alcohol', 'PolC','Agriculture', 'space', 'phylogeny'), 
                        x = c(1,2,1.5,1.5,1.5), y = c(1,1,-1,2,3))
plot_model_set(m,manual_layout = positions, edge_width = 0.5)

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


# separate models (only useful because of conditional effects plotting function)
model40<- brm (m40, data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
              iter=4000,control = list(adapt_delta = 0.99)) 
summary(m40)

model41<- brm (m41, data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
               iter=4000,control = list(adapt_delta = 0.99)) 
summary(model41)


# use one FULL MODEL
model4<- brm (m40 + m41 + set_rescor(FALSE), data=thedata, cores=4, chains=2, data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
              iter=4000,control = list(adapt_delta = 0.99), prior=c(prior(exponential(1), class = sd, resp = Alcohol),
                                                                    prior(exponential(1), class = sd, resp = PolC))) 
summary(model4)
  # intercepts = related to cutpoints
  # SimplexParameters 1 - 4 = δ (StatRe) or ζ (brms)
  # expected difference between categories i and (i-1) = proportion of the overall difference between lowest and highest categories
  # moAgr = b/beta = average change between two adjacent categories of Agr = expected average difference between two adjacent categories of x
  
 save(model4, file=paste0(outfol,"models/model4R.rds"))
 save(model40, file=paste0(outfol,"models/model40R.rds"))
 save(model41, file=paste0(outfol,"models/model41R.rds"))

load(file=paste0(outfol,"models/model4R.rds"))
load(file=paste0(outfol,"models/model40R.rds"))
load(file=paste0(outfol,"models/model41R.rds"))

#inspect priors
get_prior(model4)


# total effect of Agr = moAgriculture * 4 - this is when alcohol is kept constant
res <- as_draws_df(model4)
res$bE <- res$bsp_PolC_moAgriculture*4
median_qi(res$bE, .width = .89)[1:3]

res$bE <- res$bsp_Alcohol_moAgriculture*4
median_qi(res$bE, .width = .89)[1:3]

# conditional_effects(model3) # Predictions are treated as continuous variables in 'conditional_effects' by default which is likely invalid for ordinal families. Please set 'categorical' to TRUE. 
# conditional_effects(model3, "Agriculture", categorical = TRUE) # Error: Argument 'categorical' may only be used for categorical or ordinal models.
  # because model model41 is binary - build this plot on model40 only

conditional_effects(model40, "Alcohol", categorical = TRUE)
  # mean estimate of the probability of values in each PolC category for each Alcohol category

res <- as_draws_df(model4)
hist(res$b_PolC_Alcohol1, xlab="", main="Alcohol --> PolC, with agr, space & phylogeny", breaks=40, col="#eda31d", border="#e9b4a3")
median(res$b_PolC_Alcohol1)
sum(res$b_PolC_Alcohol1>=0)/length(res$b_PolC_Alcohol1)
hist(res$bsp_Alcohol_moAgriculture, xlab="", main="", breaks=40, col="#eda31d", border="#e9b4a3")
hist(res$bsp_PolC_moAgriculture, xlab="", main="", breaks=40, col="#eda31d", border="#e9b4a3")

# # EXPORT
# capture.output(summary(model4), file=paste0(outfol,"models/model4_summary.txt"))

 

#############
#### model ## political complexity ~ agriculture alone
#############

library(brms)
prior <- c(#set_prior("student_t(3, 0, 2.5)", class = "b"),
  set_prior("exponential(1)", class = "sd"))

m42 <- brms:: bf(PolC ~ mo(Agriculture) + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), family = cumulative())
model42<- brm (m42, data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
               iter=7000,control = list(adapt_delta = 0.99)) 
summary(model42)

 save(model42, file=paste0(outfol,"models/model42R.rds"))
 capture.output(summary(model42), file=paste0(outfol,"models/model42_summary.txt"))

load(file=paste0(outfol,"models/model42R.rds"))
# total effect of Agr = moAgriculture * 4 
res <- as_draws_df(model42)
res$bE <- res$bsp_moAgriculture*4
median_qi(res$bE, .width = .89)[1:3]


# end_time <- Sys.time()
# time <- end_time - start_time
# cat(
#   paste0("I'm done with script ", script_name, ".\n",
#          "It took ", round(time[[1]], 2), " ",  units(time)),
#   ".\n"
# )
# }


########################
#### different spaces ## no difference in results
########################

# updated varcovar fxn
library(tidyverse)
source("scripts/HS_varcov.spatial.3D.R")

coords <- thedata[, colnames(thedata)%in% c("lat","lon")]
rdist.earth_dists <- fields::rdist.earth(x1 = coords, miles = FALSE)
dimnames(rdist.earth_dists) <- list(thedata$Language_ID2, thedata$Language_ID2)

# gb parameters
par(mfrow=c(1,3))
kappa_vals = c(2,2,2.5,5)
phi_vals = c(1.15,2,3,5)

for(i in 1:length(kappa_vals)){
  spmat = varcov.spatial.3D(thedata[,c("lon", "lat")], cov.pars = c(1,phi_vals[i]), kappa = kappa_vals[i])$varcov
  dimnames(spmat) = list(thedata$Language_ID2, thedata$Language_ID2)
  #plot(spmat ~ rdist.earth_dists, pch=19, col=scales::alpha("black", 0.5),cex=0.5)
  #mtext(paste0("kappa = ", kappa_vals[i], " & phi = ", phi_vals[i]))
  
  m40 <- brms:: bf(PolC ~ Alcohol + mo(Agriculture) + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spmat)),
                   family = cumulative())
  m41 <- brms:: bf(Alcohol ~ mo(Agriculture) + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spmat)), 
                   family = bernoulli()) 
  
  m4_3d_gb<- brm (m40 + m41 + set_rescor(FALSE), data=thedata, cores=4, chains=2, data2 = list(phylo_covar_mat = phylo_covar_mat, spmat = spmat),
                  iter=5000,control = list(adapt_delta = 0.99), prior=c(prior(exponential(1), class = sd, resp = Alcohol),
                                                                        prior(exponential(1), class = sd, resp = PolC))) 
  
  # save model and extract alcohol parameters - mean effect
  save(m4_3d_gb, file=paste0(outfol,"models_GBspace/R4model_gb_",kappa_vals[i],"_",phi_vals[i],".rds"))
  capture.output(summary(m4_3d_gb), file=paste0(outfol,"models_GBspace/R4model_gb_",kappa_vals[i],"_",phi_vals[i],"summary.txt"))
}



#########################
######## alternatives ### category specific effects 
#########################

m_seq40 <-  brm(formula = PolC ~ Alcohol + mo(Agriculture) + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                family = sratio(), data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                iter=4000,control = list(adapt_delta = 0.99), save_pars = save_pars(all = TRUE))
m_seq40_cs <- brm(formula = PolC ~ cs(Alcohol) + mo(Agriculture) + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                  family = sratio(), data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                  iter=7000,control = list(adapt_delta = 0.99), save_pars = save_pars(all = TRUE))

#save(m_seq40, file=paste0(outfol, "models_otherCAT/m_seq40.rds"))
#save(m_seq40_cs, file=paste0(outfol, "models_otherCAT/m_seq40_cs.rds"))
modelc <- loo(m_seq40, m_seq40_cs, moment_match = T)
# save(modelc, file=paste0(outfol, "models_otherCAT/model4_comp.rds"))
# capture.output(print(modelc), file=paste0(outfol,"models_otherCAT/model4_comp_summary.txt"))

# check for pareto > 0.7
loo_model <- loo(m_seq40, moment_match = T)
pareto_k_values <- loo_model$diagnostics$pareto_k
print(pareto_k_values)
which(pareto_k_values>0.7) # one 0.7838113
