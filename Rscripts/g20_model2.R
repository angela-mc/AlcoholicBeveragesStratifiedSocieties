source(file = "scripts/00_set_up.R")

# output_files <- c(paste0(outfol,"models/model2_summary.txt"),
#                   paste0(outfol,"models/model2R.rds"))
# 
# if(all(file.exists(output_files))){
#   cat(paste0("The output files for this model run already exists.\n", 
#              output_files), ".\n")
#   
# }else{
#   script_name <- "g20_model2.R"
#   start_time <- Sys.time()
  
  

################
######## data ## political complexity ~ alcohol + phylogeny + space
################


# wide
read.csv(file=paste0(outfol,"workingdata/thedata_wide.csv"), stringsAsFactors = F)->thedata
colnames(thedata)[colnames(thedata)%in%"samplewide"]<-"Alcohol"
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
prior <- c(#set_prior("student_t(3, 0, 2.5)", class = "b"),
           set_prior("exponential(1)", class = "sd"))


#############
#### model ## ordered response ~ binary predictor
#############

library(phylopath)
library(tidyr)
library(tidyverse)

m <- define_model_set( model2 = c(PolC ~ Alcohol + space + phylogeny, Alcohol ~ space + phylogeny))
positions <- data.frame(name = c('Alcohol', 'PolC', 'space', 'phylogeny'), 
                        x = c(1,2,1.5,1.5), y = c(1,1,2,3))
plot_model_set(m,manual_layout = positions, edge_width = 0.5)

m2 <- brms:: bf(PolC ~ Alcohol + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                family = cumulative()) 
  # PC = ordered response so fam = cumulative - we can add probit or no
  # Alcohol = binary

model2<- brm (m2, data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
              iter=4000,control = list(adapt_delta = 0.99)) 
summary(model2)
  # intercepts = related to cutpoints
  # Alcohol1 = effect of alcohol
 save(model2, file=paste0(outfol,"models/model2R.rds"))

library(ggplot2)
load(file=paste0(outfol,"models/model2R.rds"))
ce_data<-conditional_effects(model2, "Alcohol", categorical = TRUE, theme=theme_classic())
  # probability of values in each PolC category for Alcohol absent (0) or present (1) 

res <- as_draws_df(model2)
hist(res$b_Alcohol1, xlab="", main="Alcohol --> PolC, with space & phylogeny", breaks=40, col="#eda31d", border="#e9b4a3")
median(res$b_Alcohol1)

# EXPORT
capture.output(summary(model2), file=paste0(outfol,"models/model2_summary.txt"))

# inspect priors
get_prior(model2)

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
  
  m2_3d_gb <- brm(formula = PolC ~ Alcohol + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spmat)), 
                  family = cumulative(), data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spmat = spmat),
                  iter=5000,control = list(adapt_delta = 0.99))
  # save model and extract alcohol parameters - mean effect
  save(m2_3d_gb, file=paste0(outfol,"models_GBspace/R2model_gb_",kappa_vals[i],"_",phi_vals[i],".rds"))
  capture.output(summary(m2_3d_gb), file=paste0(outfol,"models_GBspace/R2model_gb_",kappa_vals[i],"_",phi_vals[i],"summary.txt"))
}



#########################
######## alternatives ### category specific effects 
#########################

# save pars to be able to compare
m_seq2 <- brm(formula = PolC ~ Alcohol + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
              family = sratio(), data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
              iter=4000,control = list(adapt_delta = 0.99), save_pars = save_pars(all = TRUE))
m_seq2_cs <- brm(formula = PolC ~ cs(Alcohol) + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                 family = sratio(), data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                 iter=4000,control = list(adapt_delta = 0.99), save_pars = save_pars(all = TRUE))
#save(m_seq2, file=paste0(outfol, "models_otherCAT/m_seq2.rds"))
#save(m_seq2_cs, file=paste0(outfol, "models_otherCAT/m_seq2_cs.rds"))

# check for pareto > 0.7
loo_model <- loo(m_seq2, moment_match = T)
pareto_k_values <- loo_model$diagnostics$pareto_k
print(pareto_k_values)
which(pareto_k_values>0.7) # one observation but it is 0.71 

loo_model <- loo(m_seq2_cs, moment_match = T)
pareto_k_values <- loo_model$diagnostics$pareto_k
print(pareto_k_values)
which(pareto_k_values>0.7) # 2 observations only, 0.76 and 0.87  

modelc <- loo(m_seq2, m_seq2_cs,moment_match = T) # use moment_match for comparison to alleviate pareto problems
# save(modelc, file=paste0(outfol, "models_otherCAT/model2_comp.rds"))
# capture.output(print(modelc), file=paste0(outfol,"models_otherCAT/model2_comp_summary.txt"))
