source(file = "scripts/00_set_up.R")

# output_files <- c(paste0(outfol,"models/model30_envPolC.rds",
#                          paste0(outfol,"models/model30_envAlc.rds"),
#                          paste0(outfol,"models/model32_pc1.rds"),
#                          paste0(outfol,"models/model32_pc1_summary.txt"),
#                          paste0(outfol,"models/model3_pc1R.rds"),
#                          paste0(outfol,"models/model3_pc1_summary.txt"))
#   
# 
# )
# 
# if(all(file.exists(output_files))){
#   cat(paste0("The output files for this model run already exists.\n", 
#              output_files), ".\n")
#   
# }else{
# 
#   script_name <- "g30_model3.R"
#   start_time <- Sys.time()
  


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


#################
###### model 3 ## political complexity ~ alcohol + space + phy + env (pc1)
#################

library(brms)
prior <- c(#set_prior("student_t(3, 0, 2.5)", class = "b"),
  set_prior("exponential(1)", class = "sd"))

library(phylopath)
library(tidyr)
library(tidyverse)
 
m <- define_model_set( model3 = c(PolC ~ Alcohol + space + phylogeny +PC1 ,
                                  Alcohol ~ space + phylogeny +PC1 ))
positions <- data.frame(name = c('Alcohol', 'PolC', 'space', 'phylogeny', "PC1"),
                        x = c(1,2,1.5,1.5, 1), y = c(1,1,2,3,-1))
plot_model_set(m,manual_layout = positions, edge_width = 0.5)

m33 <- brms:: bf(PolC ~ Alcohol + PC1 + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                family = cumulative()) 
m34 <-  brms:: bf(Alcohol ~ PC1 + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                  family = bernoulli()) 

model3_pc1<- brm (m33 + m34 + set_rescor(FALSE), data=thedata, cores=4, chains=2, 
              data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
              iter=7000,control = list(adapt_delta = 0.99),
              prior=c(prior(exponential(1), class = sd, resp = Alcohol),prior(exponential(1), class = sd, resp = PolC))) 
summary(model3_pc1)

 save(model3_pc1, file=paste0(outfol,"models/model3_pc1R.rds"))
 capture.output(summary(model3_pc1), file=paste0(outfol,"models/model3_pc1_summary.txt"))

model32_pc1<- brm (m33, data=thedata, cores=4, chains=2,
                   data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                   iter=7000,control = list(adapt_delta = 0.99),
                   prior=prior)
 save(model32_pc1, file=paste0(outfol,"models/model32_pc1.rds")) # for cond effects plots
 capture.output(summary(model32_pc1), file=paste0(outfol,"models/model32_pc1_summary.txt"))

 #inspect priors
 get_prior(model32_pc1)

###################################
###### environment and variables ##
###################################

mbiv1 <- brm (Alcohol ~ PC1, data=thedata, cores=4, chains=2,
              data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
              control = list(adapt_delta = 0.99), family=bernoulli)
 save(mbiv1, file=paste0(outfol,"models/model30_envAlc.rds"))

mbiv1 <- brm (PolC ~ PC1, data=thedata, cores=4, chains=2,
              data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
              control = list(adapt_delta = 0.99), family=cumulative)
 save(mbiv1, file=paste0(outfol,"models/model30_envPolC.rds"))
 
 
#  end_time <- Sys.time()
#  time <- end_time - start_time
#  cat(
#    paste0("I'm done with script ", script_name, ".\n",
#           "It took ", round(time[[1]], 2), " ",  units(time)),
#    ".\n"
#  )
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
   
   m33 <- brms:: bf(PolC ~ Alcohol + PC1 + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spmat)), 
                    family = cumulative()) 
   m34 <-  brms:: bf(Alcohol ~ PC1 + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spmat)), 
                     family = bernoulli()) 
   
   m3_3d_gb<- brm (m33 + m34 + set_rescor(FALSE), data=thedata, cores=4, chains=2, 
                     data2 = list(phylo_covar_mat = phylo_covar_mat, spmat = spmat),
                     iter=7000,control = list(adapt_delta = 0.99),
                     prior=c(prior(exponential(1), class = sd, resp = Alcohol),prior(exponential(1), class = sd, resp = PolC))) 
   summary(m3_3d_gb)
   # save model and extract alcohol parameters - mean effect
   save(m3_3d_gb, file=paste0(outfol,"models_GBspace/R3model_gb_",kappa_vals[i],"_",phi_vals[i],".rds"))
   capture.output(summary(m3_3d_gb), file=paste0(outfol,"models_GBspace/R3model_gb_",kappa_vals[i],"_",phi_vals[i],"summary.txt"))
 }
 
 

 # save(m3_3d_gb, file=paste0(outfol,"models_GBspace/model3_gb.rds"))
 # capture.output(summary(m3_3d_gb), file=paste0(outfol,"models_GBspace/model3_gb_summary.txt"))
 
 
 #########################
 ######## alternatives ### category specific effects 
 #########################
 
 m_seq32_pc1 <- brm(formula = PolC ~ Alcohol + PC1 + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                    family = sratio(), data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                    iter=7000,control = list(adapt_delta = 0.99), save_pars = save_pars(all = TRUE))
 m_seq32_pc1_cs <- brm(formula = PolC ~ cs(Alcohol) + PC1 + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                       family = sratio(), data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                       iter=7000,control = list(adapt_delta = 0.99), save_pars = save_pars(all = TRUE))
 #save(m_seq32_pc1, file=paste0(outfol, "models_otherCAT/m_seq32_pc1.rds"))
 #save(m_seq32_pc1_cs, file=paste0(outfol, "models_otherCAT/m_seq32_pc1_cs.rds"))
 
 modelc <- loo(m_seq32_pc1, m_seq32_pc1_cs, moment_match = T) 
 # capture.output(print(modelc), file=paste0(outfol,"models_otherCAT/model3_comp_summary.txt"))
 # save(modelc, file=paste0(outfol, "models_otherCAT/model3_comp.rds"))
 
 # check for pareto > 0.7
 loo_model <- loo(m_seq32_pc1, moment_match = T)
 pareto_k_values <- loo_model$diagnostics$pareto_k
 print(pareto_k_values)
 which(pareto_k_values>0.7) # none 
 
 # check for pareto > 0.7
 loo_model <- loo(m_seq32_pc1_cs, moment_match = T)
 pareto_k_values <- loo_model$diagnostics$pareto_k
 print(pareto_k_values)
 which(pareto_k_values>0.7) # one 0.7714494
 
 