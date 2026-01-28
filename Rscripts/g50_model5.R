source(file = "scripts/00_set_up.R")


# output_files <- c(
# paste0(outfol,"models/model51R_pc1.rds"),
# paste0(outfol,"models/model51_pc1_summary.txt"),
# paste0(outfol,"models/model5_pc1_R.rds"),
# paste0(outfol,"models/model5_pc1_summary.txt"),
# paste0(outfol,"models/model56_pc1.rds"),
# paste0(outfol,"models/model56_pc1_summary.txt"),
# paste0(outfol,"models/model50_envAgr.rds")
#   )
# if(all(file.exists(output_files))){
#   cat(paste0("The output files for this model run already exists.\n", 
#              output_files), ".\n")
# }else{
# script_name <- "g50_model5.R"
# start_time <- Sys.time()



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
#### model ## Political Complexity ~ agriculture + space + phylogeny + env(no alcohol)
#############

library(phylopath)
library(tidyr)
library(tidyverse)
library(tidybayes)

m <- define_model_set( model5 = c(PolC ~ Agriculture + space + phylogeny +PC1,
                                  Agriculture ~ space + phylogeny +PC1))
positions <- data.frame(name = c('Agriculture', 'PolC', 'space', 'phylogeny', "PC1"),
                        x = c(1,2,1.5,1.5, 1), y = c(1,1,2,3,-1))
plot_model_set(m,manual_layout = positions, edge_width = 0.5)

m51_pc1 <- brms:: bf(PolC ~ PC1 + mo(Agriculture) + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                 family = cumulative()) 
m52_pc1 <- brms:: bf(Agriculture ~ PC1 + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                 family = cumulative()) 
model51_pc1<- brm (m51_pc1 + m52_pc1 + set_rescor(FALSE), data=thedata, cores=4, chains=2,   iter=8000,control = list(adapt_delta = 0.99),
               data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
               prior=c(prior(exponential(1), class = sd, resp = Agriculture),
                       prior(exponential(1), class = sd, resp = PolC))) 
summary(model51_pc1)
 save(model51_pc1, file=paste0(outfol,"models/model51R_pc1.rds"))
 capture.output(summary(model51_pc1), file=paste0(outfol,"models/model51_pc1_summary.txt"))


#############
#### model ## Political Complexity ~ alcohol + agriculture + space + phylogeny + env - pc1
#############

m <- define_model_set(model5 = c(PolC ~ Agriculture + Alcohol+ space + phylogeny +PC1,
                                  Agriculture ~ space + phylogeny +PC1,
                                  Alcohol ~ space + phylogeny + PC1,
                                  Alcohol ~ Agriculture))
positions <- data.frame(name = c('Alcohol', 'PolC' ,'Agriculture','space', 'phylogeny', "PC1"),
                        x = c(1,2,2,1.5,1.5, 1), y = c(1,1,0.25,2,3,-1))
plot_model_set(m,manual_layout = positions, edge_width = 0.5)


m56 <- brms:: bf(PolC ~ Alcohol + PC1 + mo(Agriculture) + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                 family = cumulative()) 
m57 <- brms:: bf(Agriculture ~ PC1  + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                 family = cumulative()) 
m58 <- brms:: bf(Alcohol ~ PC1 + mo(Agriculture) + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                 family = bernoulli()) 

model5_pc1<- brm (m56 + m57 + m58+ set_rescor(FALSE), data=thedata, cores=4, chains=2,   iter=8000,control = list(adapt_delta = 0.99),
              data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
              prior=c(prior(exponential(1), class = sd, resp = Agriculture),
                      prior(exponential(1), class = sd, resp = PolC),
                      prior(exponential(1), class = sd, resp = Alcohol))) 
summary(model5_pc1)

 save(model5_pc1, file=paste0(outfol,"models/model5_pc1_R.rds"))
 capture.output(summary(model5_pc1), file=paste0(outfol,"models/model5_pc1_summary.txt"))

res <- as_draws_df(model5_pc1)
hist(res$b_PolC_Alcohol1, xlab="", main="PolC ~ Alcohol, model5_pc1", breaks=40, col="#eda31d", border="#e9b4a3")
median(res$b_PolC_Alcohol1)
sum(res$b_PolC_Alcohol1>=0)/length(res$b_PolC_Alcohol1) 

# for cond effects plots
model56_pc1<- brm (m56, data=thedata, cores=4, chains=2,   iter=8000,control = list(adapt_delta = 0.99),
                   data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                   prior=prior) 
 save(model56_pc1, file=paste0(outfol,"models/model56_pc1.rds"))
 capture.output(summary(model56_pc1), file=paste0(outfol,"models/model56_pc1_summary.txt"))

 #inspect priors
 get_prior(model56_pc1)
 

#############
#### model ## environment and agriculture
#############

mbiv1 <- brm (Agriculture ~ PC1, data=thedata, cores=4, chains=2,
              data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
              control = list(adapt_delta = 0.99), family=cumulative)
 save(mbiv1, file=paste0(outfol,"models/model50_envAgr.rds"))
 
 
 end_time <- Sys.time()
 
 time <- end_time - start_time
 
 cat(
   paste0("I'm done with script ", script_name, ".\n",
          "It took ", round(time[[1]], 2), " ",  units(time)),
   ".\n"
 )
 
 
 ########################
 #### different spaces ## no difference in results - UPDATED 
 ########################
 
 # updated varcovar fxn
 library(tidyverse)
 source("scripts/HS_varcov.spatial.3D.R")
 
 prior <- c(#set_prior("student_t(3, 0, 2.5)", class = "b"),
   set_prior("exponential(1)", class = "sd"))
 
 # gb parameters
 par(mfrow=c(1,3))
 kappa_vals = c(2,2,2.5,5)
 phi_vals = c(1.15,2,3,5)
 
 coords <- thedata[, colnames(thedata)%in% c("lat","lon")]
 rdist.earth_dists <- fields::rdist.earth(x1 = coords, miles = FALSE)
 dimnames(rdist.earth_dists) <- list(thedata$Language_ID2, thedata$Language_ID2)
 
 for(i in 1:length(kappa_vals)){
   spmat = varcov.spatial.3D(thedata[,c("lon", "lat")], cov.pars = c(1,phi_vals[i]), kappa = kappa_vals[i])$varcov
   dimnames(spmat) = list(thedata$Language_ID2, thedata$Language_ID2)
   #plot(spmat ~ rdist.earth_dists, pch=19, col=scales::alpha("black", 0.5),cex=0.5)
   #mtext(paste0("kappa = ", kappa_vals[i], " & phi = ", phi_vals[i]))
   
   m56 <- brms:: bf(PolC ~ Alcohol + PC1 + mo(Agriculture) + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spmat)), 
                    family = cumulative()) 
   m57 <- brms:: bf(Agriculture ~ PC1  + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spmat)), 
                    family = cumulative()) 
   m58 <- brms:: bf(Alcohol ~ PC1 + mo(Agriculture) + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spmat)), 
                    family = bernoulli()) 
   
   m5_3d_gb<- brm (m56 + m57 + m58+ set_rescor(FALSE), data=thedata, cores=4, chains=2,   iter=8000,control = list(adapt_delta = 0.99),
        data2 = list(phylo_covar_mat = phylo_covar_mat, spmat = spmat),
        prior=c(prior(exponential(1), class = sd, resp = Agriculture),
                prior(exponential(1), class = sd, resp = PolC),
                prior(exponential(1), class = sd, resp = Alcohol))) 
   # save model and extract alcohol parameters - mean effect
   save(m5_3d_gb, file=paste0(outfol,"models_GBspace/R5model_gb_",kappa_vals[i],"_",phi_vals[i],".rds"))
   capture.output(summary(m5_3d_gb), file=paste0(outfol,"models_GBspace/R5model_gb_",kappa_vals[i],"_",phi_vals[i],"summary.txt"))
 }
 
 
 
 #########################
 ######## alternatives ### category specific effects 
 #########################
 
 m_seq56_pc1 <-  brm(formula = PolC ~ Alcohol + PC1 + mo(Agriculture) + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                     family = sratio(), data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                     iter=4000,control = list(adapt_delta = 0.99),save_pars = save_pars(all = TRUE))
 m_seq56_pc1_cs <- brm(formula = PolC ~ cs(Alcohol) + PC1 + mo(Agriculture) + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                       family = sratio(), data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                       iter=7000,control = list(adapt_delta = 0.99), save_pars = save_pars(all = TRUE))
 
 #save(m_seq56_pc1, file=paste0(outfol, "models_otherCAT/m_seq56_pc1.rds"))
 #save(m_seq56_pc1_cs, file=paste0(outfol, "models_otherCAT/m_seq56_pc1_cs.rds"))

 modelc <- loo(m_seq56_pc1, m_seq56_pc1_cs,moment_match = T)
 # capture.output(print(modelc), file=paste0(outfol,"models_otherCAT/model5_comp_summary.txt"))
 # save(modelc, file=paste0(outfol, "models_otherCAT/model5_comp.rds"))
 