source(file = "scripts/00_set_up.R")

################
######## data ## political complexity ~ alcohol + phylogeny + space + env
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
#### model ## Political Complexity ~ agriculture + space + phylogeny + env + alcohol
#############

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

 save(model5_pc1, file=paste0(outfol,"models_cons/cmodel5_pc1_R.rds"))
 capture.output(summary(model5_pc1), file=paste0(outfol,"models_cons/cmodel5_pc1_summary.txt"))

res <- as_draws_df(model5_pc1)
hist(res$b_PolC_Alcohol1, xlab="", main="PolC ~ Alcohol, model5_pc1", breaks=40, col="#eda31d", border="#e9b4a3")
median(res$b_PolC_Alcohol1)
sum(res$b_PolC_Alcohol1>=0)/length(res$b_PolC_Alcohol1) 


#############
#### model ## Political Complexity ~ agriculture + space + phylogeny + env 
#############

m51_pc1 <- brms:: bf(PolC ~ PC1 + mo(Agriculture) + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                     family = cumulative()) 
m52_pc1 <- brms:: bf(Agriculture ~ PC1 + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                     family = cumulative()) 

model51_pc1<- brm (m51_pc1 + m52_pc1 + set_rescor(FALSE), data=thedata, cores=4, chains=2,   iter=8000,control = list(adapt_delta = 0.99),
                   data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                   prior=c(prior(exponential(1), class = sd, resp = Agriculture),
                           prior(exponential(1), class = sd, resp = PolC))) 
summary(model51_pc1)
 save(model51_pc1, file=paste0(outfol,"models_cons/cmodel51R_pc1.rds"))
 capture.output(summary(model51_pc1), file=paste0(outfol,"models_cons/cmodel51_pc1_summary.txt"))
