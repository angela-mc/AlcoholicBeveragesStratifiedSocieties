
rm(list=ls())
source(file = "scripts/00_set_up.R")
dir <-  paste0(outfol, "models_otheri")
if(!dir.exists(dir)) dir.create(dir)

library(brms)
library(readxl)

##############
## read data #
##############

otheri<-read_excel(paste0(projfol,"Dataset/OTHER_INTOX_dataset.xlsx")) # Indigenous_intoxicant
table(otheri$Indigenous_intoxicant)

# read the rest of the data
# wide
read.csv(file=paste0(outfol,"workingdata/thedata_wide.csv"), stringsAsFactors = F)->thedata
colnames(thedata)[colnames(thedata)%in%"samplewide"]<-"Alcohol"
thedata$Alcohol<-as.factor(thedata$Alcohol) # binary response
thedata$PolC<-factor(thedata$PolC, levels=levels(as.factor(thedata$PolC)), ordered=T) # ordered reponse

otheri <- otheri[otheri$SCCS_ID%in%thedata$SCCS_ID,]
otheri<-otheri[match(thedata$SCCS_ID, otheri$SCCS_ID),]
thedata$Indigenous_intoxicant <- otheri$Indigenous_intoxicant

thedata <-thedata[!thedata$Indigenous_intoxicant%in%"NA",]
thedata$Indigenous_intoxicant <- as.factor(thedata$Indigenous_intoxicant)


# #############
# #### model ## w/o phy & space
# #############
# 
# # other i -> pol c
# m1 <- brms:: bf(PolC ~ Indigenous_intoxicant, family = cumulative()) 
# model1<- brm (m1, data=thedata, cores=4, chains=2) # no prior
# summary(model1)
# 
# # other i -> alcohol
# m2 <- brms:: bf(Alcohol ~ Indigenous_intoxicant, family = bernoulli()) 
# model2<- brm (m2, data=thedata, cores=4, chains=2) # no prior
# summary(model2)
# 
# # pol c <- alcohol + other i
# m3 <- brms:: bf(PolC ~ Indigenous_intoxicant + Alcohol, family = cumulative())
# model3 <-  brm (m3, data=thedata, cores=4, chains=2) # no prior
# summary(model3)
# 
# model4 <- brm(m3 +m2  + set_rescor(FALSE), data=thedata, cores=4, chains=2)
# summary(model4)


#############
#### model ## with space and phy
#############

thedata$Language_ID <- thedata$SCCS_ID
thedata$Language_ID2 <- thedata$SCCS_ID

#checking for duplicated locations and jittering them
duplicate_coords = thedata[duplicated(thedata[,c("lon", "lat")]) | duplicated(thedata[,c("lon", "lat")], fromLast = TRUE),"Language_ID"]
duplicate_rowid = thedata$Language_ID %in% duplicate_coords
thedata$lat[duplicate_rowid] = jitter(thedata$lat[duplicate_rowid], factor = 1)
thedata$lon[duplicate_rowid] = jitter(thedata$lon[duplicate_rowid], factor = 1)

library(ape)
thetree <- read.nexus(file=paste0(projfol,"EDGEtree/edgetree_sccs.nex"))
dtips <- thetree$tip.label[!thetree$tip.label%in%thedata$SCCS_ID]
thetree <- drop.tip(thetree, tip=dtips)
phylo_covar_mat<- ape::vcv.phylo(thetree, corr=TRUE)
source("scripts/space_correl.R")
spatial_covar_mat <- space_mat(thedata)
prior <- c(#set_prior("student_t(3, 0, 2.5)", class = "b"),
  set_prior("exponential(1)", class = "sd"))


# OTHER I -> POL C
m1_sp <- brms:: bf(PolC ~ Indigenous_intoxicant + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), family = cumulative()) 
model1_sp<- brm (m1_sp, data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                 iter=4000,control = list(adapt_delta = 0.99)) 
summary(model1_sp) # Indigenous_intoxicant1     0.66      0.64    -0.49
res <- as_draws_df(model1_sp)
hist(res$b_Indigenous_intoxicant1, xlab="", main="Other i --> PolC, with space & phylogeny", breaks=40, col="#eda31d", border="#e9b4a3")
median(res$b_Indigenous_intoxicant1)
sum(res$b_Indigenous_intoxicant1>=0)/length(res$b_Indigenous_intoxicant1) # 0.86225
outl <- list(summary(model1_sp),sum(res$b_Indigenous_intoxicant1>=0)/length(res$b_Indigenous_intoxicant1))
names(outl) <- c("model summary", "proportion of posterior coefficient >= 0") 
capture.output(outl, file=paste0(outfol,"models_otheri/modelA_summary.txt"))


# other i -> alcohol
m2_sp <- brms:: bf(Alcohol ~ Indigenous_intoxicant + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)),  family = bernoulli()) 
model2_sp<- brm (m2_sp, data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                 iter=4000,control = list(adapt_delta = 0.99)) 
summary(model2_sp) # Indigenous_intoxicant1     3.01      1.39     0.92 
res <- as_draws_df(model2_sp)
hist(res$b_Indigenous_intoxicant1, xlab="", main="Other i --> Alcohol, with space & phylogeny", breaks=40, col="#eda31d", border="#e9b4a3")
median(res$b_Indigenous_intoxicant1)
sum(res$b_Indigenous_intoxicant1>=0)/length(res$b_Indigenous_intoxicant1) # 0.999
outl <- list(summary(model2_sp),sum(res$b_Indigenous_intoxicant1>=0)/length(res$b_Indigenous_intoxicant1))
names(outl) <- c("model summary", "proportion of posterior coefficient >= 0") 
capture.output(outl, file=paste0(outfol,"models_otheri/modelB_summary.txt"))

# pol c <- alcohol + other i
m3_sp <- brms:: bf(PolC ~ Alcohol + Indigenous_intoxicant + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)),
                 family = cumulative())
model3_sp<- brm (m40, data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
               iter=4000,control = list(adapt_delta = 0.99)) 
summary(model3_sp)

# the dag
model4_sp <- brm(m3_sp +m2_sp  + set_rescor(FALSE), data=thedata,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                 iter=4000,control = list(adapt_delta = 0.99), prior=c(prior(exponential(1), class = sd, resp = Alcohol),
                                                                       prior(exponential(1), class = sd, resp = PolC)))
summary(model4_sp)
res <- as_draws_df(model4_sp)
outl <- list(summary(model4_sp),sum(res$b_PolC_Alcohol1>=0)/length(res$b_PolC_Alcohol1))
names(outl) <- c("model summary", "proportion of posterior coefficient PolC_Alcohol1 >= 0") 
capture.output(outl, file=paste0(outfol,"models_otheri/modelC_summary.txt"))


# alcohol political complexity in this sample
malc <- brms:: bf(PolC ~ Alcohol + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), family = cumulative()) 
model_alc<- brm (malc, data=thedata, cores=4, chains=2, prior = prior,data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                 iter=4000,control = list(adapt_delta = 0.99)) 
summary(model_alc) # 1.29      0.56     0.28 
res <- as_draws_df(model_alc)
hist(res$b_Alcohol1, xlab="", main="Alcohol--> PolC, with space & phylogeny", breaks=40, col="#eda31d", border="#e9b4a3")
median(res$b_Alcohol1)
sum(res$b_Alcohol1>=0)/length(res$b_Alcohol1) 
outl <- list(summary(model_alc),sum(res$b_Alcohol1>=0)/length(res$b_Alcohol1))
names(outl) <- c("model summary", "proportion of posterior coefficient >= 0") 
capture.output(outl, file=paste0(outfol,"models_otheri/model_alc_summary.txt"))


# save all models
save(model1_sp, file=paste0(outfol,"models_otheri/modelA.rds"))
save(model2_sp, file=paste0(outfol,"models_otheri/modelB.rds"))
save(model4_sp, file=paste0(outfol,"models_otheri/modelC.rds"))
save(model_alc, file=paste0(outfol,"models_otheri/model_alc.rds"))


# make figure for paper
source("scripts/h0_fxn.R")
load(file=paste0(outfol,"models_otheri/modelA.rds"))
load(file=paste0(outfol,"models_otheri/modelB.rds"))
load(file=paste0(outfol,"models_otheri/modelC.rds"))
load(file=paste0(outfol,"models_otheri/model_alc.rds"))

pdf(paste0(outfol,"models_otheri/OtherI.pdf"), width=9, height=8)
par(mfrow=c(2,2), mar=c(3,3,5,2),xpd=F)
res <- as_draws_df(model2_sp)
hist(res$b_Indigenous_intoxicant1, xlab="", breaks=60, col="#eda31d", border="#e9b4a3", main="") # Alcohol ~ Indigenous_intoxicant
mediqr (res$b_Indigenous_intoxicant1)
  # median 2.79236, range: 1.548622 -> 4.036099
sum(res$b_Indigenous_intoxicant1>=0)/length(res$b_Indigenous_intoxicant1) # 99%
mtext(side=3,adj=0,text=paste0("a Alcohol ~ Indigenous intoxicant"), font=1)
abline(v=0, lty=2, col="gray70")

res <- as_draws_df(model1_sp)
hist(res$b_Indigenous_intoxicant1, xlab="", breaks=40, col="#eda31d", border="#e9b4a3", main="") # PolC ~ Indigenous_intoxicant
mediqr (res$b_Indigenous_intoxicant1)
  # median 0.6272978, range: 0.02474569 -> 1.22985
sum(res$b_Indigenous_intoxicant1>=0)/length(res$b_Indigenous_intoxicant1) # 0.86225
mtext(side=3,adj=0,text=paste0("b Political complexity ~ Indigenous intoxicant"), font=1)
abline(v=0, lty=2, col="gray70")

crange <- range(c(as_draws_df(model_alc)$b_Alcohol1,as_draws_df(model4_sp)$b_PolC_Alcohol1))
res <- as_draws_df(model_alc)
mediqr (res$b_Alcohol1)
  # median 1.269638, range: 0.7209058 -> 1.818369
sum(res$b_Alcohol1>=0)/length(res$b_Alcohol1) # 0.99525
hist(res$b_Alcohol1, xlab="", main="", breaks=40, col="#eda31d", border="#e9b4a3", xlim=crange) # PolC ~ Alcohol
mtext(side=3,adj=0,text=paste0("c Political complexity ~ Alcohol"), font=1)
abline(v=0, lty=2, col="gray70")

res <- as_draws_df(model4_sp)
hist(res$b_PolC_Alcohol1, xlab="", main="", breaks=60, col="#eda31d", border="#e9b4a3", xlim=crange) # PolC ~ Alcohol + Indigenous_intoxicant
mediqr (res$b_PolC_Alcohol1)
# median 1.255053, range: 0.6812246 -> 1.828882
sum(res$b_PolC_Alcohol1>=0)/length(res$b_PolC_Alcohol1) # 0.989625
mtext(side=3,adj=0,text=paste0("d Political complexity ~ Alcohol", "\n", "                                      + Indigenous intoxicant"), font=1)
abline(v=0, lty=2, col="gray70")
dev.off()

