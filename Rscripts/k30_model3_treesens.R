source(file = "scripts/00_set_up.R")

################
######## data ## political complexity ~ alcohol + phylogeny + space
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

library(brms)
library(parallel)
detectCores()
ncores=15
thedata$Language_ID <- thedata$SCCS_ID
thedata$Language_ID2 <- thedata$SCCS_ID

#checking for duplicated locations and jittering them
duplicate_coords = thedata[duplicated(thedata[,c("lon", "lat")]) | duplicated(thedata[,c("lon", "lat")], fromLast = TRUE),"Language_ID"]
duplicate_rowid = thedata$Language_ID %in% duplicate_coords
thedata$lat[duplicate_rowid] = jitter(thedata$lat[duplicate_rowid], factor = 1)
thedata$lon[duplicate_rowid] = jitter(thedata$lon[duplicate_rowid], factor = 1)

source("scripts/space_correl.R")
spatial_covar_mat <- space_mat(thedata)


#############
#### model ## 3
#############

load(file=paste0(outfol,"treesensitivity/edge_sample50.rds"))
ntrees<-length(strees)

# r<-mclapply(1:ntrees, function(i){
for(i in 1:ntrees){
  thetree <- strees[[i]]
  read.csv(file=paste0(projfol,"EDGEtree/treetips.csv"), stringsAsFactors = F)->treetips
  # prune to matched tips
  tokeep<-thetree$tip.label[thetree$tip.label%in%treetips$edgetips]
  pruned_tree <- ape::keep.tip(thetree, tokeep)
  treetips<-treetips[treetips$edgetips%in%pruned_tree$tip.label,]
  matchtips<-treetips[treetips$edgetips%in%pruned_tree$tip.label,]
  matchtips <- matchtips[match(pruned_tree$tip.label, matchtips$edgetips),]
  pruned_tree$tip.label <- matchtips$alcohol.SCCS_ID
  # pruned to data
  dtips <- pruned_tree$tip.label[!pruned_tree$tip.label%in%thedata$SCCS_ID]
  pruned_tree <- drop.tip(pruned_tree, tip=dtips)
  phylo_covar_mat<- ape::vcv.phylo(pruned_tree, corr=TRUE)
  
  prior <- c(#set_prior("student_t(3, 0, 2.5)", class = "b"),
    set_prior("exponential(1)", class = "sd"))
  
  m3 <- brms:: bf(PolC ~ Alcohol + PC1 + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                  family = cumulative()) 
  m32 <-  brms:: bf(Alcohol ~ PC1+ (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                    family = bernoulli()) 
  
  model3_pc1<- brm (m3 + m32 + set_rescor(FALSE), data=thedata, cores=4, chains=2, 
                data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                iter=7000,control = list(adapt_delta = 0.99),
                prior=c(prior(exponential(1), class = sd, resp = Alcohol),prior(exponential(1), class = sd, resp = PolC))) 
  
  save(model3_pc1, file=paste0(outfol,"treesensitivity/model3/model3R_pc1_",i,".rds"))
  capture.output(summary(model3_pc1), file=paste0(outfol,"treesensitivity/model3/model3_pc1_summary_",i,".txt"))
}
#}, mc.cores=ncores)
