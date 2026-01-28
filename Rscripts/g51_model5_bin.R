source(file = "scripts/00_set_up.R")


# output_files <- c(
#   paste0(outfol,"models/binaryAgriculture/model5_pc1_binR.rds"),
#   paste0(outfol,"models/binaryAgriculture/model5_pc1_bin_summary.txt")
# )
# if(all(file.exists(output_files))){
#   cat(paste0("The output files for this model run already exists.\n", 
#              output_files), ".\n")
# }else{
# script_name <- "g51_model5.R"
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

thedata$Agriculture[thedata$Agriculture==1]<-0 # binarize agriculture
thedata$Agriculture[thedata$Agriculture>0]<-1
thedata$Agriculture<-factor(thedata$Agriculture) # binary predictor


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
#### model ## Political Complexity ~ alcohol + agriculture + space + phylogeny + env
#############

m56 <- brms:: bf(PolC ~ Alcohol + PC1 + Agriculture + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                 family = cumulative()) 
m57 <- brms:: bf(Agriculture ~ PC1  + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                 family = bernoulli()) 
m58 <- brms:: bf(Alcohol ~ PC1 + Agriculture + (1|gr(Language_ID, cov = phylo_covar_mat)) + (1|gr(Language_ID2, cov=spatial_covar_mat)), 
                 family = bernoulli()) 

model5_pc1_bin<- brm (m56 + m57 + m58+ set_rescor(FALSE), data=thedata, cores=4, chains=2,   iter=8000,control = list(adapt_delta = 0.99),
                  data2 = list(phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat),
                  prior=c(prior(exponential(1), class = sd, resp = Agriculture),
                          prior(exponential(1), class = sd, resp = PolC),
                          prior(exponential(1), class = sd, resp = Alcohol))) 
summary(model5_pc1_bin)

 save(model5_pc1_bin, file=paste0(outfol,"models/binaryAgriculture/model5_pc1_binR.rds"))
 capture.output(summary(model5_pc1_bin), file=paste0(outfol,"models/binaryAgriculture/model5_pc1_bin_summary.txt"))

res <- as_draws_df(model5_pc1_bin)
hist(res$b_PolC_Alcohol1, xlab="", main="PolC ~ Alcohol, model5_pc1", breaks=40, col="#eda31d", border="#e9b4a3")
median(res$b_PolC_Alcohol1)
sum(res$b_PolC_Alcohol1>=0)/length(res$b_PolC_Alcohol1) 


end_time <- Sys.time()

time <- end_time - start_time

cat(
  paste0("I'm done with script ", script_name, ".\n",
         "It took ", round(time[[1]], 2), " ",  units(time)),
  ".\n"
)


}
