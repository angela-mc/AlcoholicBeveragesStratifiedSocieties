source(file = "scripts/00_set_up.R")
library(bayestestR)

############################
###### check convergence ###
############################

# model 2
allf<-list.files(paste0(outfol,"treesensitivity/model2"))
allf<-allf[grep(".rds", allf)]


pdf(paste0(outfol,"treesensitivity/model2/CONV_rhat.pdf"), height=20, width=20)
par(mfrow=c(5,5), mar=c(3,3,6,3))
for(j in 1:length(allf))
  {
    load(paste0(outfol,"treesensitivity/model2/",allf[j]))
   
    # print(summary(model2))
    rstan::check_hmc_diagnostics(model2$fit)
    
    rhat_vals = rhat(model2)
    hist(rhat_vals, main=paste("Model", j))
    essj<-effective_sample(model2)
    col <- ifelse(any(essj$ESS<200), "red","black") 
    mtext(text=paste(essj$ESS, collapse = "\n"),col=col,side=3,adj=1)
    
    print(paste("Model", j))
   
  }
dev.off()
j=49 # There were 1 divergent transitions after warmup. 
j=20 # ESS values < 200


# model 3
allf<-list.files(paste0(outfol,"treesensitivity/model3"))
allf<-allf[grep(".rds", allf)]


pdf(paste0(outfol,"treesensitivity/model3/CONV_rhat.pdf"), height=20, width=20)
par(mfrow=c(5,5), mar=c(3,3,6,3))
for(j in 1:length(allf))
{
  load(paste0(outfol,"treesensitivity/model3/",allf[j]))
  
  # print(summary(model2))
  rstan::check_hmc_diagnostics(model3_pc1$fit)
  
  rhat_vals = rhat(model3_pc1)
  hist(rhat_vals, main=paste("Model", j))
  essj<-effective_sample(model3_pc1)
  col <- ifelse(any(essj$ESS<200), "red","black") 
  mtext(text=paste(essj$ESS, collapse = "\n"),col=col,side=3,adj=1)
  
  print(paste("Model", j))
  
}
dev.off()
j=34 # ESS values < 200 (see the colour red, even if all ESS values did not fit the screen, if any < 200 - col = red)



# model 4
allf<-list.files(paste0(outfol,"treesensitivity/model4"))
allf<-allf[grep(".rds", allf)]


pdf(paste0(outfol,"treesensitivity/model4/CONV_rhat.pdf"), height=20, width=20)
par(mfrow=c(5,5), mar=c(3,3,6,3))
for(j in 1:length(allf))
{
  load(paste0(outfol,"treesensitivity/model4/",allf[j]))
  
  # print(summary(model2))
  rstan::check_hmc_diagnostics(model4$fit)
  
  rhat_vals = rhat(model4)
  hist(rhat_vals, main=paste("Model", j))
  essj<-effective_sample(model4)
  col <- ifelse(any(essj$ESS<200), "red","black") 
  mtext(text=paste(essj$ESS, collapse = "\n"),col=col,side=3,adj=1)
  
  print(paste("Model", j))
  
}
dev.off()
# all j are ok


# model 5
allf<-list.files(paste0(outfol,"treesensitivity/model5"))
allf<-allf[grep(".rds", allf)]


pdf(paste0(outfol,"treesensitivity/model5/CONV_rhat.pdf"), height=20, width=20)
par(mfrow=c(5,5), mar=c(3,3,6,3))
for(j in 1:length(allf))
{
  load(paste0(outfol,"treesensitivity/model5/",allf[j]))
  
  # print(summary(model2))
  rstan::check_hmc_diagnostics(model5_pc1$fit)
  
  rhat_vals = rhat(model5_pc1)
  hist(rhat_vals, main=paste("Model", j))
  essj<-effective_sample(model5_pc1)
  col <- ifelse(any(essj$ESS<200), "red","black") 
  mtext(text=paste(essj$ESS, collapse = "\n"),col=col,side=3,adj=1)
  
  print(paste("Model", j))
  
}
dev.off()
# all j are ok


##################
###### summary ###
##################

# model 2
allf<-list.files(paste0(outfol,"treesensitivity/model2"))
allf<-allf[grep(".rds", allf)]

b_alcoholm2<-list()
for(i in 1:50)
  if(!i%in%c(49,20)) # see above convergence checks
  {
  load(paste0(outfol,"treesensitivity/model2/",allf[i]))
  res <- as_draws_df(model2)
  median(res$b_Alcohol1)
  b_alcoholm2<-c(b_alcoholm2,res$b_Alcohol1)
  print(i)
}
b_alcoholm2<-unlist(b_alcoholm2)
#save(b_alcoholm2, file=paste0(outfol, "treesensitivity/model2_summary.rds"))


# model 3
allf<-list.files(paste0(outfol,"treesensitivity/model3"))
allf<-allf[grep(".rds", allf)]

b_alcoholm3<-list()
for(i in 1:50)
  if(!i%in%c(34)) # see above convergence checks
  {
  load(paste0(outfol,"treesensitivity/model3/",allf[i]))
  res <- as_draws_df(model3_pc1)
  median(res$b_PolC_Alcohol1)
  b_alcoholm3<-c(b_alcoholm3,res$b_PolC_Alcohol1)
  print(i)
}
b_alcoholm3<-unlist(b_alcoholm3)
#save(b_alcoholm3, file=paste0(outfol, "treesensitivity/model3_summary.rds"))


# model 4
allf<-list.files(paste0(outfol,"treesensitivity/model4"))
allf<-allf[grep(".rds", allf)]

b_alcoholm4<-list()
for(i in 1:50){
  load(paste0(outfol,"treesensitivity/model4/",allf[i]))
  res <- as_draws_df(model4)
  median(res$b_PolC_Alcohol1)
  b_alcoholm4<-c(b_alcoholm4,res$b_PolC_Alcohol1)
  print(i)
}
b_alcoholm4<-unlist(b_alcoholm4)
#save(b_alcoholm4, file=paste0(outfol, "treesensitivity/model4_summary.rds"))


# model 5
allf<-list.files(paste0(outfol,"treesensitivity/model5"))
allf<-allf[grep(".rds", allf)]

b_alcoholm5<-list()
for(i in 1:50){
  load(paste0(outfol,"treesensitivity/model5/",allf[i]))
  res <- as_draws_df(model5_pc1)
  median(res$b_PolC_Alcohol1)
  b_alcoholm5<-c(b_alcoholm5,res$b_PolC_Alcohol1)
  print(i)
}
b_alcoholm5<-unlist(b_alcoholm5)
#save(b_alcoholm5, file=paste0(outfol, "treesensitivity/model5_summary.rds"))



#################
###### figure ###
#################

library(ape)
library(brms)
library(tidybayes)

load(file=paste0(outfol,"models/model2R.rds"))
load(file=paste0(outfol,"models/model3_pc1R.rds"))
load(file=paste0(outfol,"models/model4R.rds"))
load(file=paste0(outfol,"models/model5_pc1_R.rds"))

load(file=paste0(outfol, "treesensitivity/model5_summary.rds"))
load(file=paste0(outfol, "treesensitivity/model4_summary.rds"))
load(file=paste0(outfol, "treesensitivity/model3_summary.rds"))
load(file=paste0(outfol, "treesensitivity/model2_summary.rds"))

# each panel separately with a-d - this is only about comparing across 50 trees vs in the mcc
pdf(paste0(outfol,"plots/tree_sensitivity.pdf"), height=6, width=8)
par(mfrow=c(2,2), mar=c(3,2,3,6))
hist(b_alcoholm2, xlab="", main="",breaks=60, col=scales::alpha("gold",1), border=NA,freq = F, yaxt="n", ylab="")
res <- as_draws_df(model2)
hist(res$b_Alcohol1, add=T, col=scales::alpha("royalblue",0.35), border=NA, breaks=40, freq=F)
mtext(side=3,adj=0,text=paste0("a) Model 2"), line=0, font=1,cex = 1.5)

hist(b_alcoholm3, xlab="", main="",breaks=70, col=scales::alpha("gold",1), border=NA,freq = F, yaxt="n", ylab="")
res <- as_draws_df(model3_pc1)
hist(res$b_PolC_Alcohol1, add=T, col=scales::alpha("royalblue",0.35), border=NA, breaks=50, freq=F)
mtext(side=3,adj=0,text=paste0("b) Model 3"), line=0, font=1,cex = 1.5)

hist(b_alcoholm4, xlab="", main="",breaks=60, col=scales::alpha("gold",1), border=NA,freq = F, yaxt="n", ylab="")
res <- as_draws_df(model4)
hist(res$b_PolC_Alcohol1, add=T, col=scales::alpha("royalblue",0.35), border=NA, breaks=40, freq=F)
mtext(side=3,adj=0,text=paste0("c) Model 4"), line=0, font=1,cex = 1.5)

hist(b_alcoholm5, xlab="", main="",breaks=60, col=scales::alpha("gold",1), border=NA,freq = F, yaxt="n", ylab="")
res <- as_draws_df(model5_pc1)
hist(res$b_PolC_Alcohol1, add=T, col=scales::alpha("royalblue",0.35), border=NA, breaks=40, freq=F)
mtext(side=3,adj=0,text=paste0("d) Model 5"), line=0, font=1,cex = 1.5)

par(xpd=T)
legend(x=3,y=0.8, legend=c("Posterior", "MCC"), fill = c("royalblue","gold"),col = c("royalblue","gold"), bty='n',cex=1.5)
dev.off()
