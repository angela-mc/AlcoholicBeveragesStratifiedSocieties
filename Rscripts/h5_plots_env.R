source(file = "scripts/00_set_up.R")

source("scripts/h0_fxn.R")

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
table(thedata$PolC)


#########################################
######## environment and our variables ## for extended data 6
#########################################

load(file=paste0(outfol,"models/model30_envAlc.rds")) ; mbiv1 -> mod_alc; rm(mbiv1)
load(file=paste0(outfol,"models/model30_envPolC.rds")) ; mbiv1 -> mod_polc; rm(mbiv1)
load(file=paste0(outfol,"models/model50_envAgr.rds")) ; mbiv1 -> mod_agr; rm(mbiv1)

allrange1<-min(c(as_draws_df(mod_alc)$b_PC1, as_draws_df(mod_polc)$b_PC1,as_draws_df(mod_agr)$b_PC1))
allrange2<-max(c(as_draws_df(mod_alc)$b_PC1, as_draws_df(mod_polc)$b_PC1,as_draws_df(mod_agr)$b_PC1))
tylim<-c(allrange1-0.001, allrange2+0.001)

#pdf(file=paste0(outfol, "plots/env_variables.pdf"), width=14, height=5)
par(mfrow=c(1,3), mar=c(3,5,3,1))
res <- as_draws_df(mod_alc)
hist(res$b_PC1, xlab="", main="", breaks=40, col="#e3556e", border="#d64f6d", xlim=tylim, cex.lab=2,cex.axis=2)
mtext("a) Alcohol", adj=0,cex=2)
median(res$b_PC1)
abline(v=0, lty=2, col="gray70")
sum(res$b_PC1>=0)/length(res$b_PC1)

res <- as_draws_df(mod_agr)
hist(res$b_PC1, xlab="", main="", breaks=40, col="#e3556e", border="#d64f6d", xlim=tylim, cex.lab=2,cex.axis=2)
mtext("b) Agriculture ", adj=0,cex=2)
median(res$b_PC1)
abline(v=0, lty=2, col="gray70")
sum(res$b_PC1>=0)/length(res$b_PC1)

res <- as_draws_df(mod_polc)
hist(res$b_PC1, xlab="", main="", breaks=40, col="#e3556e", border="#d64f6d", xlim=tylim, cex.lab=2,cex.axis=2)
mtext("c) Political Complexity", adj=0,cex=2)
median(res$b_PC1)
abline(v=0, lty=2, col="gray70")
sum(res$b_PC1<=0)/length(res$b_PC1)

dev.off()



