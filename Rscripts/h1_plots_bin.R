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

thedata$Agriculture[thedata$Agriculture==1]<-0 # binarize agriculture
thedata$Agriculture[thedata$Agriculture>0]<-1
thedata$Agriculture<-factor(thedata$Agriculture) # binary predictor



#########################
#### alcohol comparison # 
#########################

load(file=paste0(outfol,"models/model4R.rds"))
load(file=paste0(outfol,"models/model5_pc1_R.rds"))

load(file=paste0(outfol,"models/binaryAgriculture/model4R_bin.rds"))
load(file=paste0(outfol,"models/binaryAgriculture/model5_pc1_binR.rds"))


pdf(paste0(outfol,"plots/alcohol_onPC_binaryAgricultureComparison.pdf"), height=8, width=7)
par(mfrow=c(2,1), mar=c(3,3,1.5,0),xpd=T)

allrange1<-min(c(as_draws_df(model4)$b_PolC_Alcohol1, as_draws_df(model5_pc1)$b_PolC_Alcohol1,as_draws_df(model4_bin)$b_PolC_Alcohol1, as_draws_df(model5_pc1_bin)$b_PolC_Alcohol1))
allrange2<-max(c(as_draws_df(model4)$b_PolC_Alcohol1, as_draws_df(model5_pc1)$b_PolC_Alcohol1,as_draws_df(model4_bin)$b_PolC_Alcohol1, as_draws_df(model5_pc1_bin)$b_PolC_Alcohol1))
tylim<-c(allrange1-2, allrange2+1)

res <- as_draws_df(model4)
mstats1 <- mediqr (res$b_PolC_Alcohol1)
hist(res$b_PolC_Alcohol1, xlab="", main="",breaks=40, col=scales::alpha("#eda31d",0.5), border=scales::alpha("#eda31d",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim)

res <- as_draws_df(model4_bin)
mstats2 <- mediqr (res$b_PolC_Alcohol1)
hist(res$b_PolC_Alcohol1, xlab="", main="",breaks=60,col=scales::alpha("#618b44",0.5), border=scales::alpha("#618b44",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim,add=T)
par(xpd=F)
abline(v=0, lty=2, col="gray70")
par(xpd=T)

points(x=c(mstats2[2], mstats2[3]),y=c(-0.06,-0.06), type='l', lwd=4, col=scales::alpha("#618b44",0.75))
points(x=mstats2[1],y=c(-0.06), pch=18,  col="#618b44",cex=2)
points(x=c(mstats1[2], mstats1[3]),y=c(-0.03,-0.03), type='l', lwd=4, col=scales::alpha("#eda31d",0.75))
points(x=mstats1[1],y=c(-0.03), pch=18,  col="#eda31d",cex=2)
legend(x=4,y=1, fill=c("#eda31d","#618b44"), legend=c("ordinal", "binary"), border=NA, bty = "n")
#mtext(side=3,adj=0,text=paste0("Model4", "\n", "PC ~ Alcohol + Agri", "\n", "phylo-space"), line=-2, font=2)
mtext(side=3,adj=0,text=paste0("a) Model 4 - wide"), line=-1, font=2)
mtext(side=3,adj=0,text=paste0("Political Complexity", " ~", "\n","              ", "Alcohol", "\n","              ","Agriculture","\n","              ", "phylo-space" ), line=-5, font=1)

res <- as_draws_df(model5_pc1)
mstats3 <- mediqr (res$b_PolC_Alcohol1)
hist(res$b_PolC_Alcohol1, xlab="", main="",breaks=60,col=scales::alpha("#eda31d",0.5), border=scales::alpha("#eda31d",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim)
par(xpd=F)
abline(v=0, lty=2, col="gray70")
par(xpd=T)

res <- as_draws_df(model5_pc1_bin)
mstats4 <- mediqr (res$b_PolC_Alcohol1)
hist(res$b_PolC_Alcohol1, xlab="", main="",breaks=70,col=scales::alpha("#618b44",0.5), border=scales::alpha("#618b44",0.25),freq=F, yaxt="n", ylab="", xaxt="n", xlim=tylim,add=T)

points(x=c(mstats4[2], mstats4[3]),y=c(-0.06,-0.06), type='l', lwd=4, col=scales::alpha("#618b44",0.75))
points(x=mstats4[1],y=c(-0.06), pch=18,  col="#618b44",cex=2)
points(x=c(mstats3[2], mstats3[3]),y=c(-0.03,-0.03), type='l', lwd=4, col=scales::alpha("#eda31d",0.75))
points(x=mstats3[1],y=c(-0.03), pch=18,  col="#eda31d",cex=2)
mtext(side=3,adj=0,text=paste0("b) Model 5 - wide"), line=-0.6, font=2)
mtext(side=3,adj=0,text=paste0("Political Complexity", " ~", "\n","              ", "Alcohol", "\n","              ","Agriculture","\n","              ","Env-productivity","\n","              ","phylo-space" ), line=-6, font=1)

par(xpd=F)
axis(side=1, tck=0.002, lwd.ticks = 2, line = +0.9, lwd=1.5, cex.axis=1.5,xlim=tylim, at=c(-1,0,1,2,3,4,5,6))

dev.off()


#########################
#### alcohol comparison # conservative EXT DATA fig s5
#########################

load(file=paste0(outfol,"models_cons/cmodel4R.rds"))
load(file=paste0(outfol,"models_cons/cmodel5_pc1_R.rds"))

load(file=paste0(outfol,"models_cons/binaryAgriculture/cmodel4R_bin.rds"))
load(file=paste0(outfol,"models_cons/binaryAgriculture/cmodel5_pc1_binR.rds"))


#pdf(paste0(outfol,"plots/alcohol_onPC_binaryAgricultureComparison_conservative.pdf"), height=8, width=7)
par(mfrow=c(2,1), mar=c(3,3,1.5,0),xpd=T)

allrange1<-min(c(as_draws_df(model4)$b_PolC_Alcohol1, as_draws_df(model5_pc1)$b_PolC_Alcohol1,as_draws_df(model4_bin)$b_PolC_Alcohol1, as_draws_df(model5_pc1_bin)$b_PolC_Alcohol1))
allrange2<-max(c(as_draws_df(model4)$b_PolC_Alcohol1, as_draws_df(model5_pc1)$b_PolC_Alcohol1,as_draws_df(model4_bin)$b_PolC_Alcohol1, as_draws_df(model5_pc1_bin)$b_PolC_Alcohol1))
tylim<-c(allrange1-2, allrange2+1)

res <- as_draws_df(model4)
mstats1 <- mediqr (res$b_PolC_Alcohol1)
hist(res$b_PolC_Alcohol1, xlab="", main="",breaks=40, col=scales::alpha("#eda31d",0.5), border=scales::alpha("#eda31d",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim)

res <- as_draws_df(model4_bin)
mstats2 <- mediqr (res$b_PolC_Alcohol1)
hist(res$b_PolC_Alcohol1, xlab="", main="",breaks=60,col=scales::alpha("#618b44",0.5), border=scales::alpha("#618b44",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim,add=T)
par(xpd=F)
abline(v=0, lty=2, col="gray70")
par(xpd=T)

points(x=c(mstats2[2], mstats2[3]),y=c(-0.06,-0.06), type='l', lwd=4, col=scales::alpha("#618b44",0.75))
points(x=mstats2[1],y=c(-0.06), pch=18,  col="#618b44",cex=2)
points(x=c(mstats1[2], mstats1[3]),y=c(-0.03,-0.03), type='l', lwd=4, col=scales::alpha("#eda31d",0.75))
points(x=mstats1[1],y=c(-0.03), pch=18,  col="#eda31d",cex=2)
legend(x=4,y=1, fill=c("#eda31d","#618b44"), legend=c("ordinal", "binary"), border=NA, bty = "n")
#mtext(side=3,adj=0,text=paste0("Model4", "\n", "PC ~ Alcohol + Agri", "\n", "phylo-space"), line=-2, font=2)
mtext(side=3,adj=0,text=paste0("c) Model 4 - conservative"), line=-1, font=2)
mtext(side=3,adj=0,text=paste0("Political Complexity", " ~", "\n","              ", "Alcohol", "\n","              ","Agriculture","\n","              ", "phylo-space" ), line=-5, font=1)

res <- as_draws_df(model5_pc1)
mstats3 <- mediqr (res$b_PolC_Alcohol1)
hist(res$b_PolC_Alcohol1, xlab="", main="",breaks=60,col=scales::alpha("#eda31d",0.5), border=scales::alpha("#eda31d",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim)
par(xpd=F)
abline(v=0, lty=2, col="gray70")
par(xpd=T)

res <- as_draws_df(model5_pc1_bin)
mstats4 <- mediqr (res$b_PolC_Alcohol1)
hist(res$b_PolC_Alcohol1, xlab="", main="",breaks=70,col=scales::alpha("#618b44",0.5), border=scales::alpha("#618b44",0.25),freq=F, yaxt="n", ylab="", xaxt="n", xlim=tylim,add=T)

points(x=c(mstats4[2], mstats4[3]),y=c(-0.06,-0.06), type='l', lwd=4, col=scales::alpha("#618b44",0.75))
points(x=mstats4[1],y=c(-0.06), pch=18,  col="#618b44",cex=2)
points(x=c(mstats3[2], mstats3[3]),y=c(-0.03,-0.03), type='l', lwd=4, col=scales::alpha("#eda31d",0.75))
points(x=mstats3[1],y=c(-0.03), pch=18,  col="#eda31d",cex=2)
mtext(side=3,adj=0,text=paste0("d) Model 5 - conservative"), line=-0.6, font=2)
mtext(side=3,adj=0,text=paste0("Political Complexity", " ~", "\n","              ", "Alcohol", "\n","              ","Agriculture","\n","              ","Env-productivity","\n","              ","phylo-space" ), line=-6, font=1)

par(xpd=F)
axis(side=1, tck=0.002, lwd.ticks = 2, line = +0.9, lwd=1.5, cex.axis=1.5,xlim=tylim, at=c(-1,0,1,2,3,4,5,6))

x <- dev.off()


