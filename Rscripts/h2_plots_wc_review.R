source(file = "scripts/00_set_up.R")

source("scripts/h0_fxn.R")


#########################
#### all models alcohol # compare conservative and wide data - for main figure 4
#########################

load(file=paste0(outfol,"models/model1R.rds")) ; model1 -> wmodel1 ; rm(model1)
load(file=paste0(outfol,"models/model2R.rds")) ; model2 -> wmodel2 ; rm(model2)
load(file=paste0(outfol,"models/model3_pc1R.rds")) ; model3_pc1 -> wmodel3_pc1 ; rm(model3_pc1) 
load(file=paste0(outfol,"models/model4R.rds")) ; model4 -> wmodel4 ; rm(model4)
load(file=paste0(outfol,"models/model5_pc1_R.rds")) ; model5_pc1 -> wmodel5_pc1 ; rm(model5_pc1)

load(file=paste0(outfol,"models_cons/cmodel1R.rds")) ; model1 -> cmodel1 ; rm(model1)
load(file=paste0(outfol,"models_cons/cmodel2R.rds")) ; model2 -> cmodel2 ; rm(model2)
load(file=paste0(outfol,"models_cons/cmodel3_pc1R.rds")) ; model3_pc1 -> cmodel3_pc1 ; rm(model3_pc1)
load(file=paste0(outfol,"models_cons/cmodel4R.rds")) ; model4 -> cmodel4 ; rm(model4)
load(file=paste0(outfol,"models_cons/cmodel5_pc1_R.rds")) ; model5_pc1 -> cmodel5_pc1 ; rm(model5_pc1)

#pdf(paste0(outfol,"plots_review/figure3_rev.pdf"), height=10, width=9)
par(mfrow=c(5,1), mar=c(4.5,3,1.5,1))

allrange1<- min(c(as_draws_df(cmodel1)$b_Alcohol1, as_draws_df(cmodel2)$b_Alcohol1,as_draws_df(cmodel3_pc1)$b_PolC_Alcohol1, as_draws_df(cmodel4)$b_PolC_Alcohol1, as_draws_df(cmodel5_pc1)$b_PolC_Alcohol1,
                  as_draws_df(wmodel1)$b_Alcohol1, as_draws_df(wmodel2)$b_Alcohol1,as_draws_df(wmodel3_pc1)$b_PolC_Alcohol1, as_draws_df(wmodel4)$b_PolC_Alcohol1, as_draws_df(wmodel5_pc1)$b_PolC_Alcohol1))
allrange2<- max(c(as_draws_df(cmodel1)$b_Alcohol1, as_draws_df(cmodel2)$b_Alcohol1,as_draws_df(cmodel3_pc1)$b_PolC_Alcohol1, as_draws_df(cmodel4)$b_PolC_Alcohol1, as_draws_df(cmodel5_pc1)$b_PolC_Alcohol1,
                  as_draws_df(wmodel1)$b_Alcohol1, as_draws_df(wmodel2)$b_Alcohol1,as_draws_df(wmodel3_pc1)$b_PolC_Alcohol1, as_draws_df(wmodel4)$b_PolC_Alcohol1, as_draws_df(wmodel5_pc1)$b_PolC_Alcohol1))
tylim<-c(allrange1-4, allrange2+2.5)

#m1
res <- as_draws_df(wmodel1)
mstats2 <- mediqr2 (res$b_Alcohol1)
hist(res$b_Alcohol1, xlab="", main="",breaks=60,col=scales::alpha("#eda31d",0.5), border=scales::alpha("#eda31d",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim)
res <- as_draws_df(cmodel1)
mstats1 <- mediqr2 (res$b_Alcohol1)
hist(res$b_Alcohol1, xlab="", main="",breaks=40, col=scales::alpha("#009dc8",0.5), border=scales::alpha("#009dc8",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim,add=T)
par(xpd=F)
abline(v=0, lty=2, col="gray70")
par(xpd=T)

points(x=c(mstats2$lower66, mstats2$upper66),y=c(-0.07,-0.07),  type='l', lwd=4, col=scales::alpha("#eda31d",0.75))
points(x=c(mstats2$lower95, mstats2$upper95),y=c(-0.07,-0.07),  type='l', lwd=2, col=scales::alpha("#eda31d",0.75))
points(x=mstats2$mean,y=c(-0.07), pch=18,  col="#eda31d",cex=2)
points(x=c(mstats1$lower66, mstats1$upper66),y=c(-0.2,-0.2),type='l', lwd=4, col=scales::alpha("#009dc8",0.75))
points(x=c(mstats1$lower95, mstats1$upper95),y=c(-0.2,-0.2),type='l', lwd=2, col=scales::alpha("#009dc8",0.75))
points(x=mstats1$mean,y=c(-0.2), pch=18,  col="#009dc8",cex=2)
legend(x=3,y=1.7, fill=c("#eda31d","#009dc8"), legend=c("wide", "conservative"), border=NA, bty = "n", cex=1.5)
# mtext(side=3,adj=0,text=paste0("Model 1", "\n", "Political Complexity ~ Alcohol"), line=-2, font=2)
mtext(side=3,adj=0,text=paste0("Model 1"), line=-2, font=2)
mtext(side=3,adj=0,text=paste0("Political Complexity", " ~", "\n","              ", "Alcohol"), line=-5, font=1)


#m2
res <- as_draws_df(wmodel2)
mstats2 <- mediqr2 (res$b_Alcohol1)
hist(res$b_Alcohol1, xlab="", main="",breaks=70,col=scales::alpha("#eda31d",0.5), border=scales::alpha("#eda31d",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim)
res <- as_draws_df(cmodel2)
mstats1 <- mediqr2 (res$b_Alcohol1)
hist(res$b_Alcohol1, xlab="", main="",breaks=90, col=scales::alpha("#009dc8",0.5), border=scales::alpha("#009dc8",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim,add=T)
par(xpd=F)
abline(v=0, lty=2, col="gray70")
par(xpd=T)

points(x=c(mstats2$lower66, mstats2$upper66),y=c(-0.05,-0.05), type='l', lwd=4, col=scales::alpha("#eda31d",0.75))
points(x=c(mstats2$lower95, mstats2$upper95),y=c(-0.05,-0.05), type='l', lwd=2, col=scales::alpha("#eda31d",0.75))
points(x=mstats2$mean,y=c(-0.05), pch=18,  col="#eda31d",cex=2)
points(x=c(mstats1$lower66, mstats1$upper66),y=c(-0.143,-0.143), type='l', lwd=4, col=scales::alpha("#009dc8",0.75))
points(x=c(mstats1$lower95, mstats1$upper95),y=c(-0.143,-0.143), type='l', lwd=2, col=scales::alpha("#009dc8",0.75))
points(x=mstats1$mean,y=c(-0.143), pch=18,  col="#009dc8",cex=2)
#mtext(side=3,adj=0,text=paste0("Model 2", "\n", "Political complexity ~ Alcohol", "\n", "phylo-space"), line=-1, font=2)
mtext(side=3,adj=0,text=paste0("Model 2"), line=-1.5, font=2)
mtext(side=3,adj=0,text=paste0("Political Complexity", " ~", "\n","              ", "Alcohol", "\n","              ", "phylo-space"), line=-6, font=1)

# m3
res <- as_draws_df(wmodel3_pc1)
mstats2 <- mediqr2 (res$b_PolC_Alcohol1)
hist(res$b_PolC_Alcohol1, xlab="", main="",breaks=100,col=scales::alpha("#eda31d",0.5), border=scales::alpha("#eda31d",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim)
res <- as_draws_df(cmodel3_pc1)
mstats1 <- mediqr2 (res$b_PolC_Alcohol1)
hist(res$b_PolC_Alcohol1, xlab="", main="",breaks=140, col=scales::alpha("#009dc8",0.5), border=scales::alpha("#009dc8",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim,add=T)
par(xpd=F)
abline(v=0, lty=2, col="gray70")
par(xpd=T)

points(x=c(mstats2$lower66, mstats2$upper66),y=c(-0.045,-0.045), type='l', lwd=4, col=scales::alpha("#eda31d",0.75))
points(x=c(mstats2$lower95, mstats2$upper95),y=c(-0.045,-0.045), type='l', lwd=2, col=scales::alpha("#eda31d",0.75))
points(x=mstats2$mean,y=c(-0.045), pch=18,  col="#eda31d",cex=2)
points(x=c(mstats1$lower66, mstats1$upper66),y=c(-0.125,-0.125), type='l', lwd=4, col=scales::alpha("#009dc8",0.75))
points(x=c(mstats1$lower95, mstats1$upper95),y=c(-0.125,-0.125), type='l', lwd=2, col=scales::alpha("#009dc8",0.75))
points(x=mstats1$mean,y=c(-0.125), pch=18,  col="#009dc8",cex=2)
mtext(side=3,adj=0,text=paste0("Model 3"), line=-1, font=2)
mtext(side=3,adj=0,text=paste0("Political Complexity", " ~", "\n","              ", "Alcohol", "\n","              ","Env-productivity","\n","              ", "phylo-space" ), line=-7, font=1)

#m4
res <- as_draws_df(wmodel4)
mstats2 <- mediqr2 (res$b_PolC_Alcohol1)
hist(res$b_PolC_Alcohol1, xlab="", main="",breaks=60,col=scales::alpha("#eda31d",0.5), border=scales::alpha("#eda31d",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim)
res <- as_draws_df(cmodel4)
mstats1 <- mediqr2 (res$b_PolC_Alcohol1)
hist(res$b_PolC_Alcohol1, xlab="", main="",breaks=120, col=scales::alpha("#009dc8",0.5), border=scales::alpha("#009dc8",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim,add=T)
par(xpd=F)
abline(v=0, lty=2, col="gray70")
par(xpd=T)

points(x=c(mstats2$lower66, mstats2$upper66),y=c(-0.045,-0.045), type='l', lwd=4, col=scales::alpha("#eda31d",0.75))
points(x=c(mstats2$lower95, mstats2$upper95),y=c(-0.045,-0.045), type='l', lwd=2, col=scales::alpha("#eda31d",0.75))
points(x=mstats2$mean,y=c(-0.045), pch=18,  col="#eda31d",cex=2)
points(x=c(mstats1$lower66, mstats1$upper66),y=c(-0.13,-0.13), type='l', lwd=4, col=scales::alpha("#009dc8",0.75))
points(x=c(mstats1$lower95, mstats1$upper95),y=c(-0.13,-0.13), type='l', lwd=2, col=scales::alpha("#009dc8",0.75))
points(x=mstats1$mean,y=c(-0.13), pch=18,  col="#009dc8",cex=2)
mtext(side=3,adj=0,text=paste0("Model 4"), line=-1, font=2)
mtext(side=3,adj=0,text=paste0("Political Complexity", " ~", "\n","              ", "Alcohol", "\n","              ","Agriculture","\n","              ", "phylo-space" ), line=-7, font=1)

#m5
res <- as_draws_df(wmodel5_pc1)
mstats2 <- mediqr2 (res$b_PolC_Alcohol1)
hist(res$b_PolC_Alcohol1, xlab="", main="",breaks=100, col=scales::alpha("#eda31d",0.5), border=scales::alpha("#eda31d",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim,ylim=c(0,1))
res <- as_draws_df(cmodel5_pc1)
mstats1 <- mediqr2 (res$b_PolC_Alcohol1)
hist(res$b_PolC_Alcohol1, xlab="", main="",breaks=100,col=scales::alpha("#009dc8",0.5), border=scales::alpha("#009dc8",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim,add=T)
par(xpd=F)
abline(v=0, lty=2, col="gray70")
par(xpd=T)

points(x=c(mstats2$lower66, mstats2$upper66),y=c(-0.06,-0.06), type='l', lwd=4, col=scales::alpha("#eda31d",0.75))
points(x=c(mstats2$lower95, mstats2$upper95),y=c(-0.06,-0.06), type='l', lwd=2, col=scales::alpha("#eda31d",0.75))
points(x=mstats2$mean,y=c(-0.06), pch=18,  col="#eda31d",cex=2)
points(x=c(mstats1$lower66, mstats1$upper66),y=c(-0.15,-0.15), type='l', lwd=4, col=scales::alpha("#009dc8",0.75))
points(x=c(mstats1$lower95, mstats1$upper95),y=c(-0.15,-0.15), type='l', lwd=2, col=scales::alpha("#009dc8",0.75))
points(x=mstats1$mean,y=c(-0.15), pch=18,  col="#009dc8",cex=2)
mtext(side=3,adj=0,text=paste0("Model 5"), line=-0.6, font=2)
mtext(side=3,adj=0,text=paste0("Political Complexity", " ~", "\n","              ", "Alcohol", "\n","              ","Agriculture","\n","              ","Env-productivity","\n","              ","phylo-space" ), line=-8, font=1)

axis(side=1, tck=0.002, lwd.ticks = 2, line = +1.4, lwd=1.5, cex.axis=1.5,xlim=tylim, at=c(-2, 0,2,4,6,8,10))
dev.off()

res <- as_draws_df(wmodel2)
median(res$b_Alcohol1)
sum(res$b_Alcohol1>0)/length(res$b_Alcohol1)

res <- as_draws_df(wmodel4)
median(res$b_PolC_Alcohol1)
sum(res$b_PolC_Alcohol1>0)/length(res$b_PolC_Alcohol1) # 84

res <- as_draws_df(wmodel3_pc1)
median(res$b_PolC_Alcohol1)
sum(res$b_PolC_Alcohol1>0)/length(res$b_PolC_Alcohol1) # 99

res <- as_draws_df(wmodel5_pc1)
median(res$b_PolC_Alcohol1)
sum(res$b_PolC_Alcohol1>0)/length(res$b_PolC_Alcohol1) 

res <- as_draws_df(cmodel5_pc1)
median(res$b_PolC_Alcohol1)
sum(res$b_PolC_Alcohol1>0)/length(res$b_PolC_Alcohol1) # 96%



##################
#### agriculture # # on political complexity with and without environment - compare conservative and wide data (no alcohol) - for EXT DATA fig 2
##################
# figure s6

load(file=paste0(outfol,"models/model42R.rds")) ; model42 -> wmodel42 ; rm(model42)
load(file=paste0(outfol,"models/model51R_pc1.rds")) ; model51_pc1 -> wmodel51 ; rm(model51_pc1)
load(file=paste0(outfol,"models_cons/cmodel42R.rds")) ; model42 -> cmodel42 ; rm(model42)
load(file=paste0(outfol,"models_cons/cmodel51R_pc1.rds")) ; model51_pc1 -> cmodel51 ; rm(model51_pc1)

#pdf(paste0(outfol,"plots_review/figureS6_review.pdf"), height=8, width=7)
par(mfrow=c(2,1), mar=c(3,3,2,2))

res42_c<- as_draws_df(cmodel42)
res42_c$bE <- res42_c$bsp_moAgriculture*4
res42_w<- as_draws_df(wmodel42)
res42_w$bE <- res42_w$bsp_moAgriculture*4

res51_w<- as_draws_df(wmodel51)
res51_w$bE <- res51_w$bsp_PolC_moAgriculture*4
res51_c<- as_draws_df(cmodel51)
res51_c$bE <- res51_c$bsp_PolC_moAgriculture*4

allrange1 <- min(c(res42_c$bE, res42_w$bE, res51_w$bE, res51_c$bE))
allrange2 <- max(c(res42_c$bE, res42_w$bE, res51_w$bE, res51_c$bE))
tylim<-c(allrange1-5.5, allrange2+0.01)

mstats2 <- mediqr2 (res42_w$bE)
hist(res42_w$bE, xlab="", main="",breaks=80,col=scales::alpha("#eda31d",0.5), border=scales::alpha("#eda31d",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim)
mstats1 <- mediqr2 (res42_c$bE)
hist(res42_c$bE, xlab="", main="",breaks=80, col=scales::alpha("#009dc8",0.5), border=scales::alpha("#009dc8",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim,add=T)
par(xpd=F)
abline(v=0, lty=2, col="gray70")
par(xpd=T)

points(x=c(mstats2$lower66, mstats2$upper66),y=c(-0.03,-0.03), type='l', lwd=4, col=scales::alpha("#eda31d",0.75))
points(x=c(mstats2$lower95, mstats2$upper95),y=c(-0.03,-0.03), type='l', lwd=2, col=scales::alpha("#eda31d",0.75))
points(x=mstats2$mean,y=c(-0.03), pch=18,  col="#eda31d",cex=2)
points(x=c(mstats1$lower66, mstats1$upper66),y=c(-0.06,-0.06), type='l', lwd=4, col=scales::alpha("#009dc8",0.75))
points(x=c(mstats1$lower95, mstats1$upper95),y=c(-0.06,-0.06), type='l', lwd=2, col=scales::alpha("#009dc8",0.75))
points(x=mstats1$mean,y=c(-0.06), pch=18,  col="#009dc8",cex=2)
mtext(side=3,adj=0,text=paste0("a)"), line=-1, font=1)
mtext(side=3,adj=0,text=paste0("Political Complexity", " ~", "\n","              ", " ","Agriculture","\n","              ", " phylo-space" ), line=-4, font=1)
legend(x=6,y=0.7, fill=c("#eda31d","#009dc8"), legend=c("wide", "conservative"), border=NA, bty = "n")

mstats2 <- mediqr2 (res51_w$bE)
hist(res51_w$bE, xlab="", main="",breaks=80,col=scales::alpha("#eda31d",0.5), border=scales::alpha("#eda31d",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim)
mstats1 <- mediqr2 (res51_c$bE)
hist(res51_c$bE, xlab="", main="",breaks=100, col=scales::alpha("#009dc8",0.5), border=scales::alpha("#009dc8",0.25),freq = F, yaxt="n", ylab="", xaxt="n", xlim=tylim,add=T)
par(xpd=F)
abline(v=0, lty=2, col="gray70")
par(xpd=T)

points(x=c(mstats2$lower66, mstats2$upper66),y=c(-0.02,-0.02), type='l', lwd=4, col=scales::alpha("#eda31d",0.75))
points(x=c(mstats2$lower95, mstats2$upper95),y=c(-0.02,-0.02), type='l', lwd=2, col=scales::alpha("#eda31d",0.75))
points(x=mstats2$mean,y=c(-0.02), pch=18,  col="#eda31d",cex=2)
points(x=c(mstats1$lower66, mstats1$upper66),y=c(-0.045,-0.045), type='l', lwd=4, col=scales::alpha("#009dc8",0.75))
points(x=c(mstats1$lower95, mstats1$upper95),y=c(-0.045,-0.045), type='l', lwd=2, col=scales::alpha("#009dc8",0.75))
points(x=mstats1$mean,y=c(-0.045), pch=18,  col="#009dc8",cex=2)
mtext(side=3,adj=0,text=paste0("b)"), line=-1, font=1)
mtext(side=3,adj=0,text=paste0("Political Complexity", " ~", "\n","              ", " ","Agriculture","\n","              ", " phylo-space", "\n","              "," Env-productivity" ), line=-5, font=1)

par(xpd=F)
axis(side=1, tck=0.002, lwd.ticks = 2, line = +1, lwd=1.5, cex.axis=1.5,xlim=tylim, at=c(0,2,4,6,8,10,12))

dev.off()



