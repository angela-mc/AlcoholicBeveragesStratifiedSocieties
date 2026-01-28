source(file = "scripts/00_set_up.R")
source("scripts/h0_fxn.R")


dir <-  paste0(outfol, "plots_review/")
if(!dir.exists(dir)){
  dir.create(dir)
}


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


###################################
##### conditional effects plots ### additional inspection of results
###################################

# the estimated probabilities of the ordinal response variable falling into or below each category given specific values of the predictor variable.
  # For binary predictors, it calculates these probabilities separately for the cases where the predictor equals 0 and where it equals 1.

pdf(paste0(outfol, "plots/cond_effects_alcPolC_allmodel_panel.pdf"), height=7, width=12)

par(mfrow=c(2,3), mar=c(2.5,6,3.5,1))
load(file=paste0(outfol,"models/model1R.rds"))
ce_data<-conditional_effects(model1, "Alcohol", categorical = TRUE)
ce_data<-ce_data[[1]]
condeff_plot(ce_data = ce_data, title=paste0("Model 1"))
mtext("a)",adj=0,cex=1.1)

load(file=paste0(outfol,"models/model2R.rds"))
ce_data<-conditional_effects(model2, "Alcohol", categorical = TRUE)
ce_data<-ce_data[[1]]
condeff_plot(ce_data = ce_data, title=paste0("Model 2"))
mtext("b)",adj=0,cex=1.1)

load(file=paste0(outfol,"models/model32_pc1.rds"))
ce_data<-conditional_effects(model32_pc1, "Alcohol", categorical = TRUE)
ce_data<-ce_data[[1]]
condeff_plot(ce_data = ce_data, title = paste0("Model 3"))
mtext("c)",adj=0,cex=1.1)

load(file=paste0(outfol,"models/model40R.rds"))
ce_data<-conditional_effects(model40, "Alcohol", categorical = TRUE)
ce_data<-ce_data[[1]]
condeff_plot(ce_data = ce_data, title=paste0("Model 4"))
mtext("d)",adj=0,cex=1.1)

load(file=paste0(outfol,"models/model56_pc1.rds")) 
ce_data<-conditional_effects(model56_pc1, "Alcohol", categorical = TRUE)
ce_data<-ce_data[[1]]
condeff_plot(ce_data = ce_data, title=paste0("Model 5"))
mtext("e)",adj=0,cex=1.1)

tylim<-c(0,0.85)
a<-1:5
b<- seq(from=0, to=0.85, length=5)
plot(b~a,xaxt = 'n', yaxt = 'n', bty='n', xlab="", pch='', ylim=tylim,ylab="", cex.lab=2, main="", cex.main=2)
cols<-c("darkolivegreen4","cornflowerblue","gold","darkorchid","coral1")
legend(x=1,y=0.8, legend=c("none", "local com", "local + 1","local + 2", "local + 3"), fill=cols,border=cols,title = "Political complexity", bty="n",cex=1.3)
legend(x=1,y=0.4, legend=c("alcohol absent", "alcohol present"), pch=c(19,17), bty="n",cex=1.3)

x <- dev.off()
