source(file = "scripts/00_set_up.R")

dir <-  paste0(outfol, "cultdata/")
if(!dir.exists(dir)){
  dir.create(dir)
}

##############################
######## cultural variables ##
##############################

# SCCS157 Political Integration - Measurement of cultural complexity
# SCCS151 Agriculture - Agriculture Intensity

alcohol<-read.csv(paste0(projfol,"Dataset/ALCOHOL_dataset.csv"),stringsAsFactors = F) 
sccs_folder<-paste0("../dplace-data-master/datasets/SCCS/")

read.csv(paste0(sccs_folder,"codes.csv"), stringsAsFactors = F)-> codes
read.csv(paste0(sccs_folder,"data.csv"), stringsAsFactors = F)-> data

d1<-data[data$var_id%in%"SCCS157",]
d1<-d1[match(alcohol$SCCS_ID, d1$soc_id),]
d1<-d1[,c("soc_id","year","var_id", "code")]
capture.output(codes[codes$var_id%in%"SCCS157",], file=paste0(outfol,"cultdata/sccs157.txt"))

d2<-data[data$var_id%in%"SCCS151",]
d2<-d2[match(alcohol$SCCS_ID, d2$soc_id),]
d2<-d2[,c("soc_id","year","var_id", "code")]
capture.output(codes[codes$var_id%in%"SCCS151",], file=paste0(outfol,"cultdata/sccs151.txt"))

cdata<-data.frame(alcohol$SCCS_ID)
cdata$PolC<-d1$code
cdata$Agriculture<-d2$code
write.csv(cdata, file=paste0(outfol,"cultdata/cdata.csv"))
 