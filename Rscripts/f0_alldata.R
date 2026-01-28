source(file = "scripts/00_set_up.R")

################
######## data ##
################

read.csv(file=paste0(outfol,"cultdata/cdata.csv"),stringsAsFactors = F)->cdata
read.csv(file=paste0(outfol,"envpca/envdata.csv"), stringsAsFactors = F)->envdata
read.nexus(file=paste0(projfol,"EDGEtree/edgetree_sccs.nex"))->tree

alcohol<-read.csv(paste0(projfol,"Dataset/ALCOHOL_dataset.csv"),stringsAsFactors = F) 
alcohol<-alcohol[complete.cases(alcohol$samplewide) & alcohol$SCCS_ID%in%tree$tip.label,]
sw<-alcohol$SCCS_ID
sw<-sw[sw%in%envdata$alcohol.SCCS_ID & sw%in%cdata$alcohol.SCCS_ID]
 save(sw, file=paste0(outfol,"workingdata/samplewide_tree_alldata.rds"))

alcohol<-read.csv(paste0(projfol,"Dataset/ALCOHOL_dataset.csv"),stringsAsFactors = F) 
alcohol<-alcohol[complete.cases(alcohol$samplecons) & alcohol$SCCS_ID%in%tree$tip.label,]
sc<-alcohol$SCCS_ID
sc<-sc[sc%in%envdata$alcohol.SCCS_ID & sc%in%cdata$alcohol.SCCS_ID]
save(sc, file=paste0(outfol,"workingdata/samplecons_tree_alldata.rds"))


################
######## data ## sample wide
################

read.csv(file=paste0(outfol,"cultdata/cdata.csv"),stringsAsFactors = F)->cdata
read.csv(file=paste0(outfol,"envpca/envdata.csv"), stringsAsFactors = F)->envdata
alcohol<-read.csv(paste0(projfol,"Dataset/ALCOHOL_dataset.csv"),stringsAsFactors = F) 
load(file=paste0(outfol,"workingdata/samplewide_tree_alldata.rds"))

thedata<-alcohol[,c("SCCS_ID","Society_name","glottocode","samplewide")]
thedata<-thedata[thedata$SCCS_ID%in%sw,]

cdata<-cdata[cdata$alcohol.SCCS_ID%in%thedata$SCCS_ID,]
cdata<-cdata[match(thedata$SCCS_ID, cdata$alcohol.SCCS_ID),]
thedata$PolC<-cdata$PolC
thedata$Agriculture<-cdata$Agriculture

envdata<-envdata[envdata$alcohol.SCCS_ID%in%thedata$SCCS_ID,]
envdata<-envdata[match(thedata$SCCS_ID, envdata$alcohol.SCCS_ID),]
thedata$PC1<-envdata$PC1
thedata$PC2<-envdata$PC2
thedata$PC3<-envdata$PC3

# add longitude and latitude
sccs_folder<-paste0(projfol,"dplace-data-master/datasets/SCCS/")
read.csv(paste0(sccs_folder,"societies.csv"), stringsAsFactors = F)->societies
societies<-societies[societies$id%in%thedata$SCCS_ID,]
societies<-societies[match(thedata$SCCS_ID, societies$id),]
thedata$lat<-societies$Lat
thedata$lon<-societies$Long

 write.csv(thedata, file=paste0(outfol,"workingdata/thedata_wide.csv"))


################
######## data ## sample cons
################
 
read.csv(file=paste0(outfol,"cultdata/cdata.csv"),stringsAsFactors = F)->cdata
read.csv(file=paste0(outfol,"envpca/envdata.csv"), stringsAsFactors = F)->envdata
alcohol<-read.csv(paste0(projfol,"Dataset/ALCOHOL_dataset.csv"),stringsAsFactors = F) 
load(file=paste0(outfol,"workingdata/samplecons_tree_alldata.rds"))

thedata<-alcohol[,c("SCCS_ID","Society_name","glottocode","samplecons")]
thedata<-thedata[thedata$SCCS_ID%in%sc,]

cdata<-cdata[cdata$alcohol.SCCS_ID%in%thedata$SCCS_ID,]
cdata<-cdata[match(thedata$SCCS_ID, cdata$alcohol.SCCS_ID),]
thedata$PolC<-cdata$PolC
thedata$Agriculture<-cdata$Agriculture

envdata<-envdata[envdata$alcohol.SCCS_ID%in%thedata$SCCS_ID,]
envdata<-envdata[match(thedata$SCCS_ID, envdata$alcohol.SCCS_ID),]
thedata$PC1<-envdata$PC1
thedata$PC2<-envdata$PC2
thedata$PC3<-envdata$PC3

# add longitude and latitude
sccs_folder<-paste0(projfol,"dplace-data-master/datasets/SCCS/")
read.csv(paste0(sccs_folder,"societies.csv"), stringsAsFactors = F)->societies
societies<-societies[societies$id%in%thedata$SCCS_ID,]
societies<-societies[match(thedata$SCCS_ID, societies$id),]
thedata$lat<-societies$Lat
thedata$lon<-societies$Long

write.csv(thedata, file=paste0(outfol,"workingdata/thedata_cons.csv"))


