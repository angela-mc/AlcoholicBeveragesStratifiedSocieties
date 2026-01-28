source(file = "scripts/00_set_up.R")

######################
######## raw data #### 
######################

library(readxl)
alcohol<-read_excel(paste0(projfol,"Dataset/ALCOHOL_dataset.xlsx")) # Alcohol_WB = wine & beers = variable of interest
# alcohol<-read_excel(paste0(projfol,"Dataset/ALCOHOL_dataset.xlsx"),na = "NA") # Alcohol_WB = wine & beers = variable of interest
table(alcohol$Alcohol_all)
table(alcohol$Alcohol_WB)

scaffold<-read.csv(paste0(projfol,"output/workingdata/societies.csv"),stringsAsFactors = F) # dplace societies file
scaffold<-scaffold[match(alcohol$SCCS_ID, scaffold$id),]
alcohol$glottocode<-scaffold$glottocode

  # Sample 1: Present = 7,8; Absent = 0,1,2,3,4,5
  # Sample 2: Present = 7,8; Absent = 0,2,3

# sample 1 (wide)
sum(alcohol$Alcohol_WB%in%c(7,8))
sum(alcohol$Alcohol_WB%in%c(0,1,2,3,4,5))
sum(alcohol$Alcohol_WB%in%c(6,NA))
sum(alcohol$Alcohol_WB%in%c(6))
sum(alcohol$Alcohol_WB%in%c(NA))

# sample 2 (conservative)
sum(alcohol$Alcohol_WB%in%c(7,8))
sum(alcohol$Alcohol_WB%in%c(0,2,3))
sum(alcohol$Alcohol_WB%in%c(1,4,5,6,NA))
sum(alcohol$Alcohol_WB%in%NA)

alcohol$samplewide<--1
alcohol[alcohol$Alcohol_WB%in%c(7,8),]$samplewide<-1
alcohol[alcohol$Alcohol_WB%in%c(0,1,2,3,4,5),]$samplewide<-0
alcohol[alcohol$Alcohol_WB%in%c(6),]$samplewide<-NA #the Value 6 represents alcohol being present, but origins unsure. This is treated same as missing data for the purposes of this study.
alcohol[alcohol$Alcohol_WB%in%NA | alcohol$Alcohol_WB%in%"NA",]$samplewide<-NA
table(alcohol$samplewide)

alcohol$samplecons<--1
alcohol[alcohol$Alcohol_WB%in%c(8),]$samplecons<-1
alcohol[alcohol$Alcohol_WB%in%c(0,2,3),]$samplecons<-0
alcohol[alcohol$Alcohol_WB%in%c(1,4,5,6,7),]$samplecons<-NA
alcohol[alcohol$Alcohol_WB%in%NA | alcohol$Alcohol_WB%in%"NA",]$samplecons<-NA
table(alcohol$samplecons)


# write.csv(alcohol, file=paste0(projfol,"Dataset/ALCOHOL_dataset.csv")) # FROM NOW ON - THIS IS THE working DATASET
alcohol<-read.csv(paste0(projfol,"Dataset/ALCOHOL_dataset.csv")) # Alcohol_WB = wine & beers = variable of interest