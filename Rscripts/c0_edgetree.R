source(file = "scripts/00_set_up.R")

######################
######## raw data ####
######################

read.nexus(file=paste0(projfol,"EDGEtree/globaltrees_May/mcc_cah_full.nex"))-> edget 
  # mcc containing all languages (mcc was done in TreeAnnotator, with common ancestor heights)
  # mcc was done by taking the Bouckaert posterior trees, upload in TreeAnnotator, setting common ancestor heights
  # Bouckaert, R. et al. [Pre-print] Global language diversification is linked to socio-ecology and threat status.  (2022).
alcohol<-read.csv(paste0(projfol,"Dataset/ALCOHOL_dataset.csv"),stringsAsFactors = F) 
treetips<-data.frame(alcohol$glottocode, alcohol$Society_name, alcohol$SCCS_ID, edgetips="", note="")

t <- edget$tip.label[edget$tip.label%in%treetips$alcohol.glottocode]
for(i in 1:length(t)) treetips$edgetips[treetips$alcohol.glottocode%in%t[i]]<-t[i]
treetips[!treetips$edgetips%in%"",]$note<-"glotto_InstantMatch" # match in the edge tree based on glottocode


# some manual matching; use as information: 
  # glottolog: https://glottolog.org/glottolog/language
  # dplace: https://d-place.org/societysets#1/29/153

  # info on edge tips from Quentin - misses one taxa added in the latest (May) version
  # dftips<-read.csv(file=paste0(projfol,"EDGEtree/edgetree_labels.csv")) 

treetips[treetips$alcohol.glottocode%in%"nama1265",]$edgetips<-"nama1264" # parent
treetips[treetips$alcohol.glottocode%in%"town1238",]$edgetips<-"bemb1257" # parent
treetips[treetips$alcohol.glottocode%in%"nyak1261",]$edgetips<-"nyak1260" # parent
treetips[treetips$alcohol.glottocode%in%"west1508",]$edgetips<-"basq1248" # grand - parent
treetips[treetips$alcohol.glottocode%in%"bauu1243",]$edgetips<-"fiji1243" # great-grand-parent
treetips[treetips$alcohol.glottocode%in%"east2533",]$edgetips<-"aleu1260" # parent
treetips[treetips$alcohol.glottocode%in%"alba1270",]$edgetips<-"nort2961" # parent 
treetips[treetips$alcohol.glottocode%in%"east2545",]$edgetips<-"kash1280" # or "sout2984" - both closely-related taxa "niece"
treetips[treetips$alcohol.glottocode%in%"vall1251",]$edgetips<-"yoku1256" # parent
treetips[treetips$alcohol.glottocode%in%"nort1551",]$edgetips<-"nort2954" # parent
treetips[treetips$alcohol.glottocode%in%"omah1248",]$edgetips<-"omah1247" # parent
treetips[treetips$alcohol.glottocode%in%"hava1249",]$edgetips<-"hava1248" # parent
treetips[treetips$alcohol.glottocode%in%"toho1246",]$edgetips<-"toho1245" # parent

treetips[treetips$alcohol.glottocode%in%"nkun1238",]$edgetips<-"mong1338" # parent
treetips[treetips$alcohol.glottocode%in%"asan1239",]$edgetips<-"akan1250" # grand-parent
treetips[treetips$alcohol.glottocode%in%"taln1239",]$edgetips<-"fare1241" # parent
treetips[treetips$alcohol.glottocode%in%"vedd1240",]$edgetips<-"sinh1246" # sister

treetips[treetips$alcohol.glottocode%in%"tana1285",]$edgetips<-"plat1254" # parent
treetips[treetips$alcohol.glottocode%in%"copp1244",]$edgetips<-"west2618" # parent
treetips[treetips$alcohol.glottocode%in%"twan1247",]$edgetips<-"stra1244" # closely-related taxa "niece"
treetips[treetips$alcohol.glottocode%in%"gros1243",]$edgetips<-"arap1274" # sister 
treetips[treetips$alcohol.glottocode%in%"clas1250",]$edgetips<-"nort2957" # sister - sister to Aztec
treetips[treetips$alcohol.glottocode%in%"isla1278",]$edgetips<-"gari1256" # sister
treetips[treetips$alcohol.glottocode%in%"mura1272",]$edgetips<-"gali1262" # parent
treetips[treetips$alcohol.glottocode%in%"tupi1273",]$edgetips<-"nhen1239" # sister
treetips[treetips$alcohol.glottocode%in%"nort2971",]$edgetips<-"leng1262" # parent
treetips[treetips$alcohol.glottocode%in%"abip1241",]$edgetips<-"moco1246" # closely related taxa "niece"

treetips[treetips$alcohol.glottocode%in%"bana1287",]$edgetips<-"gilb1244" # parent
treetips[treetips$alcohol.glottocode%in%"boro1274",]$edgetips<-"nige1253" # parent
treetips[treetips$alcohol.glottocode%in%"eyak1241",]$edgetips<-"" 
treetips[treetips$alcohol.glottocode%in%"hogg1238",]$edgetips<-"taha1241" # parent
treetips[treetips$alcohol.glottocode%in%"nort3051",]$edgetips<-"soma1255" # parent
treetips[treetips$alcohol.glottocode%in%"rali1241",]$edgetips<-"mars1254" # parent
treetips[treetips$alcohol.glottocode%in%"bass1257",]$edgetips<-"west2369" # great-grand-parent
treetips[treetips$alcohol.glottocode%in%"oldk1249",]$edgetips<-"nort2684" # sister

treetips[treetips$alcohol.glottocode%in%c("nama1265","town1238","nyak1261","east2533","alba1270",
                                           "vall1251","nort1551","omah1248","hava1249","toho1246",
                                           "taln1239","tana1285","copp1244","mura1272",
                                           "bana1287","boro1274","hogg1238","nort3051","rali1241","nkun1238"),]$note<-"parent"
treetips[treetips$alcohol.glottocode%in%c("vedd1240","gros1243","clas1250","isla1278","oldk1249","tupi1273"),]$note<-"sister"
treetips[treetips$alcohol.glottocode%in%c("west1508","bauu1243","kash1280","bass1257","east2545",
                                           "twan1247","abip1241","akan1250","asan1239"),]$note<-"other-closely-related-taxa"
treetips[treetips$alcohol.glottocode%in%c("nort2971"),]$note<-"parent"

# QUESTION MARK - no match
treetips[treetips$edgetips%in%"",]
  # anci1244 (Ancient Hebrews) match to hebr1245_ModernHebrew - decided against
  # akka1240 (Babylonians) match to # - https://glottolog.org/resource/languoid/id/assy1241 - decided against
  # lati1261 (Romans) match to -- 

  # hait1244 match to -- (creole)
  # sara1340 match to -- (creole)
  
  # akab1249 match to --
  # klam1254 match to -- isolate?
  # natc1249 match to -- isolate?
  # yuro1248 match to --
  # eyak1241 match to --

treetips[which(duplicated(treetips$edgetips)),] 
table(treetips$note)
treetips[treetips$note%in%"other-closely-related-taxa",]

################
# PRUNE TREE  ##
################

tokeep<-edget$tip.label[edget$tip.label%in%treetips$edgetips]
pruned_tree <- ape::keep.tip(edget, tokeep)
matchtips<-treetips[treetips$edgetips%in%pruned_tree$tip.label,]
matchtips <- matchtips[match(pruned_tree$tip.label, matchtips$edgetips),]
pruned_tree$tip.label <- matchtips$alcohol.SCCS_ID
write.nexus(pruned_tree, file=paste0(projfol,"EDGEtree/edgetree_sccs.nex"))
write.csv(treetips, file=paste0(projfol,"EDGEtree/treetips.csv"))


