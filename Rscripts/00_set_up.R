
rm(list=ls())
projfol<-"../" # overall folder
outfol<-"../output/"

library(ape)
library(phylopath)
library(brms)
library(tidyverse)
library(tidybayes)

verbose <- FALSE


dir <-  paste0(outfol, "workingdata")
if(!dir.exists(dir)){
  dir.create(dir)
}

dir <-  paste0(outfol, "models/")
if(!dir.exists(dir)){
  dir.create(dir)
}

dir <-  paste0(outfol, "models/binaryAgriculture")
if(!dir.exists(dir)){
  dir.create(dir)
}


dir <-  paste0(outfol, "models_cons/")
if(!dir.exists(dir)){
  dir.create(dir)
}

dir <-  paste0(outfol, "models_cons/binaryAgriculture/")
if(!dir.exists(dir)){
  dir.create(dir)
}

dir <-  paste0(outfol, "plots/")
if(!dir.exists(dir)){
  dir.create(dir)
}

dir <-  paste0(outfol, "treesensitivity/model2/")
if(!dir.exists(dir)){
  dir.create(dir)
}

dir <-  paste0(outfol, "treesensitivity/model3/")
if(!dir.exists(dir)){
  dir.create(dir)
}

dir <-  paste0(outfol, "treesensitivity/model4/")
if(!dir.exists(dir)){
  dir.create(dir)
}

dir <-  paste0(outfol, "/treesensitivity/model5/")
if(!dir.exists(dir)){
  dir.create(dir)
}



