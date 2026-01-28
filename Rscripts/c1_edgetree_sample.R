source(file = "scripts/00_set_up.R")

dir <-  paste0(outfol, "treesensitivity/")
if(!dir.exists(dir)){
  dir.create(dir)
}

read.nexus(file=paste0(projfol,"EDGEtree/globaltrees_May/edge6636-March-2023-no-metadata.trees"))-> edget 
  # these trees can be found from Bouckaert et al:Bouckaert, R. et al. [Pre-print] Global language diversification is linked to socio-ecology and threat status.  (2022).

sampletrees<-sample(x=1:1000, size=50)
strees <- edget[sampletrees]
save(strees, file=paste0(outfol,"treesensitivity/edge_sample50.rds"))
save(sampletrees,file=paste0(outfol,"treesensitivity/sample50.rds"))
