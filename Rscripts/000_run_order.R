#start_time_all_scripts <- Sys.time()

source("scripts/00_set_up.R")
source("scripts/a0_rawdata.R")
source("scripts/c0_edgetree.R")
source("scripts/c1_edgetree_sample.R")
source("scripts/d0_envPCA.R")
source("scripts/e0_cultData.R")
source("scripts/f0_alldata.R")
source("scripts/g_fxns.R")
source("scripts/g10_model1.R")
source("scripts/g20_model2.R")
source("scripts/g30_model3.R")
source("scripts/g40_model4.R")
source("scripts/g41_model4_bin.R")
source("scripts/g50_model5.R")
source("scripts/g51_model5_bin.R")
source("scripts/i10_model1.R")
source("scripts/i20_model2.R")
source("scripts/i30_model3.R")
source("scripts/i40_model4.R")
source("scripts/i41_model4_bin.R")
source("scripts/i50_model5.R")
source("scripts/i51_model5_bin.R")

source("scripts/h0_fxn.R")
source("scripts/h0_plots.R")
source("scripts/h1_plots_bin.R")

source("scripts/n_otherintox.R")
source("scripts/o_paragraphs.R")

# end_time_all_scripts <- Sys.time()
# 
# time <- end_time_all_scripts - start_time_all_scripts
# 
# cat(
#   paste0("I'm done with all scripts.\n",
#          "It took ", round(time[[1]], 2), " ",  units(time)),
#   ".\n"
# )


#tree sensitivty (50 sample of posterior)
source("scripts/k20_model2_treesens.R")
source("scripts/k30_model3_treesens.R")
source("scripts/k40_model4_treesens.R")
source("scripts/k50_model5_treesens.R")
source("scripts/l_treesensitvity_summary.R")