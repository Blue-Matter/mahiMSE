## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(dpi=85, fig.width=7.5, fig.height = 5)
#library(RPC)
#library(openMSE)
#library(Inverts)
#library(Inverts.GD)
#library(Inverts.MC)

## ----install_openMSE_devtools,eval=F------------------------------------------
# install.packages('openMSE')

## ----install_github_packages,eval=F-------------------------------------------
# remotes::install_github("blue-matter/mahiMSE")

## ----load_library, eval=F-----------------------------------------------------
# library(mahiMSE)
# vignette('mahiMSE')

## ----QS_load, eval=F----------------------------------------------------------
# library(mahiMSE)

## ----QS_hist, eval=F----------------------------------------------------------
# Hist = SimulateMOM(miniOM_1)

## ----QS_stat_quo_mse, eval=F, out.width="50%", out.height="50%", fig.height = 5----
# StatQuo = HalfEff = mahiMP
# formals(HalfEff)$Effort = 1/2
# 
# myMSE = ProjectMOM(Hist, c('StatQuo','HalfEff'))               # constant current catch levels

## ----QS_stat_quo_mse_plot, eval=T, out.width="50%", out.height="50%", fig.height = 5----
#mahiplot(myMSE)

## ----help_MPMC, eval=F--------------------------------------------------------
# ?mahiMP

## ----args_MP_MC_datarich, eval=F----------------------------------------------
# 

## ----Eff_DR_Fdisc, eval=F-----------------------------------------------------
# 

## ----echo=FALSE, eval=FALSE---------------------------------------------------
# sfStop()

