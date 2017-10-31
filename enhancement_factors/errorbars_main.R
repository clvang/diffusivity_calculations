
rm(list=ls(all=TRUE))   #remove all variables in  workspace
graphics.off()  #close all graphics windows


source('/Users/changlvang/mygitFiles/diffusivity_calculations/enhancement_factors/mc_analysis.R')
source('/Users/changlvang/mygitFiles/diffusivity_calculations/enhancement_factors/enhancement_factors.R')
source('/Users/changlvang/mygitFiles/diffusivity_calculations/enhancement_factors/tdelay_tdecay.R')

P <- 0.95
N <- 1000000
# call function to calculate error bars for effective diffusivities
D_eff <- mc_analysis( N=N, P=P )

Deffective <- D_eff$D_N95
do_mcN95 <- D_eff$do_mcN95
K_mcN95 <- D_eff$K_mcN95
Yo_mcN95 <- D_eff$Yo_mcN95

# call function to calculate error bars for enhancement factors
enhancement_factors(N=N,P=P)

# call function to calculate error bars for tdelay/tdecay ratio
tdelay_tdecay(do_mcN95=do_mcN95, K_mcN95=K_mcN95, Yo_mcN95=Yo_mcN95)










