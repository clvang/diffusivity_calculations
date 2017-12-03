
rm(list=ls(all=TRUE))   #remove all variables in  workspace
library(readxl)
library(moments)

source('/Users/changlvang/mygitFiles/diffusivity_calculations/errorbar_calculations/mc_analysis.R')
source('/Users/changlvang/mygitFiles/diffusivity_calculations/errorbar_calculations/enhancement_factors.R')
source('/Users/changlvang/mygitFiles/diffusivity_calculations/errorbar_calculations/tdelay_tdecay.R')
source('/Users/changlvang/mygitFiles/diffusivity_calculations/errorbar_calculations/generate_fcprops_from_xlsx.R')
source('/Users/changlvang/mygitFiles/diffusivity_calculations/errorbar_calculations/functionPlot.R')

P <- 0.95
N <- 1000000

importedData <- read_excel("experiment_data_key.xlsx",sheet="HepHexPropGlyCSV",
					col_names=TRUE)
data_table <- data.frame( matrix(ncol = 23, nrow = nrow(importedData)) )
data_table <- setNames(data_table, c("FLEX_ID","DmcLow_N95","Dmc_AV_N95",
									"DmcHigh_N95","Dmc_MP_N95",
									"DUS_lower", "DUS_bar", "DUS_upper",
									"tdtv_low95", "tdtv_avg", "tdtv_high95", 
									"EF_low95","EF_avg","EF_high95",
									"tdecay_low95","tdecay_avg","tdecay_high95,",
									"tdelay_low95","tdelay_avg","tdelay_high95",
									"eps_lower", "eps_avg", "eps_upper") )
interestExp <- c(2,6,14,53)  #index of experiments with very large/small U_{epsilon}

for (i in 1:nrow(importedData)){
# for (i in 47:47){	
	graphics.off()  #close all graphics windows

	expname <- importedData$FLEX_ID[i]

	print(paste0("============ START data for ", expname," ============") )
	print( paste0("i= ",i) )
	#generate required fcprops.txt file from .xlsx data
	generate_fcprops_from_xlsx(i=i, importedData=importedData)

	# call function to calculate error bars for effective diffusivities
	D_eff <- mc_analysis( N=N, P=P , expname=expname)

	# the folllowing values are needed as inputs
	# to the function tdelay_tdecay
	Deffective <- D_eff$D_N95
	do_mcN95 <- D_eff$do_mcN95
	K_mcN95 <- D_eff$K_mcN95
	Yo_mcN95 <- D_eff$Yo_mcN95
	eps_N95 <- D_eff$eps_N95
	LHS_mcN95 <- D_eff$LHS_mcN95
	tau_mcN95 <- D_eff$tau_mcN95

	# # functional evaluatiosn to check solution "space"
	# if(i %in% interestExp){
	# 	functionPlot(tau_vector=tau_mcN95, LHS_vector=LHS_mcN95, 
	# 		eps_vector=eps_N95, expname=expname)		
	# }

	# call function to calculate error bars for enhancement factors
	EF <- enhancement_factors(N=N,P=P, Deffective = Deffective, expname=expname)


	# call function to calculate error bars for time_delay to
	# time_decay ratio
	delaydecayTimes <- tdelay_tdecay(do_mcN95=do_mcN95, 
						K_mcN95=K_mcN95, Yo_mcN95=Yo_mcN95, expname=expname)

	# add current row of data to overall dataframe
	row_data <- c(importedData$FLEX_ID[i],
				D_eff$DN95_lower, D_eff$DN95_bar,D_eff$DN95_upper, D_eff$DN95_most_probable,
				D_eff$DUS_lower, D_eff$DUS_bar, D_eff$DUS_upper, 
				delaydecayTimes$tdtv_lower, delaydecayTimes$tdtv_bar, delaydecayTimes$tdtv_upper,
				EF$EF_lower, EF$EF_bar, EF$EF_upper,
				delaydecayTimes$tv_lower, delaydecayTimes$tv_bar, delaydecayTimes$tv_upper,
				delaydecayTimes$td_lower, delaydecayTimes$td_bar, delaydecayTimes$td_upper,
				D_eff$eps_lower, D_eff$eps_bar, D_eff$eps_upper)

	data_table[i, ] <- row_data

	copy_command <- paste0("cp fcprops.txt"," ",importedData$FLEX_ID[i],"_fcprops.txt")
	system(copy_command)

	print(paste0("============ END data for ",importedData$FLEX_ID[i]," ============") )
	print(" ")
	# readline('--press enter--')
}

write.csv(data_table, file="data_table_exported.csv")

#create directory to store output plots and text files (if applicable)
createdir_command <- paste0("mkdir output_folder")
system(createdir_command)

#move all output files to output_folder directory
move_command <- paste0("mv *.pdf *.txt *.csv *.png output_folder/")
system(move_command)

graphics.off()


