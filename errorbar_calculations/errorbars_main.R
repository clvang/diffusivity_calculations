
rm(list=ls(all=TRUE))   #remove all variables in  workspace
library(readxl)

source('/Users/changlvang/mygitFiles/diffusivity_calculations/errorbar_calculations/mc_analysis.R')
source('/Users/changlvang/mygitFiles/diffusivity_calculations/errorbar_calculations/enhancement_factors.R')
source('/Users/changlvang/mygitFiles/diffusivity_calculations/errorbar_calculations/tdelay_tdecay.R')
source('/Users/changlvang/mygitFiles/diffusivity_calculations/errorbar_calculations/generate_fcprops_from_xlsx.R')

P <- 0.95
N <- 1000000

importedData <- read_excel("FLEX_DATA_MASTER_Rev.05_test.xlsx",sheet="HepHexPropGlyCSV",
					col_names=TRUE)
data_table <- data.frame( matrix(ncol = 21, nrow = nrow(importedData)) )
data_table <- setNames(data_table, c("FLEX_ID","DmcLow_N95","Dmc_AV_N95",
									"DmcHigh_N95","Dmc_MP_N95",
									"DNS_bar","DNS_SD", "DUS_bar",
									"DUS_SD", "tdtv_low95", "tdtv_avg",
									"tdtv_high95", "EF_low95","EF_avg","EF_high95",
									"tdecay_low95","tdecay_avg","tdecay_high95,",
									"tdelay_low95","tdelay_avg","tdelay_high95") )

for (i in 1:nrow(importedData)){
# for (i in 53:53){	
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


	# call function to calculate error bars for enhancement factors
	EF <- enhancement_factors(N=N,P=P, Deffective = Deffective, expname=expname)


	# call function to calculate error bars for time_delay to
	# time_decay ratio
	delaydecayTimes <- tdelay_tdecay(do_mcN95=do_mcN95, 
						K_mcN95=K_mcN95, Yo_mcN95=Yo_mcN95, expname=expname)

	td_mcN95 <- delaydecayTimes$td_mcN95
	d <- density(td_mcN95)
	td_most_probable <- d$x[which(d$y==max(d$y))]
	td_lower <- quantile(td_mcN95,probs=c((1-P)/2,(1+P)/2))[[1]]
	td_upper <- quantile(td_mcN95,probs=c((1-P)/2,(1+P)/2))[[2]] 
	td_bar <- mean(td_mcN95 )

	dev.new()
	pdf(paste0(expname,"_tdelaydist.pdf") )
	hist(td_mcN95,prob=TRUE,n=100,  
		main=paste0(expname,": Distribution of t_delay"),
		xlab=expression("t"[delay]), col="lightgreen")
	d <-density(td_mcN95)
	lines(d,col="black",lwd=2)
	abline(v=td_upper,col='red',lwd=2.3,lty="dashed")
	abline(v=td_lower,col='red',lwd=2.3,lty="dotted")
	abline(v=td_bar,col='red',lwd=2.9)
	abline(v=td_most_probable,col='black',lwd=2.9)
	legend("topright", c("MC"), col=c("red"), lwd=2)


	tv <- delaydecayTimes$tv
	d <- density(tv)
	tv_most_probable <- d$x[which(d$y==max(d$y))]
	tv_lower <- quantile(tv,probs=c((1-P)/2,(1+P)/2))[[1]]
	tv_upper <- quantile(tv,probs=c((1-P)/2,(1+P)/2))[[2]] 
	tv_bar <- mean(tv )

	dev.new()
	pdf(paste0(expname,"_tdecaydist.pdf") )	
	hist(tv,prob=TRUE,n=100,  
		main=paste0(expname,": Distribution of t_decay"),
		xlab=expression("t"[decay]), col="lightgreen")
	d <-density(tv)
	lines(d,col="black",lwd=2)
	abline(v=tv_upper,col='red',lwd=2.3,lty="dashed")
	abline(v=tv_lower,col='red',lwd=2.3,lty="dotted")
	abline(v=tv_bar,col='red',lwd=2.9)
	abline(v=tv_most_probable,col='black',lwd=2.9)
	legend("topright", c("MC"), col=c("red"), lwd=2)

	row_data <- c(importedData$FLEX_ID[i],
				D_eff$DN95_lower, D_eff$DN95_bar,D_eff$DN95_upper, D_eff$DN95_most_probable,
				D_eff$DNS_bar, D_eff$sigma_DNS,
				D_eff$DUS_bar, D_eff$sigma_DUS,
				delaydecayTimes$tdtv_lower, delaydecayTimes$tdtv_bar, delaydecayTimes$tdtv_upper,
				EF$EF_lower, EF$EF_bar, EF$EF_upper,
				tv_lower,tv_bar,tv_upper,
				td_lower,td_bar,td_upper)

	data_table[i, ] <- row_data

	# copy_command <- paste0("cp fcprops.txt"," ",importedData$FLEX_ID[i],"_fcprops.txt")
	# system(copy_command)

	print(paste0("============ END data for ",importedData$FLEX_ID[i]," ============") )
	print(" ")

}


write.csv(data_table, file="data_table_exported.csv")

