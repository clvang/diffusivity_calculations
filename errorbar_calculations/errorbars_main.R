
rm(list=ls(all=TRUE))   #remove all variables in  workspace
graphics.off()  #close all graphics windows


source('/Users/changlvang/mygitFiles/diffusivity_calculations/enhancement_factors/mc_analysis.R')
source('/Users/changlvang/mygitFiles/diffusivity_calculations/enhancement_factors/enhancement_factors.R')
source('/Users/changlvang/mygitFiles/diffusivity_calculations/enhancement_factors/tdelay_tdecay.R')

P <- 0.95
N <- 1000000
# call function to calculate error bars for effective diffusivities
D_eff <- mc_analysis( N=N, P=P )

# the folllowing values are needed as inputs
# to the function tdelay_tdecay
Deffective <- D_eff$D_N95
do_mcN95 <- D_eff$do_mcN95
K_mcN95 <- D_eff$K_mcN95
Yo_mcN95 <- D_eff$Yo_mcN95

# call function to calculate error bars for enhancement factors
EF <- enhancement_factors(N=N,P=P)$EF

# call function to calculate error bars for tdelay/tdecay ratio
delaydecayTimes <- tdelay_tdecay(do_mcN95=do_mcN95, 
					K_mcN95=K_mcN95, Yo_mcN95=Yo_mcN95)

td_mcN95 <- delaydecayTimes$td_mcN95
d <- density(td_mcN95)
td_most_probable <- d$x[which(d$y==max(d$y))]
td_lower <- quantile(td_mcN95,probs=c((1-P)/2,(1+P)/2))[[1]]
td_upper <- quantile(td_mcN95,probs=c((1-P)/2,(1+P)/2))[[2]] 
td_bar <- mean(td_mcN95 )

dev.new()
hist(td_mcN95,prob=TRUE,n=100,  
	main=paste0("Distribution of t_delay"),
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
td_most_probable <- d$x[which(d$y==max(d$y))]
td_lower <- quantile(tv,probs=c((1-P)/2,(1+P)/2))[[1]]
td_upper <- quantile(tv,probs=c((1-P)/2,(1+P)/2))[[2]] 
td_bar <- mean(tv )

dev.new()
hist(tv,prob=TRUE,n=100,  
	main=paste0("Distribution of t_decay"),
	xlab=expression("t"[decay]), col="lightgreen")
d <-density(tv)
lines(d,col="black",lwd=2)
abline(v=td_upper,col='red',lwd=2.3,lty="dashed")
abline(v=td_lower,col='red',lwd=2.3,lty="dotted")
abline(v=td_bar,col='red',lwd=2.9)
abline(v=td_most_probable,col='black',lwd=2.9)
legend("topright", c("MC"), col=c("red"), lwd=2)



