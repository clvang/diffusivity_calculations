#fucntion to calculate uncertainties for infinite dilution
#molecular diffusivities, DAB or DBA, using Monte-Carlo method


tynCalus_MC <- function(N,P,T,uT,sigmaA,MW_A,rhoA,u_sigmaA,u_mwA,u_rhoA,
						sigmaB,MW_B,rhoB,etaB,u_sigmaB,u_mwB,u_rhoB,u_etaB, pltTitle)
{
	cf <- 2.0
	set.seed(5)
	#NOTE: convert all uncertainties to standard deviations by 
	#dividing by 2.0 for RNs from normal distribution
	T <- rnorm(n=N,mean=T,sd=uT/cf)
	MW_A <- rnorm(n=N,mean=MW_A,sd=u_mwA/cf)
	MW_B <- rnorm(n=N,mean=MW_B,sd=u_mwB/cf)
	rhoA <- rnorm(n=N,mean=rhoA,sd=u_rhoA/cf)
	rhoB <- rnorm(n=N,mean=rhoB,sd=u_rhoB/cf)
	etaB <- rnorm(n=N,mean=etaB,sd=u_etaB/cf)

	# if surface tension data is available, then
	# go ahead and generate random variables
	if ( sigmaA != 0 ){  
		sigmaA <- rnorm(n=N,mean=sigmaA,sd=u_sigmaA/cf)
		sigmaB <- rnorm(n=N,mean=sigmaB,sd=u_sigmaB/cf)		
	} 

	#use Tyn-Calus equation to evalute D_{AB} or D_{BA}
	if ( length(sigmaA) == N ){  # if surface tension info is available (Heptane/Hexadecane experiments)
		C1 <- 8.93E-8
		DAB <- C1 * ((MW_B*1E3/rhoB)^0.267)/((MW_A*1E3/rhoA)^0.433) * 
				(T/etaB) * (sigmaB/sigmaA)^(0.15) / 1E4 			# m^2/s
	}else{   # if surface tension data is NOT available (Propanol/Glycerol experiment)

		# this form of the Tyn-Calus equation assumes that
		# surface tension of the solvent and solute are essentiall
		# equal (See Reid, Poling, Prausnitz, Sec. 11.25)
		C1 <- 8.93E-8
		DAB <- C1 * ((MW_B*1E3/rhoB)^0.267)/((MW_A*1E3/rhoA)^0.433) * 
				(T/etaB) * (1.0)^(0.15) / 1E4 			# m^2/s using second form of Tyn-Calus
														# since surface tension info is not available
														# for n-Propanol at BP :(	
	}


	# calculate lower and upper 95% confidence bands, and average value
	DAB_lower <- quantile(DAB,probs=c((1-P)/2,(1+P)/2))[[1]]
	DAB_upper <- quantile(DAB,probs=c((1-P)/2,(1+P)/2))[[2]]
	DAB_bar <- mean(DAB)
	d <- density(DAB)
	DAB_most_probable <- d$x[which(d$y==max(d$y))]

	print(paste0("MC 95% lower D limit: ", DAB_lower) )
	print(paste0("MC most probable D: ", DAB_most_probable))
	print(paste0("MC average value D: ", DAB_bar))
	print(paste0("MC 95% upper D limit: ", DAB_upper))

	# dev.new()
	# # pdf(paste0(pltTitle,"_MCHistogram.pdf"))
	# hist(DAB,prob=TRUE,n=100,xlim=c(0,3.5)*DAB_most_probable,
	# 	main=paste0("Histogram of ",pltTitle),
	# 	xlab=pltTitle)
	# d <-density(DAB)
	# lines(d,col="black",lwd=2)
	# abline(v=DAB_upper,col='green',lwd=1.5,lty="dashed")
	# abline(v=DAB_lower,col='green',lwd=1.5,lty="dotted")
	# abline(v=DAB_most_probable,col='green',lwd=2.3)
	# abline(v=DAB_bar,col='green',lwd=3.2)
	# legend("topright", c("MC", "TSM"), col=c("green", "red"), lwd=2)

	list(Dmolecular = DAB)

}