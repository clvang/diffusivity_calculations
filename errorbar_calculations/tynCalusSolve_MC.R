#fucntion to calculate uncertainties for infinite dilution
#molecular diffusivities, DAB or DBA, using Monte-Carlo method

	# D <- tynCalusSolve_MC(N = N, 
	# 				P = P, 					# P-value
	# 				T = T, 					# K
	# 				uT = uT,				# p/m K
	# 				sigmaA_T = sigmaA_T,	# dynes/cm or erg/cm^2
	# 				MW_A = MW_A,			# kg/kmol or g/gmol
	# 				rhoA_T = rhoA_T,		# kg/m^3
	# 				rhoA_nbp = rhoA_nbp,    # kg/m^3
	# 				u_sigmaA_T = u_sigmaA_T,	# dynes/cm or erg/cm^2
	# 				u_mwA = u_mwA,			# kg/kmol or g/gmol
	# 				u_rhoA_T = u_rhoA_T,		# kg/m^3
	# 				u_rhoA_nbp = u_rhoA_nbp,	# kg/m^3
	# 				sigmaB_T = sigmaB_T,		# dynes/cm or erg/cm^2
	# 				MW_B = MW_B,				# kg/kmol or g/gmol
	# 				rhoB_T = rhoB_T,			# kg/m^3
	# 				rhoB_nbp = rhoB_nbp,		# kg/m^3
	# 				etaB_T = etaB_T,			# cP
	# 				u_sigmaB_T = u_sigmaB_T,	# dynes/cm or erg/cm^2
	# 				u_mwB = u_mwB,				# kg/kmol or g/gmol
	# 				u_rhoB_T = u_rhoB_T,		# kg/m^3
	# 				u_rhoB_nbp = u_rhoB_nbp,	# kg/m^3
	# 				u_etaB_T = u_etaB_T,		# cP
	# 				pltTitle = pltTitle,		
	# 				expname = expname )			

tynCalusSolve_MC <- function(N,			
							P,
							T,
							uT,
							sigmaA_T,
							MW_A,
							rhoA_T,
							rhoA_nbp,
							u_sigmaA_T,
							u_mwA,
							u_rhoA_T,
							u_rhoA_nbp,
							sigmaB_T,
							MW_B,
							rhoB_T,
							rhoB_nbp,
							etaB_T,
							u_sigmaB_T,
							u_mwB,
							u_rhoB_T,
							u_rhoB_nbp,
							u_etaB_T, 
							pltTitle,
							expname )
{
	cf <- 2.0
	set.seed(5)
	#NOTE: convert all uncertainties to standard deviations by 
	#dividing by 2.0 for RNs from normal distribution
	T <- rnorm(n=N,mean=T,sd=uT/cf)
	MW_A <- rnorm(n=N,mean=MW_A,sd=u_mwA/cf)
	MW_B <- rnorm(n=N,mean=MW_B,sd=u_mwB/cf)
	rhoA_T <- rnorm(n=N,mean=rhoA_T,sd=u_rhoA_T/cf)
	rhoB_T <- rnorm(n=N,mean=rhoB_T,sd=u_rhoB_T/cf)
	rhoA_nbp <- rnorm(n=N,mean=rhoA_nbp,sd=u_rhoA_nbp/cf)
	rhoB_nbp <- rnorm(n=N,mean=rhoB_nbp,sd=u_rhoB_nbp/cf)	
	etaB_T <- rnorm(n=N,mean=etaB_T,sd=u_etaB_T/cf)

	sigmaA_T <- rnorm(n=N,mean=sigmaA_T,sd=u_sigmaA_T/cf) 
	sigmaB_T <- rnorm(n=N,mean=sigmaB_T,sd=u_sigmaB_T/cf)		

	# calculate Parachors of each component
	P_A <- (MW_A/ (rhoA_T*1e-3) )*(sigmaA_T)^(1/4)   # ( cm^3-g^(1/4) ) / ( s^(1/2)-mol )
	P_B <- (MW_B/ (rhoB_T*1e-3) )*(sigmaB_T)^(1/4)   # ( cm^3-g^(1/4) ) / ( s^(1/2)-mol )

	# calculate molar volumes of each component at their respective
	# boiling point tempertures
	V_A <- MW_A / (rhoA_nbp*1e-3) 	# cm^3/mol
	V_B <- MW_B / (rhoB_nbp*1e-3)   # cm^3/mol

	#use Tyn-Calus equation to evalute D_{AB} or D_{BA}
	C1 <- 8.93E-8
	DAB <- C1 * ( (V_A^(1/6))/ (V_B^(1/3)) ) * ( (P_B / P_A)^(0.6) ) * (T/etaB_T) # cm^2/s
	DAB <- DAB/1e4  	# m^2/s

	# DAB <- C1 * ((MW_B*1E3/rhoB_T)^0.267)/((MW_A*1E3/rhoA_T)^0.433) * 
	# 		(T/etaB) * (sigmaB_T/sigmaA_T)^(0.15) / 1E4 			# m^2/s

	# calculate lower and upper 95% confidence bands, and average value
	DAB_lower <- quantile(DAB,probs=c((1-P)/2,(1+P)/2), type=1)[[1]]
	DAB_upper <- quantile(DAB,probs=c((1-P)/2,(1+P)/2), type=1)[[2]]
	DAB_bar <- mean(DAB)
	d <- density(DAB)
	DAB_most_probable <- d$x[which(d$y==max(d$y))]

	print(paste0("MC 95% lower D limit: ", DAB_lower) )
	print(paste0("MC most probable D: ", DAB_most_probable))
	print(paste0("MC average value D: ", DAB_bar))
	print(paste0("MC 95% upper D limit: ", DAB_upper))

	dev.new()
	pdf(paste0(expname,"_",pltTitle,"dist.pdf"))
	hist(DAB,prob=TRUE,n=100,
		main=paste0("Histogram of ",pltTitle," :",expname),
		xlab=pltTitle)
	d <-density(DAB)
	lines(d,col="black",lwd=2)
	abline(v=DAB_upper,col='green',lwd=1.5,lty="dashed")
	abline(v=DAB_lower,col='green',lwd=1.5,lty="dotted")
	abline(v=DAB_most_probable,col='green',lwd=2.3)
	abline(v=DAB_bar,col='green',lwd=3.2)
	legend("topright", c("MC", "TSM"), col=c("green", "red"), lwd=2)

	list(Dmolecular = DAB)

}