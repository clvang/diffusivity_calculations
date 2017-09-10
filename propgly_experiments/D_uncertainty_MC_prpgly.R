

D_uncertainty_MC <- function(N,P,T,uT,MW_A,rhoA,u_mwA,u_rhoA,
						MW_B,rhoB,etaB,u_mwB,u_rhoB,u_etaB,pltTitle)
{

	#Calculate D_AB (or D_BA) using Monte Carlo
	T <- rnorm(n=N,mean=T,sd=uT)
	MW_A <- rnorm(n=N,mean=MW_A,sd=u_mwA)
	MW_B <- rnorm(n=N,mean=MW_B,sd=u_mwB)
	rhoA <- rnorm(n=N,mean=rhoA,sd=u_rhoA)
	rhoB <- rnorm(n=N,mean=rhoB,sd=u_rhoB)
	etaB <- rnorm(n=N,mean=etaB,sd=u_etaB)

	C1 <- 8.93E-8
	DAB <- C1 * ((MW_B*1E3/rhoB)^0.267)/((MW_A*1E3/rhoA)^0.433) * 
			(T/etaB) * (1.0)^(0.15) / 1E4 			# m^2/s using second form of Tyn-Calus
													# since surface tension info is not available
													# for n-Propanol at BP :(
	DAB_lower <- quantile(DAB,probs=c((1-P)/2,(1+P)/2))[[1]]
	DAB_upper <- quantile(DAB,probs=c((1-P)/2,(1+P)/2))[[2]]
	DAB_bar <- mean(DAB)
	d <- density(DAB)
	DAB_most_probable <- d$x[which(d$y==max(d$y))]

	print(paste0("MC lower D limit: ", DAB_lower) )
	print(paste0("MC most probable D: ", DAB_most_probable))
	print(paste0("MC upper D limit: ", DAB_upper))

	dev.new()
	pdf(paste0(pltTitle,"_MCHistogram.pdf"))
	hist(DAB,prob=TRUE,n=100,xlim=c(0,3.5)*DAB_most_probable,
		main=paste0("Histogram of ",pltTitle),
		xlab=pltTitle)
	d <-density(DAB)
	lines(d,col="black",lwd=2)
	abline(v=DAB_upper,col='green',lwd=1.5,lty="dashed")
	abline(v=DAB_lower,col='green',lwd=1.5,lty="dotted")
	abline(v=DAB_most_probable,col='green',lwd=2.3)
	legend("topright", c("MC"), col=c("green", "red"), lwd=2)
	dev.off()

}