

 enhancement_factors <- function(N, P, Deffective, expname){

	source('/Users/changvang/mygitFiles/diffusivity_calculations/errorbar_calculations/tynCalus.R')


	# call function to solve for infinite dilution molecular diffusivities
	# for Heptane/Hexadecane experiments
	Dmolecular <- tynCalus( N, P , expname)
	DAB <- Dmolecular$DAB 
	DBA <- Dmolecular$DBA

	# calculate enhancement factors
	EF <- (Deffective)/(0.5*(DAB + DBA) )


	EF_lower <- quantile(EF,probs=c((1-P)/2,(1+P)/2), type=1)[[1]]
	EF_upper <- quantile(EF,probs=c((1-P)/2,(1+P)/2), type=1)[[2]] 
	EF_bar <- mean(EF )
	d <- density(EF)
	EF_most_probable <- d$x[which(d$y==max(d$y))]


	print("********** Enhancement Factors Monte Carlo Results *****************")
	print("--- results assume independent variables are from NORMAL distribution ---")
	print("MC 95% uncertainty LOWER EF limit [m^2/s]: ")
	print(EF_lower)
	print("MC 95% uncertainty MEAN value EF [m^2/s]: ")
	print(EF_bar)
	print("MC 95% uncertainty UPPER EF limit [m^2/s]: ")
	print(EF_upper)
	print("MC 95% uncertainty MOST PROBABLE EF [m^2/s]: ")
	print(EF_most_probable)

	dev.new()
	pdf(paste0(expname,"_EFdist.pdf") )
	hist(EF,prob=TRUE,n=100,  
		main=paste0(expname, ": Distribution of EF"),
		xlab="EF", col="lightgreen")
	d <-density(EF)
	lines(d,col="black",lwd=2)
	abline(v=EF_upper,col='red',lwd=2.3,lty="dashed")
	abline(v=EF_lower,col='red',lwd=2.3,lty="dotted")
	abline(v=EF_bar,col='red',lwd=2.9)
	legend("topright", c("MC"), col=c("red"), lwd=2)

	list(EF = EF,
		EF_lower = EF_lower,
		EF_bar = EF_bar,
		EF_upper = EF_upper)
}
