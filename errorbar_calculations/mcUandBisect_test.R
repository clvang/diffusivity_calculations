

rm(list=ls(all=TRUE))   #remove all variables in  workspace

	calcRoot <- function(tau_vector, LHS_vector, K_vector, err_tol, N){
	# SUBROUTINE mc_uncertainty( tau_vector, LHS_vector, err_tol, N, D_vector, eps_vector )
		retvals <- .Fortran("mc_uncertainty",
					tau_vector=as.double(tau_vector),
					LHS_vector=as.double(LHS_vector),
					K_vector = as.double(K_vector),
					err_tol = as.double(err_tol),
					N = as.integer(N),
					D_vector = double(N),
					eps_vector = double(N) )
		#there are the returned values
		list(dc_mcVals=retvals$D_vector,
		eps_mcVals = retvals$eps_vector)         #values for D
	}

	if (!is.loaded('mc_uncertainty')){
		dyn.load("mc_uncertainty.so")
	}

	library(triangle)
	library(moments)
	source('/Users/changvang/mygitFiles/diffusivity_calculations/errorbar_calculations/read_props.R')
	source('/Users/changvang/mygitFiles/diffusivity_calculations/errorbar_calculations/YoCorrection_Function.R')
	source('/Users/changvang/mygitFiles/diffusivity_calculations/errorbar_calculations/functionPlot.R')

	N <- 1000000
	P <- 0.05
	expname <- "testername"
	# # read data from fcprops.txt file
	experimental_parameters <- read_props()
	p <- experimental_parameters$p 
	dc_sq <- experimental_parameters$dc_sq 
	K <- experimental_parameters$K 
	do_measured <- experimental_parameters$do_measured 
	yo <- experimental_parameters$yo 
	err_tol <- experimental_parameters$err_tol
	Uk_sq <- experimental_parameters$Uk_sq
	Udo_sq <- experimental_parameters$Udo_sq
	Udc_sq <- experimental_parameters$Udc_sq
	UYo_sq <- experimental_parameters$UYo_sq
	sol_id <- experimental_parameters$sol_id
	dod1_ratio <- experimental_parameters$dod1_ratio

	y_ofc <- 1.0  # (defined) mass fraction of low volatility comp at onset of fc

	if (sol_id == 1){

		if ( yo == 0.05){
			# heptane95-hexadecane5 d_o corrections
			if (p == 1) {
			    percent_increase <- 0.03  
			} else {
			    percent_increase <- 0.06
			}	
		}
		if ( yo == 0.20){
			# heptane80-hexadecane20 d_o corrections
			if (p == 1) {
			    percent_increase <- 0.03  
			} else {
			    percent_increase <- 0.06
			}	
		}

	}

	if (sol_id == 2){
		# propanol-glycerol d_o corrections
		if (p == 1) {
		    percent_increase <- 0.03  #
		} else {
		    percent_increase <- 0.05  #at 3 atm
		}
	}


	d_o <- do_measured * (1.0 + percent_increase) #do accounting for droplet swelling
	tau_o <- log(d_o / sqrt(dc_sq) )

	# NOTE:
	# 	tau = LOG(d_o / d_c)
	# 	LHS = y_ofc / yo

	print("----- BEGIN random number generation----- ")
	#---- generate NORMAL random variables for 95% uncertianties --- ###
	## convert 95% uncertainties to standard deviations for normal distributions
	cf <- 2.0
	sigma_do <- sqrt(Udo_sq) / cf
	sigma_dc <- sqrt(Udc_sq) / cf
	sigma_K <- sqrt(Uk_sq) / cf
	sigma_Y <- sqrt(UYo_sq) / cf	

	# THIS IS WHERE YOU WOULD PASS IN dod1_ratio and yo into a function
	# which solves for Yo_mcN95 (the corrected mass fraction of the low
	# volatility component when the igniters are turned on)
	# readline(prompt="Generating RNs for Yo: ")
	Yo_mcN952 <- rnorm( n=N, mean=yo, sd=sigma_Y )
	print('Correcting for changes in Yo...')
	Yo_mcN95 <- YoCorrection_Function(N=N, cf=cf, sol_id=sol_id,
									  dod1_ratio=dod1_ratio,
									  yo=yo, UYo_sq=UYo_sq,
									  Udo_sq=Udo_sq,
									  do_measured=do_measured)
	print("completed correction for Yo...")

	hist(Yo_mcN95,prob=T,n=100)
	quartz()
	hist(Yo_mcN952,prob=T,n=100)
	readline(prompt="--pause--")
	# readline(prompt="Completed Genearting RNs for Yo: ")

	LHS_mcN95 <- y_ofc / Yo_mcN95

	## sample larger than needed values for dc and do
	## then throw out-of-range values away, element-wise
	factor <- 2  #this just specifies factor to increase number of random numbers genrated
	#set.seed(5)
	do_temp <- rnorm( n=N*factor, mean=d_o, sd= sigma_do )
	dc_mcN95 <- 1
	i <- 5
	while(length(dc_mcN95) < N){
		#set.seed(i)
		dc_mcN95 <- rnorm(n=N*factor, mean = sqrt(dc_sq), sd=sigma_dc )
		index <- which(do_temp > dc_mcN95)
		dc_mcN95 <- dc_mcN95[index]
		i <- i + 1
	}
	dc_mcN95seed <- i-1   #store seed for future debug/reproducibility purposes

	dc_mcN95 <- dc_mcN95[1:N]   #select only N values from the within range values
	do_mcN95 <- do_temp[index]
	do_mcN95 <- do_mcN95[1:N]
	tau_mcN95 <- log(do_mcN95 / dc_mcN95) # tau > 0

	K_mcN95 <- 1
	i <- 5
	while( length(K_mcN95) < N ){
		#set.seed(i)	
		K_mcN95 <- rnorm(n=N*factor, mean= K, sd= sigma_K)
		K_mcN95 <- K_mcN95[which( K_mcN95 > 0)]
		i <- i + 1	
	}
	K_mcN95seed <- i-1
	K_mcN95 <- K_mcN95[1:N]


	print("----- END random number generation----- ")

	#--------- calculate D values by solving asymptotic equation directly ------- ###
	# -------- output is just a single value is NOT an MC analysis value  ------- ###
	print("Solving for single point results...")
	results_nonMC <- calcRoot(tau_vector = tau_o, 
						LHS_vector = y_ofc/yo, 
						K_vector = K,
						err_tol = err_tol, 
						N = 1)
	eps_nonMC <- results_nonMC$eps_mcVals
	D_nonMC <- (results_nonMC$dc_mcVals)* (1/1000)^2  #convert to m^2/s
	print("Completed solving for single point result...")

	#---------------- calculate normal 95% uncertainties ------------------------ ###
	print("Solving for normal 95% results...")	
	# readline(prompt="Solving for normal 95% results...")
	results_N95 <- calcRoot(tau_vector=tau_mcN95, 
						LHS_vector=LHS_mcN95, 
						K_vector=K_mcN95,
						err_tol =err_tol, 
						N = N)

	D_N95 <- results_N95$dc_mcVals
	# readline(prompt="Completed solving for normal 95% results...")		
	print("Completed solving for normal 95% results...")		
	readline('---paused---')


	functionPlot(tau_vector=tau_mcN95, LHS_vector=LHS_mcN95, expname=expname)

	# DN95_skewness <- skewness(D_N95)
	# if( abs(DN95_skewness) > 0.6 ){
	# 	print("Distribution for D is too skewed...")
	# 	print("generating RNs for K from triangular distribution...")
	# 	# generate K from triangular distribution
	# 	K_low <- K - sqrt(Uk_sq)
	# 	K_upper <- K + sqrt(Uk_sq)
	# 	if( K_low < 0){
	# 		print("min of K RNs is negative...")
	# 		print("setting lower RN limit = 0...")
	# 		K_low <- 0
	# 	}

	# 	K_mcN95 <- rtriangle(n=N, a=K_low, b=K_upper, c=K)

	# 	results_N95 <- calcRoot(tau_vector=tau_mcN95, 
	# 						LHS_vector=LHS_mcN95, 
	# 						K_vector=K_mcN95,
	# 						err_tol =err_tol, 
	# 						N = N)		
	# 	D_N95 <- results_N95$dc_mcVals
	# 	DN95_skewness <- skewness(D_N95)		
	# }

	# D_N95 <- D_N95* (1/1000)^2  #convert to m^2/s
	# DN95_lower <- quantile(D_N95,probs=c( (1-P)/2,(1+P)/2), type=1)[[1]]
	# DN95_upper <- quantile(D_N95,probs=c( (1-P)/2,(1+P)/2), type=1)[[2]] 
	# DN95_bar <- mean( D_N95 )
	# d <- density(D_N95)
	# DN95_most_probable <- d$x[which(d$y==max(d$y))]


	# print("********** Effective Diffusivity Monte Carlo Results *****************")
	# print("--- results assume independent variables are from NORMAL distribution ---")
	# print("MC 95% uncertainty LOWER D limit [m^2/s]: ")
	# print(DN95_lower)
	# print("MC 95% uncertainty MEAN value D [m^2/s]: ")
	# print(DN95_bar)
	# print("MC 95% uncertainty UPPER D limit [m^2/s]: ")
	# print(DN95_upper)
	# print("MC 95% uncertainty MOST PROBABLE D [m^2/s]: ")
	# print(DN95_most_probable)

	# # plot distribution in D for normal RNs
	# dev.new()
	# pdf(paste0(expname,"_Ddist.pdf") )
	# hist(D_N95,prob=TRUE,n=100,  
	# 	main=paste0(expname,": Distribution of D, skw=",signif(DN95_skewness,digits=3) ),
	# 	xlab="D [m^2/s]", col="lightgreen")	
	# d <- density(D_N95)
	# lines(d,col="black",lwd=2)
	# abline(v=DN95_upper,col='red',lwd=2.3,lty="dashed")
	# abline(v=DN95_lower,col='red',lwd=2.3,lty="dotted")
	# abline(v=DN95_bar,col='red',lwd=2.9)
	# legend("topright", c("MC"), col=c("red"), lwd=2)

	# # plot distribution in tau 
	# tau_lower <- quantile(tau_mcN95,probs=c((1-P)/2,(1+P)/2), type=1)[[1]]
	# tau_upper <- quantile(tau_mcN95,probs=c((1-P)/2,(1+P)/2), type=1)[[2]] 
	# tau_bar <- mean(tau_mcN95)
	# d <- density(tau_mcN95)
	# tau_most_probable <- d$x[which(d$y==max(d$y))]	

	# dev.new()
	# pdf(paste0(expname,"_Taudist.pdf") )
	# hist(tau_mcN95,prob=TRUE,n=100,  
	# 	main=paste0(expname,": Distribution of tau"),
	# 	xlab="tau", col="lightgreen")	
	# lines(d,col="black",lwd=2)
	# abline(v=tau_upper,col='red',lwd=2.3,lty="dashed")
	# abline(v=tau_lower,col='red',lwd=2.3,lty="dotted")
	# abline(v=tau_bar,col='red',lwd=2.9)
	# legend("topright", c("MC"), col=c("red"), lwd=2)	

	# # plot distribution in K 
	# K_lower <- quantile(K_mcN95,probs=c((1-P)/2,(1+P)/2), type=1)[[1]]
	# K_upper <- quantile(K_mcN95,probs=c((1-P)/2,(1+P)/2), type=1)[[2]] 
	# K_bar <- mean(K_mcN95)
	# d <- density(K_mcN95)
	# K_most_probable <- d$x[which(d$y==max(d$y))]	

	# dev.new()
	# pdf(paste0(expname,"_Kdist.pdf") )
	# hist(K_mcN95,prob=TRUE,n=100,  
	# 	main=paste0(expname,": Distribution of K"),
	# 	xlab="K", col="lightgreen")	
	# lines(d,col="black",lwd=2)
	# abline(v=K_upper,col='red',lwd=2.3,lty="dashed")
	# abline(v=K_lower,col='red',lwd=2.3,lty="dotted")
	# abline(v=K_bar,col='red',lwd=2.9)
	# legend("topright", c("MC"), col=c("red"), lwd=2)	

	# # plot distribution in epsilon values
	# eps_N95 <- results_N95$eps_mcVals
	# eps_lower <- quantile(eps_N95,probs=c((1-P)/2,(1+P)/2), type=1)[[1]]
	# eps_upper <- quantile(eps_N95,probs=c((1-P)/2,(1+P)/2), type=1)[[2]] 
	# eps_bar <- mean(eps_N95)
	# d <- density(eps_N95)
	# eps_most_probable <- d$x[which(d$y==max(d$y))]	

	# dev.new()
	# pdf(paste0(expname,"_epsdist.pdf") )
	# hist(eps_N95,prob=TRUE,n=100,  
	# 	main=paste0(expname,": Distribution of Epsilon"),
	# 	xlab="Epsilon", col="lightgreen")	
	# lines(d,col="black",lwd=2)
	# abline(v=eps_upper,col='red',lwd=2.3,lty="dashed")
	# abline(v=eps_lower,col='red',lwd=2.3,lty="dotted")
	# abline(v=eps_bar,col='red',lwd=2.9)
	# legend("topright", c("MC"), col=c("red"), lwd=2)	


graphics.off()

