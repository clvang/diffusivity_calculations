
mc_analysis <- function(N, P, expname)
{
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

	#set.seed(5)
	Yo_mcN95 <- rnorm( n=N, mean=yo, sd=sigma_Y )

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

	# -----------------
	# U_tau <- 2.0*sqrt( (1/d_o)^2*Udo_sq + (1/sqrt(dc_sq))^2*Udc_sq )
	# tau_L <- tau_o - U_tau 
	# tau_R <- tau_o + U_tau
	# if( tau_L < 0 ){
	# 	tau_L <- 0
	# }
	# tau_mcN95 <- rtriangle(n=N, a=tau_L,b=tau_R, c=tau_o)
	# tau_mcN95 <- rnorm(n=N, mean=tau_o, sd=U_tau/2) # can run into tau<0 values when doing this

	# sample larger than needed values for K
	# then throw out-of-range values away
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


	#-------------------------------------------------------------------- ###

	#---- generate UNIFORM random variables for STANDARD uncertianties --- ###
	## convert 95% uncertainties to standard deviations for uniform distributions
	do_max <- d_o + sqrt(Udo_sq)
	do_min <- d_o - sqrt(Udo_sq)
	sigma_do <- (do_max - do_min) / sqrt(12)

	dc_max <- sqrt(dc_sq) + sqrt(Udc_sq)
	dc_min <- sqrt(dc_sq) - sqrt(Udc_sq)
	sigma_dc <- (dc_max - dc_min) / sqrt(12)

	K_max <- K + sqrt(Uk_sq)
	K_min <- K - sqrt(Uk_sq)
	sigma_K <- (K_max - K_min) / sqrt(12)

	Y_max <- yo + sqrt(UYo_sq)
	Y_min <- yo - sqrt(UYo_sq)
	sigma_Y <- (Y_max - Y_min) / sqrt(12)

	#set.seed(5)
	do_temp <- runif(n=N*factor, min=do_min , max=do_max )

	#set.seed(5)
	Yo_mcUS <- runif( n=N, min=Y_min, max=Y_max )
	LHS_mcUS <- y_ofc / Yo_mcUS
	## sample larger than needed values for dc
	## then throw out-of-range values away, element-wise
	dc_mcUS <- 1
	i <- 5
	while( length(dc_mcUS) < N){
		#set.seed(i)
		dc_mcUS <- runif(n=N*factor,min=dc_min,max=dc_max )
		index <- which(do_temp > dc_mcUS)
		dc_mcUS <- dc_mcUS[index]
		i <- i + 1
	}
	dc_mcUSseed <- i-1
	dc_mcUS <- dc_mcUS[1:N]
	do_mcUS <- do_temp[index]
	do_mcUS <- do_mcUS[1:N]
	tau_mcUS <- log(do_mcUS / dc_mcUS)
	## sample K s.t. K > 0
	if (K_min < 0){
		K_min <- 0
	}
	#set.seed(5)
	K_mcUS <- runif(n=N, min=K_min, max=K_max)
	print("----- END random number generation----- ")

	#--------- calculate D values by solving asymptotic equation directly ------- ###
	# -------- output is just a single value is NOT an MC analysis value  ------- ###
	results_nonMC <- calcRoot(tau_vector = tau_o, 
						LHS_vector = y_ofc/yo, 
						K_vector = K,
						err_tol = err_tol, 
						N = 1)
	eps_nonMC <- results_nonMC$eps_mcVals
	D_nonMC <- (results_nonMC$dc_mcVals)* (1/1000)^2  #convert to m^2/s


	#---------------- calculate normal 95% uncertainties ------------------------ ###
	results_N95 <- calcRoot(tau_vector=tau_mcN95, 
						LHS_vector=LHS_mcN95, 
						K_vector=K_mcN95,
						err_tol =err_tol, 
						N = N)

	D_N95 <- results_N95$dc_mcVals

	DN95_skewness <- skewness(D_N95)
	if( abs(DN95_skewness) > 0.6 ){
		print("Distribution for D is too skewed...")
		print("generating RNs for K from triangular distribution...")
		# generate K from triangular distribution
		K_low <- K - sqrt(Uk_sq)
		K_upper <- K + sqrt(Uk_sq)
		if( K_low < 0){
			print("min of K RNs is negative...")
			print("setting lower RN limit = 0...")
			K_low <- 0
		}

		K_mcN95 <- rtriangle(n=N, a=K_low, b=K_upper, c=K)

		results_N95 <- calcRoot(tau_vector=tau_mcN95, 
							LHS_vector=LHS_mcN95, 
							K_vector=K_mcN95,
							err_tol =err_tol, 
							N = N)		
		D_N95 <- results_N95$dc_mcVals
		DN95_skewness <- skewness(D_N95)		
	}

	D_N95 <- D_N95* (1/1000)^2  #convert to m^2/s
	DN95_lower <- quantile(D_N95,probs=c( (1-P)/2,(1+P)/2), type=1)[[1]]
	DN95_upper <- quantile(D_N95,probs=c( (1-P)/2,(1+P)/2), type=1)[[2]] 
	DN95_bar <- mean( D_N95 )
	d <- density(D_N95)
	DN95_most_probable <- d$x[which(d$y==max(d$y))]


	print("********** Effective Diffusivity Monte Carlo Results *****************")
	print("--- results assume independent variables are from NORMAL distribution ---")
	print("MC 95% uncertainty LOWER D limit [m^2/s]: ")
	print(DN95_lower)
	print("MC 95% uncertainty MEAN value D [m^2/s]: ")
	print(DN95_bar)
	print("MC 95% uncertainty UPPER D limit [m^2/s]: ")
	print(DN95_upper)
	print("MC 95% uncertainty MOST PROBABLE D [m^2/s]: ")
	print(DN95_most_probable)

	# plot distribution in D for normal RNs
	dev.new()
	pdf(paste0(expname,"_Ddist.pdf") )
	hist(D_N95,prob=TRUE,n=100,  
		main=paste0(expname,": Distribution of D, skw=",signif(DN95_skewness,digits=3) ),
		xlab="D [m^2/s]", col="lightgreen")	
	d <- density(D_N95)
	lines(d,col="black",lwd=2)
	abline(v=DN95_upper,col='red',lwd=2.3,lty="dashed")
	abline(v=DN95_lower,col='red',lwd=2.3,lty="dotted")
	abline(v=DN95_bar,col='red',lwd=2.9)
	legend("topright", c("MC"), col=c("red"), lwd=2)

	# plot distribution in tau 
	tau_lower <- quantile(tau_mcN95,probs=c((1-P)/2,(1+P)/2), type=1)[[1]]
	tau_upper <- quantile(tau_mcN95,probs=c((1-P)/2,(1+P)/2), type=1)[[2]] 
	tau_bar <- mean(tau_mcN95)
	d <- density(tau_mcN95)
	tau_most_probable <- d$x[which(d$y==max(d$y))]	

	dev.new()
	pdf(paste0(expname,"_Taudist.pdf") )
	hist(tau_mcN95,prob=TRUE,n=100,  
		main=paste0(expname,": Distribution of tau"),
		xlab="tau", col="lightgreen")	
	lines(d,col="black",lwd=2)
	abline(v=tau_upper,col='red',lwd=2.3,lty="dashed")
	abline(v=tau_lower,col='red',lwd=2.3,lty="dotted")
	abline(v=tau_bar,col='red',lwd=2.9)
	legend("topright", c("MC"), col=c("red"), lwd=2)	

	# plot distribution in K 
	K_lower <- quantile(K_mcN95,probs=c((1-P)/2,(1+P)/2), type=1)[[1]]
	K_upper <- quantile(K_mcN95,probs=c((1-P)/2,(1+P)/2), type=1)[[2]] 
	K_bar <- mean(K_mcN95)
	d <- density(K_mcN95)
	K_most_probable <- d$x[which(d$y==max(d$y))]	

	dev.new()
	pdf(paste0(expname,"_Kdist.pdf") )
	hist(K_mcN95,prob=TRUE,n=100,  
		main=paste0(expname,": Distribution of K"),
		xlab="K", col="lightgreen")	
	lines(d,col="black",lwd=2)
	abline(v=K_upper,col='red',lwd=2.3,lty="dashed")
	abline(v=K_lower,col='red',lwd=2.3,lty="dotted")
	abline(v=K_bar,col='red',lwd=2.9)
	legend("topright", c("MC"), col=c("red"), lwd=2)	

	# plot distribution in epsilon values
	eps_N95 <- results_N95$eps_mcVals
	eps_lower <- quantile(eps_N95,probs=c((1-P)/2,(1+P)/2), type=1)[[1]]
	eps_upper <- quantile(eps_N95,probs=c((1-P)/2,(1+P)/2), type=1)[[2]] 
	eps_bar <- mean(eps_N95)
	d <- density(eps_N95)
	eps_most_probable <- d$x[which(d$y==max(d$y))]	

	dev.new()
	pdf(paste0(expname,"_epsdist.pdf") )
	hist(eps_N95,prob=TRUE,n=100,  
		main=paste0(expname,": Distribution of Epsilon"),
		xlab="Epsilon", col="lightgreen")	
	lines(d,col="black",lwd=2)
	abline(v=eps_upper,col='red',lwd=2.3,lty="dashed")
	abline(v=eps_lower,col='red',lwd=2.3,lty="dotted")
	abline(v=eps_bar,col='red',lwd=2.9)
	legend("topright", c("MC"), col=c("red"), lwd=2)	


	#---------------- calculate uniform standard error ------------------------ ###
	D_US <- calcRoot(tau_vector=tau_mcUS, 
						LHS_vector=LHS_mcUS, 
						K_vector=K_mcUS,
						err_tol =err_tol, 
						N = N)$dc_mcVals
	D_US <- D_US * (1/1000)^2  #convert to m^2/s
	DUS_bar <- mean(D_US)
	d <- density(D_US)
	DUS_upper <- quantile(D_US,probs=c((1-P)/2,(1+P)/2), type=1)[[1]] #DUS_bar + sigma_DUS #(sigma_DUS/sqrt(N))
	DUS_lower <- quantile(D_US,probs=c((1-P)/2,(1+P)/2), type=1)[[2]] #DUS_bar - sigma_DUS #(sigma_DUS/sqrt(N))

	print("********************************************************************")
	print("--- results assume independent variables are from UNIFORM distribution ---")
	print("MC lower bound D limit [m^2/s]: ")
	print(DUS_lower)
	print("MC MEAN value D [m^2/s]: ")
	print(DUS_bar )
	print("MC upper bound D limit [m^2/s]: ")
	print(DUS_upper)


	list(eps_nonMC = eps_nonMC,
		D_nonMC = D_nonMC,
		DN95_lower = DN95_lower, #---- required to output to screen ---- #
		DN95_bar = DN95_bar, 	  #					|
		DN95_upper = DN95_upper,  #					|
		DN95_most_probable = DN95_most_probable,#   |
		DUS_lower = DUS_lower,
		DUS_bar = DUS_bar,		  #					|
		DUS_upper = DUS_upper,
		D_N95 = D_N95,  		 #--- required for calculations of tdtv ratio ---#
		do_mcN95= do_mcN95,      #             |
		K_mcN95=K_mcN95,  		 # 			   |
		Yo_mcN95=Yo_mcN95,		 #-----------------------------------------------
		eps_N95 = eps_N95,
		eps_lower = eps_lower,
		eps_bar = eps_bar,
		eps_upper = eps_upper,
		LHS_mcN95 = LHS_mcN95,
		tau_mcN95 = tau_mcN95)       

}


