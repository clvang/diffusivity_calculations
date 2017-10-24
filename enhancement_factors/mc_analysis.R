
mc_analysis <- function(N, P)
{
	calcRoot <- function(tau_vector, LHS_vector, K_vector, err_tol, N){
	# SUBROUTINE mc_uncertainty( tau_vector, LHS_vector, err_tol, N, D_vector )
		retvals <- .Fortran("mc_uncertainty",
					tau_vector=as.double(tau_vector),
					LHS_vector=as.double(LHS_vector),
					K_vector = as.double(K_vector),
					err_tol = as.double(err_tol),
					N = as.integer(N),
					D_vector = double(N) )
		#there are the returned values
		list(dc_mcVals=retvals$D_vector)         #values for D
	}

	if (!is.loaded('mc_uncertainty')){
		dyn.load("mc_uncertainty.so")
	}

	# read data from fcprops.txt file
	data_import <- read.table('fcprops.txt',skip=4,nrows=11, sep="!")
	data_numeric <- as.numeric(as.character(data_import$V1[seq(1,10)]))

	p 			<- data_numeric[1]    # chamber pressure [atm]
	dc_sq 		<- data_numeric[2]    # drop diameter @ onset of fc squared [mm^2]
	K 			<- data_numeric[3]    # burning rate constant prior to onset of fc [mm^2/s]
	do_measured <- data_numeric[4]    # initial drop diameter measured [mm]
	yo 			<- data_numeric[5]    # initial mass frac of low volitiliy component
	err_tol 	<- data_numeric[6]    # error tolerrance for bisection method
	Uk_sq       <- data_numeric[7]^2  # uncertainty in K squared (U_k ^2) [mm^2/s]^2
	Udo_sq      <- data_numeric[8]^2  # uncertainty in do squared (U_do^2) [mm^2]
	Udc_sq      <- data_numeric[9]^2  # uncertainty in dc squared (U_dc^2) [mm^2]
	UYo_sq      <- (data_numeric[10]*yo)^2  #uncertainty in Yo squared
	sol_id      <- data_import$V1[11] #solvent id (1) - Heptane (2)- Propanol


	y_ofc <- 1.0  # (defined) mass fraction of low volatility comp at onset of fc

	if (sol_id == 1){

		if ( yo == 0.05){
			# heptane95-hexadecane5 d_o corrections
			if (p == 1) {
			    percent_increase <- 0.055  #so far these numbers are only valid for hep-hex exp's
			} else {
			    percent_increase <- 0.086
			}	
		}
		if ( yo == 0.20){
			# heptane80-hexadecane20 d_o corrections
			if (p == 1) {
			    percent_increase <- 0.047  
			} else {
			    percent_increase <- 0.076
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
	factor <- 2  #this just specifies factor to increase number of random numbers genrated

	set.seed(5)
	do_temp <- rnorm( n=N*factor, mean=d_o, sd= sqrt(Udo_sq) )

	set.seed(5)
	Yo_mcN95 <- rnorm( n=N, mean=yo, sd=sqrt(UYo_sq) )

	LHS_mcN95 <- y_ofc / Yo_mcN95
	## sample larger than needed values for dc
	## then throw out-of-range values away, element-wise
	dc_mcN95 <- 1
	i <- 5
	while(length(dc_mcN95) < N){
		set.seed(i)
		dc_mcN95 <- rnorm(n=N*factor, mean = sqrt(dc_sq), sd=sqrt(Udc_sq) )
		index <- which(do_temp > dc_mcN95)
		dc_mcN95 <- dc_mcN95[index]
		i <- i + 1
	}
	dc_mcN95seed <- i-1   #store seed for future debug/reproducibility purposes

	dc_mcN95 <- dc_mcN95[1:N]   #select only N values from the within range values
	do_mcN95 <- do_temp[index]
	do_mcN95 <- do_mcN95[1:N]
	tau_mcN95 <- log(do_mcN95 / dc_mcN95) # tau > 0
	## sample larger than needed values for K
	## then throw out-of-range values away
	set.seed(5)
	K_temp <- rnorm(n=N*factor, mean= K, sd= sqrt(Uk_sq))
	K_mcN95 <- 1
	i <- 5
	while( length(K_mcN95) < N ){
		set.seed(i)	
		K_mcN95 <- rnorm(n=N*factor, mean= K, sd= sqrt(Uk_sq))
		K_mcN95 <- K_mcN95[which( K_mcN95 > 0)]
		i <- i + 1	
	}
	K_mcN95seed <- i-1
	K_mcN95 <- K_mcN95[1:N]

	print("----- END random number generation----- ")

	#---------------- calculate normal 95% uncertainties ------------------------ ###
	D_N95 <- calcRoot(tau_vector=tau_mcN95, 
						LHS_vector=LHS_mcN95, 
						K_vector=K_mcN95,
						err_tol =err_tol, 
						N = N)$dc_mcVals
	D_N95 <- D_N95* (1/1000)^2  #convert to m^2/s
	DN95_lower <- quantile(D_N95,probs=c((1-P)/2,(1+P)/2))[[1]]
	DN95_upper <- quantile(D_N95,probs=c((1-P)/2,(1+P)/2))[[2]] 
	DN95_bar <- mean(D_N95 )
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

	# pdf('D_distribution.pdf',width=8, height=5)
	hist(D_N95,prob=TRUE,n=80,  
		main=paste0("Distribution of D"),
		xlab="D [m^2/s]", xlim=c(-3.2E-10, 4.0E-09), col="lightgreen")
	d <-density(D_N95)
	lines(d,col="black",lwd=2)
	abline(v=DN95_upper,col='red',lwd=2.3,lty="dashed")
	abline(v=DN95_lower,col='red',lwd=2.3,lty="dotted")
	abline(v=DN95_bar,col='red',lwd=2.9)
	legend("topright", c("MC"), col=c("red"), lwd=2)
	# dev.off()

	#print out values to screen in a format that's easy to cut and paste into excel
	print(" ")
	print(paste(DN95_lower," ",DN95_bar," ",DN95_upper," ",DN95_most_probable) )
	print(" ")

	list(D_N95 = D_N95)

}


