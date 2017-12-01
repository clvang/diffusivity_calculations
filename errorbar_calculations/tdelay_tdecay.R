
#------------- calculating error bars for viscous decay time ----------#

tdelay_tdecay <- function(do_mcN95, K_mcN95, Yo_mcN95, expname)
{
	library(triangle)
	source('/Users/changlvang/mygitFiles/diffusivity_calculations/errorbar_calculations/read_props.R')
	source('/Users/changlvang/mygitFiles/diffusivity_calculations/errorbar_calculations/staticProps.R')

	#read data from fcfprops.txt file
	experimental_parameters <- read_props()
	Uo <- experimental_parameters$Uo 
	UUo_sq <- experimental_parameters$UUo_sq
	td <- experimental_parameters$td 
	Utd_sq <- experimental_parameters$Utd_sq
	sol_id <- experimental_parameters$sol_id

	# call staticProps function to read in the appropriate
	# static properties for the component densities,
	# viscosities, and molecualr weights
	static_properties <- staticProps(sol_id)
	mu_Ao <- static_properties$mu_Ao
	rho_Ao <- static_properties$rho_Ao
	MW_A <- static_properties$MW_A
	mu_Bo <- static_properties$mu_Bo
	rho_Bo <- static_properties$rho_Bo
	MW_B <- static_properties$MW_B

	# errors for static viscosities and densitities
	mu_relative_error <- 0.15   			#viscosity 10% relative error for 95% confidence 
	rho_relative_error <- 0.01 				#density 10% relative error for 95% confidnece
	u_muAo <- mu_Ao*mu_relative_error 		#relative error U_mu10/mu_10 = 15% for 95% confidence
	u_rhoAo <- rho_Ao*rho_relative_error    #relative error U_rho/rho = 15% for 95% confidence
	u_muBo <- mu_Bo*mu_relative_error 
	u_rhoBo <- rho_Bo*rho_relative_error

	#uncertainties in molecular weights
	u_mwA <- 0.01*MW_A		# uncertainty MW_A, p/m [kg/kmol]
	u_mwB <- 0.01*MW_B      # uncertainty MW_B, p/m [kg/kmol]

	#convert 95% uncertainties to standard deviations form
	#generating values from a normal distribution
	cf <- 2.0
	sigma_mwA <- u_mwA / cf 
	sigma_mwB <- u_mwB / cf 
	sigma_muAo <- u_muAo / cf 
	sigma_muBo <- u_muBo / cf 
	sigma_rhoAo <- u_rhoAo / cf 
	sigma_rhoBo <- u_rhoBo / cf 
	sigma_Uo <- sqrt(UUo_sq) / cf 

	#generate random variables and compute viscosity of liquid mixture
	set.seed(5)
	MW_solute <- rnorm(n=N,mean=MW_A,sd=sigma_mwA)  #hexadecane
	MW_solvent <- rnorm(n=N,mean=MW_B,sd=sigma_mwB)  #heptane

	X_A <- (Yo_mcN95/MW_solute)/( (Yo_mcN95/MW_solute) + (1- Yo_mcN95)/MW_solvent )
	mu_A <- rnorm(n=N,mean=mu_Ao,sd=sigma_muAo)
	mu_B <- rnorm(n=N,mean=mu_Bo,sd=sigma_muBo)

	mu_mix <- exp( X_A*log(mu_A) + (1-X_A)*log(mu_B) + 
				X_A*(1-X_A)*1.08*(1.343-X_A*0.685) )    #mixture static viscosity N-s/m^2 at 1 atm, 298 K

	rhoA_mcN95 <- rnorm(n=N,mean=rho_Ao,sd=sigma_rhoAo)
	rhoB_mcN95 <- rnorm(n=N,mean=rho_Bo,sd=sigma_rhoBo)
	rho_mix <- 1/( (Yo_mcN95/rhoA_mcN95) + (1-Yo_mcN95)/rhoB_mcN95 )  #mixture density kg/m^3 at 1 atm, 298 K

	nu_mix <- (mu_mix / rho_mix ) * 1e6  #kinematic visocity of liquid mixture (mm^2/s)
	# print( paste0("avg kinematic vis.: ",mean(nu_mix)," mm^2/s") )
	# readline("--press enter--")

	# generate rand numbers for Uo 
	Uo_mcN95 <- rnorm( n=N, mean=Uo, sd=sigma_Uo )

	#calculate viscous decay time
	beta <- 0.1	
	tv <- (do_mcN95^2/(4*nu_mix) )*log(Uo_mcN95*do_mcN95/(2*K_mcN95*beta) )

	# generate random numbers for time delay
	numFrame <- 2  # 95% uncertainty for measurement of t_delay
	L_a <- td - (numFrame/30)
	L_b <- td + (numFrame/30)

	if ( L_a < 0){
		td_mcN95 <- rtriangle(N, a=1e-3, b=L_b, c=td )
	}else{
		td_mcN95 <- rtriangle(N, a=L_a, b=L_b, c=td )
	}


	td_tv <- td_mcN95 / tv 

	tdtv_lower <- quantile(td_tv,probs=c((1-P)/2,(1+P)/2), type=1)[[1]]
	tdtv_upper <- quantile(td_tv,probs=c((1-P)/2,(1+P)/2), type=1)[[2]] 
	tdtv_bar <- mean(td_tv )
	d <- density(td_tv)
	tdtv_most_probable <- d$x[which(d$y==max(d$y))]

	print("******** Viscous Decay Time Monte Carlo Results ***************")
	print("--- results assume independent variables are from NORMAL distribution ---")
	print("MC 95% uncertainty LOWER td/tv limit [s]: ")
	print(tdtv_lower)
	print("MC 95% uncertainty MEAN value td/tv [s]: ")
	print(tdtv_bar)
	print("MC 95% uncertainty UPPER td/tv limit [s]: ")
	print(tdtv_upper)
	print("MC 95% uncertainty MOST PROBABLE td/tv [s]: ")
	print(tdtv_most_probable)

	# plot distirubtion of td/tv ratio
	dev.new()
	pdf(paste0(expname,"_tdtvdist.pdf") )	
	hist(td_tv,prob=TRUE,n=100,  
		main=paste0(expname,": Distribution of td/tv"),
		xlab=expression("t"[d]*"/t"[v]), col="lightgreen")
	# d <-density(td_tv)
	lines(d,col="black",lwd=2)
	abline(v=tdtv_upper,col='red',lwd=2.3,lty="dashed")
	abline(v=tdtv_lower,col='red',lwd=2.3,lty="dotted")
	abline(v=tdtv_bar,col='red',lwd=2.9)
	legend("topright", c("MC"), col=c("red"), lwd=2)	

	# also plot distributions of t_delay for
	# troubleshooting purposes later
	d <- density(td_mcN95)
	td_most_probable <- d$x[which(d$y==max(d$y))]
	td_lower <- quantile(td_mcN95,probs=c((1-P)/2,(1+P)/2), type=1)[[1]]
	td_upper <- quantile(td_mcN95,probs=c((1-P)/2,(1+P)/2), type=1)[[2]] 
	td_bar <- mean(td_mcN95 )

	dev.new()
	pdf(paste0(expname,"_tdelaydist.pdf") )
	hist(td_mcN95,prob=TRUE,n=100,  
		main=paste0(expname,": Distribution of t_delay"),
		xlab=expression("t"[delay]), col="lightgreen")
	# d <-density(td_mcN95)
	lines(d,col="black",lwd=2)
	abline(v=td_upper,col='red',lwd=2.3,lty="dashed")
	abline(v=td_lower,col='red',lwd=2.3,lty="dotted")
	abline(v=td_bar,col='red',lwd=2.9)
	abline(v=td_most_probable,col='black',lwd=2.9)
	legend("topright", c("MC"), col=c("red"), lwd=2)

	# plot distribution of t_decay for troubleshooting
	# purposes later
	d <- density(tv)
	tv_most_probable <- d$x[which(d$y==max(d$y))]
	tv_lower <- quantile(tv,probs=c((1-P)/2,(1+P)/2), type=1)[[1]]
	tv_upper <- quantile(tv,probs=c((1-P)/2,(1+P)/2), type=1)[[2]] 
	tv_bar <- mean( tv )

	dev.new()
	pdf(paste0(expname,"_tdecaydist.pdf") )	
	hist(tv,prob=TRUE,n=100,  
		main=paste0(expname,": Distribution of t_decay"),
		xlab=expression("t"[decay]), col="lightgreen")
	# d <-density(tv)
	lines(d,col="black",lwd=2)
	abline(v=tv_upper,col='red',lwd=2.3,lty="dashed")
	abline(v=tv_lower,col='red',lwd=2.3,lty="dotted")
	abline(v=tv_bar,col='red',lwd=2.9)
	abline(v=tv_most_probable,col='black',lwd=2.9)
	legend("topright", c("MC"), col=c("red"), lwd=2)	

	# output variables 
	list( tdtv_lower = tdtv_lower,
		tdtv_bar = tdtv_bar,
		tdtv_upper = tdtv_upper,
		td_lower =td_lower,		
		td_bar = td_bar, 
		td_upper = td_upper,
		tv_lower = tv_lower,
		tv_bar = tv_bar,
		tv_upper = tv_upper)
}


