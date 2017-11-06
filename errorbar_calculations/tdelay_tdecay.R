
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
	mu_relative_error <- 0.15   			#viscosity 15% relative error for 95% confidence 
	rho_relative_error <- 0.15 				#density 15% relative error for 95% confidnece
	u_muAo <- mu_Ao*mu_relative_error 		#relative error U_mu10/mu_10 = 15% for 95% confidence
	u_rhoAo <- rho_Ao*rho_relative_error    #relative error U_rho/rho = 15% for 95% confidence
	u_muBo <- mu_Bo*mu_relative_error 
	u_rhoBo <- rho_Bo*rho_relative_error

	#uncertainties in molecular weights
	u_mwA <- 0.01*MW_A		# uncertainty MW_A, p/m [kg/kmol]
	u_mwB <- 0.01*MW_B      # uncertainty MW_B, p/m [kg/kmol]

	#generate random variables and compute viscosity of liquid mixture
	set.seed(5)
	MW_solute <- rnorm(n=N,mean=MW_A,sd=u_mwA)  #hexadecane
	MW_solvent <- rnorm(n=N,mean=MW_B,sd=u_mwB)  #heptane

	X_A <- (Yo_mcN95/MW_solute)/( (Yo_mcN95/MW_solute) + (1- Yo_mcN95)/MW_solvent )
	mu_A <- rnorm(n=N,mean=mu_Ao,sd=u_muAo)
	mu_B <- rnorm(n=N,mean=mu_Bo,sd=u_muBo)

	mu_mix <- exp( X_A*log(mu_A) + (1-X_A)*log(mu_B) + 
				X_A*(1-X_A)*1.08*(1.343-X_A*0.685) )    #mixture static viscosity N-s/m^2 at 1 atm, 298 K

	rhoA_mcN95 <- rnorm(n=N,mean=rho_Ao,sd=u_rhoAo)
	rhoB_mcN95 <- rnorm(n=N,mean=rho_Bo,sd=u_rhoBo)
	rho_mix <- 1/( (Yo_mcN95/rhoA_mcN95) + (1-Yo_mcN95)/rhoB_mcN95 )  #mixture density kg/m^3 at 1 atm, 298 K

	nu_mix <- (mu_mix / rho_mix ) * 1e6  #kinematic visocity of liquid mixture (mm^2/s)

	# generate rand numbers for Uo 
	Uo_mcN95 <- rnorm(n=N,mean=Uo,sd=sqrt(UUo_sq) )

	#calculate viscous decay time
	beta <- 0.1	
	tv <- (do_mcN95^2/(4*nu_mix) )*log(Uo_mcN95*do_mcN95/(2*K_mcN95*beta) )

	#generate random variables for td and calculate ratio td/tv
	td_mcN95 <- rnorm(n=N, mean= td, sd= sqrt(Utd_sq))
	#if it turns out that negative values are generated for td_mcN95,
	#then generate random variables from triangular distribution instead
	if ( min(td_mcN95) < 0 ){
		numFrame <- 2
		L_a <- td - (numFrame/30)
		L_b <- td + (numFrame/30)

		if ( L_a < 0){
			td_mcN95 <- rtriangle(N, a=1e-3, b=L_b, c=td )
		}else{
			td_mcN95 <- rtriangle(N, a=L_a, b=L_b, c=td )
		}

	}


	td_tv <- td_mcN95 / tv 

	tdtv_lower <- quantile(td_tv,probs=c((1-P)/2,(1+P)/2))[[1]]
	tdtv_upper <- quantile(td_tv,probs=c((1-P)/2,(1+P)/2))[[2]] 
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

	dev.new()
	pdf(paste0(expname,"_tdtvdist.pdf") )	
	hist(td_tv,prob=TRUE,n=100,  
		main=paste0(expname,": Distribution of td/tv"),
		xlab=expression("t"[d]*"/t"[v]), col="lightgreen")
	d <-density(td_tv)
	lines(d,col="black",lwd=2)
	abline(v=tdtv_upper,col='red',lwd=2.3,lty="dashed")
	abline(v=tdtv_lower,col='red',lwd=2.3,lty="dotted")
	abline(v=tdtv_bar,col='red',lwd=2.9)
	legend("topright", c("MC"), col=c("red"), lwd=2)	

	list( tdtv_lower = tdtv_lower,
		tdtv_bar = tdtv_bar,
		tdtv_upper = tdtv_upper,
		td_mcN95=td_mcN95, 
		tv=tv)
}