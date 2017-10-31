
#------------- calculating error bars for viscous decay time ----------#

tdelay_tdecay <- function(do_mcN95, K_mcN95, Yo_mcN95)
{

	source('/Users/changlvang/mygitFiles/diffusivity_calculations/enhancement_factors/read_props.R')

	#read data from fcfprops.txt file
	experimental_parameters <- read_props()
	Uo <- experimental_parameters$Uo 
	UUo_sq <- experimental_parameters$UUo_sq
	td <- experimental_parameters$td 
	Utd_sq <- experimental_parameters$Utd_sq

	# script to calculate viscous delay time along with its error bars

	#first, define component viscosities and densities
	#see the code decayTimes.m for the origin of these numbers
	mu_relative_error <- 0.15   #15% relative error for 95% confidence 
	rho_relative_error <- 0.15
	mu_Ao <- 2.891e-3  			#static viscosity of Hexadecane (N-s/m^2) at 1 atm, 298 K
	rho_Ao <- 7.702e2  			#density of Hexadecane (kg/m^3) at 1 atm, 298 K
	u_muAo <- mu_Ao*mu_relative_error 		#relative error U_mu10/mu_10 = 15% 
											#for 95% confidence
	u_rhoAo <- rho_Ao*rho_relative_error    #relative error U_rho/rho = 15% for 95% confidence

	mu_Bo <- 3.871e-4  #static viscosity of Heptane at 1 atm, 298K (kg/m^3)
	rho_Bo <- 6.821e2  #density of Heptane at 1 atm, 298 K (kg/m^3)
	u_muBo <- mu_Bo*mu_relative_error 
	u_rhoBo <- rho_Bo*rho_relative_error

	#define molecular weights of each species
	MW_A <- 2.264460E+02 	# molecualr weight of solute, Hexadecane (kg/kmol)
	u_mwA <- 0.01*MW_A		# uncertainty MW_A, p/m [kg/kmol]
	MW_B <- 1.002040E+02  	# molecular weight of solvent, Heptane (kg/kmol)
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
	beta <- 0.1

	#calculate viscous delay time
	tv <- (do_mcN95^2/(4*nu_mix) )*log(Uo_mcN95*do_mcN95/(2*K_mcN95*beta) )

	#calculate ratio td/tv
	factor <- 2
	td_mcN95 <- 1
	i <- 5
	while( length(td_mcN95) < N ){
		set.seed(i)	
		td_mcN95 <- rnorm(n=N*factor, mean= td, sd= sqrt(Utd_sq))
		td_mcN95 <- td_mcN95[which( td_mcN95 > 0)]
		i <- i + 1	
	}
	td_mcN95seed <- i-1
	td_mcN95 <- td_mcN95[1:N]

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
	hist(td_tv,prob=TRUE,n=100,  
		main=paste0("Distribution of td/tv"),
		xlab=expression("t"[d]*"/t"[v]), col="lightgreen")
	d <-density(td_tv)
	lines(d,col="black",lwd=2)
	abline(v=tdtv_upper,col='red',lwd=2.3,lty="dashed")
	abline(v=tdtv_lower,col='red',lwd=2.3,lty="dotted")
	abline(v=tdtv_bar,col='red',lwd=2.9)
	legend("topright", c("MC"), col=c("red"), lwd=2)	

	print(mean(td_mcN95))
	print(td)

}