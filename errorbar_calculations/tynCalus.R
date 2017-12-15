# function to calculate infinite dilution molecular diffusivities
# returns D_AB and D_BA

tynCalus <- function(N, P, expname)
{
	source('/Users/changvang/mygitFiles/diffusivity_calculations/errorbar_calculations/tynCalusSolve_MC.R')
	source('/Users/changvang/mygitFiles/diffusivity_calculations/errorbar_calculations/read_props.R')
	source('/Users/changvang/mygitFiles/diffusivity_calculations/errorbar_calculations/hephex_1atmProps.R')	
	source('/Users/changvang/mygitFiles/diffusivity_calculations/errorbar_calculations/propgly_1atmProps.R')
	source('/Users/changvang/mygitFiles/diffusivity_calculations/errorbar_calculations/propgly_3atmProps.R')	

	# # read data from fcprops.txt file
	experimental_parameters <- read_props()
	p <- experimental_parameters$p 
	sol_id <- experimental_parameters$sol_id  #solvent id (1) - Heptane (2) - Propanol

	if ( sol_id == 1 ){     # Heptane/Hexadecane experiments

		# all Heptane/Hexadecane are at 1 atm only
		if ( p == 1){  

			# read in thermal phyiscal props for heptane-hexadecane
			hephexprops_1atm <- hephex_1atmProps()  
			T <- hephexprops_1atm$T

			rhoB_nbp <- hephexprops_1atm$rhoB_nbp
			rhoB_T <- hephexprops_1atm$rhoB_T 
			sigmaB_T <- hephexprops_1atm$sigmaB_T
			etaB_T <- hephexprops_1atm$etaB_T 
			MW_B <- hephexprops_1atm$MW_B

			rhoA_nbp <- hephexprops_1atm$rhoA_nbp
			rhoA_T <- hephexprops_1atm$rhoA_T
			sigmaA_T <- hephexprops_1atm$sigmaA_T
			etaA_T <- hephexprops_1atm$etaA_T
			MW_A <- hephexprops_1atm$MW_A

			# --------------------------------------------------------------------#
		}

	}else{  			     # Propanol/Glycerol experiments

		if (p == 1){
		#-------------------- Propanol/Glycerol 1 atm properties --------------------#

			# read in thermal phyiscal props for propanol-glycerol @ 1atm
			propglyprops_1atm <- propgly_1atmProps()  
			T <- propglyprops_1atm$T

			rhoB_nbp <- propglyprops_1atm$rhoB_nbp
			rhoB_T <- propglyprops_1atm$rhoB_T 
			sigmaB_T <- propglyprops_1atm$sigmaB_T
			etaB_T <- propglyprops_1atm$etaB_T 
			MW_B <- propglyprops_1atm$MW_B

			rhoA_nbp <- propglyprops_1atm$rhoA_nbp
			rhoA_T <- propglyprops_1atm$rhoA_T
			sigmaA_T <- propglyprops_1atm$sigmaA_T
			etaA_T <- propglyprops_1atm$etaA_T
			MW_A <- propglyprops_1atm$MW_A

		#------------------------------------------------------------------------#
		}
		if (p == 3){
		#-------------------- Propanol/Glycerol 3 atm properties --------------------#
			# read in thermal phyiscal props for propanol-glycerol @ 3atm
			propglyprops_3atm <- propgly_3atmProps()  
			T <- propglyprops_3atm$T

			rhoB_nbp <- propglyprops_3atm$rhoB_nbp
			rhoB_T <- propglyprops_3atm$rhoB_T 
			sigmaB_T <- propglyprops_3atm$sigmaB_T
			etaB_T <- propglyprops_3atm$etaB_T 
			MW_B <- propglyprops_3atm$MW_B

			rhoA_nbp <- propglyprops_3atm$rhoA_nbp
			rhoA_T <- propglyprops_3atm$rhoA_T
			sigmaA_T <- propglyprops_3atm$sigmaA_T
			etaA_T <- propglyprops_3atm$etaA_T
			MW_A <- propglyprops_3atm$MW_A
		#------------------------------------------------------------------------#
		}
	} # if sold_id

	uT <- 15.0    			      # uncertainty T, p/m 15 [K]  
	u_sigmaA_T <- 0.15*sigmaA_T   # uncertainty sigmaA at T, p/m [dynes/cm] changed to 15% since at boiling point
	u_sigmaB_T <- 0.15*sigmaB_T   # uncertainty sigmaB at T, p/m [dynes/cm] changed to 15% since at boiling point
	u_mwA <- 0.01*MW_A		      # uncertainty MW_A, p/m [kg/kmol]
	u_mwB <- 0.01*MW_B            # uncertainty MW_B, p/m [kg/kmol]
	u_rhoA_T <- 0.01*rhoA_T       # uncertainty rhoA at T, p/m [kg/m^3]
	u_rhoB_T <- 0.01*rhoB_T       # uncertainty rhoB at T, p/m [kg/m^3]
	u_rhoA_nbp <- 0.01*rhoA_nbp   # uncertainty rhoA at NBP of component A, p/m [kg/m^3]
	u_rhoB_nbp <- 0.01*rhoB_nbp   # uncertainty rhoB at NPT of component B, p/m [kg/m^3]

	u_etaB_T <- 0.15*etaB_T       # uncertainty etaB, p/m [cP] changed to 15% since at boiling point
	u_etaA_T <- 0.15*etaA_T       # uncertainty etaA, p/m [cP] changed to 15% since at boiling point


	#Calculate D_AB using Monte Carlo
	pltTitle <- "DAB"
	print("********* Molecular Diffusivities Monte Carlo Results ****************")	
	print("D_AB uncertainties.........................")
	D <- tynCalusSolve_MC(N = N, 
					P = P, 					# P-value
					T = T, 					# K
					uT = uT,				# p/m K
					sigmaA_T = sigmaA_T,	# dynes/cm or erg/cm^2
					MW_A = MW_A,			# kg/kmol or g/gmol
					rhoA_T = rhoA_T,		# kg/m^3
					rhoA_nbp = rhoA_nbp,    # kg/m^3
					u_sigmaA_T = u_sigmaA_T,	# dynes/cm or erg/cm^2
					u_mwA = u_mwA,			# kg/kmol or g/gmol
					u_rhoA_T = u_rhoA_T,		# kg/m^3
					u_rhoA_nbp = u_rhoA_nbp,	# kg/m^3
					sigmaB_T = sigmaB_T,		# dynes/cm or erg/cm^2
					MW_B = MW_B,				# kg/kmol or g/gmol
					rhoB_T = rhoB_T,			# kg/m^3
					rhoB_nbp = rhoB_nbp,		# kg/m^3
					etaB_T = etaB_T,			# cP
					u_sigmaB_T = u_sigmaB_T,	# dynes/cm or erg/cm^2
					u_mwB = u_mwB,				# kg/kmol or g/gmol
					u_rhoB_T = u_rhoB_T,		# kg/m^3
					u_rhoB_nbp = u_rhoB_nbp,	# kg/m^3
					u_etaB_T = u_etaB_T,		# cP
					pltTitle = pltTitle,		
					expname = expname )		
	DAB <- D$Dmolecular

	print(" ")
	print("D_BA uncertainties.........................")
	#Calculate D_BA using Monte Carlo
	pltTitle <- "DBA"
	D <- tynCalusSolve_MC(N = N, 
					P = P, 					# P-value
					T = T, 					# K
					uT = uT,				# p/m K
					sigmaA_T = sigmaB_T,	# dynes/cm or erg/cm^2
					MW_A = MW_B,			# kg/kmol or g/gmol
					rhoA_T = rhoB_T,		# kg/m^3
					rhoA_nbp = rhoB_nbp,    # kg/m^3
					u_sigmaA_T = u_sigmaB_T,	# dynes/cm or erg/cm^2
					u_mwA = u_mwB,			# kg/kmol or g/gmol
					u_rhoA_T = u_rhoB_T,		# kg/m^3
					u_rhoA_nbp = u_rhoB_nbp,	# kg/m^3
					sigmaB_T = sigmaA_T,		# dynes/cm or erg/cm^2
					MW_B = MW_A,				# kg/kmol or g/gmol
					rhoB_T = rhoA_T,			# kg/m^3
					rhoB_nbp = rhoA_nbp,		# kg/m^3
					etaB_T = etaA_T,			# cP
					u_sigmaB_T = u_sigmaA_T,	# dynes/cm or erg/cm^2
					u_mwB = u_mwA,				# kg/kmol or g/gmol
					u_rhoB_T = u_rhoA_T,		# kg/m^3
					u_rhoB_nbp = u_rhoA_nbp,	# kg/m^3
					u_etaB_T = u_etaA_T,		# cP
					pltTitle = pltTitle,		
					expname = expname )	

	DBA <- D$Dmolecular

	list(DAB=DAB, DBA=DBA)

}
