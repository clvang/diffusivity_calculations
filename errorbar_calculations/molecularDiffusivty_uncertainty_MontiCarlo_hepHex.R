# function to calculate infinite dilution molecular diffusivities
# returns D_AB and D_BA

molecularDiffusivty_uncertainty_MontiCarlo_hepHex <- function(N, P)
{
	source('/Users/changlvang/mygitFiles/diffusivity_calculations/errorbar_calculations/tynCalus_MC.R')
	source('/Users/changlvang/mygitFiles/diffusivity_calculations/errorbar_calculations/read_props.R')

	# # read data from fcprops.txt file
	experimental_parameters <- read_props()
	p <- experimental_parameters$p 
	sol_id <- experimental_parameters$sol_id  #solvent id (1) - Heptane (2) - Propanol

	if ( sol_id == 1 ){     # Heptane/Hexadecane experiments

		# all Heptane/Hexadecane are at 1 atm only
		if ( p == 1){  
			# ------------- Heptane/Hexadecane 1 atm properties -------------------#
			# NOTE: These values are from imperical relations provided
			# by CHERIC, YAW, etc.  See the contents of the folder
			# "hephex_experiments" for the origins of these numbers

			# B ==> Solvent (Heptane)
			rhoB <- 6.139854E+02  	# kg/m^3 @ solvent NBP
			MW_B <- 1.002040E+02  	# kg/kmol
			sigmaB <- 12.6191    	# solvant surface tension @ NBP [erg/cm^2]
			etaB <- 1.983129E-01   # solvant liquid viscosity at NBP [cP]

			# A ==> Solute (Hexadecane)
			rhoA <- 7.252842E+02 	# kg/m^3 @ solvent NBP
			MW_A <- 2.264460E+02 	# kg/kmol
			sigmaA <- 21.03   		# solute surface tension @ BP [erg/cm^2]
			etaA <- 9.338927E-01    # solute liquid viscosity at BP [cP]

			T <- 3.716167E+02 		# boiling point of solvent @ 1 atm [K]
			# --------------------------------------------------------------------#
		}

	}else{  			     # Propanol/Glycerol experiments

		if (p == 1){
			#---------------- Propanol/Glycerol 1 atm properties --------------------#
			# B ==> Solvent (Propanol) - 1 atm
			rhoB <- 733.2       	# kg/m^3 @ BP: 1atm => 97.35C=370.5 K (check)
			MW_B <- 60.1        	# kg/kmol (check)
			etaB <- 0.467     		# solvant liquid viscosity at BP [cP] CHERIC (check)
			sigmaB <- 0    			# place holder only. surface tenion data not available

			# A ==> Solute (Glycerol) - 1 atm
			rhoA <- 1212.7   	    # kg/m^3 @ BP: 1atm => 97.3C=370.5 K   (check)
			MW_A <- 92.1        	# kg/kmol (check)
			etaA <- 14.385          # solute liquid viscosity at BP [cP] Vargaftik
			sigmaA <- 0 			# place holder

			T <- 370.5 				# boiling point of solvent @ 1 atm [K] - 1 atm
			#------------------------------------------------------------------------#
		}
		if (p == 3){
			#---------------- Propanol/Glycerol 3 atm properties --------------------#
			# B ==> Solvent (Propanol) - 3 atm
			rhoB <- 697.4       	# kg/m^3 @ BP: 3atm => 129.85C=403K (check)
			MW_B <- 60.1        	# kg/kmol (check)
			etaB <- 0.313     		# solvant liquid viscosity at BP [cP] CHERIC (check)
			sigmaB <- 0    			# place holder only. surface tenion data not available			

			# A ==> Solute (Glycerol) - 3 atm
			rhoA <- 1191.0  	    # kg/m^3 @ BP: 3atm => 129.85C=403K   (check)
			MW_A <- 92.1        	# kg/kmol (check)
			etaA <- 3.5255          # solute liquid viscosity at BP [cP] Vargaftik
			sigmaA <- 0 			# place holder			

			T <- 403.0 				# boiling point of solvent @ 3 atm [K] - 1 atm
			#------------------------------------------------------------------------#
		}
	} # if sold_id

	uT <- 15.0    			# uncertainty T, p/m 15 [K]  
	u_sigmaA <- 0.15*sigmaA # uncertainty sigmaA, p/m [erg/cm^2]
	u_sigmaB <- 0.15*sigmaB # uncertainty sigmaB, p/m [erg/cm^2]
	u_mwA <- 0.01*MW_A		# uncertainty MW_A, p/m [kg/kmol]
	u_mwB <- 0.01*MW_B      # uncertainty MW_B, p/m [kg/kmol]
	u_rhoA <- 0.10*rhoA     # uncertainty rhoA, p/m [kg/m^3]
	u_rhoB <- 0.10*rhoB     # uncertainty rhoB, p/m [kg/m^3]
	u_etaB <- 0.15*etaB     # uncertainty etaB, p/m [cP]
	u_etaA <- 0.15*etaA     # uncertainty etaA, p/m [cP]


	#Calculate D_AB using Monte Carlo
	pltTitle <- "DAB"
	print("********* Molecular Diffusivities Monte Carlo Results ****************")	
	print("D_AB uncertainties.........................")
	D <- tynCalus_MC(N = N, 
					P = P, 
					T = T, 
					uT = uT,
					sigmaA = sigmaA,
					MW_A = MW_A,
					rhoA = rhoA,
					u_sigmaA = u_sigmaA,
					u_mwA = u_mwA,
					u_rhoA = u_rhoA,
					sigmaB = sigmaB,
					MW_B = MW_B,
					rhoB = rhoB,
					etaB = etaB,
					u_sigmaB = u_sigmaB,
					u_mwB = u_mwB,
					u_rhoB = u_rhoB,
					u_etaB = u_etaB,
					pltTitle = pltTitle )
	DAB <- D$Dmolecular

	print(" ")
	print("D_BA uncertainties.........................")
	#Calculate D_BA using Monte Carlo
	pltTitle <- "DBA"
	D <- tynCalus_MC(N = N,
					 P = P,
					 T = T,
					 uT = uT,
					 sigmaA = sigmaB,
					 MW_A = MW_B,
					 rhoA = rhoB,
					 u_sigmaA = u_sigmaB,
					 u_mwA = u_mwB,
					 u_rhoA = u_rhoB,
					 sigmaB = sigmaA,
					 MW_B = MW_A,
					 rhoB = rhoA,
					 etaB = etaA,
					 u_sigmaB = u_sigmaA,
					 u_mwB = u_mwA,
					 u_rhoB = u_rhoA,
					 u_etaB = u_etaA,
					 pltTitle = pltTitle)
	DBA <- D$Dmolecular

	list(DAB=DAB, DBA=DBA)

}
