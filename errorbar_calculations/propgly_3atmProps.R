# function to defines the thermo-physical properties for heptane and hexadecane
# these include, densities, surface tensions, and molecular weight, required
# as inputs to the Tyn-Calus equation

propgly_3atmProps <- function(){

	T <- 403.0 				# boiling point of solvent @ 3 atm [K] (129.85 C)

	# B ==============> Solvent (Propanol) - 3 atm =======================>

	# linear interp to find density at the normal boiling point of the solvent 
	# This value willbe used to compute molar volume of solvent, V_{A}
	Tnb_prop <- 97.8 + 273.15 # normal boiling point of propanol, [C] --Vargaftik
	T_1 <- 90 				  #Temperature 1, [C] 			--Vargaftik Table Value
	rho_1 <- 0.7425e3		  #density at T_1, [kg/m^3]  	--Vargaftik Table Value
	T_2 <- 100 				  #Temperature 2, [C]			--Vargaftik Table Value
	rho_2 <- 0.7325e3 		  #density at T_2, [kg/m^3]		--Vargaftik Table Value
	rhoB_nbp <- rho_1 + ( (Tnb_prop-273.15) -T_1)*(rho_2-rho_1)/(T_2-T_1)  # kg/m^3

	# define the density of the solvent at the temperture of the droplet surface
	# this value will be used to compute parachor of the Solvent, P_{A}
	T_1 <- 120 				  #Temperature 1, [C] 			--Vargaftik Table Value
	rho_1 <- 0.7110e3		  #density at T_1, [kg/m^3]  	--Vargaftik Table Value
	T_2 <- 130 				  #Temperature 2, [C]			--Vargaftik Table Value
	rho_2 <- 0.6995e3 		  #density at T_2, [kg/m^3]		--Vargaftik Table Value
	rhoB_T <- rho_1 + ((T-273.15)-T_1)*(rho_2-rho_1)/(T_2-T_1)  # kg/m^3

	# find surface tension of solvent at the temperature of the drplet surface
	# (assumed to be equal to the normal boiling point of the solvent)
	# this quantity will be used to calculate Parachor of the solvent, P_{A}
	A_sigma <- 44.605
	B_sigma <- 536.78
	n_sigma <- 0.8
	# Equation from Yaw, p. 725  181.83K <= T <= 652.19 K
	sigmaB_T <- A_sigma*(1-T/B_sigma)^n_sigma  #[dynes/cm] 

	# linear interp to find liquid viscosity of the solvent at the 
	# droplet surface (assumed to be equal to the normal boiling point of the solvent)
	T_1 <- 120				#Temp. 1, [C] 					--Vargaftik Table Value
	eta_1 <- 0.337e-3		#Viscosity at T_1 [N-s/m^2] 	--Vargaftik Table Value
	T_2 <- 130				#Temp. 2, [C] 					--Vargaftik Table Value
	eta_2 <- 0.291e-3		#Viscosity at T_2 [N-s/m^2]		--Vargaftik Table Value
	etaB_T <- ( eta_1 + ( (T-273.15)-T_1 )*(eta_2-eta_1)/(T_2-T_1) )*1e3  # [cP]

	# molecular weight of solvent
	MW_B <- 60.1        	  # kg/kmol (check)			


	# A ===============> Solute (Glycerol) - 3 atm =========================>

	# linear interp to find density of the solute at the normal boilng 
	# point of the solute. This value will be used to compute molar volume
	# of the solute, V_{B}
	Tnb_gly <- 290+273.15     # normal bp of glycerol, [K]  --Vargaftik
	A_rho <- 0.3488			  # --Yaw Table Value
	B_rho <- 0.25408		  # --Yaw Table Value
	C_rho <- 850.00			  # --Yaw Table Value
	n_rho <- 0.15410		  # --Yaw Table Value
	# eqn. taken from Yaw (Tmin = 291.33, Tmax= 850.00) for Glycerol
	rhoA_nbp <- A_rho*B_rho^(-(1-Tnb_gly/C_rho)^n_rho)  #g/ml
	rhoA_nbp <- rhoA_nbp*(1/1000)*(1000/0.001)          #kg/m^3

	# define density of the solute at the temperature of the droplet surface, T,
	# assumed to be approximate the boiling point of the solvent. This value
	# will be used to calculate parachor, P_{B}
	T_1 <- 120 				  #Temperature 1, [C] 			--Vargaftik Table Value
	rho_1 <- 1188		      #density at T_1, [kg/m^3]  	--Vargaftik Table Value
	T_2 <- 140 				  #Temperature 2, [C]			--Vargaftik Table Value
	rho_2 <- 1167 		      #density at T_2, [kg/m^3]		--Vargaftik Table Value
	rhoA_T <- rho_1 + ( (T-273.15)-T_1 )*(rho_2-rho_1)/(T_2-T_1)  # kg/m^3

	# compute surface tension of solute at the droplet surface temperture, T. Required
	# for computing parachor of solute, P_{B}
	A_sigma <- 143.408
	B_sigma <- 850.00
	n_sigma <- 1.22222
	# eqn from Yaw, p. 726 291.33 K <= T <= 850.00
	sigmaA_T <- A_sigma*(1-T/B_sigma)^n_sigma  #[dynes/cm] 

	# solute liquid viscosity at the droplet surface temperature, T. Rquired
	# for computing parachor of solut, P_{B}
	T_1 <- 120			#Temp. 1, [C] 					--Vargaftik Table Value
	eta_1 <- 5.2e-3		#Viscosity at T_1 [N-s/m^2] 	--Vargaftik Table Value
	T_2 <- 140			#Temp. 2, [C] 					--Vargaftik Table Value
	eta_2 <- 1.8e-3		#Viscosity at T_2 [N-s/m^2]		--Vargaftik Table Value
	etaA_T <- ( eta_1 + ( (T-273.15)-T_1 )*(eta_2-eta_1)/(T_2-T_1) )*1e3  # [cP]

	# specify molecular weight of solute
	MW_A <- 92.1        	# kg/kmol (check)

	# #---------------- Propanol/Glycerol 3 atm properties --------------------#

	# # B ==> Solvent (Propanol) - 3 atm
	# rhoB <- 697.4       	# kg/m^3 @ BP: 3atm => 129.85C=403K (check)
	# MW_B <- 60.1        	# kg/kmol (check)
	# etaB <- 0.313     		# solvant liquid viscosity at BP [cP] CHERIC (check)
	# sigmaB <- 0    			# place holder only. surface tenion data not available			

	# # A ==> Solute (Glycerol) - 3 atm
	# rhoA <- 1191.0  	    # kg/m^3 @ BP: 3atm => 129.85C=403K   (check)
	# MW_A <- 92.1        	# kg/kmol (check)
	# etaA <- 3.5255          # solute liquid viscosity at BP [cP] Vargaftik
	# sigmaA <- 0 			# place holder			

	# #------------------------------------------------------------------------#

	list(T=T, rhoB_nbp=rhoB_nbp, rhoB_T=rhoB_T, sigmaB_T=sigmaB_T, etaB_T=etaB_T, 
		MW_B=MW_B, rhoA_nbp=rhoA_nbp, rhoA_T=rhoA_T, sigmaA_T=sigmaA_T, etaA_T=etaA_T, MW_A=MW_A)
}
