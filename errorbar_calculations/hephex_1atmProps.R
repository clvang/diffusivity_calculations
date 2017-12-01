# function to defines the thermo-physical properties for heptane and hexadecane
# these include, densities, surface tensions, and molecular weight, required
# as inputs to the Tyn-Calus equation

hephex_1atmProps <- function(){

	T <- 3.716167E+02 		# boiling point of solvent @ 1 atm [K] (CHERIC)	(98.47 deg. C)	

	# ------------- Heptane/Hexadecane 1 atm properties -------------------#
	# NOTE: These values are from imperical relations provided
	# by CHERIC, YAW, etc.  See the contents of the folder
	# "hephex_experiments" for the origins of these numbers
	# B ==============> Solvent (Heptane) - 1 atm =======================>

	# linear interp to find DENSITY at the normal boiling point of the solvent 
	# This value willbe used to compute molar volume of solvent, V_{A}
	Tnb_hep <- 3.716167E+02   #normal boiling point of propanol, [C] --Vargaftik
	T_1 <- 90 				   #Temperature 1, [C] 			--Vargaftik Table Value
	rho_1 <- 0.6218e3		   #density at T_1, [kg/m^3]  	--Vargaftik Table Value
	T_2 <- 100 				   #Temperature 2, [C]			--Vargaftik Table Value
	rho_2 <- 0.6124e3 		   #density at T_2, [kg/m^3]		--Vargaftik Table Value
	rhoB_nbp <- rho_1 + ( (Tnb_hep-273.15)-T_1 )*(rho_2-rho_1)/(T_2-T_1)  # kg/m^3

	# define the DENSITY of the solvent at the temperture of the droplet surface
	# this value will be used to compute parachor of the Solvent, P_{A}
	rhoB_T <- rhoB_nbp 		# kg/m^3

	# find surface tension of solvent at the temperature of the drplet surface
	# (assumed to be equal to the normal boiling point of the solvent)
	# this quantity will be used to calculate Parachor of the solvent, P_{A}
	A_sigma <- 53.649
	B_sigma <- 540.20
	n_sigma <- 1.24310
	# Equation from Yaw, p. 725  182.55K <= T <= 540.20 K
	sigmaB_T <- A_sigma*(1-T/B_sigma)^n_sigma  #[dynes/cm] or erg/cm^2

	# linear interp to find liquid VISCOSITY of the solvent at the 
	# droplet surface (assumed to be equal to the normal boiling point of the solvent)
	T_1 <- 90				#Temp. 1, [C] 					--Vargaftik Table Value
	eta_1 <- 2170e-7		#Viscosity at T_1 [N-s/m^2] 	--Vargaftik Table Value
	T_2 <- 100				#Temp. 2, [C] 					--Vargaftik Table Value
	eta_2 <- 1980e-7		#Viscosity at T_2 [N-s/m^2]		--Vargaftik Table Value
	etaB_T <- ( eta_1 + ( (T-273.15)-T_1 )*(eta_2-eta_1)/(T_2-T_1) )*1e3  # [cP]

	# molecular weight of solvent
	MW_B <- 1.002040E+02  	# kg/kmol or g/gmol 


	# A ===============> Solute (Hexadecane) - 1 atm =========================>

	# linear interp to find density of the solute at the normal boilng 
	# point of the solute. This value will be used to compute molar volume
	# of the solute, V_{B}
	Tnb_hex <- 560.2          # normal bp of glycerol, [K]  --Vargaftik
	A_rho <- 0.2471			  # --Yaw Table Value
	B_rho <- 0.26626		  # --Yaw Table Value
	C_rho <- 723.00			  # --Yaw Table Value
	n_rho <- 0.28571		  # --Yaw Table Value
	# eqn. taken from Yaw (Tmin = 291.31, Tmax= 723.00) for hexadecane
	rhoA_nbp <- A_rho*B_rho^(-(1-Tnb_hex/C_rho)^n_rho)  #g/ml
	rhoA_nbp <- rhoA_nbp*(1/1000)*(1000/0.001)          #kg/m^3

	# define density of the solute at the temperature of the droplet surface, T,
	# assumed to be approximate the boiling point of the solvent. This value
	# will be used to calculate parachor, P_{B}
	T_1 <- 90 			    	  #Temperature 1, [C] 			--Vargaftik Table Value
	rho_1 <- 0.7249e3		      #density at T_1, [kg/m^3]  	--Vargaftik Table Value
	T_2 <- 100 		     		  #Temperature 2, [C]			--Vargaftik Table Value
	rho_2 <- 0.7179e3 		      #density at T_2, [kg/m^3]		--Vargaftik Table Value
	rhoA_T <- rho_1 + ( (T-273.15)-T_1 )*(rho_2-rho_1)/(T_2-T_1)  # kg/m^3

	# compute surface tension of solute at the droplet surface temperture, T. Required
	# for computing parachor of solute, P_{B}
	A_sigma <- 56.805
	B_sigma <- 723.00
	n_sigma <- 1.39290
	# eqn from Yaw, p. 812 291.31 K <= T <= 723.00
	sigmaA_T <- A_sigma*(1-T/B_sigma)^n_sigma  #[dynes/cm] or erg/cm^2

	# solute liquid viscosity at the droplet surface temperature, T. Rquired
	# for computing parachor of solut, P_{B}
	T_1 <- 90			    #Temp. 1, [C] 					--Vargaftik Table Value
	eta_1 <- 1.014e-3		#Viscosity at T_1 [N-s/m^2] 	--Vargaftik Table Value
	T_2 <- 100			    #Temp. 2, [C] 					--Vargaftik Table Value
	eta_2 <- 0.892e-3		#Viscosity at T_2 [N-s/m^2]		--Vargaftik Table Value
	etaA_T <- ( eta_1 + ( (T-273.15)-T_1 )*(eta_2-eta_1)/(T_2-T_1) )*1e3  # [cP]

	# specify molecular weight of solute
	MW_A <- 2.264460E+02 	# kg/kmol or g/gmol


	# # B ==> Solvent (Heptane)
	# # Tnb_hep <-  			# normal boiling point of heptane
	# rhoB <- 6.139854E+02  	# kg/m^3 @ solvent NBP
	# MW_B <- 1.002040E+02  	# kg/kmol
	# sigmaB <- 12.6191    	# solvant surface tension @ NBP [erg/cm^2]
	# etaB <- 1.983129E-01   # solvant liquid viscosity at NBP [cP]

	# # A ==> Solute (Hexadecane)
	# rhoA <- 7.252842E+02 	# kg/m^3 @ solvent NBP
	# MW_A <- 2.264460E+02 	# kg/kmol
	# sigmaA <- 21.03   		# solute surface tension @ BP [erg/cm^2]
	# etaA <- 9.338927E-01    # solute liquid viscosity at BP [cP]

	# --------------------------------------------------------------------#

	list(T=T, rhoB_nbp=rhoB_nbp, rhoB_T=rhoB_T, sigmaB_T=sigmaB_T, etaB_T=etaB_T, 
		MW_B=MW_B, rhoA_nbp=rhoA_nbp, rhoA_T=rhoA_T, sigmaA_T=sigmaA_T, etaA_T=etaA_T, MW_A=MW_A)
}
