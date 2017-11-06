
# Change properties accordingly for the pressure!!!!

rm(list=ls(all=TRUE))   #remove all variables in  workspace


source('/Users/changlvang/mygitFiles/diffusivity_calculations/propgly_experiments/D_uncertainty_MC_prpgly.R')

#See text by Shaw, p. 147 for a similar example
# #--------------------------- 1 atm propertis ----------------------------#
# # B ==> Solvent (Propanol) - 1 atm
# rhoB <- 733.2       	# kg/m^3 @ BP: 1atm => 97.35C=370.5 K (check)
# MW_B <- 60.1        	# kg/kmol (check)
# etaB <- 0.467     		# solvant liquid viscosity at BP [cP] CHERIC (check)

# # A ==> Solute (Glycerol) - 1 atm
# rhoA <- 1212.7   	    # kg/m^3 @ BP: 1atm => 97.3C=370.5 K   (check)
# MW_A <- 92.1        	# kg/kmol (check)
# etaA <- 14.385          # solute liquid viscosity at BP [cP] Vargaftik

# T <- 370.5 				# boiling point of solvent @ 1 atm [K] - 1 atm
# #------------------------------------------------------------------------#

#--------------------------- 3 atm propertis ----------------------------#
# B ==> Solvent (Propanol) - 3 atm
rhoB <- 697.4       	# kg/m^3 @ BP: 3atm => 129.85C=403K (check)
MW_B <- 60.1        	# kg/kmol (check)
etaB <- 0.313     		# solvant liquid viscosity at BP [cP] CHERIC (check)

# A ==> Solute (Glycerol) - 3 atm
rhoA <- 1191.0  	    # kg/m^3 @ BP: 3atm => 129.85C=403K   (check)
MW_A <- 92.1        	# kg/kmol (check)
etaA <- 3.5255          # solute liquid viscosity at BP [cP] Vargaftik

T <- 403.0 				# boiling point of solvent @ 3 atm [K] - 1 atm
#------------------------------------------------------------------------#


uT <- 15.0    			# uncertainty T, p/m [K]
u_mwA <- 0.01*MW_A		# uncertainty MW_A, p/m [kg/kmol]
u_mwB <- 0.01*MW_B      # uncertainty MW_B, p/m [kg/kmol]
u_rhoA <- 0.10*rhoA     # uncertainty rhoA, p/m [kg/m^3]
u_rhoB <- 0.10*rhoB     # uncertainty rhoB, p/m [kg/m^3]
u_etaB <- 0.15*etaB     # uncertainty etaB, p/m [cP]
u_etaA <- 0.15*etaA     # uncertainty etaA, p/m [cP]

#these uncertaintes calculated using TSM -- Fortran code
# DAB_upper_TSM <- 6.364718E-09   
# DAB_lower_TSM <- 3.682655E-09
# DAB_TSM <- 5.023686E-09

# DBA_upper_TSM <- 2.757853E-09
# DBA_lower_TSM <- 1.158306E-09
# DBA_TSM <- 1.958079E-09

#Calculate D_AB using Monte Carlo
N <- 1000000
P <- 0.95
pltTitle <- "DAB"
print("D_AB uncertainties.........................")
D_uncertainty_MC(N,P,T,uT,MW_A,rhoA,u_mwA,u_rhoA,
						MW_B,rhoB,etaB,u_mwB,u_rhoB,u_etaB, pltTitle)
print(" ")
print("D_BA uncertainties.........................")
#Calculate D_BA using Monte Carlo
pltTitle <- "DBA"
D_uncertainty_MC(N,P,T,uT,MW_B,rhoB,u_mwB,u_rhoB,
						MW_A,rhoA,etaA,u_mwA,u_rhoA,u_etaA,pltTitle)






