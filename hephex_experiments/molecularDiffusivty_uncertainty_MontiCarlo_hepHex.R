
rm(list=ls(all=TRUE))   #remove all variables in  workspace


source('/Users/changvang/mygitFiles/diffusivity_calculations/hephex_experiments/tynCalus_MC.R')

#See text by Shaw, p. 147 for a similar example

# B ==> Solvent (Heptane)
rhoB <- 6.139854E+02  	# kg/m^3
MW_B <- 1.002040E+02  	# kg/kmol
sigmaB <- 12.6191    	# solvant surface tension @ BP [erg/cm^2]
etaB <- 1.983129E-01   # solvant liquid viscosity at BP [cP]

# A ==> Solute (Hexadecane)
rhoA <- 7.252842E+02 	# kg/m^3
MW_A <- 2.264460E+02 	# kg/kmol
sigmaA <- 21.03   		# solute surface tension @ BP [erg/cm^2]
etaA <- 9.338927E-01    # solute liquid viscosity at BP [cP]

T <- 3.716167E+02 		# boiling point of solvent @ 1 atm [K]

uT <- 15.0    			# uncertainty T, p/m [K]
u_sigmaA <- 0.15*sigmaA # uncertainty sigmaA, p/m [erg/cm^2]
u_sigmaB <- 0.15*sigmaB # uncertainty sigmaB, p/m [erg/cm^2]
u_mwA <- 0.01*MW_A		# uncertainty MW_A, p/m [kg/kmol]
u_mwB <- 0.01*MW_B      # uncertainty MW_B, p/m [kg/kmol]
u_rhoA <- 0.10*rhoA     # uncertainty rhoA, p/m [kg/m^3]
u_rhoB <- 0.10*rhoB     # uncertainty rhoB, p/m [kg/m^3]
u_etaB <- 0.15*etaB     # uncertainty etaB, p/m [cP]
u_etaA <- 0.15*etaA     # uncertainty etaA, p/m [cP]

#these uncertaintes calculated using TSM -- Fortran code
DAB_upper_TSM <- 6.364718E-09   
DAB_lower_TSM <- 3.682655E-09
DAB_TSM <- 5.023686E-09

DBA_upper_TSM <- 2.757853E-09
DBA_lower_TSM <- 1.158306E-09
DBA_TSM <- 1.958079E-09

#Calculate D_AB using Monte Carlo
N <- 1000000
P <- 0.95
pltTitle <- "DAB"
print("D_AB uncertainties.........................")
tynCalus_MC(N,P,T,uT,sigmaA,MW_A,rhoA,u_sigmaA,u_mwA,u_rhoA,
						sigmaB,MW_B,rhoB,etaB,u_sigmaB,u_mwB,u_rhoB,u_etaB,
						DAB_upper_TSM,DAB_lower_TSM,DAB_TSM,pltTitle)
print(" ")
print("D_BA uncertainties.........................")
#Calculate D_BA using Monte Carlo
pltTitle <- "DBA"
tynCalus_MC(N,P,T,uT,sigmaB,MW_B,rhoB,u_sigmaB,u_mwB,u_rhoB,
						sigmaA,MW_A,rhoA,etaA,u_sigmaA,u_mwA,u_rhoA,u_etaA,
						DBA_upper_TSM,DBA_lower_TSM,DBA_TSM,pltTitle)






