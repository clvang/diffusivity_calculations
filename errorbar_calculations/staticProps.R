# script to define static properties of components at 1 atm and 298 K
# These information are required in the calculation of viscous decay times

staticProps <- function( sol_id ){

	# define static props for Heptane/Hexadecane
	# CHERIC, YAW imperical relations
	if (sol_id == 1){
		# Define solute (hexadecane) static viscosity and density
		mu_Ao <- 3.1025e-3 #(Vargaftik) 2.891e-3  			#static viscosity of Hexadecane (N-s/m^2) at 1 atm, 298 K
		rho_Ao <- 7.702e2  			#density of Hexadecane (kg/m^3) at 1 atm, 298 K

		# Define solvent (Heptane) static viscosity and density
		mu_Bo <- 3.883e-4 #(Vargaftik) 3.871e-4  			#static viscosity of Heptane at 1 atm, 298K (N-s/m^2)
		rho_Bo <- 5.2935e2 #(Vargaftik) 6.821e2  			#density of Heptane at 1 atm, 298 K (kg/m^3)

		#define molecular weights of each species
		MW_A <- 2.264460E+02 	# molecualr weight of solute, Hexadecane (kg/kmol)
		MW_B <- 1.002040E+02  	# molecular weight of solvent, Heptane (kg/kmol)
	}

	# define static props for Propanol/Glycerol
	if (sol_id == 2){
		# Define solute (Glycerol) static viscosity and density
		# Vargaftik tabulated values
		mu_Ao <- 1.040 #1053.2e-3  # static viscosity of Glycerol (N-s/m^2) at 1 atm, 298 K
		rho_Ao <- 1247.5 #1257.575 	# density of Glycerol (kg/m^3) at 1 atm, 298 K  (I THINK I MESSED UP ENTERING THIS VALUE HERE BEFORE)

		# Define solvent (Propanol) static viscosity and density
		# Vargaftik tabulated values
		mu_Bo <- 1.96e-3 #1.9672e-3  # static viscosity of Propanol at 1 atm, 298K (kg/m^3)
		rho_Bo <- 0.79995e3 #0.8e3 	# density of Propanol at 1 atm, 298 K (kg/m^3)

		#define molecular weights of each species
		MW_A <- 92.1  	# molecualr weight of solute, Glycerol (kg/kmol)
		MW_B <- 60.1  	# molecular weight of solvent, Propanol (kg/kmol)	
	}

	list(MW_A=MW_A, MW_B=MW_B, mu_Ao=mu_Ao, mu_Bo=mu_Bo, rho_Ao=rho_Ao, rho_Bo=rho_Bo)
}

