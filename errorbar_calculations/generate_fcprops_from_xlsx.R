
generate_fcprops_from_xlsx <- function(i, importedData){

	if (importedData$P[i]==1){
		pressure <- "1.0"
	}
	if (importedData$P[i]==3){
		pressure <- "3.0"
	}
	if (importedData$fuel[i] == "Prop95/Gly5"){
		Yo <- "0.05"
		sol_id <- "2"
	}else if( importedData$fuel[i] == "Hep95/Hex5" ){
		Yo <- "0.05"
		sol_id <- "1"
	}else{
		Yo <- "0.20"
		sol_id <- "1"
	}

	fileConn <- file("fcprops.txt")
	writeLines(c("...................................................................................",
		".    Header for input file containing experimental",
		".	 measurements for dc, do, yo, p, etc...",
		"...................................................................................",
		paste(pressure,"                   !chamber pressure [atm]"),
		paste(as.numeric(importedData$d_c[i])^2,"      !dc_sq, drop diameter squared at onset of flame contraction [mm^2]"),
		paste(importedData$K[i],"                 !K, burning rate on onset of flame contraction [mm^2/s]"), 
		paste(importedData$d0[i],"                 !do_measured = initial droplet diameter measured [mm]"),
		paste(Yo,"                  !yo=initial low volatility mass fraction (either Y of hexadecane or glycerol)"),
		paste("1.0E-12","               !err_tol = error tolerance for norms"),
		paste(importedData$Uncertainty_K[i],"    !95% Uncertainty in K (mm^2/s)"),
		paste("0.1","                   !95% Uncertainty in measurement of d_o (mm)"),
		paste("0.1","                   !95% Uncertainty in measurement of d_c (mm)"),
		paste("0.02","                  !95% relative uncertainty in Yo e.g. U_Yo/Yo (percent)"),
		paste("800.0","                 !Characteristic velocity, Uo (mm/s)"),
		paste("150.0","                 !95% Uncertainty in measurement of Uo (mm/s)"),
		paste(importedData$t_delay[i],"               !time delay [s]"),
		paste("0.67","                  !95% uncertainty in measurement of time delay [s] (p/m two frames)"),
		paste(sol_id,"                  !solvent (1)-Heptane (2)-Propanol"),
		paste(importedData$dod1_ratio[i],"				!diameter ratio at needle retraction to igniter on, do/d1") ),
		fileConn)

	close(fileConn)

}




