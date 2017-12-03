
read_props <- function()
{
		# read data from fcprops.txt file
	data_import <- read.table('fcprops.txt',skip=4,nrows=15, sep="!")
	data_numeric <- as.numeric(as.character(data_import$V1[seq(1,14)]))

	p 			<- data_numeric[1]    # chamber pressure [atm]
	dc_sq 		<- data_numeric[2]    # drop diameter @ onset of fc squared [mm^2]
	K 			<- data_numeric[3]    # burning rate constant prior to onset of fc [mm^2/s]
	do_measured <- data_numeric[4]    # initial drop diameter measured [mm]
	yo 			<- data_numeric[5]    # initial mass frac of low volitiliy component
	err_tol 	<- data_numeric[6]    # error tolerrance for bisection method
	Uk_sq       <- data_numeric[7]^2  # uncertainty in K squared (U_k ^2) [mm^2/s]^2
	Udo_sq      <- data_numeric[8]^2  # uncertainty in do squared (U_do^2) [mm^2]
	Udc_sq      <- data_numeric[9]^2  # uncertainty in dc squared (U_dc^2) [mm^2]
	UYo_sq      <- (data_numeric[10]*yo)^2  #uncertainty in Yo squared
	Uo          <- data_numeric[11]
	UUo_sq 		<- data_numeric[12]^2 # uncertainty in Uo (the characteristic liquid velocity) squared [mm/s]
	td          <- data_numeric[13]   # time delay [s]
	Utd_sq      <- data_numeric[14]^2   # uncertainty in time delay [s]
	sol_id      <- data_import$V1[15] #solvent id (1) - Heptane (2)- Propanol

	list(p=p, dc_sq=dc_sq,
		K=K, do_measured=do_measured,
		yo=yo, err_tol=err_tol,
		Uk_sq=Uk_sq, Udo_sq=Udo_sq,
		Udc_sq=Udc_sq, UYo_sq=UYo_sq,
		Uo=Uo, UUo_sq=UUo_sq,
		td=td, Utd_sq=Utd_sq, sol_id=sol_id)
}
