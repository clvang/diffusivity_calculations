
! 	Po - pressure in kPa

SUBROUTINE DABuncertainty(MWsolvent,rhoSatSolventNBP, MWsolute,rhoSoluteAtSolventTsatNBP, &
						Po, AantoinneSolvent, BantoinneSolvent, CantoinneSolvent, DantoinneSolvent, &
						TminAntoinneSolvent, TmaxAntoinneSolvent, &
						ArhoSolvent,BrhoSolvent,CrhoSolvent,nrhoSolvent, &
						ArhoSolute,BrhoSolute,CrhoSolute,nrhoSolute, &
						AandradeSolvent, BandradeSolvent,sigmaSolvent, sigmaSolute, DAB, & 
						U_T, U_MW_relative, U_sigma_relative,U_Dab_95, T)
	IMPLICIT NONE
	REAL, INTENT(IN) :: MWsolvent,rhoSatSolventNBP, MWsolute,rhoSoluteAtSolventTsatNBP, &
						Po, AantoinneSolvent, BantoinneSolvent, CantoinneSolvent, DantoinneSolvent, &
						TminAntoinneSolvent, TmaxAntoinneSolvent, &
						ArhoSolvent,BrhoSolvent,CrhoSolvent,nrhoSolvent, &
						ArhoSolute,BrhoSolute,CrhoSolute,nrhoSolute, &
						AandradeSolvent, BandradeSolvent,sigmaSolvent, sigmaSolute, DAB, &
						U_T, U_MW_relative, U_sigma_relative, T
	REAL, INTENT(OUT) :: U_Dab_95
	REAL, DIMENSION(10) :: dP, deltaDP, P_plus, P_minus, T_plus, T_minus, dTdP 
	INTEGER :: i, N 
	REAL :: Va, Vb, U_Pv_sq, U_ToT_sq, U_Arho_sq, U_Brho_sq, U_rhoAOrhoA_sq, &
			U_VaOVa_sq, U_MWA_sq, U_rhoBOrhoB_sq, U_VbOVb_sq, U_MWB_sq, U_AetaB_sq, &
			U_BetaB_sq, U_etaBOetaB_sq, U_sigmaAosigmaA_sq, &
			U_sigmaBosigmaB_sq, U_DabODab_sq, U_Dab, U_Crho_sq, U_nrho_sq

	!compute molar volume [cm^3/mol] at NBP solvent
	Vb = (MWsolvent/rhoSatSolventNBP)*(1E3)	

	!compute molar volume[cm^3/mol] at NBP of solvent
	Va = (MWsolute/rhoSoluteAtSolventTsatNBP)*(1E3)

	!Calculate relative uncertainty in T squared, (U_T/T)^2
	N = 10
	deltaDP(1) = 0.01*Po
	!estimate sensitivity coefficent dT/dP
	DO i=1,N  
		deltaDP(i+1) = 0.01*Po / ( (i+1)*2.0 ) 
		P_plus(i) = Po + deltaDP(i)
		P_minus(i) = Po - deltaDP(i)

		CALL TsatAntoinne(P_plus(i), AantoinneSolvent,&
			BantoinneSolvent,CantoinneSolvent,DantoinneSolvent, &
			TminAntoinneSolvent, TmaxAntoinneSolvent,T_plus(i))

		CALL TsatAntoinne(P_minus(i), AantoinneSolvent,&
			BantoinneSolvent,CantoinneSolvent,DantoinneSolvent, &
			TminAntoinneSolvent, TmaxAntoinneSolvent,T_minus(i))

		dTdP(i) = ( T_plus(i) - T_minus(i) ) / ( 2.0*deltaDP(i)  )
	END DO
	! CALL TsatAntoinne(Po, AantoinneSolvent,&
	! 	BantoinneSolvent,CantoinneSolvent,DantoinneSolvent, &
	! 	TminAntoinneSolvent, TmaxAntoinneSolvent,T)

	U_Pv_sq = (0.1*Po)**2   					! U_{P_v}^2, Assume U_{P_v} = 10 percent  (U_Pv_sq not really required)

	!U_ToT_sq = (dTdP(N) / T)**2 * U_Pv_sq     ! (U_T/T)^2 
	U_ToT_sq = (U_T/T)**2  					    ! (U_T/T)^2 assume U_T = 15 K

	!Calculate relative uncertainty in density of solute (A) squared, (U_rho_A/ rho_A)^2
	U_Arho_sq = 0.01**2 		! (U_A/A)^2
	U_Brho_sq = 0.01**2			! (U_B/B)^2
	U_Crho_sq = 0.01**2 	    ! (U_C/C)^2
	U_nrho_sq = 0.01**2 		! (U_n/n)^2
	U_rhoAOrhoA_sq = ( 1. / ArhoSolute)**2 * U_Arho_sq + &
		( (1. - (T/CrhoSolute))**(nrhoSolute) )**2 * U_Brho_sq + &
		( nrhoSolute*(1- (T/CrhoSolute))**( nrhoSolute - 1) * LOG(BrhoSolute) / CrhoSolute )**2 * (U_ToT_sq * T**2) + &
		( nrhoSolute*T*(1- (T/CrhoSolute))**( nrhoSolute - 1) * LOG(BrhoSolute) / CrhoSolute )**2 * U_Crho_sq + &
		( (1- (T/CrhoSolute))**( nrhoSolute ) * LOG(BrhoSolute) * LOG(1- (T/CrhoSolute))  )**2 * U_nrho_sq

	!Calculate relative uncertainty in density of solvant (B) squared: (U_rho_B/ rho_B)^2
	U_rhoBOrhoB_sq = ( 1. / ArhoSolvent)**2 * U_Arho_sq + &
		( (1. - (T/CrhoSolvent))**(nrhoSolvent) / BrhoSolvent )**2 * U_Brho_sq + &
		( nrhoSolvent*(1- (T/CrhoSolvent))**( nrhoSolvent - 1) * LOG(BrhoSolvent) / CrhoSolvent )**2 * (U_ToT_sq * T**2) + &
		( nrhoSolvent*T*(1- (T/CrhoSolvent))**( nrhoSolvent - 1) * LOG(BrhoSolvent) / CrhoSolvent )**2 * U_Crho_sq + &
		( (1- (T/CrhoSolvent))**( nrhoSolvent ) * LOG(BrhoSolvent) * LOG(1- (T/CrhoSolvent))  )**2 * U_nrho_sq

	! Calculate uncertainty in molar volume of solute (A), Va: (U_{V_A}/V_A)^2
	U_MWA_sq = ( U_MW_relative*MWsolute )**2 		! U_{MW_A}^2
	U_VaOVa_sq = (U_MWA_sq/(MWsolute**2)) + U_rhoAOrhoA_sq   

	! Calculate uncertainty in molar volume of solvant (B) squared: (U_{V_B}/V_B)^2
	U_MWB_sq = ( U_MW_relative*MWsolvent )**2 		! U_{MW_A}^2
	U_VbOVb_sq = ( U_MWB_sq/(MWsolvent**2) ) + U_rhoBOrhoB_sq

	! Calculate uncertainty in viscosity of solvant (B) squared: ( U_{\eta_B} / \eta_{B} )^2
	U_AetaB_sq = (0.01*AandradeSolvent)**2 			! U_{A}^2 = ( 0.01*A )^2
	U_BetaB_sq = (0.01*BandradeSolvent)**2			! U_{B}^2 = ( 0.01*B )^2
	U_etaBOetaB_sq = U_AetaB_sq + ( U_BetaB_sq / (T**2) ) + ( BandradeSolvent/T )**2 * (U_ToT_sq)

	! Define uncertainties in surface tension
	U_sigmaAosigmaA_sq = (U_sigma_relative)**2  	! ( U_sigmaA / sigmaA )**2
	U_sigmaBosigmaB_sq = (U_sigma_relative)**2 		! ( U_sigmaB / sigmaB )**2. 
	
	! Calculate relative uncertainty in D_{AB}: ( U_{D_{AB}}/D_{AB} )^2
	U_DabODab_sq = U_ToT_sq + 0.187489*U_VaOVa_sq + 0.071289*U_VbOVb_sq + U_etaBOetaB_sq + &
					0.0225*U_sigmaAosigmaA_sq + 0.0225*U_sigmaBosigmaB_sq

	! Calculate uncertainty in D_{AB}: U_{D_{AB}}
	U_Dab = DAB*SQRT( U_DabODab_sq )

	! Calculate 95% expanded uncertainty in D_{AB}:
	! See Coleman, Steele Eqn. 3.19
	U_Dab_95 = 2.0*U_Dab   

	! WRITE(*,100) (U_Dab_95 / DAB )*100
	! 100 FORMAT("U_Dab_95 / DAB ...........", ES14.6,' % (95 percent confidence)')	


END SUBROUTINE DABuncertainty




