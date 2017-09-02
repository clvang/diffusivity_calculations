
! 	Po - pressure in kPa

SUBROUTINE DABuncertainty(MWsolvent,rhoSatSolventNBP, MWsolute,rhoSoluteAtSolventTsatNBP, &
						Po, AantoinneSolvent, BantoinneSolvent, CantoinneSolvent, DantoinneSolvent, &
						TminAntoinneSolvent, TmaxAntoinneSolvent, &
						ArhoSolvent,BrhoSolvent,CrhoSolvent,nrhoSolvent, &
						ArhoSolute,BrhoSolute,CrhoSolute,nrhoSolute, &
						BandradeSolvent,sigmaSolvent, sigmaSolute, DAB, & 
						U_T_relative, U_MW_relative, U_sigma_relative )
	IMPLICIT NONE
	REAL, INTENT(IN) :: MWsolvent,rhoSatSolventNBP, MWsolute,rhoSoluteAtSolventTsatNBP, &
						Po, AantoinneSolvent, BantoinneSolvent, CantoinneSolvent, DantoinneSolvent, &
						TminAntoinneSolvent, TmaxAntoinneSolvent, &
						ArhoSolvent,BrhoSolvent,CrhoSolvent,nrhoSolvent, &
						ArhoSolute,BrhoSolute,CrhoSolute,nrhoSolute, &
						BandradeSolvent,sigmaSolvent, sigmaSolute, DAB, &
						U_T_relative, U_MW_relative, U_sigma_relative
	REAL, DIMENSION(10) :: dP, deltaDP, P_plus, P_minus, T_plus, T_minus, dTdP 
	INTEGER :: i, N 
	REAL :: Va, Vb, U_Pv_sq, U_ToT_sq, U_Arho_sq, U_Brho_sq, T, U_rhoAOrhoA_sq, &
			U_VaOVa_sq, U_MWA_sq, U_rhoBOrhoB_sq, U_VbOVb_sq, U_MWB_sq, U_AetaB_sq, &
			U_BetaB_sq, U_etaBOetaB_sq, U_sigmaAosigmaA_sq, &
			U_sigmaBosigmaB_sq, U_DabODab_sq, U_Dab, U_Dab_95

	!compute molar volume [cm^3/mol] at NBP solvent
	Vb = (MWsolvent/rhoSatSolventNBP)*(1E3)	

	!compute molar volume[cm^3/mol] at NBP of solvent
	Va = (MWsolute/rhoSoluteAtSolventTsatNBP)*(1E3)

	!Calculate relative uncertainty in T squared, (U_T/T)^2
	N = 10
	deltaDP(1) = 0.01*Po
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
	CALL TsatAntoinne(Po, AantoinneSolvent,&
		BantoinneSolvent,CantoinneSolvent,DantoinneSolvent, &
		TminAntoinneSolvent, TmaxAntoinneSolvent,T)

	U_Pv_sq = (0.1*Po)**2   					! U_{P_v}^2, Assume U_{P_v} = 10 percent  (U_Pv_sq not really required)

	! U_ToT_sq = (dTdP(N) / T)**2 * U_Pv_sq     ! (U_T/T)^2 
	U_ToT_sq = (U_T_relative/T)**2  					! (U_T/T)^2 assume U_T = 15 K
	WRITE(*,10) U_ToT_sq
	10 FORMAT("U_ToT_sq..................", ES14.6)

	!Calculate relative uncertainty in density of solute (A) squared, (U_rho_A/ rho_A)^2
	U_Arho_sq = 0.0 	! uncertainties in A coefficient in equation = 0
	U_Brho_sq = 0.0		! uncertainties in B coefficient in equation = 0
	U_rhoAOrhoA_sq = ( 1. / ArhoSolute)**2 * U_Arho_sq + &
		( (1. - (T/CrhoSolute))**(nrhoSolute) / BrhoSolute )**2 * U_Brho_sq + &
		( nrhoSolute*(1- (T/CrhoSolute))**( nrhoSolute - 1) * LOG(BrhoSolute) / CrhoSolute )**2 * (U_ToT_sq * T**2)
	WRITE(*,20) U_rhoAOrhoA_sq 					! (U_rho_A/ rho_A)^2
	20 FORMAT("U_rhoAOrhoA_sq............", ES14.6)

	!Calculate relative uncertainty in density of solvant (B) squared: (U_rho_B/ rho_B)^2
	U_Arho_sq = 0.0
	U_Brho_sq = 0.0
	U_rhoBOrhoB_sq = ( 1. / ArhoSolvent)**2 * U_Arho_sq + &
		( (1. - (T/CrhoSolvent))**(nrhoSolvent) / BrhoSolvent )**2 * U_Brho_sq + &
		( nrhoSolvent*(1- (T/CrhoSolvent))**( nrhoSolvent - 1) * LOG(BrhoSolvent) / CrhoSolvent )**2 * (U_ToT_sq * T**2)
	WRITE(*,40) U_rhoBOrhoB_sq 					! (U_rho_B/ rho_B)^2
	40 FORMAT("U_rhoBOrhoB_sq............", ES14.6)

	! Calculate uncerntainty in molar volume of solute (A), Va: (U_{V_A}/V_A)^2
	U_MWA_sq = ( U_MW_relative*MWsolute )**2 		!U_{MW_A}^2. Assume U_{MW_A} = 1 percent
	U_VaOVa_sq = (U_MWA_sq/(MWsolute**2)) + U_rhoAOrhoA_sq   
	WRITE(*,30) U_VaOVa_sq 						! (U_{V_A}/V_A)^2
	30 FORMAT("U_VaOVa_sq................", ES14.6)

	! Calculate uncerntainty in molar volume of solvant (B) squared: (U_{V_B}/V_B)^2
	U_MWB_sq = ( U_MW_relative*MWsolvent )**2 			! U_{MW_A}^2. Assume U_{MW_A} = 1 percent
	U_VbOVb_sq = ( U_MWB_sq/(MWsolvent**2) ) + U_rhoBOrhoB_sq
	WRITE(*,50) U_VbOVb_sq 						! (U_{V_B}/V_B)^2
	50 FORMAT("U_VbOVb_sq................", ES14.6)


	! Calculate uncertainty in viscosity of solvant (B) squared: ( U_{\eta_B} / \eta_{B} )^2
	U_AetaB_sq = 0.0 							! U_{A}^2
	U_BetaB_sq = 0.0 							! U_{B}^2
	U_etaBOetaB_sq = U_AetaB_sq + ( U_BetaB_sq / (T**2) ) + ( BandradeSolvent/T )**2 * (U_ToT_sq)
	WRITE(*,60) U_etaBOetaB_sq 					! (U_{V_B}/V_B)^2
	60 FORMAT("U_etaBOetaB_sq............", ES14.6)

	! Calculate relative uncertainty in D_{AB}: ( U_{D_{AB}}/D_{AB} )^2
	U_sigmaAosigmaA_sq = (U_sigma_relative)**2  		! ( U_sigmaA / sigmaA )**2. Assume relative uncertainty in sigmaA = 15 %
	U_sigmaBosigmaB_sq = (U_sigma_relative)**2 			! ( U_sigmaB / sigmaB )**2. Assume relative uncertainty in sigmaB = 15 %
	U_DabODab_sq = U_ToT_sq + 0.187489*U_VaOVa_sq + 0.071289*U_VbOVb_sq + U_etaBOetaB_sq + &
					0.0225*U_sigmaAosigmaA_sq + 0.0225*U_sigmaBosigmaB_sq
	WRITE(*,70) U_DabODab_sq 						
	70 FORMAT("U_DabODab_sq..............", ES14.6)

	! Calculate uncertainty in D_{AB}: U_{D_{AB}}
	U_Dab = DAB*SQRT( U_DabODab_sq )
	WRITE(*,80) U_Dab 						
	80 FORMAT("U_Dab.....................", ES14.6,' m^2/s')

	! Calculate 95% expanded uncertainty in D_{AB}:
	! See Coleman, Steele Eqn. 3.19
	U_Dab_95 = 2*U_Dab
	WRITE(*,90) U_Dab_95						
	90 FORMAT("U_Dab_95..................", ES14.6,' m^2/s')	

	WRITE(*,100) (U_Dab_95 / DAB )*100
	100 FORMAT("U_Dab_95 / DAB ...........", ES14.6,' %')	


END SUBROUTINE DABuncertainty




