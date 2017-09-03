PROGRAM molecularDiff

!  Purpose:
!    This program calculates the molecular species diffusivities
!	 of a liquid mixture with two components.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!     10/11/16    C. Vang             Original code
!	  10/14/16    C. Vang			  Added temperature range for Tsat
!									  Removed calculation of TsatSoluteNBP
!
IMPLICIT NONE


CHARACTER(len=30) :: filename, filenameOUT, filenameOUText, chemFormula
CHARACTER(len=100) :: junk
INTEGER :: status, i, j, N
REAL, DIMENSION(20) :: soluteProps, solventProps
REAL, DIMENSION(3) :: uncertaintyInputs
REAL, DIMENSION(300) :: Tvector
REAL :: TminRhoSolute, TmaxRhoSolute, ArhoSolute, BrhoSolute, CrhoSolute, &
nrhoSolute, TminAntoinneSolute, TmaxAntoinneSolute, AantoinneSolute, &
BantoinneSolute, CantoinneSolute, DantoinneSolute, TminAndradeSolute, &
TmaxAndradeSolute, AandradeSolute, BandradeSolute, PcriticalSolute, &
TcriticalSolute, MWsolute
REAL :: TminRhoSolvent, TmaxRhoSolvent, ArhoSolvent, BrhoSolvent, CrhoSolvent, &
nrhoSolvent, TminAntoinneSolvent, TmaxAntoinneSolvent, AantoinneSolvent, &
BantoinneSolvent, CantoinneSolvent, DantoinneSolvent, TminAndradeSolvent, &
TmaxAndradeSolvent, AandradeSolvent, BandradeSolvent, PcriticalSolvent, &
TcriticalSolvent, MWsolvent
REAL :: Pchamber, TsatSolventNBP, rhoSatSolventNBP, TsatSolventChamber, &
VSLsolventChamber, sigmaSolvent
REAL :: TsatSoluteNBP, rhoSoluteAtSolventTsatNBP, TsatSoluteChamber, &
VSLsoluteChamber, sigmaSolute
REAL :: DAB, DBA, maxTMIN, minTMAX 						 
REAL :: U_T_relative, U_MW_relative, U_sigma_relative, U_Dab_95, U_Dba_95    !quantities for uncertainty calculations


filename = "inputData.txt"

!Open the file and check for errors on Open 
OPEN(UNIT=1, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=status)

openif : IF (status == 0) THEN
	WRITE(*,*) 'OPENING OF FILE IS SUCCESSFUL!'

	readjunkLines: DO i=1,4 			!read lines 1-4
	  READ(1,*,IOSTAT=status) junk
	END DO readjunkLines

	READ(1,*,IOSTAT = status) Pchamber  !read line 5
	WRITE(*,5) Pchamber
	5 FORMAT(/,/,"Chamber Pressure [atm]: ", ES7.1)

	READ(1,*,IOSTAT=status) junk !read line 6 (headerline)

	!Read Solvent Properties (The component in excess)
	readSolventProps: DO i=1,20  !read lines 7-25
	  READ(1,*,IOSTAT=status) solventProps(i)
	  IF (status /= 0) EXIT
	END DO readSolventProps

	TminRhoSolvent = solventProps(1)
	TmaxRhoSolvent = solventProps(2)
	ArhoSolvent    = solventProps(3)
	BrhoSolvent    = solventProps(4)
	CrhoSolvent    = solventProps(5)
	nrhoSolvent    = solventProps(6)
	TminAntoinneSolvent = solventProps(7)
	TmaxAntoinneSolvent = solventProps(8)
	AantoinneSolvent    = solventProps(9)
	BantoinneSolvent    = solventProps(10)
	CantoinneSolvent    = solventProps(11)
	DantoinneSolvent    = solventProps(12)
	TminAndradeSolvent  = solventProps(13)
	TmaxAndradeSolvent  = solventProps(14)
	AandradeSolvent     = solventProps(15)
	BandradeSolvent     = solventProps(16)
	PcriticalSolvent    = solventProps(17)
	TcriticalSolvent    = solventProps(18)
	MWsolvent           = solventProps(19)
	sigmaSolvent        = solventProps(20)

	WRITE(*,*)"==============================================="

	WRITE(*,30) TminRhoSolvent, TmaxRhoSolvent, ArhoSolvent, BrhoSolvent, CrhoSolvent, &
	nrhoSolvent, TminAntoinneSolvent, TmaxAntoinneSolvent, AantoinneSolvent, &
	BantoinneSolvent, CantoinneSolvent, DantoinneSolvent, TminAndradeSolvent, &
	TmaxAndradeSolvent, AandradeSolvent, BandradeSolvent, PcriticalSolvent, &
	TcriticalSolvent, MWsolvent, sigmaSolvent

	30 FORMAT (/,"SOLVENT PROPERTIES: ",/,  &
		"Low Temp Range Density Eqn [K]......",ES14.6,/, &
		"High Temp Range Density Eqn [K].....",ES14.6,/, &
		"A coefficient Density Eqn...........",ES14.6,/, &
		"B coefficient Density Eqn...........",ES14.6,/, &
		"C coefficient Density Eqn...........",ES14.6,/, &
		"n coefficient Density Eqn...........",ES14.6,/,/, &		 
		"Low Temp Range Antoinne Eqn [K].....",ES14.6,/, &
		"High Temp Range Antoinne Eqn [K]....",ES14.6,/, &		
		"A coefficient Antoinne Eqn..........",ES14.6,/, &
		"B coefficient Antoinne Eqn..........",ES14.6,/, &
		"C coefficient Antoinne Eqn..........",ES14.6,/, &
		"D coefficient Antoinne Eqn..........",ES14.6,/,/, &
		"Low Temp Range Andrade Eqn [K]......",ES14.6,/, &
		"High Temp Range Andrade Eqn [K].....",ES14.6,/, &
		"A coefficient Andrade Eqn...........",ES14.6,/, &
		"B coefficient Andrade Eqn...........",ES14.6,/,/, &
		"Critical Pressure [bar].............",ES14.6,/, &
		"Critical Temperature [K]............",ES14.6,/, &
		"Molecular Weight [kg/kmol]..........",ES14.6,/, &
		"Surface Tension at BP [erg/cm^2]....",ES14.6)

	READ(1,*,IOSTAT=status) junk !read line 28 (another headerline)

	!Read solute properties (component NOT in excess)
	readSoluteProps: DO i=1,20   !read line 29-48
	  READ(1,*,IOSTAT=status) soluteProps(i)
	  IF (status /= 0) EXIT
	END DO readSoluteProps

	TminRhoSolute = soluteProps(1)
	TmaxRhoSolute = soluteProps(2)
	ArhoSolute    = soluteProps(3)
	BrhoSolute    = soluteProps(4)
	CrhoSolute    = soluteProps(5)
	nrhoSolute    = soluteProps(6)
	TminAntoinneSolute = soluteProps(7)
	TmaxAntoinneSolute = soluteProps(8)
	AantoinneSolute    = soluteProps(9)
	BantoinneSolute    = soluteProps(10)
	CantoinneSolute    = soluteProps(11)
	DantoinneSolute    = soluteProps(12)
	TminAndradeSolute  = soluteProps(13)
	TmaxAndradeSolute  = soluteProps(14)
	AandradeSolute     = soluteProps(15)
	BandradeSolute     = soluteProps(16)
	PcriticalSolute    = soluteProps(17)
	TcriticalSolute    = soluteProps(18)
	MWsolute           = soluteProps(19)
	sigmaSolute 	   = soluteProps(20)

	WRITE(*,20) TminRhoSolute, TmaxRhoSolute, ArhoSolute, BrhoSolute, CrhoSolute, &
	nrhoSolute, TminAntoinneSolute, TmaxAntoinneSolute, AantoinneSolute, &
	BantoinneSolute, CantoinneSolute, DantoinneSolute, TminAndradeSolute, &
	TmaxAndradeSolute, AandradeSolute, BandradeSolute, PcriticalSolute, &
	TcriticalSolute, MWsolute, sigmaSolute

	20 FORMAT (/,"SOLUTE PROPERTIES: ",/,  &
		"Low Temp Range Density Eqn [K]......",ES14.6,/, &
		"High Temp Range Density Eqn [K].....",ES14.6,/, &
		"A coefficient Density Eqn...........",ES14.6,/, &
		"B coefficient Density Eqn...........",ES14.6,/, &
		"C coefficient Density Eqn...........",ES14.6,/, &
		"n coefficient Density Eqn...........",ES14.6,/,/, &		 
		"Low Temp Range Antoinne Eqn [K].....",ES14.6,/, &
		"High Temp Range Antoinne Eqn [K]....",ES14.6,/, &		
		"A coefficient Antoinne Eqn..........",ES14.6,/, &
		"B coefficient Antoinne Eqn..........",ES14.6,/, &
		"C coefficient Antoinne Eqn..........",ES14.6,/, &
		"D coefficient Antoinne Eqn..........",ES14.6,/,/, &
		"Low Temp Range Andrade Eqn [K]......",ES14.6,/, &
		"High Temp Range Andrade Eqn [K].....",ES14.6,/, &
		"A coefficient Andrade Eqn...........",ES14.6,/, &
		"B coefficient Andrade Eqn...........",ES14.6,/,/, &
		"Critical Pressure [bar].............",ES14.6,/, &
		"Critical Temperature [K]............",ES14.6,/, &
		"Molecular Weight [kg/kmol]..........",ES14.6,/, &
		"Surface Tension at BP [erg/cm^2]....",ES14.6)

	READ(1,*, IOSTAT = status) chemFormula !read line 50
	READ(1,*,IOSTAT=status) junk !read line 52 (another headerline)

	!Read inputs required for uncertainty calculations
	DO i=1,3   !read line 53-55
	  READ(1,*,IOSTAT=status) uncertaintyInputs(i)
	  IF (status /= 0) EXIT
	END DO 
	U_T_relative = uncertaintyInputs(1)
	U_MW_relative = uncertaintyInputs(2)
	U_sigma_relative = uncertaintyInputs(3)

	CLOSE(UNIT=1)
ELSE

	WRITE(*,*) 'Open of file NOT sucessful!'

END IF openif
	WRITE(*,*)"==============================================="


!find max of Tmin 
maxTMIN = MAX(TminRhoSolvent, TminAntoinneSolvent, TminAndradeSolvent, &
	TminRhoSolute, TminAntoinneSolute, TminAndradeSolute)
!find min of Tmax
minTMAX = MIN(TmaxRhoSolvent, TmaxAntoinneSolvent, TmaxAndradeSolvent, &
	TmaxRhoSolute, TmaxAntoinneSolute, TmaxAndradeSolute)
WRITE(*,22) maxTMIN, minTMAX
22 FORMAT("Temperature Range of Validity for All Equations: [" ES14.6, " K", ES14.6, " K]" )

maxTMIN = CEILING(maxTMIN)	!round to nearest integer greater than or euqal to argument
minTMAX = FLOOR(minTMAX)	!round to nearest integer less than or equal to argument



!open file to write DAB and DBA calculated results
filenameOUText = "_OUT.txt"
filenameOUT = filename(1:9) // filenameOUText
OPEN(UNIT=4, FILE=filenameOUT, STATUS='REPLACE', ACTION='WRITE', IOSTAT = status)
WRITE(4,*) " DAB [m^2/s]     DBA [m^2/s]    Temperature [K]" !write header for output file


!!!!!!!!!!!!!!!!!!!!!!!!!!!  Calculate SOLVENT Properties !!!!!!!!!!!!!!!!!!!!!!!!!!!
!evaluate NPB tempearture of solvent using Antoinne Eqn.
CALL TsatAntoinne(101.325, AantoinneSolvent,&
	BantoinneSolvent,CantoinneSolvent,DantoinneSolvent, &
	TminAntoinneSolvent, TmaxAntoinneSolvent,TsatSolventNBP)

!create array of temperatures in which we would like to evalute
!values of DAB and DBA
! CALL linspace(Tvector, maxTMIN, minTMAX, 300)
CALL linspace(Tvector, maxTMIN, TsatSolventNBP, 300)

N = SIZE(Tvector)
loopOverTemp : DO j=1, N
	WRITE(*,24) Tvector(j), j 
	24 	FORMAT(/,/"============SPECIFIED TEMPERATURE: " ES14.6, "K, INDEX:", I3, "============") 

	WRITE(*,40) TsatSolventNBP
	40 FORMAT("NBP Tsat for Solvent is........." ES14.6, " K       @  1.0ATM")

	!evaluate density of solvent at NBP temperature
	CALL rhoCalc(ArhoSolvent,BrhoSolvent,CrhoSolvent,nrhoSolvent, &
		TminRhoSolvent, TmaxRhoSolvent, TsatSolventNBP, rhoSatSolventNBP)
	rhoSatSolventNBP = rhoSatSolventNBP*1000 !convert to [kg/m^3]
	WRITE(*,50) rhoSatSolventNBP, TsatSolventNBP
	50 FORMAT("NBP Solvent Density............." ES14.6, " kg/m^3  @", ES14.6, "K, 1.0ATM")

	!specify temperature at droplet surface, TsatSolventChamber
	TsatSolventChamber = Tvector(j)
	WRITE(*,60) TsatSolventChamber
	60 FORMAT("Solvent Surface Temperature....." ES14.6, " K       (SPECIFIED)" )

	!evaluate liquid viscosity of solvent at BP temperature 
	!corresponding to the chamber pressure
	CALL AndradeVSL(AandradeSolvent, BandradeSolvent, &
		TminAndradeSolvent, TmaxAndradeSolvent, &
		TsatSolventChamber, VSLsolventChamber)
	WRITE(*,70) VSLsolventChamber, TsatSolventChamber, Pchamber
	70 FORMAT("Solvent VSL ....................", ES14.6, " cP      @", ES14.6, "K, ",F3.1, "ATM")


	!!!!!!!!!!!!!!!!!!!!!!!!!!!  Calculate SOLUTE Properties !!!!!!!!!!!!!!!!!!!!!!!!!!!

	!evaluate solute density at NBP of solute
	CALL rhoCalc(ArhoSolute,BrhoSolute,CrhoSolute,nrhoSolute, &
		TminRhoSolute, TmaxRhoSolute, TsatSolventNBP, rhoSoluteAtSolventTsatNBP)
	rhoSoluteAtSolventTsatNBP = rhoSoluteAtSolventTsatNBP*1000 !convert to kg/m^3
	WRITE(*,90) rhoSoluteAtSolventTsatNBP, TsatSolventNBP
	90 FORMAT("Solute Density..................", ES14.6, " kg/m^3  @", ES14.6, "K, 1.0ATM")

	!evaluate liquid viscosity of solute using Andrade Eqn at
	!BP temperature corresponding to the chamber pressure of the SOLVENT
	CALL AndradeVSL(AandradeSolute, BandradeSolute, &
		TminAndradeSolute, TmaxAndradeSolute, &
		TsatSolventChamber, VSLsoluteChamber)
	WRITE(*,100) VSLsoluteChamber, TsatSolventChamber, Pchamber
	100 FORMAT("Solute VSL......................", ES14.6, " cP      @", ES14.6, "K, ", F3.1, "ATM")


	!!!!!!!!!!!!!!!!!!!!!!!!!!!  Calculate DAB & DBA !!!!!!!!!!!!!!!!!!!!!!!!!!!
	!evaluate DAB - solute (A) diffusing into solvant (B)
	CALL tynCalus(MWsolvent,rhoSatSolventNBP, &
		MWsolute,rhoSoluteAtSolventTsatNBP,VSLsolventChamber,&
		TsatSolventChamber, sigmaSolvent, sigmaSolute, DAB)
	DAB = DAB/(1E4)
	WRITE(*,110) DAB, TsatSolventChamber
	110 FORMAT("DAB............................", ES14.6, " m^2/s    @", ES14.6, "K")

	!calculate uncertainty in DAB 
	CALL DABuncertainty(MWsolvent,rhoSatSolventNBP, MWsolute,rhoSoluteAtSolventTsatNBP, &
				Pchamber*101.325, AantoinneSolvent, BantoinneSolvent, CantoinneSolvent, DantoinneSolvent, &
				TminAntoinneSolvent, TmaxAntoinneSolvent, &
				ArhoSolvent,BrhoSolvent,CrhoSolvent,nrhoSolvent, &
				ArhoSolute,BrhoSolute,CrhoSolute,nrhoSolute, &
				AandradeSolvent, BandradeSolvent, sigmaSolvent, sigmaSolute, DAB, &
				U_T_relative, U_MW_relative, U_sigma_relative, U_Dab_95,TsatSolventChamber )
	WRITE(*,112) U_Dab_95, TsatSolventChamber 
	112 FORMAT("-> 95% Uncertainty in DAB........ +/-", ES14.6, " m^2/s @", ES14.6, "K")
	WRITE(*,115) U_Dab_95/DAB * 100.0, TsatSolventChamber
	115 FORMAT("-> Relative Uncertainty in DAB...", ES14.6, " % @", ES14.6, "K")
	WRITE(*,113) DAB+U_Dab_95 
	113 FORMAT("-> DAB Upper Bound...............", ES14.6, " m^2/s ")
	WRITE(*,114) DAB-U_Dab_95
	114 FORMAT("-> DAB Lower Bound...............", ES14.6, " m^2/s")


	!evaluate DBA - solvant (B) diffusing into solute (A)
	CALL tynCalus(MWsolute, rhoSoluteAtSolventTsatNBP, &
		MWsolvent, rhoSatSolventNBP, VSLsoluteChamber, &
		TsatSolventChamber, sigmaSolute, sigmaSolvent, DBA)
	DBA = DBA/(1E4)
	WRITE(*,120) DBA,  TsatSolventChamber
	120 FORMAT("DBA............................", ES14.6, " m^2/s   @", ES14.6, "K")

	!calculate uncertainty in DBA
	CALL DABuncertainty(MWsolute,rhoSoluteAtSolventTsatNBP, MWsolvent,rhoSatSolventNBP, &
				Pchamber*101.325, AantoinneSolute, BantoinneSolute, CantoinneSolute, DantoinneSolute, &
				TminAntoinneSolute, TmaxAntoinneSolute, &
				ArhoSolute,BrhoSolute,CrhoSolute,nrhoSolute, &
				ArhoSolvent,BrhoSolvent,CrhoSolvent,nrhoSolvent, &
				AandradeSolute, BandradeSolute, sigmaSolute, sigmaSolvent, DBA, &
				U_T_relative, U_MW_relative, U_sigma_relative, U_Dba_95,TsatSolventChamber )
	WRITE(*,122) U_Dba_95, TsatSolventChamber 
	122 FORMAT("-> 95% Uncertainty in DBA........ +/-", ES14.6, " m^2/s @", ES14.6, "K")
	WRITE(*,124) U_Dba_95/DBA * 100.0, TsatSolventChamber
	124 FORMAT("-> Relative Uncertainty in DBA...", ES14.6, " % @", ES14.6, "K")
	WRITE(*,126) DBA+U_Dba_95 
	126 FORMAT("-> DBA Upper Bound...............", ES14.6, " m^2/s ")
	WRITE(*,128) DBA-U_Dba_95
	128 FORMAT("-> DBA Lower Bound...............", ES14.6, " m^2/s"/)

	!openoutput file to write DAB and DBA values
	WRITE(4,130) DAB, DBA, TsatSolventChamber
	130 FORMAT(ES14.6, " ", ES14.6, " ", ES14.6)

END DO loopOverTemp
	CLOSE(UNIT=4)


END PROGRAM molecularDiff




