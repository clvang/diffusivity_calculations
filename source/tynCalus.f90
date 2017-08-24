
SUBROUTINE tynCalus(MWsolvent,rhoSATsolvent,MWsolute,rhoSATsolute,VSLsolvent,T, D)

!subroutine to calculate infinite diffusivities using Tyn-Calus equation
!(See Poling et. al, p 11.25, Eqn. 11-9.5)

	REAL, INTENT(IN) :: MWsolvent, rhoSATsolvent, MWsolute, rhoSATsolute, &
	VSLsolvent, T 
	REAL, INTENT(OUT) :: D
	REAL :: Vb, Va
	!compute molar volume [cm^3/mol] at NBP solvent
	Vb = (MWsolvent/rhoSATsolvent)*(1E3)	

	!compute molar volume[cm^3/mol] at NBP of solvent
	Va = (MWsolute/rhoSATsolute)*(1E3)
	D = 8.93E-08*(Vb**(0.267)/Va**(0.433) )*(T/VSLsolvent)
	
END SUBROUTINE tynCalus