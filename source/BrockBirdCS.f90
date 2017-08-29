! This funciton calculates the surface tension of a non-polard liquid using 
! the Brock-Bird corresponding states method (See Eqn. 1-9.6 in Poling et al.
! The Properties of Gases and Liquids).
! INPUT:  
!    Pc - Critical pressure of pure liquid [bar]
!    Tc - Critical temperature of pure liquid [Kelvin]
!    Tb - boiling point of liquid [Kelvin]
! OUTPUT:
!    sigma - surface tension of liquid [erg/cm^2]
! This funciton calculates the surface tension of a non-polard liquid using 
! the Brock-Bird corresponding states method (See Eqn. 1-9.6 in Poling et al.
! The Properties of Gases and Liquids).
! INPUT:  
!    Pc - Critical pressure of pure liquid [bar]
!    Tc - Critical temperature of pure liquid [Kelvin]
!    Tb - boiling point of liquid [Kelvin]
! OUTPUT:
!    sigma - surface tension of liquid [erg/cm^2]


SUBROUTINE BrockBirdCS(Pc, Tc, Tb, sigma)

	REAL, INTENT(IN) :: Pc, Tc, Tb
	REAL, INTENT(OUT) :: sigma
	REAL :: Tbr, alpha

	Tbr = Tb/Tc;

	alpha = 0.9076*( 1 + ( Tbr*LOG(Pc / 1.013) / (1 - Tbr) ) );

	sigma = Pc**(2./3.) * Tc**(1/3) * (0.132*alpha - 0.279) * (1 - Tbr)**(11./9.);

END SUBROUTINE BrockBirdCS