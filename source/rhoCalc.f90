!subroutine to calculate liquid density at a temperature (See Yaw's)
SUBROUTINE rhoCalc(A, B, C, n, Tmin, Tmax, T, rho)
	IMPLICIT NONE
	REAL, INTENT(IN) :: A, B, C, n, Tmin, Tmax, T
	REAL, INTENT(OUT) :: rho 

	IF ( (T < Tmax).AND. (T > Tmin) )THEN
		rho = A*( B **(-(1-T/C)**n) )
	ELSE
		WRITE(*,*) "T is out of range of validity for density calculation"
	END IF
END SUBROUTINE rhoCalc