!subroutine to calculate liquid viscosity using Andrade Equation
SUBROUTINE AndradeVSL(A, B, Tmin, Tmax, T, VSL)
	IMPLICIT NONE
	REAL, INTENT(IN) :: A, B, Tmin, Tmax, T 
	REAL, INTENT(OUT) :: VSL 

	IF ( (T<Tmax).AND.(T>Tmin) ) THEN
		VSL = EXP(A + B/T )
	ELSE
		WRITE(*,*) "T is out of range in Viscosity Calculation"
	END IF
END SUBROUTINE AndradeVSL