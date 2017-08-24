SUBROUTINE linspace(varArray, lowLimit, upLimit, numPts) 
	IMPLICIT NONE

	REAL, INTENT(OUT), DIMENSION(numPts) :: varArray
	REAL, INTENT(IN) :: lowLimit, upLimit
	INTEGER, INTENT(IN) :: numPts
	INTEGER :: i
	REAL :: intervalSize

	intervalSize = (upLimit - lowLimit) / numPts
	varArray(1) = lowLimit
	DO i = 2, numPts-1
		varArray(i) = varArray(i-1) + intervalSize
	END DO
	varArray(1) = lowLimit
	varArray(numPts) = upLimit
END SUBROUTINE linspace