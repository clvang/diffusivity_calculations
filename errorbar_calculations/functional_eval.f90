
SUBROUTINE functional_eval( tau_vector, LHS_vector, eps_vector, N, num_pts, expname )

  IMPLICIT NONE

  INTEGER (KIND=4), INTENT(IN) :: N, num_pts
  CHARACTER(LEN=7), INTENT(IN) :: expname
  REAL(KIND = 8), DIMENSION(N), INTENT(IN) :: tau_vector, LHS_vector, eps_vector

  INTEGER ( KIND = 4 ) :: i, j
  REAL (KIND=8), DIMENSION(num_pts) :: eps_values
  REAL (KIND=8), DIMENSION(N) :: f_value
  REAL (KIND=8) :: eps_min, eps_max, deps
  CHARACTER(LEN=16) :: filename_out, filename_header

  eps_min = MINVAL(eps_vector)
  eps_max = MAXVAL(eps_vector)
  deps = (eps_max - eps_min)/(num_pts - 1)  
  DO i = 1,num_pts
    eps_values(i) = eps_min + (i-1)*deps
  END DO

  filename_out = expname//"_data.dat"
  ! filename_header = expname//"_head.dat"
  ! OPEN(UNIT=10, FILE=filename_header, FORM="unformatted", ACCESS="stream", STATUS="unknown")
  ! OPEN(UNIT=20, FILE=filename_out, FORM="unformatted", ACCESS="stream", STATUS="unknown")
  OPEN(UNIT=20, FILE=filename_out)

  DO i = 1,1 !num_pts
    DO j = 1,N
      CALL Fx_Eval( tau_vector(j),eps_values(i),LHS_vector(j),f_value(j) )
    END DO    
     WRITE(20,*) eps_values(i), f_value
  END DO
  CLOSE(UNIT=20)

END SUBROUTINE functional_eval 


