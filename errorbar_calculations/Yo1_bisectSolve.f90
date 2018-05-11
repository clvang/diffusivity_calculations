
! subroutine to solve for mass fraction of high volatility component when the
! igniters are turned on
SUBROUTINE Yo1_bisectSolve( rho_l_1, rho_h_1, Y_ratio_o, d_ratio, err_tol, N, Yl1_vector )

  IMPLICIT NONE

  INTEGER (KIND=4), INTENT(IN) :: N
  REAL(KIND = 8), DIMENSION(N), INTENT(IN) :: rho_l_1, rho_h_1, Y_ratio_o, d_ratio
  REAL (KIND=8), INTENT(IN) :: err_tol  
  REAL (KIND=8), DIMENSION(N), INTENT(OUT) :: Yl1_vector

  REAL(KIND=8) :: a, b, f_p, f_a, f_b, p
  INTEGER :: maxit, i, j
  ! REAL (KIND=8), DIMENSION(N) :: eps_mc


  DO j=1,N 

      a = 0.0
      b = 1.0
      f_a = -a/( ((1-a)/rho_h_1(j)) + (a/rho_l_1(j)) ) - d_ratio(j)*Y_ratio_o(j)
      f_b = -b/( ((1-b)/rho_h_1(j)) + (b/rho_l_1(j)) ) - d_ratio(j)*Y_ratio_o(j)

      IF ( f_a * f_b > 0 ) THEN 
          WRITE(*,10) a, b 
  10      FORMAT('!!!!! ROOT IS NOT BOUNDED IN: [', ES8.2, ',', ES8.2 '] !!!!!!!!!!!!')  
      ELSE 

          maxit = 300
          DO i=1,maxit
              p = a + (b - a)/2.0
              f_p = -p/( ((1-p)/rho_h_1(j)) + (p/rho_l_1(j)) ) - d_ratio(j)*Y_ratio_o(j) 

              IF ( ( ABS(f_p) < err_tol ) .OR. ( (b-a)/2. < err_tol ) ) EXIT 

              IF ( f_a*f_p > 0 ) THEN 
                  a = p 
                  f_a = f_p 
              ELSE 
                  b = p 
              END IF 
          END DO   ! DO i
  !             WRITE(*,20) i, ABS(f_p)
  ! 20          FORMAT(' ',"iteration count:" I3,  "   f(epsilon):" ES14.6)     
      END IF 

      Yl1_vector(j) = p 
  !             WRITE(*,30) i, Yh1_vector(j)
  ! 30          FORMAT(' ',"iteration count:" I3,  "   Yh1:" ES14.6)           

    END DO ! DO j

END SUBROUTINE Yo1_bisectSolve 
