!------------------------- BEGIN gckpp_blas.f90 BEGIN ------------------------
! @file gckpp_blas.f90                                                        
! @author yifanwu                                                             
! @date 2014-05-02 14:22:00.239209                                            
! @brief Basic linear algebra subprogram definitions                          
!                                                                             
! A reduced set of BLAS routines optimized for Kppa-generated solvers         
!                                                                             
! This file was generated by Kppa: http://www.paratools.com/Kppa              
!-----------------------------------------------------------------------------


MODULE gckpp_blas

  USE gckpp_parameters

  IMPLICIT NONE





  CONTAINS

!----------------------------------- WCOPY -----------------------------------
! Copies vector x to vector y: y <= x                                         
! Like the BLAS {S,D}COPY(N,X,1,Y,1)                                          
!                                                                             
! @param[in]     n Vector length                                              
! @param[in]     x Vector x                                                   
! @param[out]    y Vector y                                                   
!-----------------------------------------------------------------------------
  SUBROUTINE WCOPY(n, x, y)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: x(n)
    REAL(8), INTENT(OUT) :: y(n)

        y = x
  END SUBROUTINE WCOPY


!----------------------------------- WSCAL -----------------------------------
! Constant times a vector: x <= alpha*x                                       
! Like the BLAS {S,D}SCAL(N,alpha,X,1)                                        
!                                                                             
! @param[in]     n     Vector length                                          
! @param[in]     alpha Scalar                                                 
! @param[in,out] x     Vector x                                               
!-----------------------------------------------------------------------------
  SUBROUTINE WSCAL(n, alpha, x)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: alpha
    REAL(8), INTENT(INOUT) :: x(n)

        x = alpha * x
  END SUBROUTINE WSCAL


!----------------------------------- WAXPY -----------------------------------
! Constant times a vector plus a vector: y <= y + alpha*x                     
! Like the BLAS {S,D}AXPY(N,alpha,X,1,Y,1)                                    
!                                                                             
! @param[in]     n     Vector length                                          
! @param[in]     alpha Scalar                                                 
! @param[in]     x     Vector x                                               
! @param[in,out] y     Vector y                                               
!-----------------------------------------------------------------------------
  SUBROUTINE WAXPY(n, alpha, x, y)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: alpha
    REAL(8), INTENT(IN) :: x(n)
    REAL(8), INTENT(INOUT) :: y(n)

        y = y + alpha * x
  END SUBROUTINE WAXPY


!----------------------------------- WYMXDA ----------------------------------
! Difference of two vectors divided by a constant: z <= (y - x) / alpha       
!                                                                             
! @param[in]     n     Vector length                                          
! @param[in]     x     Vector x                                               
! @param[in]     y     Vector y                                               
! @param[in]     alpha Scalar                                                 
! @param[out]    z     Vector z                                               
!-----------------------------------------------------------------------------
  SUBROUTINE WYMXDA(n, x, y, alpha, z)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: x(n)
    REAL(8), INTENT(IN) :: y(n)
    REAL(8), INTENT(IN) :: alpha
    REAL(8), INTENT(OUT) :: z(n)

        z = (y - x)/alpha
  END SUBROUTINE WYMXDA


END MODULE gckpp_blas
!--------------------------- END gckpp_blas.f90 END --------------------------
