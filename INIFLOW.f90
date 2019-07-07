! 
!======================================================================! 
!                Generate the initial flow field.                      ! 
!     The following method generates a mean flow field.                ! 
!======================================================================! 
      SUBROUTINE INIFLOW 
      USE MODGLOB 
      IMPLICIT REAL*8(A-H,O-Z) 
! 
      DNINI = DNINL 
      SPINI = SPINL 
      VXINI = VMINL*DSQRT(GAMA*SPINL/DNINL) 
! 
      DO 11 J=1,JN-1 
      DO 11 I=1,IN-1 
        DENS(I,J) = DNINI 
        VELX(I,J) = VXINI 
        VELY(I,J) = 0.D0 
        PRES(I,J) = SPINI 
11    CONTINUE 
! 
      RETURN 
      END 
!=========================== End of INIFLOW ============================ 
