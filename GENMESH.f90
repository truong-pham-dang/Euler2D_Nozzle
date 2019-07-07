! 
!======================================================================! 
!           Generate the computational mesh of the nozzle.             ! 
!======================================================================! 
      SUBROUTINE GENMESH 
      USE MODGLOB 
      IMPLICIT REAL*8(A-H,O-Z) 
! 
!------------------------- Generate the mesh  -------------------------- 
      XINL = 0.D0 
      XOUT = 4.D0 
      DXTB = 1.D0 
      DYTB = 0.042D0 
! 
      RDTB = 0.125D0*DXTB/DYTB  + 0.5D0*DYTB 
      XCEN = 1.5D0 
      YCEN = DYTB - RDTB 
! 
      DO 11 I=1,IN 
        XNOW = XINL + REAL(I-1)/REAL(IN-1)*(XOUT-XINL) 
        IF(XNOW<=1.D0 .OR. XNOW>=2.D0) THEN 
          YWAL = 0.D0 
        ELSE 
          YWAL = YCEN + DSQRT( RDTB**2 - (XCEN-XNOW)**2 ) 
        END IF 
        DO 12 J=1,JN 
          YNOW = YWAL + REAL(J-1)/REAL(JN-1)*(1.D0-YWAL) 
          X(I,J) = XNOW 
          Y(I,J) = YNOW 
12      CONTINUE 
11    CONTINUE 
! 
      OPEN(71,FILE='pltgrid.plt') 
      WRITE(71,*) 'TITLE="NONAME"' 
      WRITE(71,*) 'VARIABLES="X","Y"' 
      WRITE(71,*) 'ZONE T="NONAME", I=',IN,',J=',JN,', F=POINT' 
      DO J=1,JN 
        DO I=1,IN 
          WRITE(71,*) X(I,J), Y(I,J) 
        END DO 
      END DO 
      CLOSE(71) 
! 
!------------ Calculate the surface length vectors and ----------------- 
!------------ the areas of 2D control volumes divided  ----------------- 
!------------ by mesh lines.                           ----------------- 
      DO 21 J=1,JN-1 
      DO 21 I=1,IN 
        XLNI(I,J) = Y(I,J+1) - Y(I,J  ) 
        YLNI(I,J) = X(I,J  ) - X(I,J+1) 
        SLNI(I,J) = DSQRT( XLNI(I,J)*XLNI(I,J) + YLNI(I,J)*YLNI(I,J) ) 
21    CONTINUE 
! 
      DO 22 J=1,JN 
      DO 22 I=1,IN-1 
        XLNJ(I,J) = Y(I  ,J) - Y(I+1,J) 
        YLNJ(I,J) = X(I+1,J) - X(I  ,J) 
        SLNJ(I,J) = DSQRT( XLNJ(I,J)*XLNJ(I,J) + YLNJ(I,J)*YLNJ(I,J) ) 
22    CONTINUE 
! 
      DO 23 J=1,JN-1 
      DO 23 I=1,IN-1 
        I1 = I + 1 
        J1 = J + 1 
        SARC(I,J) = 0.5D0*( DABS( ( X(I ,J ) - X(I1,J ) )*Y(I1,J1)& 
                                 +( X(I1,J ) - X(I1,J1) )*Y(I ,J )& 
                                 +( X(I1,J1) - X(I ,J ) )*Y(I1,J ) )& 
                          + DABS( ( X(I ,J ) - X(I1,J1) )*Y(I ,J1)& 
                                 +( X(I1,J1) - X(I ,J1) )*Y(I ,J )& 
                                 +( X(I ,J1) - X(I ,J ) )*Y(I1,J1) ) ) 
23    CONTINUE 
! 
      RETURN 
      END 
!=========================== End of GENMESH ============================ 