! 
!======================================================================! 
!     Update the conservative variables at time level n+1 using        ! 
!     a simple one-step time marching scheme.                          ! 
!======================================================================! 
      SUBROUTINE UPDATEQ 
      USE MODGLOB 
      IMPLICIT REAL*8(A-H,O-Z) 
      DIMENSION RESAVE(4) 
      ALLOCATABLE :: QCON(:,:,:) 
      SAVE QCON 
! 
      IF(.NOT.ALLOCATED(QCON)) THEN 
        ALLOCATE( QCON(4,IN-1,JN-1) ) 
        DO 11 J=1,JN-1 
        DO 11 I=1,IN-1 
          QCON(1,I,J) = DENS(I,J) 
          QCON(2,I,J) = DENS(I,J)*VELX(I,J) 
          QCON(3,I,J) = DENS(I,J)*VELY(I,J) 
          QCON(4,I,J) = PRES(I,J)/GAMAM + 0.5D0*DENS(I,J)*(VELX(I,J)*VELX(I,J)+VELY(I,J)*VELY(I,J)) 
11      CONTINUE 
        OPEN(70,FILE='monitor.dat') 
      END IF 
! 
      RESM1 = 0.D0 
      RESAVE = 0.D0 
      IMAX = 0 
      JMAX = 0 
      DO 21 J=1,JN-1 
      DO 21 I=1,IN-1 
! 
!---Calculate the local time step: 
        TXLNI = 0.5D0*( XLNI(I,J) + XLNI(I+1,J) ) 
        TYLNI = 0.5D0*( YLNI(I,J) + YLNI(I+1,J) ) 
        TSLNI = 0.5D0*( SLNI(I,J) + SLNI(I+1,J) ) 
        TXLNJ = 0.5D0*( XLNJ(I,J) + XLNJ(I,J+1) ) 
        TYLNJ = 0.5D0*( YLNJ(I,J) + YLNJ(I,J+1) ) 
        TSLNJ = 0.5D0*( SLNJ(I,J) + SLNJ(I,J+1) ) 
        VNMCI = TXLNI*VELX(I,J) + TYLNI*VELY(I,J) 
        VNMCJ = TXLNJ*VELX(I,J) + TYLNJ*VELY(I,J) 
        SONIC = DSQRT(GAMA*PRES(I,J)/DENS(I,J)) 
        CHVEL = DABS(VNMCI) + DABS(VNMCJ) + SONIC*(TSLNI+TSLNJ) 
        DTCON = CFL/CHVEL 
! 
!---Update the flow variables: 
        QCON(1:4,I,J) = QCON(1:4,I,J) - DTCON*FLUX(1:4,I,J)       ! THE MAJOR FUNCTION IS JUST A SENTENCE 
        DENS(I,J) = QCON(1,I,J) 
        VELX(I,J) = QCON(2,I,J)/QCON(1,I,J) 
        VELY(I,J) = QCON(3,I,J)/QCON(1,I,J) 
        PRES(I,J) = GAMAM*( QCON(4,I,J) - 0.5D0*DENS(I,J)*(VELX(I,J)*VELX(I,J)+VELY(I,J)*VELY(I,J)) ) 
! 
        TRES1 = DABS(FLUX(1,I,J)) 
        IF(TRES1>RESM1) THEN 
          RESM1 = TRES1 
          IMAX = I 
          JMAX = J 
        END IF 
        RESAVE(1:4) = RESAVE(1:4) + DABS(FLUX(1:4,I,J)) 
21    CONTINUE 
      RESAVE = RESAVE/REAL(IN-1)/REAL(JN-1) 
! 
      IF(MOD(NPASS,10)==0) WRITE(70,7001) NPASS, IMAX, JMAX, RESAVE(1:4) 
      IF(MOD(NPASS,10)==0) WRITE(* ,7001) NPASS, IMAX, JMAX, RESAVE(1:4) 
7001  FORMAT(3(I6),10(1X,E10.4)) 
      IF(NPASS==NPASSM) CLOSE(70) 
! 
      RETURN 
      END 
!=========================== End of UPDATEQ ============================ 