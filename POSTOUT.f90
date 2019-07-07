! 
!======================================================================! 
!     Output the flow field to the data file "flwout.dat" which can    ! 
!     be opened by Tecplot.                                            ! 
!======================================================================! 
      SUBROUTINE POSTOUT(NSTEP) 
      USE MODGLOB 
      IMPLICIT REAL*8(A-H,O-Z) 
      ALLOCATABLE :: XC(:,:), YC(:,:) 
	  INTEGER     :: NSTEP
	  CHARACTER(7):: CSTEP
! 
      ALLOCATE( XC(0:IN,0:JN), YC(0:IN,0:JN) )
	  WRITE(CSTEP,'(I7.7)') NSTEP 
! 
      DO 11 J=1,JN-1 
      DO 11 I=1,IN-1 
        XC(I,J) = 0.25D0*( X(I,J) + X(I+1,J) + X(I+1,J+1) + X(I,J+1) ) 
        YC(I,J) = 0.25D0*( Y(I,J) + Y(I+1,J) + Y(I+1,J+1) + Y(I,J+1) ) 
11    CONTINUE 
! 
      DO 12 J=1,JN-1 
        XC(0 ,J) = 0.5D0*( X(1 ,J) + X(1 ,J+1) ) 
        YC(0 ,J) = 0.5D0*( Y(1 ,J) + Y(1 ,J+1) ) 
        XC(IN,J) = 0.5D0*( X(IN,J) + X(IN,J+1) ) 
        YC(IN,J) = 0.5D0*( Y(IN,J) + Y(IN,J+1) ) 
12    CONTINUE 
! 
      DO 13 I=1,IN-1 
        XC(I,0 ) = 0.5D0*( X(I,1 ) + X(I+1,1 ) ) 
        YC(I,0 ) = 0.5D0*( Y(I,1 ) + Y(I+1,1 ) ) 
        XC(I,JN) = 0.5D0*( X(I,JN) + X(I+1,JN) ) 
        YC(I,JN) = 0.5D0*( Y(I,JN) + Y(I+1,JN) ) 
13    CONTINUE 
! 
      XC(0 ,0 ) = X(1 ,1 ) 
      YC(0 ,0 ) = Y(1 ,1 ) 
      XC(IN,0 ) = X(IN,1 ) 
      YC(IN,0 ) = Y(IN,1 ) 
      XC(IN,JN) = X(IN,JN) 
      YC(IN,JN) = Y(IN,JN) 
      XC(0 ,JN) = X(1 ,JN) 
      YC(0 ,JN) = Y(1 ,JN) 
! 
      DO 21 J=1,JN-1 
        DENS(0 ,J) = 0.5D0*( DENS(0 ,J) + DENS(1   ,J) ) 
        VELX(0 ,J) = 0.5D0*( VELX(0 ,J) + VELX(1   ,J) ) 
        VELY(0 ,J) = 0.5D0*( VELY(0 ,J) + VELY(1   ,J) ) 
        PRES(0 ,J) = 0.5D0*( PRES(0 ,J) + PRES(1   ,J) ) 
        DENS(IN,J) = 0.5D0*( DENS(IN,J) + DENS(IN-1,J) ) 
        VELX(IN,J) = 0.5D0*( VELX(IN,J) + VELX(IN-1,J) ) 
        VELY(IN,J) = 0.5D0*( VELY(IN,J) + VELY(IN-1,J) ) 
        PRES(IN,J) = 0.5D0*( PRES(IN,J) + PRES(IN-1,J) ) 
21    CONTINUE 
! 
      DO 22 I=1,IN-1 
        DENS(I,0 ) = 0.5D0*( DENS(I,0 ) + DENS(I,1   ) ) 
        VELX(I,0 ) = 0.5D0*( VELX(I,0 ) + VELX(I,1   ) ) 
        VELY(I,0 ) = 0.5D0*( VELY(I,0 ) + VELY(I,1   ) ) 
        PRES(I,0 ) = 0.5D0*( PRES(I,0 ) + PRES(I,1   ) ) 
        DENS(I,JN) = 0.5D0*( DENS(I,JN) + DENS(I,JN-1) ) 
        VELX(I,JN) = 0.5D0*( VELX(I,JN) + VELX(I,JN-1) ) 
        VELY(I,JN) = 0.5D0*( VELY(I,JN) + VELY(I,JN-1) ) 
        PRES(I,JN) = 0.5D0*( PRES(I,JN) + PRES(I,JN-1) ) 
22    CONTINUE 
! 
      DENS(0 ,0 ) = DENS(1   ,0 ) 
      VELX(0 ,0 ) = VELX(1   ,0 ) 
      VELY(0 ,0 ) = VELY(1   ,0 ) 
      PRES(0 ,0 ) = PRES(1   ,0 ) 
      DENS(0 ,JN) = DENS(1   ,JN) 
      VELX(0 ,JN) = VELX(1   ,JN) 
      VELY(0 ,JN) = VELY(1   ,JN) 
      PRES(0 ,JN) = PRES(1   ,JN) 
      DENS(IN,0 ) = DENS(IN-1,0 ) 
      VELX(IN,0 ) = VELX(IN-1,0 ) 
      VELY(IN,0 ) = VELY(IN-1,0 ) 
      PRES(IN,0 ) = PRES(IN-1,0 ) 
      DENS(IN,JN) = DENS(IN-1,JN) 
      VELX(IN,JN) = VELX(IN-1,JN) 
      VELY(IN,JN) = VELY(IN-1,JN) 
      PRES(IN,JN) = PRES(IN-1,JN) 
! 
      OPEN(71,FILE='pltflow_'//TRIM(CSTEP)//'.plt') 
      WRITE(71,*) 'TITLE="EULER_FLOW"' 
      WRITE(71,*) 'VARIABLES="X","Y","DENS","VELX","VELY","SPRE","TPRE","TTEM","MACH"' 
      WRITE(71,*) 'ZONE T="EULER_FLOW", I=',IN+1,', J=',JN+1,', F=POINT' 
      DO 31 J=0,JN 
      DO 31 I=0,IN 
        STEM = PRES(I,J)/(RGAS*DENS(I,J)) 
        VMAC = DSQRT((VELX(I,J)*VELX(I,J)+VELY(I,J)*VELY(I,J))/(GAMA*PRES(I,J)/DENS(I,J))) 
        TEMP = 1.D0 + GAMMH*VMAC*VMAC 
        TTEM = STEM*TEMP 
        TPRE = PRES(I,J)*TEMP**GAMSM 
        WRITE(71,7101) XC(I,J),YC(I,J),DENS(I,J),VELX(I,J),VELY(I,J),PRES(I,J),TPRE,TTEM,VMAC 
31    CONTINUE 
      CLOSE(71) 
7101  FORMAT(20(1X,F18.8)) 
! 
      RETURN 
      END 
!=========================== End of POSTOUT ============================ 

