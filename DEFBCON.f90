! 
!======================================================================! 
!                Define the boundary conditions.                       ! 
!======================================================================! 
      SUBROUTINE DEFBCON 
      USE MODGLOB 
      IMPLICIT REAL*8(A-H,O-Z) 
! 
!---Define the B.C. at inlet: 
      TEMP = 1.D0 + GAMMH*VMINL**2 
      TPINL = SPINL*TEMP**GAMSM 
      STINL = SPINL/DNINL/RGAS 
      TTINL = STINL*TEMP 
      VXINL = VMINL*DSQRT(GAMA*RGAS*STINL) 
      DO 11 J=1,JN-1 
        VMCHB = DSQRT(VELX(1,J)**2+VELY(1,J)**2)/DSQRT(GAMA*PRES(1,J)/DENS(1,J)) 
        IF(VMCHB>=1.D0) THEN 
          DENS(0,J) = DNINL 
          VELX(0,J) = VXINL 
          VELY(0,J) = 0.D0 
          PRES(0,J) = SPINL 
        ELSE 
          SPI = DMIN1( PRES(1,J), TPINL ) 
          TMP = (TPINL/SPI)**(1.D0/GAMSM) 
          VMI = DSQRT((TMP-1.D0)/GAMMH) 
          STI = TTINL/TMP 
          DENS(0,J) = SPI/STI/RGAS 
          VELX(0,J) = VMI*DSQRT(GAMA*RGAS*STI) 
          VELY(0,J) = 0.D0 
          PRES(0,J) = SPI 
        END IF 
          DENS(-1,J) = 2.D0*DENS(0,J) - DENS(1,J) 
          VELX(-1,J) = 2.D0*VELX(0,J) - VELX(1,J) 
          VELY(-1,J) = 2.D0*VELY(0,J) - VELY(1,J) 
          PRES(-1,J) = 2.D0*PRES(0,J) - PRES(1,J) 
11    CONTINUE 
! 
!---Define the B.C. at outlet: 
      DO 21 J=1,JN-1 
        VMCHB = VELX(IN-1,J)/DSQRT(GAMA*PRES(IN-1,J)/DENS(IN-1,J)) 
        IF(VMCHB>=1.D0) THEN 
          PRES(IN,J) = PRES(IN-1,J) 
        ELSE 
          PRES(IN,J) = SPOUT 
        END IF 
        DENS(IN,J) = DENS(IN-1,J) 
        VELX(IN,J) = VELX(IN-1,J) 
        VELY(IN,J) = VELY(IN-1,J) 
        DENS(IN+1,J) = 2.D0*DENS(IN,J) - DENS(IN-1,J) 
        VELX(IN+1,J) = 2.D0*VELX(IN,J) - VELX(IN-1,J) 
        VELY(IN+1,J) = 2.D0*VELY(IN,J) - VELY(IN-1,J) 
        PRES(IN+1,J) = 2.D0*PRES(IN,J) - PRES(IN-1,J) 
21    CONTINUE 
! 
!---Define the wall B.C. on the upper side of the computational domain: 
      DO 31 I=1,IN-1 
        VNORM = VELX(I,JN-1)*XLNJ(I,JN) + VELY(I,JN-1)*YLNJ(I,JN) 
        VTEMP = 2.D0*VNORM/SLNJ(I,JN)/SLNJ(I,JN) 
        DENS(I,JN) = DENS(I,JN-1) 
        VELX(I,JN) = VELX(I,JN-1) - VTEMP*XLNJ(I,JN) 
        VELY(I,JN) = VELY(I,JN-1) - VTEMP*YLNJ(I,JN) 
        PRES(I,JN) = PRES(I,JN-1) 
        DENS(I,JN+1) = 2.D0*DENS(I,JN) - DENS(I,JN-1) 
        VELX(I,JN+1) = 2.D0*VELX(I,JN) - VELX(I,JN-1) 
        VELY(I,JN+1) = 2.D0*VELY(I,JN) - VELY(I,JN-1) 
        PRES(I,JN+1) = 2.D0*PRES(I,JN) - PRES(I,JN-1) 
31    CONTINUE 
! 
!---Define the symmetrical B.C. at the lower side of the domain: 
      DO 41 I=1,IN-1 
        VNORM = VELX(I,1)*XLNJ(I,1) + VELY(I,1)*YLNJ(I,1) 
        VTEMP = 2.D0*VNORM/SLNJ(I,1)/SLNJ(I,1) 
        DENS(I,0) = DENS(I,1) 
        VELX(I,0) = VELX(I,1) - VTEMP*XLNJ(I,1) 
        VELY(I,0) = VELY(I,1) - VTEMP*YLNJ(I,1) 
        PRES(I,0) = PRES(I,1) 
        DENS(I,-1) = 2.D0*DENS(I,0) - DENS(I,1) 
        VELX(I,-1) = 2.D0*VELX(I,0) - VELX(I,1) 
        VELY(I,-1) = 2.D0*VELY(I,0) - VELY(I,1) 
        PRES(I,-1) = 2.D0*PRES(I,0) - PRES(I,1) 
41    CONTINUE 
! 
      RETURN 
      END 
!=========================== End of DEFBCON ============================ 
