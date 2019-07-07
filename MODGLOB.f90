! 
! 
!======================================================================! 
!   The module definition in which contains the global variables and   ! 
!   arrays. These variables and arrays can be "seen" in any subroutine ! 
!   when the module is included at the begining of the subroutine.     ! 
!======================================================================! 
      MODULE MODGLOB 
! 
!---The dimensions of the computational mesh; IN is the number of mesh 
!---points along streamwise; JN is that on the normal direction. 
      INTEGER IN, JN 
! 
!---NPASS is the index of the current time-marching step; NPASSM is the 
!---maximum time-marching loops. 
      INTEGER NPASS, NPASSM 
! 
!---The following parameters are in turn the static pressure at inlet, 
!---the density at inlet, the Mach number at inlet, the static pressure 
!---at outlet, and CFL number. 
      REAL*8 SPINL, DNINL, VMINL, SPOUT, CFL 
! 
!---The following are thermodynamic parameters and their derivations, see 
!---their definitions in the code. 
      REAL*8 GAMA, RGAS, GAMAM, GAMSM, GAMMH 
! 
!---X and Y are mesh coordinates. 
!---XLNI,YLNI,SLNI are length vectors and the magnitudes of the length of the I-direction cell surfaces. 
!---XLNJ,YLNJ,SLNJ are length vectors and the magnitudes of the length of the J-direction cell surfaces. 
!---SARC are cell areas of the control volumes divided by mesh lines. 
      REAL*8, ALLOCATABLE :: & 
        X(:,:), Y(:,:),& 
        XLNI(:,:), YLNI(:,:), SLNI(:,:), XLNJ(:,:), YLNJ(:,:), SLNJ(:,:), SARC(:,:) 
! 
!---DENS,VELX,VELY,PRES are density, Vx, Vy, and static pressure. For 
!---cell-centered finite volume method as applied in this code, these 
!---flow variables are saved at the cell centers. Additionally, the pseudo 
!---cells are used to apply the boundary conditions, so that these arrays 
!---include the values at centers of pseudo cells. 
!---FLUX is the residual (sum of the advective fluxes). 
      REAL*8, ALLOCATABLE :: & 
        DENS(:,:), VELX(:,:), VELY(:,:), PRES(:,:), FLUX(:,:,:) 
! 
      END MODULE 
!====================== End of Module definition ======================= 
