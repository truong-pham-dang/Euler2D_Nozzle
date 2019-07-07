! 
!======================================================================! 
!      Calculate the residual which is the sum of numerical flux.      ! 
!======================================================================! 
      SUBROUTINE RESIDUAL 
      USE MODGLOB 
      IMPLICIT NONE 
      REAL*8::TEMAL,TEMAR,TEMBL,TEMBR,TEMDL,TEMDR,VELXL,VELXR,VELYL,VELYR,& 
              PRESL,PRESR,DENSL,DENSR,TEU1L,TEU1R,TEU2L,TEU2R,TEU3L,TEU3R,& 
              TEU4L,TEU4R,AVDENS,AVPRES,AVVELX,AVVELY,AVENTH,ENTHL,ENTHR,& 
              AVERQ2,AVSONI,AVVELT,VECTNX,VECTNY,DETVX,DETVY,DETVE,DELTA,& 
              VELTL,VELTR,DETF1,DETF2,DETF3,DETF4 
      REAL*8::TEMP1,TEMP2,TEMP3,TEMP4 
      REAL*8,PARAMETER::EPSL = 1.D0,PARA1 = 0.5D0 
	  REAL*8,ALLOCATABLE::ALUM(:),FLUXI(:,:,:),FLUXJ(:,:,:) 
      INTEGER::I,J,K 
      ALLOCATE(ALUM(3),FLUXI(4,IN,JN),FLUXJ(4,IN,JN)) 
! 
      DO 11 J = 1,JN-1 
      DO 11 I = 1,IN 
        IF(VMINL > 10.D0)THEN 
! 
!---MUSCL interpolation of conservative variables with limiter function when supersonic 
!!---U1 
          TEMP1 = DENS(I-2,J) 
          TEMP2 = DENS(I-1,J) 
          TEMP3 = DENS(I  ,J) 
          TEMP4 = DENS(I+1,J) 
          TEMAL = TEMP3 - TEMP2 
          TEMAR = TEMP4 - TEMP3 
          TEMBL = TEMP2 - TEMP1 
          TEMBR = TEMP3 - TEMP2 
          TEMDL = (TEMAL*(TEMBL**2 + EPSL) + TEMBL*(TEMAL**2 + EPSL))/& 
                  (TEMAL**2 + TEMBL**2 + 2.D0*EPSL) 
          TEMDR = (TEMAR*(TEMBR**2 + EPSL) + TEMBR*(TEMAR**2 + EPSL))/& 
                  (TEMAR**2 + TEMBR**2 + 2.D0*EPSL) 
          TEU1L = TEMP2 + 0.5D0*TEMDL 
          TEU1R = TEMP3 - 0.5D0*TEMDR 
!!---U2 
          TEMP1 = DENS(I-2,J)*VELX(I-2,J) 
          TEMP2 = DENS(I-1,J)*VELX(I-1,J) 
          TEMP3 = DENS(I  ,J)*VELX(I  ,J) 
          TEMP4 = DENS(I+1,J)*VELX(I+1,J) 
          TEMAL = TEMP3 - TEMP2 
          TEMAR = TEMP4 - TEMP3 
          TEMBL = TEMP2 - TEMP1 
          TEMBR = TEMP3 - TEMP2 
          TEMDL = (TEMAL*(TEMBL**2 + EPSL) + TEMBL*(TEMAL**2 + EPSL))/& 
                  (TEMAL**2 + TEMBL**2 + 2.D0*EPSL) 
          TEMDR = (TEMAR*(TEMBR**2 + EPSL) + TEMBR*(TEMAR**2 + EPSL))/& 
                  (TEMAR**2 + TEMBR**2 + 2.D0*EPSL) 
          TEU2L = TEMP2 + 0.5D0*TEMDL 
          TEU2R = TEMP3 - 0.5D0*TEMDR 
!!---U3 
          TEMP1 = DENS(I-2,J)*VELY(I-2,J) 
          TEMP2 = DENS(I-1,J)*VELY(I-1,J) 
          TEMP3 = DENS(I  ,J)*VELY(I  ,J) 
          TEMP4 = DENS(I+1,J)*VELY(I+1,J) 
          TEMAL = TEMP3 - TEMP2 
          TEMAR = TEMP4 - TEMP3 
          TEMBL = TEMP2 - TEMP1 
          TEMBR = TEMP3 - TEMP2 
          TEMDL = (TEMAL*(TEMBL**2 + EPSL) + TEMBL*(TEMAL**2 + EPSL))/& 
                  (TEMAL**2 + TEMBL**2 + 2.D0*EPSL) 
          TEMDR = (TEMAR*(TEMBR**2 + EPSL) + TEMBR*(TEMAR**2 + EPSL))/& 
                  (TEMAR**2 + TEMBR**2 + 2.D0*EPSL) 
          TEU3L = TEMP2 + 0.5D0*TEMDL 
          TEU3R = TEMP3 - 0.5D0*TEMDR 
!!---U4 
          TEMP1 = PRES(I-2,J)/GAMAM + 0.5D0*DENS(I-2,J)*(VELX(I-2,J)**2+VELY(I-2,J)**2) 
          TEMP2 = PRES(I-1,J)/GAMAM + 0.5D0*DENS(I-1,J)*(VELX(I-1,J)**2+VELY(I-1,J)**2) 
          TEMP3 = PRES(I  ,J)/GAMAM + 0.5D0*DENS(I  ,J)*(VELX(I  ,J)**2+VELY(I  ,J)**2) 
          TEMP4 = PRES(I+1,J)/GAMAM + 0.5D0*DENS(I+1,J)*(VELX(I+1,J)**2+VELY(I+1,J)**2) 
          TEMAL = TEMP3 - TEMP2 
          TEMAR = TEMP4 - TEMP3 
          TEMBL = TEMP2 - TEMP1 
          TEMBR = TEMP3 - TEMP2 
          TEMDL = (TEMAL*(TEMBL**2 + EPSL) + TEMBL*(TEMAL**2 + EPSL))/& 
                  (TEMAL**2 + TEMBL**2 + 2.D0*EPSL) 
          TEMDR = (TEMAR*(TEMBR**2 + EPSL) + TEMBR*(TEMAR**2 + EPSL))/& 
                  (TEMAR**2 + TEMBR**2 + 2.D0*EPSL) 
          TEU4L = TEMP2 + 0.5D0*TEMDL 
          TEU4R = TEMP3 - 0.5D0*TEMDR 
!!---change back to flow variables 
          DENSL = TEU1L 
          DENSR = TEU1R 
          VELXL = TEU2L/TEU1L 
          VELXR = TEU2R/TEU1R 
          VELYL = TEU3L/TEU1L 
          VELYR = TEU3R/TEU1R 
          PRESL = GAMAM*(TEU4L - 0.5D0*DENSL*(VELXL**2 + VELYL**2)) 
          PRESR = GAMAM*(TEU4R - 0.5D0*DENSR*(VELXR**2 + VELYR**2)) 
        ELSE 
          DENSL = DENS(I-1,J) 
          DENSR = DENS(I  ,J) 
          VELXL = VELX(I-1,J) 
          VELXR = VELX(I  ,J) 
          VELYL = VELY(I-1,J) 
          VELYR = VELY(I  ,J) 
          PRESL = PRES(I-1,J) 
          PRESR = PRES(I  ,J) 
        END IF 
!---Roe scheme: 
!!---AVERAGED variables 
!				PRINT*,DENS(I,J),I,J 
        TEMP1 = DSQRT(DENSL) 
        TEMP2 = DSQRT(DENSR) 
        TEMP3 = TEMP1 + TEMP2 
        AVDENS = TEMP1 * TEMP2 
        AVVELX = (VELXL*TEMP1 + VELXR*TEMP2)/TEMP3 
        AVVELY = (VELYL*TEMP1 + VELYR*TEMP2)/TEMP3 
        ENTHL = GAMSM*PRESL/DENSL + 0.5D0*(VELXL**2 + VELYL**2) 
        ENTHR = GAMSM*PRESR/DENSR + 0.5D0*(VELXR**2 + VELYR**2) 
        AVENTH = (ENTHL*TEMP1 + ENTHR*TEMP2)/TEMP3 
        AVERQ2 = AVVELX**2 + AVVELY**2 
        AVSONI = DSQRT(GAMAM*(AVENTH - 0.5D0*AVERQ2)) 
        VECTNX = XLNI(I,J) / SLNI(I,J) 
        VECTNY = YLNI(I,J) / SLNI(I,J) 
        AVVELT = AVVELX*VECTNX + AVVELY*VECTNY 
!!---variables for residual 
        VELTL = VELXL*VECTNX + VELYL*VECTNY 
        VELTR = VELXR*VECTNX + VELYR*VECTNY 
        DETVX = VELXR - VELXL 
        DETVY = VELYR - VELYL 
        DETVE = DSQRT(VELXR**2 + VELYR**2) - DSQRT(VELXL**2 + VELYL**2) 
        ALUM(1) = DABS(AVVELT - AVSONI) 
        ALUM(2) = DABS(AVVELT) 
        ALUM(3) = DABS(AVVELT + AVSONI) 
!!---entropy correction 
        DELTA = PARA1*DSQRT(GAMA*PRESR/DENSR) 
        DO K = 1,3 
          IF(ALUM(K) <= DELTA)THEN 
            ALUM(K) = 0.5D0*(ALUM(K)**2 + DELTA**2)/DELTA 
          END IF 
        END DO 
!!---about the Roe matrix 
        DETF1 = ALUM(1)*(0.5D0*((PRESR-PRESL) - (AVDENS*AVSONI*DETVE))/(AVSONI**2)) 
        DETF2 = ALUM(2)*((DENSR - DENSL) - (PRESR - PRESL)/(AVSONI**2)) 
        DETF3 = ALUM(2)*AVDENS 
        DETF4 = ALUM(3)*(0.5D0*((PRESR-PRESL) + (AVDENS*AVSONI*DETVE))/(AVSONI**2)) 
!!---generate the results 
        FLUXI(1,I,J) = 0.5D0*(DENSR*VELTR + DENSL*VELTL - (DETF1 + DETF2 + DETF4)) 
        FLUXI(2,I,J) = 0.5D0*(DENSR*VELXR*VELTR + DENSL*VELXL*VELTL+(PRESR+PRESL)*VECTNX - & 
                      (DETF1*(AVVELX - AVSONI*VECTNX) + DETF2*(AVVELX) +& 
                       DETF3*(DETVX  - DETVE *VECTNX) + DETF4*(AVVELX + AVSONI*VECTNX))) 
        FLUXI(3,I,J) = 0.5D0*(DENSR*VELYR*VELTR + DENSL*VELYL*VELTL+(PRESR+PRESL)*VECTNY - & 
                      (DETF1*(AVVELY - AVSONI*VECTNY) + DETF2*(AVVELY) +& 
                       DETF3*(DETVY  - DETVE *VECTNY) + DETF4*(AVVELY + AVSONI*VECTNY))) 
        FLUXI(4,I,J) = 0.5D0*(DENSR*ENTHR*VELTR + DENSL*ENTHL*VELTL - & 
                      (DETF1*(AVENTH - AVSONI*AVVELT) + DETF2*(0.5D0*AVERQ2) +& 
                       DETF3*(AVVELX*DETVX + AVVELY*DETVY - AVVELT*DETVE) +& 
                       DETF4*(AVENTH + AVSONI*AVVELT)))	 
11		CONTINUE 
! 
      DO 12 J = 1,JN 
      DO 12 I = 1,IN-1 
        IF(VMINL > 10.D0)THEN 
! 
!---MUSCL interpolation of conservative variables with limiter function when supersonic 
!!---U1 
          TEMP1 = DENS(I,J-2) 
          TEMP2 = DENS(I,J-1) 
          TEMP3 = DENS(I,J  ) 
          TEMP4 = DENS(I,J+1) 
          TEMAL = TEMP3 - TEMP2 
          TEMAR = TEMP4 - TEMP3 
          TEMBL = TEMP2 - TEMP1 
          TEMBR = TEMP3 - TEMP2 
          TEMDL = (TEMAL*(TEMBL**2 + EPSL) + TEMBL*(TEMAL**2 + EPSL))/& 
                  (TEMAL**2 + TEMBL**2 + 2.D0*EPSL) 
          TEMDR = (TEMAR*(TEMBR**2 + EPSL) + TEMBR*(TEMAR**2 + EPSL))/& 
                  (TEMAR**2 + TEMBR**2 + 2.D0*EPSL) 
          TEU1L = TEMP2 + 0.5D0*TEMDL 
          TEU1R = TEMP3 - 0.5D0*TEMDR 
!!---U2 
          TEMP1 = DENS(I,J-2)*VELX(I,J-2) 
          TEMP2 = DENS(I,J-1)*VELX(I,J-1) 
          TEMP3 = DENS(I,J  )*VELX(I,J  ) 
          TEMP4 = DENS(I,J+1)*VELX(I,J+1) 
          TEMAL = TEMP3 - TEMP2 
          TEMAR = TEMP4 - TEMP3 
          TEMBL = TEMP2 - TEMP1 
          TEMBR = TEMP3 - TEMP2 
          TEMDL = (TEMAL*(TEMBL**2 + EPSL) + TEMBL*(TEMAL**2 + EPSL))/& 
                  (TEMAL**2 + TEMBL**2 + 2.D0*EPSL) 
          TEMDR = (TEMAR*(TEMBR**2 + EPSL) + TEMBR*(TEMAR**2 + EPSL))/& 
                  (TEMAR**2 + TEMBR**2 + 2.D0*EPSL) 
          TEU2L = TEMP2 + 0.5D0*TEMDL 
          TEU2R = TEMP3 - 0.5D0*TEMDR 
!!---U3 
          TEMP1 = DENS(I,J-2)*VELY(I,J-2) 
          TEMP2 = DENS(I,J-1)*VELY(I,J-1) 
          TEMP3 = DENS(I,J  )*VELY(I,J  ) 
          TEMP4 = DENS(I,J+1)*VELY(I,J+1) 
          TEMAL = TEMP3 - TEMP2 
          TEMAR = TEMP4 - TEMP3 
          TEMBL = TEMP2 - TEMP1 
          TEMBR = TEMP3 - TEMP2 
          TEMDL = (TEMAL*(TEMBL**2 + EPSL) + TEMBL*(TEMAL**2 + EPSL))/& 
                  (TEMAL**2 + TEMBL**2 + 2.D0*EPSL) 
          TEMDR = (TEMAR*(TEMBR**2 + EPSL) + TEMBR*(TEMAR**2 + EPSL))/& 
                  (TEMAR**2 + TEMBR**2 + 2.D0*EPSL) 
          TEU3L = TEMP2 + 0.5D0*TEMDL 
          TEU3R = TEMP3 - 0.5D0*TEMDR 
!!---U4 
          TEMP1 = PRES(I,J-2)/GAMAM + 0.5D0*DENS(I,J-2)*(VELX(I,J-2)**2+VELY(I,J-2)**2) 
          TEMP2 = PRES(I,J-1)/GAMAM + 0.5D0*DENS(I,J-1)*(VELX(I,J-1)**2+VELY(I,J-1)**2) 
          TEMP3 = PRES(I,J  )/GAMAM + 0.5D0*DENS(I,J  )*(VELX(I,J  )**2+VELY(I,J  )**2) 
          TEMP4 = PRES(I,J+1)/GAMAM + 0.5D0*DENS(I,J+1)*(VELX(I,J+1)**2+VELY(I,J+1)**2) 
          TEMAL = TEMP3 - TEMP2 
          TEMAR = TEMP4 - TEMP3 
          TEMBL = TEMP2 - TEMP1 
          TEMBR = TEMP3 - TEMP2 
          TEMDL = (TEMAL*(TEMBL**2 + EPSL) + TEMBL*(TEMAL**2 + EPSL))/& 
                  (TEMAL**2 + TEMBL**2 + 2.D0*EPSL) 
          TEMDR = (TEMAR*(TEMBR**2 + EPSL) + TEMBR*(TEMAR**2 + EPSL))/& 
                  (TEMAR**2 + TEMBR**2 + 2.D0*EPSL) 
          TEU4L = TEMP2 + 0.5D0*TEMDL 
          TEU4R = TEMP3 - 0.5D0*TEMDR 
!!---change back to flow variables 
          DENSL = TEU1L 
          DENSR = TEU1R 
          VELXL = TEU2L/TEU1L 
          VELXR = TEU2R/TEU1R 
          VELYL = TEU3L/TEU1L 
          VELYR = TEU3R/TEU1R 
          PRESL = GAMAM*(TEU4L - 0.5D0*DENSL*(VELXL**2 + VELYL**2)) 
          PRESR = GAMAM*(TEU4R - 0.5D0*DENSR*(VELXR**2 + VELYR**2)) 
        ELSE 
          DENSL = DENS(I,J-1) 
          DENSR = DENS(I,J  ) 
          VELXL = VELX(I,J-1) 
          VELXR = VELX(I,J  ) 
          VELYL = VELY(I,J-1) 
          VELYR = VELY(I,J  ) 
          PRESL = PRES(I,J-1) 
          PRESR = PRES(I,J  ) 
        END IF 
! 
!---Roe scheme: 
!!---AVERAGED variables 
        TEMP1 = DSQRT(DENSL) 
        TEMP2 = DSQRT(DENSR) 
        TEMP3 = TEMP1 + TEMP2 
        AVDENS = TEMP1 * TEMP2 
        AVVELX = (VELXL*TEMP1 + VELXR*TEMP2)/TEMP3 
        AVVELY = (VELYL*TEMP1 + VELYR*TEMP2)/TEMP3 
        ENTHL = GAMSM*PRESL/DENSL + 0.5D0*(VELXL**2 + VELYL**2) 
        ENTHR = GAMSM*PRESR/DENSR + 0.5D0*(VELXR**2 + VELYR**2) 
        AVENTH = (ENTHL*TEMP1 + ENTHR*TEMP2)/TEMP3 
        AVERQ2 = AVVELX**2 + AVVELY**2 
        AVSONI = DSQRT(GAMAM*(AVENTH - 0.5D0*AVERQ2)) 
        VECTNX = XLNJ(I,J) / SLNJ(I,J) 
        VECTNY = YLNJ(I,J) / SLNJ(I,J) 
        AVVELT = AVVELX*VECTNX + AVVELY*VECTNY 
!!---variables for residual 
        VELTL = VELXL*VECTNX + VELYL*VECTNY 
        VELTR = VELXR*VECTNX + VELYR*VECTNY 
        DETVX = VELXR - VELXL 
        DETVY = VELYR - VELYL 
        DETVE = DSQRT(VELXR**2 + VELYR**2) - DSQRT(VELXL**2 + VELYL**2) 
        ALUM(1) = DABS(AVVELT - AVSONI) 
        ALUM(2) = DABS(AVVELT) 
        ALUM(3) = DABS(AVVELT + AVSONI)  
!!---entropy correction 
        DELTA = PARA1*DSQRT(GAMA*PRESR/DENSR) 
        DO K = 1,3 
          IF(ALUM(K) <= DELTA)THEN 
            ALUM(K) = 0.5D0*(ALUM(K)**2 + DELTA**2)/DELTA 
          END IF 
        END DO 
!!---about the Roe matrix 
        DETF1 = ALUM(1)*(0.5D0*((PRESR-PRESL) - (AVDENS*AVSONI*DETVE))/(AVSONI**2)) 
        DETF2 = ALUM(2)*((DENSR - DENSL) - (PRESR - PRESL)/(AVSONI**2)) 
        DETF3 = ALUM(2)*AVDENS 
        DETF4 = ALUM(3)*(0.5D0*((PRESR-PRESL) + (AVDENS*AVSONI*DETVE))/(AVSONI**2)) 
!!---generate the results 
        FLUXJ(1,I,J) = 0.5D0*(DENSR*VELTR + DENSL*VELTL - (DETF1 + DETF2 + DETF4)) 
        FLUXJ(2,I,J) = 0.5D0*(DENSR*VELXR*VELTR + DENSL*VELXL*VELTL+(PRESR+PRESL)*VECTNX - & 
                      (DETF1*(AVVELX - AVSONI*VECTNX) + DETF2*(AVVELX) +& 
                       DETF3*(DETVX  - DETVE *VECTNX) + DETF4*(AVVELX + AVSONI*VECTNX))) 
        FLUXJ(3,I,J) = 0.5D0*(DENSR*VELYR*VELTR + DENSL*VELYL*VELTL+(PRESR+PRESL)*VECTNY - & 
                      (DETF1*(AVVELY - AVSONI*VECTNY) + DETF2*(AVVELY) +& 
                       DETF3*(DETVY  - DETVE *VECTNY) + DETF4*(AVVELY + AVSONI*VECTNY))) 
        FLUXJ(4,I,J) = 0.5D0*(DENSR*ENTHR*VELTR + DENSL*ENTHL*VELTL - & 
                      (DETF1*(AVENTH - AVSONI*AVVELT) + DETF2*(0.5D0*AVERQ2) +& 
                       DETF3*(AVVELX*DETVX + AVVELY*DETVY - AVVELT*DETVE) +& 
                       DETF4*(AVENTH + AVSONI*AVVELT)))	 
12		CONTINUE 
      DO 31 J=1,JN-1 
      DO 31 I=1,IN-1 
        FLUX(1,I,J)=FLUXI(1,I+1,J)*SLNI(I+1,J)-FLUXI(1,I,J)*SLNI(I,J)& 
               +FLUXJ(1,I,J+1)*SLNJ(I,J+1)-FLUXJ(1,I,J)*SLNJ(I,J) 
        FLUX(2,I,J)=FLUXI(2,I+1,J)*SLNI(I+1,J)-FLUXI(2,I,J)*SLNI(I,J)& 
               +FLUXJ(2,I,J+1)*SLNJ(I,J+1)-FLUXJ(2,I,J)*SLNJ(I,J) 
        FLUX(3,I,J)=FLUXI(3,I+1,J)*SLNI(I+1,J)-FLUXI(3,I,J)*SLNI(I,J)& 
               +FLUXJ(3,I,J+1)*SLNJ(I,J+1)-FLUXJ(3,I,J)*SLNJ(I,J) 
        FLUX(4,I,J)=FLUXI(4,I+1,J)*SLNI(I+1,J)-FLUXI(4,I,J)*SLNI(I,J)& 
               +FLUXJ(4,I,J+1)*SLNJ(I,J+1)-FLUXJ(4,I,J)*SLNJ(I,J) 
31    CONTINUE						 
! 
      RETURN 
      END 
!=========================== End of ADVFLUX ============================ 