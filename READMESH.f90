! 
!======================================================================! 
!            Read the computational mesh of the CD nozzle.             ! 
!======================================================================! 
      SUBROUTINE READMESH 
      USE MODGLOB 
      IMPLICIT REAL*8(A-H,O-Z) 
! 
!------------------------- Generate the mesh  -------------------------- 
      OPEN ( UNIT=9, FILE='cdnoz.mesh', STATUS='OLD' )
		READ (9,*) JN, IN
		DO J = 1, JN
			DO I = 1, IN
				READ(9,*) X(I,J), Y(I,J)
			ENDDO
		ENDDO
	  CLOSE(9)
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
!=========================== End of READMESH ============================ 