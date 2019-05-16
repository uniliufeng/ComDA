      SUBROUTINE ATM_CONFIG_SET (
     O    GLON , GLAT , IMASK         ) 
!
      USE PRM , ONLY :
     1   IDIM  ,
     1   JDIM  
!
      USE ISLSCP_FILE_UNIT , ONLY :   IFILE_VEG
!
      IMPLICIT NONE
!
! OUTPUT
!
      INTEGER,INTENT(OUT)   :: IMASK            (IDIM,JDIM) 
      REAL(8),INTENT(OUT)   :: GLON  (IDIM, JDIM)
      REAL(8),INTENT(OUT)   :: GLAT  (IDIM, JDIM)
      REAL:: GLON_REAL(IDIM, JDIM)
      REAL:: GLAT_REAL(IDIM, JDIM)
!
! WORK
!
      INTEGER      :: I , J
!
! IMASK ÇÃê›íË
!
      CALL IRESET ( IMASK  , IDIM*JDIM )
      READ ( IFILE_VEG,REC=1 ) IMASK
!
      DO J=1,JDIM
      DO I=1,IDIM
        IF     ( IMASK(I,J) .EQ. 14 ) THEN 
          IMASK(I,J) = 6
        ELSEIF ( IMASK(I,J) .EQ. 15 ) THEN 
          IMASK(I,J) = 7
        ENDIF
      ENDDO
      ENDDO
!
! GLON, GLAT ÇÃê›íË
!
        READ (IFILE_VEG, REC=2) GLON_REAL
        READ (IFILE_VEG, REC=3) GLAT_REAL
        GLON(:, :) = GLON_REAL(:, :)
        GLAT(:, :) = GLAT_REAL(:, :)
!
      RETURN
      END 
