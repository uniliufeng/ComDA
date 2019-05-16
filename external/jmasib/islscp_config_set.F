      SUBROUTINE ISLSCP_CONFIG_SET (
     O    GLON , GLAT , IMASK         ) 
!
      USE PRM , ONLY :
     1   IDIM  ,
     1   JDIM  
!
      USE ISLSCP_FILE_UNIT , ONLY :   IFILE_VEG
      USE ISLSCP , ONLY : JCN_ISLSCP_GRADS 
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
      CHARACTER(4) :: CMARK
      INTEGER      :: IY , IM , ID , IH 
      INTEGER      :: I , J
!
      IF (       JCN_ISLSCP_GRADS.EQ.1 
     1     .AND. (IDIM.NE.360 .OR. JDIM.NE.180) ) THEN
        write(6,*) 'ISLSCP_CONFIG_SET ERROR' , 
     1             JCN_ISLSCP_GRADS , IDIM , JDIM 
        stop 999
      ENDIF
!
! IMASK ÇÃê›íË
!
      CALL IRESET ( IMASK  , IDIM*JDIM )
!
      IF ( JCN_ISLSCP_GRADS.EQ.0 ) THEN
        READ ( IFILE_VEG,REC=1 ) CMARK , IY , IM , ID , IH , IMASK       
        IF ( CMARK.NE.'VEG' .OR. IY.NE.1987 .OR. IM.NE.1 ) THEN
          write(6,*) 'islscp_config_set error ' , 
     1               cmark , iy , im , id , ih
          STOP 999
        ENDIF
      ELSE
        READ ( IFILE_VEG,REC=1 ) IMASK
      ENDIF
!
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
      IF( JCN_ISLSCP_GRADS .EQ.0 ) THEN 
        READ (IFILE_VEG, REC=2) CMARK , IY , IM , ID , IH , GLON_REAL
        IF ( CMARK.NE.'LON' ) THEN
          write(6,*) 'islscp_config_set error in reading GLON: ' , 
     1               cmark , iy , im , id , ih
          STOP 999
        ENDIF
!
        GLON(:, :) = GLON_REAL(:, :)
!
        READ (IFILE_VEG, REC=3) CMARK , IY , IM , ID , IH , GLAT_REAL
        IF ( CMARK.NE.'LAT' ) THEN
          write(6,*) 'islscp_config_set error in reading GLAT: ' , 
     1               cmark , iy , im , id , ih
          STOP 999
        ENDIF
        GLAT(:, :) = GLAT_REAL(:, :)
!
      ELSE
          DO J=1,JDIM
        DO I=1,IDIM
          GLON(I,J) = I - 180.5 
          GLAT(I,J) = 90.5 - J
        ENDDO
        ENDDO
      ENDIF
!
      RETURN
      END 
