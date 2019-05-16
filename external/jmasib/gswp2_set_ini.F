      SUBROUTINE GSWP2_SET_SIB_INI (
     I   IMASK              , ID_NEXT            , 
     O   TMP_CNP_NOS_PHY    , FR_WTR_CNP_NOS_PHY , FR_GLA_CNP_NOS_PHY , 
     O   TMP_CNP_SNW_PHY    , FR_WTR_CNP_SNW_PHY , FR_GLA_CNP_SNW_PHY , 
     O   TMP_GRSK_PHY       , FR_WTR_GRS_PHY     , FR_GLA_GRS_PHY     , 
     O   TMP_SNSK_PHY       , 
     O   INFO_SNOW_PHY      , 
     O   TMP_SOIL_PHY       , FR_WTR_SOIL_PHY    , FR_GLA_SOIL_PHY    , 
     O   TMP_SNSL_PHY       , FR_WTR_SNSL_PHY    , FR_GLA_SNSL_PHY    , 
     O   TMP_SNOW_PHY       , WTR_SNOW_PHY       , GLA_SNOW_PHY       , 
     O   RHO_SNOW_INV_PHY   , AGE_SNOW_PHY       , 
     O   ENG_SNOW_BUCKET_PHY, H2O_SNOW_BUCKET_PHY ) 
!
      USE PRM , ONLY :
     1   IDIM  , 
     1   JDIM  , 
     1   ISPT   
!
      USE SIBPRM , ONLY :
     1   IDP   ,
     1   ISN 
!
      USE SIBCON , ONLY : TMP_FREZ_c 
      use gswp2, only : gswp2__dataread 
!
      use calendar, only : calendar_run_getkt
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)   :: IMASK              (IDIM,JDIM)
      INTEGER,INTENT(IN)   :: ID_NEXT(5) 
!
! OUTPUT
!
      REAL(8),INTENT(OUT)  :: TMP_CNP_NOS_PHY    (IDIM,JDIM)
      REAL(8),INTENT(OUT)  :: FR_WTR_CNP_NOS_PHY (IDIM,JDIM)    
      REAL(8),INTENT(OUT)  :: FR_GLA_CNP_NOS_PHY (IDIM,JDIM)    
      REAL(8),INTENT(OUT)  :: TMP_CNP_SNW_PHY    (IDIM,JDIM)
      REAL(8),INTENT(OUT)  :: FR_WTR_CNP_SNW_PHY (IDIM,JDIM)    
      REAL(8),INTENT(OUT)  :: FR_GLA_CNP_SNW_PHY (IDIM,JDIM)    
      REAL(8),INTENT(OUT)  :: TMP_GRSK_PHY       (IDIM,JDIM)
      REAL(8),INTENT(OUT)  :: FR_WTR_GRS_PHY     (IDIM,JDIM)   
      REAL(8),INTENT(OUT)  :: FR_GLA_GRS_PHY     (IDIM,JDIM)   
      REAL(8),INTENT(OUT)  :: TMP_SNSK_PHY       (IDIM,JDIM)
!
      REAL(8),INTENT(OUT)  :: TMP_SOIL_PHY       (IDIM,JDIM,IDP)
      REAL(8),INTENT(OUT)  :: FR_WTR_SOIL_PHY    (IDIM,JDIM,IDP)
      REAL(8),INTENT(OUT)  :: FR_GLA_SOIL_PHY    (IDIM,JDIM,IDP)
      REAL(8),INTENT(OUT)  :: TMP_SNSL_PHY       (IDIM,JDIM)    
      REAL(8),INTENT(OUT)  :: FR_WTR_SNSL_PHY    (IDIM,JDIM)    
      REAL(8),INTENT(OUT)  :: FR_GLA_SNSL_PHY    (IDIM,JDIM)    
!
      REAL(8),INTENT(OUT)  :: TMP_SNOW_PHY       (IDIM,JDIM,ISN)
      REAL(8),INTENT(OUT)  :: WTR_SNOW_PHY       (IDIM,JDIM,ISN)
      REAL(8),INTENT(OUT)  :: GLA_SNOW_PHY       (IDIM,JDIM,ISN)
      REAL(8),INTENT(OUT)  :: RHO_SNOW_INV_PHY   (IDIM,JDIM,ISN)    
!
      REAL(8),INTENT(OUT)  :: AGE_SNOW_PHY       (IDIM,JDIM)
      INTEGER,INTENT(OUT)  :: INFO_SNOW_PHY      (IDIM,JDIM) 
!
      REAL(8),INTENT(OUT)  :: ENG_SNOW_BUCKET_PHY(IDIM,JDIM)
      REAL(8),INTENT(OUT)  :: H2O_SNOW_BUCKET_PHY(IDIM,JDIM)
!
      INTEGER :: I 
      INTEGER :: J 
      integer,save :: ifile_tmp = -80
      REAL(4) :: TMP4(IDIM,JDIM) 
      REAL(8) :: TMP 
      REAL(8) :: TMP_F 
      REAL(8) :: TMP_F_plus 
      REAL(8) :: TMP_F_minus
      character(20) :: C20 = 'Tair_cru' 
      integer :: idate0(5) = 0
      integer :: kt , irec 
!
! ======================
! >> ó§ñ èâä˙ílÇÃê›íË <<  
! ======================
!
      idate0(1) = id_next(1) 
      idate0(2) = 1
      idate0(3) = 1
      call calendar_run_getkt ( idate0 , id_next , 4 , kt  ) 
      irec = kt / 3 + 1 
      write(6,*) 'gswp2_set_ini : id_next= ' , id_next(:) 
      write(6,*) 'gswp2_set_ini : irec=' , irec 
      call gswp2__dataread ( ifile_tmp, C20, tmp4, id_next(1), irec) 
      close (ifile_tmp)
!
      TMP_F = TMP_FREZ_c 
      TMP_F_plus = TMP_F + 1 
      TMP_F_minus= TMP_F - 1 
      DO J=1,JDIM
      DO I=1,IDIM
      if ( imask(i,j)>= 0 ) then
        TMP   = TMP4(I,J) 
!
        TMP_CNP_NOS_PHY (I,J)  = TMP
!
        IF ( TMP .GT. TMP_F ) THEN 
          TMP_GRSK_PHY    (I,J)   = TMP + ( TMP_F_plus - TMP ) * 0.1 
          TMP_SOIL_PHY    (I,J,1) = TMP + ( TMP_F_plus - TMP ) * 0.3 
          TMP_SOIL_PHY    (I,J,2) = TMP + ( TMP_F_plus - TMP ) * 0.6 
          TMP_SOIL_PHY    (I,J,3) = TMP + ( TMP_F_plus - TMP ) * 0.9 
!
          FR_WTR_SOIL_PHY (I,J,1) = 0.5
          FR_WTR_SOIL_PHY (I,J,2) = 0.5
          FR_WTR_SOIL_PHY (I,J,3) = 0.5
          FR_GLA_SOIL_PHY (I,J,1) = 0.
          FR_GLA_SOIL_PHY (I,J,2) = 0.
          FR_GLA_SOIL_PHY (I,J,3) = 0.
        ELSE 
          TMP_GRSK_PHY    (I,J)   = TMP + ( TMP_F_minus - TMP ) * 0.1 
          TMP_SOIL_PHY    (I,J,1) = TMP + ( TMP_F_minus - TMP ) * 0.3 
          TMP_SOIL_PHY    (I,J,2) = TMP + ( TMP_F_minus - TMP ) * 0.6 
          TMP_SOIL_PHY    (I,J,3) = TMP + ( TMP_F_minus - TMP ) * 0.9 
!
          FR_WTR_SOIL_PHY (I,J,1) = 0.
          FR_WTR_SOIL_PHY (I,J,2) = 0.
          FR_WTR_SOIL_PHY (I,J,3) = 0.
          FR_GLA_SOIL_PHY (I,J,1) = 0.5
          FR_GLA_SOIL_PHY (I,J,2) = 0.5
          FR_GLA_SOIL_PHY (I,J,3) = 0.5
        ENDIF
!
        FR_WTR_CNP_NOS_PHY(I,J) = 0.
        FR_GLA_CNP_NOS_PHY(I,J) = 0. 
        FR_WTR_GRS_PHY    (I,J) = 0.
        FR_GLA_GRS_PHY    (I,J) = 0.
!
        TMP_CNP_SNW_PHY   (I,J) = 0.
        FR_WTR_CNP_SNW_PHY(I,J) = 0.
        FR_GLA_CNP_SNW_PHY(I,J) = 0.
!
        TMP_SNSK_PHY      (I,J) = 0. 
        TMP_SNSL_PHY      (I,J) = 0. 
        FR_WTR_SNSL_PHY   (I,J) = 0. 
        FR_GLA_SNSL_PHY   (I,J) = 0.
!
        AGE_SNOW_PHY      (I,J)   = 10*24*3600 
        INFO_SNOW_PHY     (I,J)   = -1 
!
        TMP_SNOW_PHY    (I,J,1) = 0.
        TMP_SNOW_PHY    (I,J,2) = 0.
        TMP_SNOW_PHY    (I,J,3) = 0.
        TMP_SNOW_PHY    (I,J,4) = 0.
        WTR_SNOW_PHY    (I,J,1) = 0.
        WTR_SNOW_PHY    (I,J,2) = 0.
        WTR_SNOW_PHY    (I,J,3) = 0.
        WTR_SNOW_PHY    (I,J,4) = 0.
        GLA_SNOW_PHY    (I,J,1) = 0.
        GLA_SNOW_PHY    (I,J,2) = 0.
        GLA_SNOW_PHY    (I,J,3) = 0.
        GLA_SNOW_PHY    (I,J,4) = 0.
        RHO_SNOW_INV_PHY(I,J,1) = 0.
        RHO_SNOW_INV_PHY(I,J,2) = 0.
        RHO_SNOW_INV_PHY(I,J,3) = 0.
        RHO_SNOW_INV_PHY(I,J,4) = 0.
!
        ENG_SNOW_BUCKET_PHY(I,J) = 0.  
        H2O_SNOW_BUCKET_PHY(I,J) = 0.  
!
      endif
      ENDDO
      ENDDO
!
      RETURN
      END
