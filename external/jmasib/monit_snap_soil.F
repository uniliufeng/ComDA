      SUBROUTINE MONIT_SNAP_SOIL  (
     I    NMAX             , NTYP             , 
     I    IMEAN            , INFO_MEAN        ,        
     I    DELT_SOIL        , DELT_CNP         ,
     I    CVR_SNW_NEXT     , 
     I    TMP_SOIL_NEXT    , TMP_SNSL_NEXT    , TMP_SOIL4 , 
     I    FR_WTR_SOIL_NEXT , FR_WTR_SNSL_NEXT , 
     I    FR_GLA_SOIL_NEXT , FR_GLA_SNSL_NEXT     
     1   )  
!
      USE SIBPRM , ONLY : 
     1    IDX , 
     1    IDP  
!
      USE SIBCON , ONLY :
     1    CAP_WTR_c  ,
     1    CAP_GLA_c  , 
     1    CAP_SOIL_BY_VOL_c , 
     1    TMP_FREZ_c , 
     1    HL_W2G_c   , 
     1    RHO_WTR_c  , 
     1    RHO_GLA_c  , 
     1    CAP_SOIL4_c  
!
!
      USE SIBVEG , ONLY :
     1           POR_SOIL_co            ,  ! �󌄗�   (ITYP)  
     1           DEP_SOIL_co            ,  ! �y��w�̐[��    (ITYP,IDP)
     1           POR_SOIL_INV_co        ,  ! �󌄗�   (ITYP)  
     1           DEP_SOIL_INV_co           ! �y��w�̐[��    (ITYP,IDP)
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB_NEXT , INA, ICN, ICS,
     1                       IMONIT_LEVEL , IMONIT_ALL ,
     1                       IMONIT_MEAN_MRI ,
     1                       IMONIT_SNAP_NPD ,
     1                       IMONIT_ALMA 
!
      USE COM_STEP_SIB0109 , ONLY : ICN_SIB0109_MONITOR
#endif
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)    :: IMEAN        ! ���ϒl���o�� 1 <-- �s�v??
      INTEGER,INTENT(IN)    :: INFO_MEAN    ! �� -> delt_mean=DELT_CNP
      INTEGER,INTENT(IN)    :: NMAX       
      INTEGER,INTENT(IN)    :: NTYP(IDX)
!
!     INTEGER,INTENT(IN)    :: NMAX_01(3)
      REAL(8),INTENT(IN)    :: CVR_SNW_NEXT    (IDX) 
      REAL(8),INTENT(IN)    :: DELT_SOIL
      REAL(8),INTENT(IN)    :: DELT_CNP 
!
      REAL(8),INTENT(IN)    :: TMP_SOIL4            (IDX) 
      REAL(8),INTENT(IN)    :: TMP_SOIL_NEXT        (IDX,IDP) 
      REAL(8),INTENT(IN)    :: TMP_SNSL_NEXT        (IDX,IDP) 
      REAL(8),INTENT(IN)    :: FR_WTR_SOIL_NEXT     (IDX,IDP) 
      REAL(8),INTENT(IN)    :: FR_WTR_SNSL_NEXT     (IDX,IDP) 
      REAL(8),INTENT(IN)    :: FR_GLA_SOIL_NEXT     (IDX,IDP) 
      REAL(8),INTENT(IN)    :: FR_GLA_SNSL_NEXT     (IDX,IDP) 
!
! WORK
!     
#ifdef MONYOS
      REAL(8),PARAMETER :: ONE  = 1. 
      REAL(8),PARAMETER :: D100 = 100.D0 
      CHARACTER(7) :: CMARK 
      REAL(8)      :: WORK  (IDX)
      REAL(8)      :: WORK_N(IDX)
      REAL(8)      :: WORK_S(IDX)
      REAL(8)      :: WTR_N(IDX)
      REAL(8)      :: WTR_S(IDX)
      REAL(8)      :: GLA_N(IDX)
      REAL(8)      :: GLA_S(IDX)
      REAL(8)      :: DELT_MEAN
      REAL(8)      :: DELT_SNAP
#endif
!
      INTEGER   ::  I
      INTEGER   ::  K
!
      LOGICAL,SAVE :: LFIRST
!
      DATA  LFIRST /.TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_SNAP_SOIL 99/10 POCHAKA'

        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
      IF ( INFO_MEAN .LT. 0 ) THEN
        DELT_MEAN = DELT_CNP 
      ELSE
        DELT_MEAN = DELT_SOIL
      ENDIF
      DELT_SNAP   = DELT_CNP
!
! ==========================
! >> ���ԕ��ρE�i�q���ϒl <<
! ==========================
!
      IF ( IMONIT_LEVEL .GE. IMONIT_MEAN_MRI ) THEN
      IF ( IMEAN .EQ. 1 ) THEN
!
      DO K=1,IDP
!
! �y�뉷�x
!
        DO I=1,NMAX
          WORK(I) =   TMP_SOIL_NEXT(I,K) * (1-CVR_SNW_NEXT(I))
     1              + TMP_SNSL_NEXT(I,K) * CVR_SNW_NEXT(I)
        ENDDO
        CMARK = 'TMPLxM' 
        WRITE(CMARK (5:5)  ,'(I1)') K 
        CALL MONIT_SIB_NEXT (CMARK,WORK,INA,DELT_MEAN,ONE)      
!
! �t���� 
!
        DO I=1,NMAX
          WORK(I) =   FR_WTR_SOIL_NEXT(I,K) * (1-CVR_SNW_NEXT(I))
     1              + FR_WTR_SNSL_NEXT(I,K) * CVR_SNW_NEXT(I)
        ENDDO
        CMARK = 'FRWLxM' 
        WRITE(CMARK (5:5)  ,'(I1)') K 
        CALL MONIT_SIB_NEXT (CMARK,WORK,INA,DELT_MEAN,ONE)      
!
! �X��
!
        DO I=1,NMAX
          WORK(I) =   FR_GLA_SOIL_NEXT(I,K) * (1-CVR_SNW_NEXT(I))
     1              + FR_GLA_SNSL_NEXT(I,K) * CVR_SNW_NEXT(I)
        ENDDO
        CMARK = 'FRGLxM' 
        WRITE(CMARK (5:5)  ,'(I1)') K 
        CALL MONIT_SIB_NEXT (CMARK,WORK,INA,DELT_MEAN,ONE)      
!
      ENDDO
!
        CMARK = 'TMPL4M' 
        CALL MONIT_SIB_NEXT (CMARK,TMP_SOIL4,INA,DELT_MEAN,ONE)      
!
      ENDIF
      ENDIF
!
! * for ALMA *********
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
!
!  A_TSLx  
!
        DO K=1,IDP
          DO I=1,NMAX
            WORK(I) =   TMP_SOIL_NEXT(I,K) * (1-CVR_SNW_NEXT(I))
     1                + TMP_SNSL_NEXT(I,K) * CVR_SNW_NEXT(I)
          ENDDO
          CMARK = 'A_TSLx' 
          WRITE(CMARK (6:6)  ,'(I1)') K 
          CALL MONIT_SIB_NEXT (CMARK,WORK,INA,DELT_MEAN,ONE)      
        ENDDO
!
!  A_MSLx  
!
        DO K=1,IDP
          WORK(:) = 0. 
          DO I=1,NMAX
          WORK(I) = WORK(I) 
     1             + (   (FR_WTR_SOIL_NEXT(I,K)+FR_GLA_SOIL_NEXT(I,K)) 
     1                    * (1-CVR_SNW_NEXT(I))
     1                 + (FR_WTR_SNSL_NEXT(I,K)+FR_GLA_SNSL_NEXT(I,K)) 
     1                    * CVR_SNW_NEXT(I) 
     1               ) * POR_SOIL_co(NTYP(I)) * DEP_SOIL_co(NTYP(I),K)
     1                 * RHO_WTR_c
          ENDDO
          CMARK = 'A_MSLx' 
          WRITE(CMARK (6:6)  ,'(I1)') K 
          CALL MONIT_SIB_NEXT (CMARK,WORK,INA,DELT_MEAN,ONE)      
        ENDDO
!
!
!  A_WETSL  
!
        WORK_S(:) = 0. 
        WORK_N(:) = 0. 
        DO K = 1,IDP
        DO I = 1,NMAX
          WORK_S(I) = WORK_S(I) 
     1             + (   (FR_WTR_SOIL_NEXT(I,K)+FR_GLA_SOIL_NEXT(I,K)) 
     1                    * (1-CVR_SNW_NEXT(I))
     1                 + (FR_WTR_SNSL_NEXT(I,K)+FR_GLA_SNSL_NEXT(I,K)) 
     1                    * CVR_SNW_NEXT(I) 
     1               ) * POR_SOIL_co(NTYP(I)) * DEP_SOIL_co(NTYP(I),K)
          WORK_N(I) = WORK_N(I) 
     1                 + POR_SOIL_co(NTYP(I)) * DEP_SOIL_co(NTYP(I),K)
        ENDDO
        ENDDO
        DO I=1,NMAX
          WORK(I) = WORK_S(I) / WORK_N(I) 
        ENDDO
        CMARK = 'A_WETSLA' 
        CALL MONIT_SIB_NEXT (CMARK,WORK,INA,DELT_MEAN,ONE)      
!
!  A_MSLAX  ( �u�Ԓl ) 
!
        WORK(:) = 0. 
        DO K=1,IDP
          DO I=1,NMAX
          WORK(I) = WORK(I) 
     1             + (   (FR_WTR_SOIL_NEXT(I,K)+FR_GLA_SOIL_NEXT(I,K)) 
     1                    * (1-CVR_SNW_NEXT(I))
     1                 + (FR_WTR_SNSL_NEXT(I,K)+FR_GLA_SNSL_NEXT(I,K)) 
     1                    * CVR_SNW_NEXT(I) 
     1               ) * POR_SOIL_co(NTYP(I)) * DEP_SOIL_co(NTYP(I),K)
     1                 * RHO_WTR_c
          ENDDO
        ENDDO
        CMARK = 'A_MSLAX' 
        CALL MONIT_SIB_NEXT (CMARK,WORK,INA,DELT_SNAP,ONE)      
!
      ENDIF

!
! ======================
! >> �u�ԁE�i�q���ϒl <<
! ======================
!
      IF (       IMONIT_LEVEL .GE. IMONIT_SNAP_NPD
     1     .AND. ICN_SIB0109_MONITOR.EQ.1
     1    ) THEN
!
      DO K=1,IDP
!
! �y�뉷�x
!
        DO I=1,NMAX
          WORK(I) =   TMP_SOIL_NEXT(I,K) * (1-CVR_SNW_NEXT(I))
     1              + TMP_SNSL_NEXT(I,K) * CVR_SNW_NEXT(I)
        ENDDO
        CMARK = 'SoilTx' 
        WRITE(CMARK (6:6)  ,'(I1)') K 
        CALL MONIT_SIB_NEXT (CMARK,WORK,INA,DELT_SNAP,ONE)      
!
! �t���� 
!
        DO I=1,NMAX
          WORK(I) =   FR_WTR_SOIL_NEXT(I,K) * (1-CVR_SNW_NEXT(I))
     1              + FR_WTR_SNSL_NEXT(I,K) * CVR_SNW_NEXT(I)
        ENDDO
        CMARK = 'SoilWx' 
        WRITE(CMARK (6:6)  ,'(I1)') K 
        CALL MONIT_SIB_NEXT (CMARK,WORK,INA,DELT_SNAP,D100)      
!
! �X��
!
        DO I=1,NMAX
          WORK(I) =   FR_GLA_SOIL_NEXT(I,K) * (1-CVR_SNW_NEXT(I))
     1              + FR_GLA_SNSL_NEXT(I,K) * CVR_SNW_NEXT(I)
        ENDDO
        CMARK = 'SoilIx' 
        WRITE(CMARK (6:6)  ,'(I1)') K 
        CALL MONIT_SIB_NEXT (CMARK,WORK,INA,DELT_SNAP,D100)      
!
      ENDDO
      ENDIF
!
! ==========================
! >> �y��A���x�`�F�b�N�p�� <<
! ==========================
!
      IF (  IMONIT_LEVEL .GE. IMONIT_ALL  ) THEN
!
! -------------------------
! > �y��A���x�`�F�b�N�p <
! -------------------------
!
      DO K=1,IDP
!
!  �t����
!
      DO I = 1,NMAX
        WTR_N(I) = FR_WTR_SOIL_NEXT(I,K) * RHO_WTR_c 
     1                * POR_SOIL_co(NTYP(I)) * DEP_SOIL_co(NTYP(I),K)
      ENDDO
!
      DO I = 1,NMAX
        WTR_S(I) = FR_WTR_SNSL_NEXT(I,K) * RHO_WTR_c 
     1                * POR_SOIL_co(NTYP(I)) * DEP_SOIL_co(NTYP(I),K)
      ENDDO
!
!  �X��
!
      DO I = 1,NMAX
        GLA_N(I) = FR_GLA_SOIL_NEXT(I,K) * RHO_WTR_c 
     1                * POR_SOIL_co(NTYP(I)) * DEP_SOIL_co(NTYP(I),K)       
      ENDDO
!
      DO I = 1,NMAX
        GLA_S(I) = FR_GLA_SNSL_NEXT(I,K) * RHO_WTR_c 
     1                * POR_SOIL_co(NTYP(I)) * DEP_SOIL_co(NTYP(I),K)
      ENDDO
!
!  �G�l���M�[
!
      DO I = 1,NMAX
        WORK_N(I) 
     1       = (
     1             CAP_SOIL_BY_VOL_c * DEP_SOIL_co(NTYP(I),K) 
     1                        * ( 1 - POR_SOIL_co(NTYP(I)) ) 
     1           + CAP_WTR_c  * WTR_N(I)
     1           + CAP_GLA_c  * GLA_N(I) )
     1         * ( TMP_SOIL_NEXT(I,K) - TMP_FREZ_c ) 
     1         - HL_W2G_c * GLA_N(I) 
      ENDDO
!
      DO I = 1,NMAX
        WORK_S(I) 
     1       = (
     1             CAP_SOIL_BY_VOL_c * DEP_SOIL_co(NTYP(I),K) 
     1                        * ( 1 - POR_SOIL_co(NTYP(I)) ) 
     1           + CAP_WTR_c  * WTR_S(I)
     1           + CAP_GLA_c  * GLA_S(I) )
     1         * ( TMP_SNSL_NEXT(I,K) - TMP_FREZ_c ) 
     1         - HL_W2G_c * GLA_S(I) 
      ENDDO
! 
!  �u�Ԓl
!
        IF ( ICN_SIB0109_MONITOR.EQ.1 ) THEN
          CMARK = 'WTRLNxX' 
          WRITE(CMARK (6:6)  ,'(I1)') K 
          CALL MONIT_SIB_NEXT (CMARK,WTR_N,ICN,DELT_SNAP,ONE)      
!
          CMARK = 'WTRLSxX' 
          WRITE(CMARK (6:6)  ,'(I1)') K 
          CALL MONIT_SIB_NEXT (CMARK,WTR_S,ICS,DELT_SNAP,ONE)      

          CMARK = 'GLALNxX' 
          WRITE(CMARK (6:6)  ,'(I1)') K 
          CALL MONIT_SIB_NEXT (CMARK,GLA_N,ICN,DELT_SNAP,ONE)      
!
          CMARK = 'GLALSxX' 
          WRITE(CMARK (6:6)  ,'(I1)') K 
          CALL MONIT_SIB_NEXT (CMARK,GLA_S,ICS,DELT_SNAP,ONE)      
!
          CMARK = 'ENGLNxX' 
          WRITE(CMARK (6:6)  ,'(I1)') K 
          CALL MONIT_SIB_NEXT (CMARK,WORK_N,ICN,DELT_SNAP,ONE)      
!
          CMARK = 'ENGLSxX' 
          WRITE(CMARK (6:6)  ,'(I1)') K 
          CALL MONIT_SIB_NEXT (CMARK,WORK_S,ICS,DELT_SNAP,ONE)      
        ENDIF
!
!  ���̑��u�Ԓl
!
        IF ( ICN_SIB0109_MONITOR.EQ.1 ) THEN
          CMARK = 'TMPLNxX' 
          WRITE(CMARK (6:6)  ,'(I1)') K 
          CALL MONIT_SIB_NEXT (CMARK,TMP_SOIL_NEXT(1,K),
     1                               ICS,DELT_SNAP,ONE) 
!
          CMARK = 'TMPLSxX' 
          WRITE(CMARK (6:6)  ,'(I1)') K 
          CALL MONIT_SIB_NEXT (CMARK,TMP_SNSL_NEXT(1,K),
     1                               ICS,DELT_SNAP,ONE) 
        ENDIF  
!
!  ���ԕ��ϒl
!
        IF ( IMEAN .EQ. 1 ) THEN
!
!    �t����
!
          CMARK = 'WTRLNxM' 
          WRITE(CMARK (6:6)  ,'(I1)') K 
          CALL MONIT_SIB_NEXT (CMARK,WTR_N,ICN,DELT_MEAN,ONE)      
          CMARK = 'WTRLSxM' 
          WRITE(CMARK (6:6)  ,'(I1)') K 
          CALL MONIT_SIB_NEXT (CMARK,WTR_S,ICS,DELT_MEAN,ONE)      
!
!  �X��
!
          CMARK = 'GLALNxM' 
          WRITE(CMARK (6:6)  ,'(I1)') K 
          CALL MONIT_SIB_NEXT (CMARK,GLA_N,ICN,DELT_MEAN,ONE)      
          CMARK = 'GLALSxM' 
          WRITE(CMARK (6:6)  ,'(I1)') K 
          CALL MONIT_SIB_NEXT (CMARK,GLA_S,ICS,DELT_MEAN,ONE)      
!
!  �G�l���M�[
!
          CMARK = 'ENGLNxM' 
          WRITE(CMARK (6:6)  ,'(I1)') K 
          CALL MONIT_SIB_NEXT (CMARK,WORK_N,ICN,DELT_MEAN,ONE)      
!
          CMARK = 'ENGLSxM' 
          WRITE(CMARK (6:6)  ,'(I1)') K 
          CALL MONIT_SIB_NEXT (CMARK,WORK_S,ICS,DELT_MEAN,ONE)      
        ENDIF
      ENDDO
!
      ENDIF
!
! ================
! >> �y���l�w <<
! ================
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
!   ���x�p�G�l���M�[�X�i�b�v
!
        IF ( ICN_SIB0109_MONITOR.EQ.1 ) THEN
         WORK(:) = TMP_SOIL4(:) - TMP_FREZ_c 
!
         CMARK = 'ENGLN4X' 
         CALL MONIT_SIB_NEXT(CMARK,WORK,ICN,DELT_SNAP,CAP_SOIL4_c)     
!       
         CMARK = 'ENGLS4X' 
         CALL MONIT_SIB_NEXT(CMARK,WORK,ICS,DELT_SNAP,CAP_SOIL4_c)      
!
         CMARK = 'ENGL4X'
         CALL MONIT_SIB_NEXT(CMARK,WORK,INA,DELT_SNAP,CAP_SOIL4_c)      
!
         CMARK = 'TMPL4X'
         CALL MONIT_SIB_NEXT(CMARK,TMP_SOIL4,INA,DELT_SNAP,ONE)      
        ENDIF        
      ENDIF
!
#if defined(MONIT_JMA)
      IF ( IMONIT_LEVEL >= IMONIT_SNAP_NPD .AND. 
     &     IMONIT_LEVEL <  IMONIT_ALL ) THEN
         CMARK = 'TMPL4X'
         CALL MONIT_SIB_NEXT(CMARK,TMP_SOIL4,INA,DELT_SNAP,ONE)      
      ENDIF        
#endif
!
! ==================
! >> �y�낻�̑��p <<
! ==================
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
      IF ( IMEAN .EQ. 1 ) THEN
!
      DO K=1,IDP
!
!  ���x
!
        CMARK = 'TMPLNxM' 
        WRITE(CMARK (6:6)  ,'(I1)') K 
        CALL MONIT_SIB_NEXT (CMARK,TMP_SOIL_NEXT(1,K),ICN,DELT_MEAN,ONE)      
        CMARK = 'TMPLSxM' 
        WRITE(CMARK (6:6)  ,'(I1)') K 
        CALL MONIT_SIB_NEXT (CMARK,TMP_SNSL_NEXT(1,K),ICS,DELT_MEAN,ONE)      
!
!  ������
!
        CMARK = 'FRWLNxM' 
        WRITE(CMARK (6:6)  ,'(I1)') K 
        CALL MONIT_SIB_NEXT (CMARK,FR_WTR_SOIL_NEXT(1,K),ICN,
     1                       DELT_MEAN,ONE)      
        CMARK = 'FRWLSxM' 
        WRITE(CMARK (6:6)  ,'(I1)') K 
        CALL MONIT_SIB_NEXT (CMARK,FR_WTR_SNSL_NEXT(1,K),ICS,
     1                       DELT_MEAN,ONE)      
!
!  �X����
!
        CMARK = 'FRGLNxM' 
        WRITE(CMARK (6:6)  ,'(I1)') K 
        CALL MONIT_SIB_NEXT (CMARK,FR_GLA_SOIL_NEXT(1,K),ICN,
     1                       DELT_MEAN,ONE)      
        CMARK = 'FRGLSxM' 
        WRITE(CMARK (6:6)  ,'(I1)') K 
        CALL MONIT_SIB_NEXT (CMARK,FR_GLA_SNSL_NEXT(1,K),ICS,
     1                       DELT_MEAN,ONE)      
!
!
      ENDDO
!
      ENDIF
      ENDIF
#endif
!
      RETURN
      END
!=========================================================
      SUBROUTINE MR_MONIT_SNAP_SOIL
!
      USE SIBPRM , ONLY :
     1    ISN , 
     1    IDP 
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : 
     1                       IMONIT_LEVEL , IMONIT_ALL ,
     1                       IMONIT_MEAN_MRI ,
     1                       IMONIT_SNAP_NPD , IMONIT_ALMA 
#endif

!
      IMPLICIT NONE
!
      CHARACTER( 7)  :: CMARK 
      CHARACTER(32)  :: CTITLE 
      CHARACTER(13)  :: CUNIT
!
      INTEGER        :: K
!
      LOGICAL,SAVE ::  LFIRST 
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_REGIST_MONIT_SNAP_SOIL 99/10 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
!
! ==========================
! >> ���ԕ��ρE�i�q���ϗ� <<
! ==========================
!
      IF ( IMONIT_LEVEL .GE. IMONIT_MEAN_MRI ) THEN
!
        DO K=1,IDP
!
!     �y�뉷�x
!
        CMARK  = 'TMPLxM' 
        CTITLE = 'TEMP OF SOIL  x LEV AREA MEAN (TMEAN)'
        WRITE(CMARK  ( 5: 5)  ,'(I1)') K 
        WRITE(CTITLE (15:15)  ,'(I1)') K 
        CUNIT  = 'K'
        CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
!     �y�됅����
!
        CMARK  = 'FRWLxM' 
        CTITLE = 'FRACTION OF WATER OF SOIL x LEV AREA MEAN (TMEAN)'
        WRITE(CMARK  ( 5: 5)  ,'(I1)') K 
        WRITE(CTITLE (27:27)  ,'(I1)') K 
        CUNIT  = '0-1'
        CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
!     �y��X��
!
        CMARK  = 'FRGLxM' 
        CTITLE = 'FRACTION OF GLA OF SOIL   x LEV AREA MEAN (TMEAN)'
        WRITE(CMARK  ( 5: 5)  ,'(I1)') K 
        WRITE(CTITLE (27:27)  ,'(I1)') K 
        CUNIT  = '0-1'
        CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
        ENDDO
      ENDIF
!
! for ALMA
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
        DO K=1,IDP
          CMARK  = 'A_TSLx' 
          WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
          CTITLE = 'ALMA SoilTemp (TMEAN)'
          CUNIT  = 'K'
          CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
        ENDDO
        DO K=1,IDP
          CMARK  = 'A_MSLx' 
          WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
          CTITLE = 'ALMA SoilMoist (TMEAN)'
          CUNIT  = 'kg/m**2'
          CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
        ENDDO
        CMARK  = 'A_WETSLA' 
        CTITLE = 'ALMA SoilWet (TMEAN)'
        CUNIT  = '0-1'
        CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'A_MSLAX' 
        CTITLE = 'ALMA Dt DelSoilMoist (SoilMoist Tsnap)'
        CUNIT  = 'kg/m**2'
        CALL MONIT_REGIST_SIB_snap ( CMARK , CTITLE , CUNIT )

      ENDIF
!
! ==========================
! >> �X�i�b�v�E�i�q���ϗ� <<
! ==========================
!
      IF ( IMONIT_LEVEL .GE. IMONIT_SNAP_NPD ) THEN 
!
        DO K=1,IDP
!
!   �y�뉷�x
!
          CMARK  = 'SoilTx'
          CTITLE = 'Temp Soil x'
          WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
          WRITE(CTITLE (11:11)  ,'(I1)') K 
          CUNIT  = 'K'
          CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
!   �y�됅��
!
          CMARK  = 'SoilWx'
          CTITLE = 'Fraction of Water in Soil x'
          WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
          WRITE(CTITLE (27:27)  ,'(I1)') K 
          CUNIT  = '%'
          CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
!   �y�됅��
!
          CMARK  = 'SoilIx'
          CTITLE = 'Fraction of Ice   in Soil x'
          WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
          WRITE(CTITLE (27:27)  ,'(I1)') K 
          CUNIT  = '%'
          CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
        ENDDO
      ENDIF
!
! ==================
! >> �̈�ʃ��j�^ <<
! ==================
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN 
!
      DO K=1,IDP
!
! ================================
! >> �y����w�A���x�`�F�b�N�p <<
! ================================
!
!  �t����
!
      CMARK  = 'WTRLNxX' 
      CTITLE = 'WATER OF SOIL x LEV NOS (SNAP)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (15:15)  ,'(I1)') K 
      CUNIT  = 'KG/M**2'
      CALL MONIT_REGIST_SIB_SNAP ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WTRLSxX' 
      CTITLE = 'WATER OF SOIL x LEV SNW (SNAP)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (15:15)  ,'(I1)') K 
      CUNIT  = 'KG/M**2'
      CALL MONIT_REGIST_SIB_SNAP ( CMARK , CTITLE , CUNIT )
!
!  �X��
!
      CMARK  = 'GLALNxX' 
      CTITLE = 'GLA OF SOIL   x LEV NOS (SNAP)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (15:15)  ,'(I1)') K 
      CUNIT  = 'KG/M**2'
      CALL MONIT_REGIST_SIB_SNAP ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'GLALSxX' 
      CTITLE = 'GLA OF SOIL   x LEV SNW (SNAP)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (15:15)  ,'(I1)') K 
      CUNIT  = 'KG/M**2'
      CALL MONIT_REGIST_SIB_SNAP ( CMARK , CTITLE , CUNIT )
!
!  �G�l���M�[
!
      CMARK  = 'ENGLNxX' 
      CTITLE = 'ENERGY OF SOIL x LEV NOS (SNAP)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (15:15)  ,'(I1)') K 
      CUNIT  = 'J/M**2'
      CALL MONIT_REGIST_SIB_SNAP ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'ENGLSxX' 
      CTITLE = 'ENERGY OF SOIL x LEV SNW (SNAP)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (15:15)  ,'(I1)') K 
      CUNIT  = 'J/M**2'
      CALL MONIT_REGIST_SIB_SNAP ( CMARK , CTITLE , CUNIT )
!
! ========================
! >> �y����w���̑��p <<
! ========================
!
!  ���x
!
      CMARK  = 'TMPLNxM' 
      CTITLE = 'TEMP OF SOIL   x LEV NOS (TMEAN)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (15:15)  ,'(I1)') K 
      CUNIT  = 'K'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'TMPLSxM' 
      CTITLE = 'TEMP OF SOIL   x LEV SNW (TMEAN)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (15:15)  ,'(I1)') K 
      CUNIT  = 'K'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
!  ���x
!
      CMARK  = 'TMPLNxX' 
      CTITLE = 'TEMP OF SOIL   x LEV NOS (SNAP)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (15:15)  ,'(I1)') K 
      CUNIT  = 'K'
      CALL MONIT_REGIST_SIB_SNAP ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'TMPLSxX' 
      CTITLE = 'TEMP OF SOIL   x LEV SNW (SNAP)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (15:15)  ,'(I1)') K 
      CUNIT  = 'K'
      CALL MONIT_REGIST_SIB_SNAP ( CMARK , CTITLE , CUNIT )
!
!  �t����
!
      CMARK  = 'WTRLNxM' 
      CTITLE = 'WATER OF SOIL  x LEV NOS (TMEAN)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (15:15)  ,'(I1)') K 
      CUNIT  = 'KG/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WTRLSxM' 
      CTITLE = 'WATER OF SOIL  x LEV SNW (TMEAN)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (15:15)  ,'(I1)') K 
      CUNIT  = 'KG/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
!  ������
!
      CMARK  = 'FRWLNxM' 
      CTITLE = 'FRACTION OF WATER OF SOIL x LEV NOS (TMEAN)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (27:27)  ,'(I1)') K 
      CUNIT  = '0-1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'FRWLSxM' 
      CTITLE = 'FRACTION OF WATER OF SOIL x LEV SNW (TMEAN)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (27:27)  ,'(I1)') K 
      CUNIT  = '0-1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
!
!  �X��
!
      CMARK  = 'GLALNxM' 
      CTITLE = 'GLA OF SOIL    x LEV NOS (TMEAN)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (15:15)  ,'(I1)') K 
      CUNIT  = 'KG/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'GLALSxM' 
      CTITLE = 'GLA OF SOIL    x LEV SNW (TMEAN)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (15:15)  ,'(I1)') K 
      CUNIT  = 'KG/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
!  �X����
!
      CMARK  = 'FRGLNxM' 
      CTITLE = 'FRACTION OF GLA OF SOIL   x LEV NOS (TMEAN)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (27:27)  ,'(I1)') K 
      CUNIT  = ' '
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'FRGLSxM' 
      CTITLE = 'FRACTION OF GLA OF SOIL   x LEV SNW (TMEAN)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (27:27)  ,'(I1)') K 
      CUNIT  = ' '
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
!
!  �G�l���M�[
!
      CMARK  = 'ENGLNxM' 
      CTITLE = 'ENERGY OF SOIL x LEV NOS (TMEAN)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (15:15)  ,'(I1)') K 
      CUNIT  = 'J/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'ENGLSxM' 
      CTITLE = 'ENERGY OF SOIL x LEV SNW (TMEAN)'
      WRITE(CMARK  ( 6: 6)  ,'(I1)') K 
      WRITE(CTITLE (15:15)  ,'(I1)') K 
      CUNIT  = 'J/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      ENDDO
!
      ENDIF
!
! --------------
! > �y���l�w <
! --------------
!
      IF ( IMONIT_LEVEL .GE. IMONIT_MEAN_MRI ) THEN
        CMARK ='TMPL4M'
        CTITLE='TEMP SOIL 4'
        CUNIT ='K'
        CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
        CMARK ='ENGLN4X'
        CTITLE='ENG SOIL 4 NOS (snap)'
        CUNIT ='J/m**2'
        CALL MONIT_REGIST_SIB_SNAP ( CMARK , CTITLE , CUNIT )
        CMARK ='ENGLS4X'
        CTITLE='ENG SOIL 4 SNW (snap)'
        CUNIT ='J/m**2'
        CALL MONIT_REGIST_SIB_SNAP ( CMARK , CTITLE , CUNIT )
        CMARK ='ENGL4X'
        CTITLE='ENG SOIL 4 (snap)'
        CUNIT ='J/m**2'
        CALL MONIT_REGIST_SIB_SNAP ( CMARK , CTITLE , CUNIT )
!
        CMARK ='TMPL4X'
        CTITLE='TEMP SOIL 4 (snap)'
        CUNIT ='K'
        CALL MONIT_REGIST_SIB_SNAP ( CMARK , CTITLE , CUNIT )
      ENDIF
!
#endif
      RETURN
      END

