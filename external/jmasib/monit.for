	SUBROUTINE MONIT(GLON, GLAT, DELT_MODEL, MON, KT0, KTP, KTSTAR,
	1 KTEND, IDATE, FSECP, IMASK)
	

	USE SIB0109 , ONLY :
     1  SIB0109_RUN_MONITOR_INITIAL
	
	USE PRM , ONLY : 
     1    IJPHY , 
     1    JLPHY , 
     1    INTRI , 
     1    INTRJ ,
     1    IDIM  , 
     1    JDIM  , 
     1    JMAX  , 
     1    ISPT  ,
     1    KMP1
     
      USE COM_RUNCONF , ONLY : 
     1    JCNIMNT 
	
	USE COM_RUNCONF_SIB0109 , ONLY : 
     1    JCN_IWL_SKIP
     1     

	USE COM_FILEUNIT , ONLY : 
     I   IMONIT_2CFL ,IMONIT_ZCFL ,IMONIT_3CFL ,
     I   IMONIT_ZETACFL ,IMONIT_3ETACFL ,
     I   IMONIT_6HRCFL ,IMONIT_DAYCFL ,IMONIT_GLOBCFL,
     1   imonit_2fl , imonit_2fl8 , imonit_sfl
	
	USE COM_JOBINFO_SIB0109 , ONLY :
     1    COM_JOBINFO_SIB0109_INI   ,
     1    IDSTAR            , 
!    1    IDATE             , 
     1    IDEND             , 
     1    IOUT_8BYTE        ,
     1    INI_MON           , 
     1    CFILE_MONIT       ,  
     1    CDIR_MONIT 
     
      USE COM_STEP_SIB0109 , only : 
     1   ICN_SIB0109_CALC_SOIL_SNOW    
     
      real(8)   :: GLON (IDIM,JDIM)    
      real(8)   :: GLAT (IDIM,JDIM)   
	real(8)   :: GLON_AVERAGE(IDIM)
      real(8)   :: GLAT_AVERAGE(JDIM)

	real(8)	  :: DELT_MODEL
	integer   :: IDATE(5)
	real(8)   :: FSECP
	INTEGER,INTENT(IN) :: IMASK     (IDIM,JDIM)

      CHARACTER(100) ::
     I  CMONIT_2_NAME, CMONIT_Z_NAME, CMONIT_3_NAME,
     I  CMONIT_ZETA_NAME, CMONIT_3ETA_NAME,
     I  CMONIT_6HR_NAME, CMONIT_DAY_NAME, CMONIT_GLOB_NAME
!
      CHARACTER(8)  :: MODEL
      CHARACTER(11) :: RESL    
      CHARACTER(4)  :: EXPR
!     CHARACTER(1)  :: EXPR1
!     CHARACTER(2)  :: EXPR2
      CHARACTER(80) :: CINF(10)
      REAL(8) :: GW          (JMAX)
      REAL(8) :: COSCLT      (JMAX)
      REAL(8) :: PA(KMP1)
      REAL(8) :: PB(KMP1)    
!
      integer  :: k 
      real(8)  :: TOTMON 
      REAL(8)  :: TOTMON_X 
      INTEGER   :: J

	!
      MODEL = 'MJSIB907'
      RESL  = '1DIMversion'
!     CHARACTER(1)  :: EXPR1
!     CHARACTER(2)  :: EXPR2
      DO J=1,JMAX
        GW(J)      = 1.
        COSCLT (J) = (JMAX-2*J+1.) / JMAX  
      ENDDO
      EXPR = 'TEST'  
!
      do 2 k=1,10
        cinf(k) = ' '
 2    continue    
!
      write(6,*) 'main : jcn_iwl_skip , jcnimnt original ' , 
     1           jcn_iwl_skip , jcnimnt   
      IF ( JCN_IWL_SKIP.LE.-1 ) THEN     ! 毎ステップ出力
        IF ( JCNIMNT .NE. -1 ) THEN 
          WRITE(6,*) 'MAIN WARNING ' , 
     1               ' JCNIMNT IS MODIFIED ' , JCNIMNT , ' TO ' , -1
        ENDIF 
        jcnimnt = -1        
      ELSE                               ! JCNIMNT 時間おき出力
        IF ( JCNIMNT .LT. 0 ) THEN 
          WRITE(6,*) 'MAIN WARNING ' , 
     1               ' JCN_IWL_SKIP IS MODIFIED ' 
          JCN_IWL_SKIP = - JCN_IWL_SKIP 
        ENDIF 
      ENDIF
      write(6,*) 'main : jcn_iwl_skip , jcnimnt ' , 
     1           jcn_iwl_skip , jcnimnt   
!
      DO 1 K=1,KMP1
        PA(K) = 0.
        PB(K) = 0.
 1    CONTINUE  
!
! 吉村モニタファイル名の設定
!
      IF ( CFILE_MONIT == ' ' ) then
        CMONIT_2_NAME = 
     1       'newsib_monit_YYYY_MM_DD_HH_' // 
     1                    'YYYY_MM_DD_HH'
        CALL REPLACE_INT(CMONIT_2_NAME, 'YYYY', IDSTAR(1))
        CALL REPLACE_INT(CMONIT_2_NAME, 'MM', IDSTAR(2))
        CALL REPLACE_INT(CMONIT_2_NAME, 'DD', IDSTAR(3))
        CALL REPLACE_INT(CMONIT_2_NAME, 'HH', IDSTAR(4))
        CALL REPLACE_INT(CMONIT_2_NAME, 'YYYY', IDEND(1))
        CALL REPLACE_INT(CMONIT_2_NAME, 'MM', IDEND(2))
        CALL REPLACE_INT(CMONIT_2_NAME, 'DD', IDEND(3))
        CALL REPLACE_INT(CMONIT_2_NAME, 'HH', IDEND(4))
      ELSE 
        CMONIT_2_NAME = CFILE_MONIT
      ENDIF
!
! 吉村モニタファイルオープン
!
      open ( imonit_2fl  , 
     1        file  = TRIM(CDIR_MONIT) // TRIM(cmonit_2_name)//'.dr',
     1        access='direct' , recl=4*IJPHY*JLPHY , 
     1        form  ='unformatted' ) 
      if ( iout_8byte.eq.1 ) then
        open ( imonit_2fl8 , 
     1        file  = TRIM(CDIR_MONIT) // TRIM(cmonit_2_name)//'_8.dr',
     1        access='direct' , recl=8*IJPHY*JLPHY , 
     1        form='unformatted' ) 
      endif
      open ( imonit_2cfl , 
     1        file  = TRIM(CDIR_MONIT) // TRIM(cmonit_2_name)//'.ctl')
!
! 吉村モニタ初期化 
!
      call monit_ini
!
      glon_average(:) = sum(glon, dim=2) / JDIM
      glat_average(:) = sum(glat, dim=1) / IDIM
      call monit_grads_ctl (
     I   IMONIT_2CFL ,IMONIT_ZCFL ,IMONIT_3CFL ,
     I   IMONIT_ZETACFL ,IMONIT_3ETACFL ,
     I   IMONIT_6HRCFL ,IMONIT_DAYCFL ,IMONIT_GLOBCFL,
     I   CMONIT_2_NAME ,CMONIT_Z_NAME ,CMONIT_3_NAME ,
     I   CMONIT_ZETA_NAME ,CMONIT_3ETA_NAME ,
     I   CMONIT_6HR_NAME ,CMONIT_DAY_NAME ,CMONIT_GLOB_NAME,
     I   cosclt, idstar, idend, KTSTAR, KTEND,         
     I   model, resl, gw, pa, pb , glon_average , glat_average)
      call monit_clear(totmon) 

	IF ( INI_MON .EQ.1 ) THEN
        write(6,*) 'initial snap shot monitoring ' , DELT_MODEL
!
        CALL SIB0109_RUN_MONITOR_INITIAL ( MON , KT0 )   
!    
        TOTMON_X = DELT_MODEL            ! 積算時間 ..  モニタをだます
!
        ICN_SIB0109_CALC_SOIL_SNOW = 1   ! 雪土壌も計算したステップであると、
                                         ! モニタをだます
!
        CALL MONIT_OUT_2
     1  ( IMONIT_2FL       , IMONIT_2FL8   , IMONIT_SFL    ,
     1    IDATE , MODEL    , RESL  , EXPR  , CINF  ,
     1    KTP   , TOTMON_X , FSECP , IMASK , GW )
!
        CALL MONIT_CLEAR  ( TOTMON_X )
!
        write(6,*) 'initial snap shot monitoring end ' 
!
      ENDIF
!



	END SUBROUTINE MONIT




	SUBROUTINE OUTPUT(IDATE, KTP, TOTMON, FSECP, IMASK)

	USE PRM , ONLY : 
     1    IDIM  , 
     1    JDIM  , 
     1    JMAX 
     
	USE COM_STEP , ONLY :
     1   ICNMNTMON                       ! mj98 com_step

	USE COM_FILEUNIT , ONLY : 
     I    IMONIT_2FL     , IMONIT_2FL8   , IMONIT_SFL

	integer   :: IDATE(5)
	real(8)   :: FSECP
	real(8)   :: TOTMON
	INTEGER,INTENT(IN) :: IMASK     (IDIM,JDIM)

	CHARACTER(8)  :: MODEL
      CHARACTER(11) :: RESL    
      CHARACTER(4)  :: EXPR
      CHARACTER(80) :: CINF(10)
      REAL(8) :: GW          (JMAX)
	integer  :: k 
      INTEGER   :: J

	MODEL = 'MJSIB907'
      RESL  = '1DIMversion'
      DO J=1,JMAX
        GW(J)      = 1.
      ENDDO
      EXPR = 'TEST'  
!
      do 2 k=1,10
        cinf(k) = ' '
 2    continue    
       

	IF ( ICNMNTMON.EQ.1 ) THEN
        CALL MONIT_OUT_2
     1    ( IMONIT_2FL     , IMONIT_2FL8   , IMONIT_SFL    ,
     1      IDATE , MODEL  , RESL  , EXPR  , CINF  ,
     1      KTP   , TOTMON , FSECP , IMASK , GW )
        CALL MONIT_CLEAR  ( TOTMON )
      ENDIF

	END SUBROUTINE OUTPUT
	
