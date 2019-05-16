      MODULE COM_JOBINFO_SIB0109
!
      INTEGER  :: IUNIT_RESTART_IN  = 10
      INTEGER  :: IDSTAR(5)
      INTEGER  :: IDEND (5)
      INTEGER  :: IDATE (5)     
      INTEGER  :: INI_MON       = 1 
!
      integer  :: nocheck_restartfile_date = 0  ! check 0 nocheck 1 
!
      INTEGER  :: IOUT_8BYTE    = 1  
      CHARACTER(80)  :: CFILE_MONIT = ' '
! 20050803
!      CHARACTER(80)  :: CDIR_RESTART= '../output/'
!      CHARACTER(80)  :: CDIR_MONIT  = '../output/'
      CHARACTER(80)  :: CDIR_RESTART= 'output/'
      CHARACTER(80)  :: CDIR_MONIT  = 'output/'

!
      NAMELIST /SIB0109_COM_JOBINFO_SIB0109/ 
     1       IUNIT_RESTART_IN  , 
     1       IDSTAR            , 
     1       IDEND             , 
     1       IDATE             ,  
     1       INI_MON           , 
     1       IOUT_8BYTE        , 
     1       NOCHECK_RESTARTFILE_DATE  , 
     1       CDIR_RESTART      , 
     1       CDIR_MONIT        , 
     1       CFILE_MONIT        
!
      CONTAINS

        SUBROUTINE COM_JOBINFO_SIB0109_INI
!
        use com_stdinout_unit , only : iunit_in , iunit_out 
!
        idate(:) = 0 

         !--- ネームリスト読み込みによる変更       
          write(iunit_out,*) 'com_jobinfo_sib0109_ini' 
          READ (iunit_in ,SIB0109_COM_JOBINFO_SIB0109)
          WRITE(iunit_out,SIB0109_COM_JOBINFO_SIB0109)

        END SUBROUTINE COM_JOBINFO_SIB0109_INI
!
      END MODULE COM_JOBINFO_SIB0109


