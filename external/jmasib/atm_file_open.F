!* ëÂãCã≠êßóÕÇÃÉtÉ@ÉCÉãÇäJÇ≠ Ç®ÇÊÇ— ï¬Ç∂ÇÈ
!      SUBROUTINE ATM_FILE_OPEN
!      SUBROUTINE ATM_FILE_CLOSE 
!* êAê∂Ç‚çÇìxÉtÉ@ÉCÉãÇäJÇ≠
!      SUBROUTINE ATM_CNFILE_OPEN
!
      SUBROUTINE ATM_FILE_OPEN (ID_NOW, ID_OPEN)
!
      USE PRM , ONLY : IDIM , JDIM 
!
      USE ATM_FILE_UNIT , ONLY :
     1    IFILE_UREF  , 
     1    IFILE_VREF  , 
     1    IFILE_TREF  , 
     1    IFILE_QREF  , 
     1    IFILE_PREF  , 
     1    IFILE_PSUR  , 
     1    IFILE_LPRC  , 
     1    IFILE_CPRC  , 
     1    IFILE_LWDN  , 
     1    IFILE_SWDN  , 
     1    IFILE_CLD   ,
     1    ID_ATM_INI  
!
      USE COM_RUNCONF_SIB0109 , ONLY :    JCN_SIB_ATM_CYCLE
!
      IMPLICIT NONE 
!
      INTEGER,INTENT(IN)      :: ID_NOW (5)
      INTEGER,INTENT(IN)      :: ID_OPEN   ! 1:í∑îg  2:ç~êÖÅEíZîg  3:ÇªÇÃëº
!
      CHARACTER(100) :: CFILE_UREF
      CHARACTER(100) :: CFILE_VREF
      CHARACTER(100) :: CFILE_TREF
      CHARACTER(100) :: CFILE_QREF
      CHARACTER(100) :: CFILE_PREF
      CHARACTER(100) :: CFILE_PSUR
      CHARACTER(100) :: CFILE_LPRC
      CHARACTER(100) :: CFILE_CPRC
      CHARACTER(100) :: CFILE_LWDN
      CHARACTER(100) :: CFILE_SWDN
      CHARACTER(100) :: CFILE_CLD 
!
      CHARACTER(4) :: CDATE
      LOGICAL,SAVE  :: LFIRST(3) = (/.FALSE., .FALSE., .FALSE./)
!
      IF(JCN_SIB_ATM_CYCLE.NE.0) THEN
         ID_ATM_INI(1,ID_OPEN) = JCN_SIB_ATM_CYCLE
      ELSE
         ID_ATM_INI(1,ID_OPEN) = ID_NOW(1)
      ENDIF
      ID_ATM_INI(2,ID_OPEN) =    1
      ID_ATM_INI(3,ID_OPEN) =    1
      ID_ATM_INI(4,ID_OPEN) =    0
      ID_ATM_INI(5,ID_OPEN) =    0
!
      CDATE = '0000'
      WRITE(CDATE(1:4), '(I4)') ID_ATM_INI(1,ID_OPEN)
!      write(6,*) 'cdate=', cdate
!
      CFILE_UREF = 'input/' // CDATE // '/uref'
      CFILE_VREF = 'input/' // CDATE // '/vref'
      CFILE_TREF = 'input/' // CDATE // '/tref' 
      CFILE_QREF = 'input/' // CDATE // '/qref'
      CFILE_PREF = 'input/' // CDATE // '/pref'
      CFILE_PSUR = 'input/' // CDATE // '/psur'
      CFILE_LPRC = 'input/' // CDATE // '/rrl' 
      CFILE_CPRC = 'input/' // CDATE // '/rrc' 
      CFILE_LWDN = 'input/' // CDATE // '/rldb'
      CFILE_SWDN = 'input/' // CDATE // '/rsdb'
      CFILE_CLD  = 'input/' // CDATE // '/cla' 
      write(6,*) 'file name=', cfile_tref
!      return
!
      IF ( ID_OPEN.LT.1 .OR. ID_OPEN.GT.3 ) THEN
          write(6,*) 'atm_file_open error id_open=', id_open
          STOP 999
      ENDIF
!
      IF ( LFIRST(ID_OPEN) ) THEN
          write(6,*) 'atm_file_close() id_close=', id_open
          CALL ATM_FILE_CLOSE(ID_OPEN) 
          LFIRST(ID_OPEN) = .FALSE.
      ENDIF

      IF(ID_OPEN.EQ.3) THEN
        OPEN ( IFILE_UREF , FILE=CFILE_UREF , FORM='UNFORMATTED' , 
     1         ACCESS='DIRECT' , RECL = 4*IDIM*JDIM   )
        OPEN ( IFILE_VREF , FILE=CFILE_VREF , FORM='UNFORMATTED' , 
     1         ACCESS='DIRECT' , RECL = 4*IDIM*JDIM   )
        OPEN ( IFILE_TREF , FILE=CFILE_TREF , FORM='UNFORMATTED' , 
     1         ACCESS='DIRECT' , RECL = 4*IDIM*JDIM   )
        OPEN ( IFILE_QREF , FILE=CFILE_QREF , FORM='UNFORMATTED' , 
     1         ACCESS='DIRECT' , RECL = 4*IDIM*JDIM   )
        OPEN ( IFILE_PREF , FILE=CFILE_PREF , FORM='UNFORMATTED' , 
     1         ACCESS='DIRECT' , RECL = 4*IDIM*JDIM   )
        OPEN ( IFILE_PSUR , FILE=CFILE_PSUR , FORM='UNFORMATTED' , 
     1         ACCESS='DIRECT' , RECL = 4*IDIM*JDIM   )
        OPEN ( IFILE_CLD  , FILE=CFILE_CLD  , FORM='UNFORMATTED' , 
     1         ACCESS='DIRECT' , RECL = 4*IDIM*JDIM   )
      ENDIF
!
      IF(ID_OPEN.EQ.2) THEN
        OPEN ( IFILE_LPRC , FILE=CFILE_LPRC , FORM='UNFORMATTED' , 
     1         ACCESS='DIRECT' , RECL = 4*IDIM*JDIM   )
        OPEN ( IFILE_CPRC , FILE=CFILE_CPRC , FORM='UNFORMATTED' , 
     1         ACCESS='DIRECT' , RECL = 4*IDIM*JDIM   )
        OPEN ( IFILE_SWDN , FILE=CFILE_SWDN , FORM='UNFORMATTED' , 
     1         ACCESS='DIRECT' , RECL = 4*IDIM*JDIM   )
      ENDIF
!
      IF(ID_OPEN.EQ.1) THEN
        OPEN ( IFILE_LWDN , FILE=CFILE_LWDN , FORM='UNFORMATTED' , 
     1       ACCESS='DIRECT' , RECL = 4*IDIM*JDIM   )
      ENDIF
!
      write(6,*) 'atm_file_open(',id_open,'):', ID_ATM_INI(1,ID_OPEN)
      LFIRST(ID_OPEN) = .TRUE.
!
      RETURN
      END
      SUBROUTINE ATM_FILE_CLOSE(ID_CLOSE) 
!
      USE ATM_FILE_UNIT , ONLY :
     1    IFILE_UREF  , 
     1    IFILE_VREF  , 
     1    IFILE_TREF  , 
     1    IFILE_QREF  , 
     1    IFILE_PREF  , 
     1    IFILE_PSUR  , 
     1    IFILE_LPRC  , 
     1    IFILE_CPRC  , 
     1    IFILE_LWDN  , 
     1    IFILE_SWDN  , 
     1    IFILE_CLD
!
      INTEGER,INTENT(IN)      :: ID_CLOSE   ! 1:í∑îg  2:ç~êÖÅEíZîg  3:ÇªÇÃëº

      IF(ID_CLOSE.EQ.3) THEN
      CLOSE ( IFILE_UREF )
      CLOSE ( IFILE_VREF )
      CLOSE ( IFILE_TREF )
      CLOSE ( IFILE_QREF )
      CLOSE ( IFILE_PREF )
      CLOSE ( IFILE_PSUR )
      CLOSE ( IFILE_CLD  )
      ENDIF
!
      IF(ID_CLOSE.EQ.2) THEN
      CLOSE ( IFILE_LPRC )
      CLOSE ( IFILE_CPRC )
      CLOSE ( IFILE_SWDN )
      ENDIF
!
      IF(ID_CLOSE.EQ.1) THEN
      CLOSE ( IFILE_LWDN )
      ENDIF
!
      RETURN
      END
!
!
      SUBROUTINE ATM_CNFILE_OPEN
!
      USE PRM , ONLY : IDIM , JDIM 
!
      USE ATM_FILE_UNIT , ONLY :
     1    IFILE_VEG  
!
      IMPLICIT NONE 
!
      CHARACTER(100) :: CFILE_VEG
!
!      CFILE_VEG  = '/data1/a1113/era15/config/typ_t63'
!      CFILE_VEG  = '/users/a1/a1113/data/typ_t106_jma'
      CFILE_VEG  = 'input/maptyp'
      write(6,*) 'config set file name=', cfile_veg,ifile_veg
!
!      OPEN ( IFILE_VEG  , FILE=CFILE_VEG  , FORM='UNFORMATTED' , 
!     1       ACCESS='DIRECT' , RECL = 2*IDIM*JDIM   )
      OPEN ( IFILE_VEG  , FILE=CFILE_VEG  , FORM='UNFORMATTED' , 
     1       ACCESS='DIRECT' , RECL = 4*IDIM*JDIM   )
!
      RETURN
      END





