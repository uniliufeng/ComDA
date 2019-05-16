!
! ISLSCP ファイルたちを開く。ファイル形式は JCN_ISLSCP_GRADS が非零の
! ときは GRADS 形式であり、そうでなければ mabiki 形式である。
!
      SUBROUTINE ISLSCP_FILE_OPEN 
!
      USE PRM , ONLY : IDIM , JDIM 
!
      USE ISLSCP_FILE_UNIT , ONLY :
     1    IFILE_MWND  , 
     1    IFILE_TEMP  , 
     1    IFILE_QREF  , 
     1    IFILE_PRSS  , 
     1    IFILE_TPRC  , 
     1    IFILE_CPRC  , 
     1    IFILE_LWDN  , 
     1    IFILE_SWDN  , 
     1    IFILE_VEG   , 
     1    IFILE_CLD
!
      USE ISLSCP , ONLY : JCN_ISLSCP_GRADS 
!
      CHARACTER(40) :: CFILE_MWND
      CHARACTER(40) :: CFILE_TEMP
      CHARACTER(40) :: CFILE_QREF
      CHARACTER(40) :: CFILE_PRSS
      CHARACTER(40) :: CFILE_TPRC
      CHARACTER(40) :: CFILE_CPRC
      CHARACTER(40) :: CFILE_LWDN
      CHARACTER(40) :: CFILE_SWDN
      CHARACTER(40) :: CFILE_VEG
      CHARACTER(40) :: CFILE_CLD 
!
      CHARACTER(11) :: CTAIL_R4
      CHARACTER(11) :: CTAIL_I4
!
      INTEGER       :: IHEAD
!
      CTAIL_R4 = '_XXX_YYY_R4'
      CTAIL_I4 = '_XXX_YYY_I4'
!
      CALL REPLACE_INT(CTAIL_R4, 'XXX', IDIM ) 
      CALL REPLACE_INT(CTAIL_I4, 'XXX', IDIM ) 
      CALL REPLACE_INT(CTAIL_R4, 'YYY', JDIM ) 
      CALL REPLACE_INT(CTAIL_I4, 'YYY', JDIM ) 
!
      IF ( JCN_ISLSCP_GRADS .EQ. 0 ) THEN
        IHEAD = 20
        CFILE_CLD  = 'input/CLD'       // CTAIL_R4 
        CFILE_MWND = 'input/MWND'      // CTAIL_R4 
        CFILE_TEMP = 'input/TEMP'      // CTAIL_R4 
        CFILE_QREF = 'input/QREF'      // CTAIL_R4 
        CFILE_PRSS = 'input/PRSS'      // CTAIL_R4 
        CFILE_TPRC = 'input/TPRC'      // CTAIL_R4 
        CFILE_CPRC = 'input/CPRC'      // CTAIL_R4 
        CFILE_LWDN = 'input/LWDN'      // CTAIL_R4 
        CFILE_SWDN = 'input/SWDN'      // CTAIL_R4 
        CFILE_VEG  = 'input/VEG_MAP'   // CTAIL_I4 
      ELSE
        IHEAD = 0
        CFILE_CLD  = '../ISLSCP_360X180/CLD'       // CTAIL_R4 
        CFILE_MWND = '../ISLSCP_360X180/MWND'      // CTAIL_R4 
        CFILE_TEMP = '../ISLSCP_360X180/TEMP'      // CTAIL_R4 
        CFILE_QREF = '../ISLSCP_360X180/QREF'      // CTAIL_R4 
        CFILE_PRSS = '../ISLSCP_360X180/PRSS'      // CTAIL_R4 
        CFILE_TPRC = '../ISLSCP_360X180/TPRC'      // CTAIL_R4 
        CFILE_CPRC = '../ISLSCP_360X180/CPRC'      // CTAIL_R4 
        CFILE_LWDN = '../ISLSCP_360X180/LWDN'      // CTAIL_R4 
        CFILE_SWDN = '../ISLSCP_360X180/SWDN'      // CTAIL_R4 
        CFILE_VEG  = '../ISLSCP_360X180/VEG_MAP'   // CTAIL_I4 
      ENDIF
!
      write(6,*) 'islscp_file_open file_names ' , CFILE_MWND 
     1 ,  CFILE_TEMP
     1 ,  CFILE_QREF
     1 ,  CFILE_PRSS
     1 ,  CFILE_TPRC
     1 ,  CFILE_CPRC
     1 ,  CFILE_LWDN
     1 ,  CFILE_SWDN
     1 ,  CFILE_VEG
     1 ,  CFILE_CLD 
!
      OPEN ( IFILE_MWND , FILE=CFILE_MWND , FORM='UNFORMATTED' , 
     1       ACCESS='DIRECT' , RECL = 4*IDIM*JDIM + IHEAD )
      OPEN ( IFILE_TEMP , FILE=CFILE_TEMP , FORM='UNFORMATTED' , 
     1       ACCESS='DIRECT' , RECL = 4*IDIM*JDIM + IHEAD )
      OPEN ( IFILE_QREF , FILE=CFILE_QREF , FORM='UNFORMATTED' , 
     1       ACCESS='DIRECT' , RECL = 4*IDIM*JDIM + IHEAD )
      OPEN ( IFILE_PRSS , FILE=CFILE_PRSS , FORM='UNFORMATTED' , 
     1       ACCESS='DIRECT' , RECL = 4*IDIM*JDIM + IHEAD )
      OPEN ( IFILE_TPRC , FILE=CFILE_TPRC , FORM='UNFORMATTED' , 
     1       ACCESS='DIRECT' , RECL = 4*IDIM*JDIM + IHEAD )
      OPEN ( IFILE_CPRC , FILE=CFILE_CPRC , FORM='UNFORMATTED' , 
     1       ACCESS='DIRECT' , RECL = 4*IDIM*JDIM + IHEAD )
      OPEN ( IFILE_LWDN , FILE=CFILE_LWDN , FORM='UNFORMATTED' , 
     1       ACCESS='DIRECT' , RECL = 4*IDIM*JDIM + IHEAD )
      OPEN ( IFILE_SWDN , FILE=CFILE_SWDN , FORM='UNFORMATTED' , 
     1       ACCESS='DIRECT' , RECL = 4*IDIM*JDIM + IHEAD )
      OPEN ( IFILE_CLD  , FILE=CFILE_CLD  , FORM='UNFORMATTED' , 
     1       ACCESS='DIRECT' , RECL = 4*IDIM*JDIM + IHEAD )
      OPEN ( IFILE_VEG  , FILE=CFILE_VEG  , FORM='UNFORMATTED' , 
     1       ACCESS='DIRECT' , RECL = 4*IDIM*JDIM + IHEAD )
!
      RETURN
      END
