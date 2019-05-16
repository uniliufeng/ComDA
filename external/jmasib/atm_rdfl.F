!==================================================================
      SUBROUTINE ATM_RDFL ( IFILE , CMARK , ID_FILE_IN , IMASK , 
     O                      DAT4 ) 
!
      USE PRM , ONLY :  IDIM , JDIM 
      USE ATM_FILE_UNIT , ONLY :   ID_ATM_INI  
      USE COM_RUNCONF_SIB0109 , ONLY :    JCN_SIB_ATM_CYCLE
!
      USE CALENDAR , ONLY : 
     1     CALENDAR_RUN_GETKT ,
     1     CALENDAR_RUN_GETID
!
      IMPLICIT NONE
!
      INTEGER     ,INTENT(IN)    :: IFILE     
      INTEGER     ,INTENT(IN)    :: IMASK(IDIM,JDIM)
      CHARACTER(4),INTENT(IN)    :: CMARK 
      INTEGER     ,INTENT(IN)    :: ID_FILE_IN (5)
      REAL(4)     ,INTENT(OUT)   :: DAT4 (IDIM*JDIM)
!
      INTEGER            :: IY , IM , ID , IH 
!     INTEGER            :: I
      INTEGER            :: KT_LAG
      CHARACTER(4)       :: CMARK_TMP
      INTEGER            :: IREC
      INTEGER            :: ID_FILE     (5)
      INTEGER            :: ID_ELEM
!
!
      IF ( CMARK .EQ. 'LWDN' ) THEN
        KT_LAG = -3
        ID_ELEM= 1
      ELSE IF ( CMARK.EQ.'LPRC' .OR. CMARK.EQ.'CPRC' .OR.
     -     CMARK.EQ.'SWDN' ) THEN
        KT_LAG = 0
        ID_ELEM= 2
      ELSE 
        KT_LAG = 0
        ID_ELEM= 3
      ENDIF
!
      CALL CALENDAR_RUN_GETID ( ID_FILE_IN , ID_FILE , KT_LAG )
      IF ( JCN_SIB_ATM_CYCLE .NE. 0 .AND.
     -     ID_FILE(1).NE.JCN_SIB_ATM_CYCLE) THEN 
          ID_FILE(1) = JCN_SIB_ATM_CYCLE
          WRITE(6,*) 'CYCLE ATM_RDFL DATA MODIFIED TO ' , ID_FILE
      ENDIF
      CALL CALENDAR_RUN_GETKT 
     1        ( ID_ATM_INI(1,ID_ELEM) , ID_FILE , 4 , KT_LAG )
      IREC = KT_LAG / 6 + 1 

      READ ( IFILE , REC=IREC )  DAT4
!
      WRITE(6,*) 'ATM_RDFL ' , CMARK , ID_FILE, IREC
!
      call atm_data_maxmin ( DAT4 , IDIM , JDIM , IMASK , CMARK )
!
      RETURN
      END
