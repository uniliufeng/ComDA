!==================================================================
      SUBROUTINE ISLSCP_RDFL ( IFILE , CMARK , ID_FILE_IN , IMASK , 
     O                         DAT4 ) 
!
      USE PRM , ONLY :  IDIM , JDIM 
      USE ISLSCP , ONLY : JCN_ISLSCP_GRADS 
      USE COM_RUNCONF_SIB0109 , ONLY :    JCN_SIB_ISLSCP_CYCLE1987
!
      USE CALENDAR , ONLY : 
     1     CALENDAR_RUN_GETKT ,
     1     CALENDAR_RUN_GETID
      USE MESSAGE, ONLY: MESSAGE_PUT
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
      INTEGER,SAVE       :: ID_ISLSCP_INI  (5)         ! 1987/01/01/00/00  
!
      DATA ID_ISLSCP_INI / 1987 , 1 , 1 , 0 , 0 /
!
      IF      ( CMARK .EQ. 'CLD ' ) THEN
        KT_LAG = 0
        CALL CALENDAR_RUN_GETID ( ID_FILE_IN , ID_FILE , KT_LAG )
!
!----modified start tok
        IF ( JCN_SIB_ISLSCP_CYCLE1987 .EQ. 1 ) THEN 
           ID_FILE(1) = 1987
#          ifdef DEBUG
             WRITE(6,*) 'CYCLE ISLSCP_RDFL DATA MODIFIED TO ' , ID_FILE
#          else
             CALL MESSAGE_PUT('ISLSCP_RDFL: CYCLE1987 YEAR CHANGE')
#          endif
        ENDIF
!----modified end tok
!
        IREC = (ID_FILE(1)-1987)*12 + ID_FILE(2)
        IF (IREC < 0 .OR. IREC > 24) THEN
          write(6,*) 'ISLSCP_RDFL CLD ERROR ' , IREC , ID_FILE
          stop 999
        ENDIF  
        ! 1986îN12åé15ì˙ÇÃóvãÅÇ…ëŒÇµÇƒÇÕ1987îN1åé15ì˙Çï‘Ç∑
        IF (IREC == 0) THEN
          IREC = 1
          ID_FILE(1:2) = (/1987, 1/)
        ENDIF
      ELSE IF ( CMARK .EQ. 'LWDN' ) THEN
!       KT_LAG = 9
        KT_LAG = 3
        CALL CALENDAR_RUN_GETID ( ID_FILE_IN , ID_FILE , KT_LAG )
!
!----modified start tok
        IF ( JCN_SIB_ISLSCP_CYCLE1987 .EQ. 1 ) THEN 
          ID_FILE(1) = 1987
#         ifdef DEBUG
            WRITE(6,*) 'CYCLE ATM_RDFL LWDN DATA MODIFIED TO ' , ID_FILE
#         else
            call message_put('ISLSCP_RDFL(LWDN): YEAR set to 1987')
#         endif
        ELSE
#         ifdef DEBUG
            WRITE(6,*) 'ISLSCP_RDFL LWDN DATA MODIFIED TO ' , ID_FILE
#         endif
        ENDIF
!----modified end tok
!
        CALL CALENDAR_RUN_GETKT ( ID_ISLSCP_INI , ID_FILE , 4 , KT_LAG )
        IREC = KT_LAG / 6 + 1 
      ELSE 
        KT_LAG = 0
        CALL CALENDAR_RUN_GETID ( ID_FILE_IN , ID_FILE , KT_LAG )
!
!----modified start tok
        IF ( JCN_SIB_ISLSCP_CYCLE1987 .EQ. 1 ) THEN 
           ID_FILE(1) = 1987
#          ifdef DEBUG
             WRITE(6,*) 'CYCLE ISLSCP_RDFL DATA MODIFIED TO ' , ID_FILE
#          else
            call message_put('ISLSCP_RDFL: YEAR set to 1987')
#          endif
        ENDIF
!----modified end tok
!
        CALL CALENDAR_RUN_GETKT(ID_ISLSCP_INI, ID_FILE, 4, KT_LAG)      
        IREC = KT_LAG / 6 + 1 
      ENDIF
!
      IF ( JCN_ISLSCP_GRADS .EQ. 0 ) THEN
        READ (IFILE, REC=IREC)  CMARK_TMP, IY, IM, ID, IH, DAT4
!
        IF (      CMARK_TMP .NE. CMARK 
     1       .OR. IY     .NE. ID_FILE(1) 
     1       .OR. IM     .NE. ID_FILE(2) 
     1       .OR. ID     .NE. ID_FILE(3) 
     1       .OR. IH     .NE. ID_FILE(4) 
     1     ) THEN
          WRITE(6,*) 'ISLSCP_RDFL ERROR. REC = ' , IREC 
          WRITE(6,*) '  MARK = ' , cmark , ' ' , cmark_tmp  
          WRITE(6,*) '  DATE = ' , IY , IM , ID, IH , ' ' , ID_FILE
          STOP 999
        ENDIF
      ELSE
        READ ( IFILE , REC=IREC )  DAT4
      ENDIF
!
#     ifdef DEBUG
      WRITE(6,*) 'ISLSCP_RDFL ' , CMARK , ID_FILE   
#     endif
!
      call islscp_data_maxmin ( DAT4 , IDIM , JDIM , IMASK , CMARK )
!
      RETURN
      END
