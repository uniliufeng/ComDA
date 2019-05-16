      SUBROUTINE monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      IMPLICIT NONE
!
      CHARACTER( 7) , INTENT(IN)  :: CMARK 
      CHARACTER(32) , INTENT(IN)  :: CTITLE 
      CHARACTER(13) , INTENT(IN)  :: CUNIT
!
      REAL(8) :: smax,smin,sl,smag
!
      smax =   100.0D0
      smin = - 100.0D0
      sl   =     0.D0
      smag =     1.D0  
!
      CALL MONIT_REGIST ( 
     1    1, '2D  ', CMARK , 'LAND', 'AVR', 'PH' ,             
     1    CTITLE   , CUNIT , 
     1    0 , smax, smin, sl, smag )              
!
      RETURN
      END            
!======================================================
      SUBROUTINE monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
      IMPLICIT NONE
!
      CHARACTER( 7) , INTENT(IN)  :: CMARK 
      CHARACTER(32) , INTENT(IN)  :: CTITLE 
      CHARACTER(13) , INTENT(IN)  :: CUNIT
!
      REAL(8) :: smax,smin,sl,smag
!
      smax =   100.0D0
      smin = - 100.0D0
      sl   =     0.D0
      smag =     1.D0  
!
      CALL MONIT_REGIST ( 
     1    1, '2D  ', CMARK , 'LAND', 'SNP', 'PH' ,             
     1    CTITLE   , CUNIT , 
     1    0 , smax, smin, sl, smag )              
!
      RETURN
      END            





