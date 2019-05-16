! sib 用 prm ファイル
! 99/05/xx 保坂征宏
!
      MODULE SIBPRM                 
!
      USE PRM , ONLY :              
!     1   IJPHY , JLPHY, NUMPE 
     1   IJPHY , JLPHY
!
      IMPLICIT NONE
!
!      INTEGER,PARAMETER :: NUM_PE = NUMPE       
!
      INTEGER,PARAMETER :: ITYO_veg = 12
      INTEGER,PARAMETER :: ITYO     = 2*ITYO_veg+1 ! 植生タイプ
      INTEGER,PARAMETER :: ITYP_veg = 40           ! 植生タイプ  12+1 
      INTEGER,PARAMETER :: ITYP     = 2*ITYP_veg+1 ! 植生タイプ
!
      INTEGER,PARAMETER :: IDX        = IJPHY      
      INTEGER,PARAMETER :: NUM_BLOCK  = JLPHY      
      INTEGER,PARAMETER :: IMON = 12      ! 12  月
      INTEGER,PARAMETER :: ICG  =  2      !  2  キャノピー・地面
      INTEGER,PARAMETER :: IWV  =  3      !  3  可視・近赤・長波
      INTEGER,PARAMETER :: IVN  =  2      !  2  可視・近赤        （短波）
      INTEGER,PARAMETER :: ILD  =  2      !  2  葉の生死
      INTEGER,PARAMETER :: IBD  =  2      !  2  直達・散乱
      INTEGER,PARAMETER :: IBS  =  2      !  2  雪なし、雪あり    
      INTEGER,PARAMETER :: IDP  =  3      !  3  土壌の層数
      INTEGER,PARAMETER :: ISN  =  4      !  4  雪の層数
      INTEGER,PARAMETER :: IWG  =  2      !  2  WTR or GLA
!
      INTEGER,PARAMETER :: IEST =  3      !  3 EVAP SUBL TRANS
!
      INTEGER,PARAMETER :: ICGB =  3      !  3 CAN,GRS,BSL
      INTEGER,PARAMETER :: ICGR =  3      !  4 CAN,GRD,REF
!
      INTEGER,PARAMETER :: ICS  =  2      !  4 CAN,SNW
      INTEGER,PARAMETER :: ICSR =  3      !  4 CAN,SNW,REF
!
      INTEGER,PARAMETER :: ISPT = 1 
!
      END MODULE SIBPRM


