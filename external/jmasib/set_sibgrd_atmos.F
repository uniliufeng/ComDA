! SiB äiéqÇ…ãlÇﬂÇÈ
! 99/03/xx ï€ç‚ê™çGÅAébíËî≈
! 99/08/04 ï€ç‚ê™çG ëÂêÚÇ≥ÇÒÇ…ìnÇ∑ÉoÅ[ÉWÉáÉìÉtÉBÉbÉNÉX

! 2000 5/10 bug fix by M.Ohizunmi

      SUBROUTINE SET_SIBGRD_ATMOS (
     I  NMAX               , 
     I  IMSKLT             ,
     I  GPDEL_PHY          , GPHAF_PHY          , GPFUL_PHY  , 
     I  GT_PHY             , GQ_PHY             , 
     I  GU_PHY             , GV_PHY             , 
     I  ZMEAN_PHY          , ZTEMP_PHY          , 
     I  PPLI_PHY           , PPCI_PHY            , 
     I  TMTX_PHY           , QMTX_PHY           , UVMTX_PHY  , 
     I  DAYTIME_1HR_PHY    , 
     O  ZMEAN_SIB          , ZTEMP_SIB          , 
     O  PPC_SIB            , PPL_SIB             ,
     O  TMP_REF_SIB        , HUM_REF_SIB        , 
     O  U_REF_SIB          , V_REF_SIB          ,
     O  P_DEL_SIB          , PSUR_SIB           , 
     O  Z_REF_SIB          , FAC_CORTMP_SIB     ,
     O  GMT_SIB            , GMQ_SIB            , GMUV_SIB  , 
     O  DAYTIME_1HR_SIB   ) 
!
      USE PRM    , ONLY : 
     1   IJPHY  , 
     1   JLPHY  , 
     1   KMAX 
!
      USE SIBPRM , ONLY : 
     1   IDX    , 
     1   IBD    ,
     1   ICG    ,
     1   IDP    ,
     1   ISN    ,
     1   ILD    ,
     1   IMON   ,
     1   ITYP   ,
     1   IWV        
!
      USE SIBCON , ONLY :
     1   GRAV_c , 
     1   R_DRY_c  , 
     1   R_BY_G_c ,
     1   CP_AIR_c , 
     1   AKAPPA_c , 
     1   fact_precipi_c 
!
      IMPLICIT NONE
!
! INPUT START
!
      INTEGER,INTENT(IN)   ::  IMSKLT   (IDX)         ! SiB äiéq -> ï®óù
!
      REAL(8),INTENT(IN)   ::  ZMEAN_PHY (IJPHY)
      REAL(8),INTENT(IN)   ::  ZTEMP_PHY (IJPHY)
      REAL(8),INTENT(IN)   ::  PPLI_PHY  (IJPHY)
      REAL(8),INTENT(IN)   ::  PPCI_PHY  (IJPHY)
      REAL(8),INTENT(IN)   ::  TMTX_PHY  (IJPHY,KMAX,2)
      REAL(8),INTENT(IN)   ::  QMTX_PHY  (IJPHY,KMAX,2)
      REAL(8),INTENT(IN)   ::  UVMTX_PHY (IJPHY,KMAX,3)
      REAL(8),INTENT(IN)   ::  DAYTIME_1HR_PHY (IJPHY)
!
! ëÂãCèÓïÒ
!
      REAL(8),INTENT(IN)   ::  GT_PHY    (IJPHY)
      REAL(8),INTENT(IN)   ::  GQ_PHY    (IJPHY)
      REAL(8),INTENT(IN)   ::  GU_PHY    (IJPHY)
      REAL(8),INTENT(IN)   ::  GV_PHY    (IJPHY)
      REAL(8),INTENT(IN)   ::  GPFUL_PHY (IJPHY)
      REAL(8),INTENT(IN)   ::  GPHAF_PHY (IJPHY)
      REAL(8),INTENT(IN)   ::  GPDEL_PHY (IJPHY)

!     INTEGER,INTENT(IN)   :: JL 
      INTEGER,INTENT(IN)   :: NMAX  
!
!  OUTPUT START
!
      REAL(8),INTENT(OUT)  ::  ZMEAN_SIB  (IDX)  ! cos(ìVí∏äp)
      REAL(8),INTENT(OUT)  ::  ZTEMP_SIB  (IDX)
      REAL(8),INTENT(OUT)  ::  PPC_SIB   (IDX)  ! ëŒó¨ê´ç~êÖ  
      REAL(8),INTENT(OUT)  ::  PPL_SIB   (IDX)  ! ëÂãKñÕã√åãê´ç~êÖ
      REAL(8),INTENT(OUT)  ::  FAC_CORTMP_SIB(IDX)  
!                                               ! (GPHAF(1)/GPFUL(1))**AKAPPA
      REAL(8),INTENT(OUT)  ::  PSUR_SIB  (IDX)  ! ínï\ãCà≥
      REAL(8),INTENT(OUT)  ::  P_DEL_SIB (IDX)  ! GPDEL(1)
      REAL(8),INTENT(OUT)  ::  DAYTIME_1HR_SIB (IDX)
!                                               
! ëÂãCèÓïÒ
!
      REAL(8),INTENT(OUT)  ::  TMP_REF_SIB (IDX)  ! ëÂãCç≈â∫ëwâ∑ìx
      REAL(8),INTENT(OUT)  ::  HUM_REF_SIB (IDX)  ! ëÂãCç≈â∫ëwÇÃî‰éº 
      REAL(8),INTENT(OUT)  ::  U_REF_SIB   (IDX)  ! ëÂãCç≈â∫ëwïóë¨
      REAL(8),INTENT(OUT)  ::  V_REF_SIB   (IDX)  ! ëÂãCç≈â∫ëwïóë¨
      REAL(8),INTENT(OUT)  ::  Z_REF_SIB   (IDX)  ! ï‚ê≥çÇÇ≥ÅH  ZM ? 
!
      REAL(8),INTENT(OUT)  ::  GMT_SIB   (IDX,2)
      REAL(8),INTENT(OUT)  ::  GMQ_SIB   (IDX,2)
      REAL(8),INTENT(OUT)  ::  GMUV_SIB  (IDX,3)

!
!-----------------------------------------------------------
!
! WORK 
      INTEGER      :: I 
      INTEGER      :: J 
!
      REAL(8),PARAMETER :: ZERO = 0.D0
!
      LOGICAL         LFIRST
! SAVE END
      DATA            LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
!       AKAPPA_c = R_DRY_c / CP_AIR_c
        WRITE(6,*) 'SET_SIBGRID_ATMOS 99/07 POCHAKA'
        LFIRST   = .FALSE. 

!------check vacant data location---------------------------

      do i=1,nmax
        J = IMSKLT(I) 
        if(gpful_phy(j).le.0.d0) then
           write(*,*) 'I=',i,'J=',j,'GPFUL=',gpful_phy(j)
        endif
      enddo  
!-----------------------------------------------------------

      END IF
!
      IF ( NMAX.LT.1 ) RETURN
!
#ifdef DEBUG
        write(*,*) 'GPHAF & GPFUL'
        write(*,*) 'NMAX=',NMAX
#endif
        
      DO 1000 I=1,NMAX 

        J = IMSKLT(I) 
!
        U_REF_SIB      (I) = GU_PHY    (J)   
        V_REF_SIB      (I) = GV_PHY    (J)   
        PSUR_SIB       (I) = GPHAF_PHY (J)   
        P_DEL_SIB      (I) = GPDEL_PHY (J)   
        FAC_CORTMP_SIB (I) = (GPHAF_PHY(J)/GPFUL_PHY(J))**AKAPPA_c 

!        write(*,*) 'J=',J,'FAC_CORTMP(',I,')=',FAC_CORTMP_SIB (I) 

        TMP_REF_SIB    (I) = GT_PHY    (J)                         
        HUM_REF_SIB    (I) = GQ_PHY    (J)                         
        Z_REF_SIB      (I) = GT_PHY    (J) * R_BY_G_c * 0.5        
     1                                     * GPDEL_PHY(J) 
     1                                     / GPHAF_PHY(J)
!
!
!        write(*,*) 'Z_REF_SIB  (',I,')=', Z_REF_SIB      (I) 
!
!                             â∫å¿èàóùí«â¡ ( 01/11/16 ) 
        PPL_SIB        (I) = MAX(PPLI_PHY(J),ZERO) * FACT_precipi_c  ! use
        PPC_SIB        (I) = MAX(PPCI_PHY(J),ZERO) * FACT_precipi_c  ! use
!
        ZMEAN_SIB      (I) = ZMEAN_PHY (J)                         ! use
        ZTEMP_SIB      (I) = ZTEMP_PHY (J)                         ! use
        DAYTIME_1HR_SIB(I) = DAYTIME_1HR_PHY (J)  
!
        GMT_SIB        (I,1) = TMTX_PHY  (J,1,1)                   
        GMT_SIB        (I,2) = TMTX_PHY  (J,1,2)                   
        GMQ_SIB        (I,1) = QMTX_PHY  (J,1,1)                   
        GMQ_SIB        (I,2) = QMTX_PHY  (J,1,2)                   
        GMUV_SIB       (I,1) = UVMTX_PHY (J,1,1)                   
        GMUV_SIB       (I,2) = UVMTX_PHY (J,1,2)                   
        GMUV_SIB       (I,3) = UVMTX_PHY (J,1,3)                   

 1000 CONTINUE
!
!------ remove below --------------
!
!      DO I=1,NMAX 
!        FAC_CORTMP_sib(I) =   1.001434440814268             ! ppp
!        z_ref_sib(I)      =  42.53747485208125              ! ppp 
!      ENDDO
!

#ifdef DEBUG
      write(*,*) 'end set_sibgrd_atmos' 
#endif
      RETURN
      END
