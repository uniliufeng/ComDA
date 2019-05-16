! 1 次元 SiB 用ツール : 大気状態などの情報を与える。
!     
! 99/07/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE TEST1D_CONFIG_SET (JL, GLON, GLAT, IMASK) 
!
      USE PRM , ONLY :
     1   IDIM  ,
     1   JDIM  ,
     1   JPHY2 ,
     1   JLPHY , 
     1   IJPHY , 
     1   INTRI , 
     1   INTRJ , 
     1   ISPT 
!
      USE SIBPRM , ONLY :
     1   IDX     
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: JL 
!
! OUTPUT
!
      INTEGER,INTENT(OUT)   :: IMASK            (ISPT*IJPHY,JLPHY) 
      REAL(8),INTENT(OUT)   :: GLON             (IDIM, JDIM)
      REAL(8),INTENT(OUT)   :: GLAT             (IDIM, JDIM)
!
      INTEGER :: IJ_SIB 
      INTEGER :: IJ_PHY
      INTEGER :: I_PHY 
!     INTEGER :: J_PHY 
!     INTEGER :: J_PHY_X 
!
      INTEGER :: I
      INTEGER :: J
!
      CALL IRESET ( IMASK           (1,JL) , IJPHY*ISPT )
!
      DO IJ_SIB = 1,IDX
        IJ_PHY = (IJ_SIB-1) / ISPT + 1 
!
        I_PHY   = MOD ( IJ_PHY - 1 , IDIM ) + 1
!       J_PHY   = ( IJ_PHY - 1 ) / IDIM + 1 
!       J_PHY_X = ( IJ_PHY - 1 ) / IDIM + 1 + (JL-1)*JPHY2  
!
        IMASK            (IJ_SIB,JL) = 6
        IF     ( I_PHY.EQ.3 ) THEN  
          IMASK          (IJ_SIB,JL) = 2
        ELSEIF ( I_PHY.EQ.4 ) THEN  
          IMASK          (IJ_SIB,JL) = 0
        ENDIF
!
      ENDDO
!
      DO I=1,IDIM
        GLON(I, :) = ( I - 1 ) * 360.D0 / IDIM 
      ENDDO
      DO J=1,JDIM
        GLAT(:, J) = 90. - ( J - 0.5 ) * 180.D0 / JDIM 
      ENDDO
!
      RETURN
      END 
!=================================================================
      SUBROUTINE TEST1D_RAD_SET (
     I   JL   ,   
     O   RVISB , RVISD , RNIRB , RNIRD , DLWB )
!
      USE PRM , ONLY :
     1   IJPHY , INTRI , INTRJ , JLPHY , IDIM , JDIM , JPHY2 
!
      IMPLICIT NONE 
! 
      INTEGER,INTENT(IN)  :: JL
!
      REAL(8),INTENT(OUT)  :: RNIRB ( IJPHY/INTRI/INTRJ, JLPHY )  
      REAL(8),INTENT(OUT)  :: RNIRD ( IJPHY/INTRI/INTRJ, JLPHY )  
      REAL(8),INTENT(OUT)  :: RVISB ( IJPHY/INTRI/INTRJ, JLPHY )  
      REAL(8),INTENT(OUT)  :: RVISD ( IJPHY/INTRI/INTRJ, JLPHY )  
      REAL(8),INTENT(OUT)  :: DLWB  ( IJPHY/INTRI/INTRJ, JLPHY )  
!
      INTEGER :: IJ 
      INTEGER :: IR 
      INTEGER :: JR 
!
      DO IJ = 1,IJPHY/INTRI/INTRJ  
!
        IR = MOD ( IJ - 1 , (IDIM/INTRI) ) + 1  
        JR = (IJ-1) / (IDIM/INTRI) + 1 + (JL-1)*(JPHY2/INTRJ)
!
        RVISB(IJ,JL) = 345.26647077356949
        RVISD(IJ,JL) =  86.050009414557906
        RNIRB(IJ,JL) = 356.34212471694264
        RNIRD(IJ,JL) =  50.879018867417585
        DLWB (IJ,JL) = 375.83957483659299
!
        IF (  IR.EQ.2 .or. IR.EQ.3 ) THEN
          RVISB(IJ,JL) =   8.11382597222277415E-5 
          RVISD(IJ,JL) =  37.058570797929633 
          RNIRB(IJ,JL) =   2.13785468910254071E-5
          RNIRD(IJ,JL) =  30.720071474259441
          DLWB (IJ,JL) = 250.96393562585149
!
          IF ( JR .EQ. 4 ) THEN                            ! 融かす
            RVISB(IJ,JL) =   8.11382597222277415E-5 + 40
            RVISD(IJ,JL) =  37.058570797929633 + 80
            RNIRB(IJ,JL) =   2.13785468910254071E-5 + 40
            RNIRD(IJ,JL) =  30.720071474259441 + 80
            DLWB (IJ,JL) = 250.96393562585149 + 500
          ENDIF
!
        ENDIF
      ENDDO
!
      RETURN
      END
!==============================================================
      SUBROUTINE TEST1D_ATMOS_SET (
     I  JL            , 
!    I  DELT_ATM      , 
     I  DELT_CNP      , 
     O  GPDEL_PHY     , GPHAF_PHY     , GPFUL_PHY      , 
     O  GT_PHY        , GQ_PHY        , 
     O  GU_PHY        , GV_PHY        , 
     O  ZMEAN_PHY     , ZTEMP_PHY     , 
     O  PPLI_PHY      , PPCI_PHY      , 
     O  GMT_PHY       , GMQ_PHY       , GMUV_PHY        )
!
      USE PRM , ONLY : 
     1    IJPHY , 
     1    JLPHY , 
     1    JPHY2 ,
     1    IDIM 
!
      INTEGER,INTENT(IN)    :: JL
      REAL(8),INTENT(IN)    :: DELT_CNP
!     REAL(8),INTENT(IN)    :: DELT_ATM
!
      REAL(8),INTENT(OUT)   :: GPDEL_PHY  (IJPHY,JLPHY) 
      REAL(8),INTENT(OUT)   :: GPHAF_PHY  (IJPHY,JLPHY)            
      REAL(8),INTENT(OUT)   :: GPFUL_PHY  (IJPHY,JLPHY)            
      REAL(8),INTENT(OUT)   :: GT_PHY     (IJPHY,JLPHY)
      REAL(8),INTENT(OUT)   :: GQ_PHY     (IJPHY,JLPHY) 
      REAL(8),INTENT(OUT)   :: GU_PHY     (IJPHY,JLPHY) 
      REAL(8),INTENT(OUT)   :: GV_PHY     (IJPHY,JLPHY) 
!
      REAL(8),INTENT(OUT)  :: ZTEMP_PHY   (IJPHY,JLPHY)
      REAL(8),INTENT(OUT)  :: ZMEAN_PHY   (IJPHY,JLPHY)
!
      REAL(8),INTENT(OUT)  :: PPLI_PHY    (IJPHY,JLPHY)
      REAL(8),INTENT(OUT)  :: PPCI_PHY    (IJPHY,JLPHY)
!
      REAL(8),INTENT(OUT)  :: GMT_PHY    (IJPHY,2)
      REAL(8),INTENT(OUT)  :: GMQ_PHY    (IJPHY,2)
      REAL(8),INTENT(OUT)  :: GMUV_PHY   (IJPHY,3)
!
      INTEGER :: I
      INTEGER :: J
      INTEGER :: IJ
!
      real(8) :: con
!
      DO IJ=1,IJPHY
!         
        I = MOD(IJ-1,IDIM)+1
        J = (IJ-1)/IDIM + 1 + (JL-1)*JPHY2  
!   
        GU_PHY   (IJ,JL)    =  -2.0313592261540601
        GV_PHY   (IJ,JL)    =    0.21227152179459852
        GT_PHY   (IJ,JL)    =  300.09161291534753
        GQ_PHY   (IJ,JL)    =   1.27222007742046894E-2
        GPHAF_PHY(IJ,JL)    =  984.65471491414473
        GPFUL_PHY(IJ,JL)    =  979.72732267533604
        GPDEL_PHY(IJ,JL)    =  9.8465377587423291  
!
        ZTEMP_PHY(IJ,JL)    =   0.71933462552619054
        ZMEAN_PHY(IJ,JL) =   0.75251968419095772
!
        PPLI_PHY(IJ,JL) =  0. * DELT_CNP / 3600
        PPCI_PHY(IJ,JL) =  0. * DELT_CNP / 3600 
!
        IF ( J.EQ.2 ) THEN
          PPLI_PHY(IJ,JL) =  5. * DELT_CNP / 3600
          PPCI_PHY(IJ,JL) =  4. * DELT_CNP / 3600 
        ENDIF
!
        IF ( I.EQ.2 .or. i.eq.3  ) THEN
          GU_PHY  (IJ,JL)    =  -2.29000474437234293E-2
          GV_PHY  (IJ,JL)    =  -0.61756492473476021
          GT_PHY  (IJ,JL)    =  261.46372842510897
          GQ_PHY  (IJ,JL)    =  1.37912566138209507E-3
!
          GPHAF_PHY(IJ,JL)   =  982.60539104214035
          GPFUL_PHY(IJ,JL)   =  977.68825399475281
          GPDEL_PHY(IJ,JL)   =  9.8260445395661691       
!
          ZTEMP_PHY(IJ,JL)   =   0.12281362226119813
          ZMEAN_PHY(IJ,JL)   =   0.16729917397673677
        ENDIF
      ENDDO
!
! 1 次元用モデルでは、
!      大気最下層の比熱相当は大きく、
!      第二層から第一層への拡散の入力は小さくしておく。 
!
      con = 1.D6                                       ! 大きい値 
      call con_set (  GMT_PHY (1,1)  , ijphy , con )
      call con_set (  GMQ_PHY (1,1)  , ijphy , con )
      call con_set (  GMUV_PHY(1,1)  , ijphy , con )
!
      con = 1.D-6 
      call con_set (  GMT_PHY (1,2)  , ijphy   , con )
      call con_set (  GMQ_PHY (1,2)  , ijphy   , con )
      call con_set (  GMUV_PHY(1,2)  , ijphy*2 , con )
!
#ifdef DDDD
      do i=1,ijphy
        GMT_PHY(i,1)      =  2.580905791152716
        GMT_PHY(i,2)      =  0.
        GMQ_PHY(i,1)      =  2.580905791152716 
        GMQ_PHY(i,2)      =  0.
      enddo
#endif 
!
      RETURN
      END








