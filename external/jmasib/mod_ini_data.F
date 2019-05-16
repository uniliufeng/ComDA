      SUBROUTINE MOD_INI_DATA (
     I      NMAX_01 , NMAX , MON  , NTYP ,            
     M      FR_WTR_CNP_NOS , FR_WTR_GRS       , FR_WTR_CNP_SNW , 
     M      FR_GLA_CNP_NOS , FR_GLA_GRS       , FR_GLA_CNP_SNW  ) 
!
      USE SIBPRM, ONLY : IDX
      USE SIBCON, ONLY : H2O_SAT_CNP_DEP_c , H2O_SAT_GRS_DEP_c  
      USE SIBVEG, ONLY : XLAI_mo
!
      INTEGER,INTENT(IN)   :: NTYP (IDX) 
      INTEGER,INTENT(IN)   :: NMAX_01(3) 
      INTEGER,INTENT(IN)   :: NMAX
      INTEGER,INTENT(IN)   :: MON
!
      REAL(8),INTENT(INOUT):: FR_WTR_CNP_NOS (IDX)
      REAL(8),INTENT(INOUT):: FR_GLA_CNP_NOS (IDX)
      REAL(8),INTENT(INOUT):: FR_WTR_CNP_SNW (IDX)
      REAL(8),INTENT(INOUT):: FR_GLA_CNP_SNW (IDX)
      REAL(8),INTENT(INOUT):: FR_WTR_GRS     (IDX)
      REAL(8),INTENT(INOUT):: FR_GLA_GRS     (IDX)
!
      INTEGER :: I
!
      REAL(8),PARAMETER :: D08 = 1.D0 
!
      DO I =1,NMAX_01(3)
        FR_WTR_CNP_NOS(I) = FR_WTR_CNP_NOS(I)                     
     1        / ( XLAI_mo   (NTYP(I),MON,1) * H2O_SAT_CNP_DEP_c ) 
        FR_WTR_GRS    (I) = FR_WTR_GRS    (I) 
     1        / ( XLAI_mo   (NTYP(I),MON,2) * H2O_SAT_GRS_DEP_c ) 
        FR_GLA_CNP_NOS(I) = FR_GLA_CNP_NOS(I) 
     1        / ( XLAI_mo   (NTYP(I),MON,1) * H2O_SAT_CNP_DEP_c ) 
        FR_GLA_GRS    (I) = FR_GLA_GRS    (I) 
     1        / ( XLAI_mo   (NTYP(I),MON,2) * H2O_SAT_GRS_DEP_c ) 
!
        FR_WTR_CNP_NOS(I) = MIN ( FR_WTR_CNP_NOS(I) , D08 )
        FR_WTR_GRS    (I) = MIN ( FR_WTR_GRS    (I) , D08 )
        FR_GLA_CNP_NOS(I) = MIN ( FR_GLA_CNP_NOS(I) , D08 )
        FR_GLA_GRS    (I) = MIN ( FR_GLA_GRS    (I) , D08 )
!
        write(6,*) 'fr_wtr_cnp mod ' , I , 
     1         FR_WTR_CNP_NOS(I) , 
     1         FR_GLA_CNP_NOS(I) , 
     1         XLAI_mo   (NTYP(I),MON,1) * H2O_SAT_CNP_DEP_c  
        write(6,*) 'fr_wtr_grs mod ' , I , 
     1         FR_WTR_GRS    (I) , 
     1         FR_GLA_GRS    (I) , 
     1         XLAI_mo   (NTYP(I),MON,2) * H2O_SAT_GRS_DEP_c  
      ENDDO
!
      DO I =NMAX_01(2),NMAX
        FR_WTR_CNP_SNW(I) = FR_WTR_CNP_SNW(I) 
     1       / ( XLAI_mo   (NTYP(I),MON,1) * H2O_SAT_CNP_DEP_c ) 
        FR_GLA_CNP_SNW(I) = FR_GLA_CNP_SNW(I) 
     1       / ( XLAI_mo   (NTYP(I),MON,1) * H2O_SAT_CNP_DEP_c ) 
!
        FR_WTR_CNP_SNW(I) = MIN ( FR_WTR_CNP_SNW(I) , D08 )
        FR_GLA_CNP_SNW(I) = MIN ( FR_GLA_CNP_SNW(I) , D08 )
!
        write(6,*) 'fr_wtr_cnp mod ' , I , 
     1              FR_WTR_CNP_SNW(I) , 
     1              FR_GLA_CNP_SNW(I) ,  
     1              XLAI_mo   (NTYP(I),MON,1) * H2O_SAT_CNP_DEP_c 
!
 
      ENDDO
!
      RETURN
      END 
