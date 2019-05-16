SUBROUTINE monit_add_2   &
( akwdi, datak, jl, delt, daysec) ! Input

  USE prm,ONLY : &
    ijmax
  USE com_step,ONLY : &
    icnmntmon
  USE com_monit,ONLY : &
    mxptr_2, &
    ikwd_2,  &
    ckwd_2,  &
    mxptr_nouse_2, &
    ikwd_nouse_2,  &
    ckwd_nouse_2,  &
    ldisk_2,  &
    ctflag_2, &
    sum_2,    &
    lsdisk_2

  IMPLICIT NONE

  CHARACTER(LEN=7),INTENT(IN) :: akwdi  !-- keyword 
  CHARACTER(LEN=8)   :: akwd_in         !-- keyword 
  integer(8)         :: ikwd_in
!
  equivalence ( akwd_in , ikwd_in )
!
  REAL(8),INTENT(IN) :: datak(ijmax)    !-- data
  INTEGER,INTENT(IN) :: jl              !-- band number
  REAL(8),INTENT(IN) :: delt            !-- time interval
  REAL(8),INTENT(IN) :: daysec          !-- * daysec

  REAL(8) :: delttt

  INTEGER :: iptr  !-- pointer for 2D monitor data.
  INTEGER :: i, ij

! write(6,*) 'monit_add_2 start' 

  akwd_in = akwdi 
!-----------      start !      -----------

 !-- find iptr which corresponds to akwdi.
  iptr = -999
  DO i = 1, mxptr_2
    IF ( ikwd_in == ikwd_2(i) ) THEN
      iptr = i
      EXIT             ! non-parallelized, non-vectorized
    END IF
  END DO
  IF ( iptr == -999 ) THEN
    DO i = 1, mxptr_nouse_2
      IF ( ikwd_in == ikwd_nouse_2(i) ) THEN
        RETURN
      END IF
    END DO
    WRITE(6,*) 'monit_add_2 : I cannot find the keyword data:',akwdi
    STOP 999
  END IF

  IF ( ldisk_2(iptr) == 0 .AND. lsdisk_2(iptr) == 0  ) THEN
   !-- not output
    RETURN
  END IF

  IF ( ctflag_2(iptr) == 'SNP' .AND. icnmntmon /= 1 ) THEN
    RETURN
  END IF

!-----    WEIGHT    -----
!HYHY  IF ( ctflag_2(iptr) /= 'SNP' ) THEN    !--'AVR'
    delttt = delt * daysec
!HYHY  ELSE
!HYHY    delttt = daysec
!HYHY  END IF

!-----    ADD   -----

! write(6,*) 'monit_add_2 ' , akwdi , delttt , datak(1)
!
  DO ij = 1,ijmax
    sum_2( ij, jl, iptr ) =   &
        sum_2( ij, jl, iptr ) + delttt*datak(ij)
  END DO


END SUBROUTINE monit_add_2
