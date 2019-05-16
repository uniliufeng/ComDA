SUBROUTINE monit_mask           &!
( akwdi, jl, delt )             !  In

  USE prm
  USE com_step,ONLY : &
    icnmntmon
  USE com_monit

  IMPLICIT NONE

  CHARACTER(LEN=7),INTENT(IN) :: akwdi
  REAL(8)         ,INTENT(IN) :: delt
  INTEGER         ,INTENT(IN) :: jl

  REAL(8) :: wmskadd
  INTEGER :: i,iptr,ka,ij

!------------       start !       ------------------

  IF ( akwdi == 'WMSK   ') THEN

    DO ka=1,kmaxm
      DO ij=1,ijmax
        wmskadd = delt*wmsk(ij,ka)
        sum_wmsk(ij,jl,ka)     = sum_wmsk(ij,jl,ka)     + wmskadd
        sum_wmsk_day(ij,jl,ka) = sum_wmsk_day(ij,jl,ka) + wmskadd
        sum_wmsk_6hr(ij,jl,ka) = sum_wmsk_6hr(ij,jl,ka) + wmskadd
      END DO
    END DO

    DO ka=1,kmaxm
      DO ij=1,ijmax
        snp_wmsk(ij,jl,ka) = delt*wmsk(ij,ka)
      END DO
    END DO

  ELSE   !-----   IF ( akwdi == 'WMSK   ')

    DO ka=1,kmaxm
      DO ij=1,ijmax
        wmskadd = delt*wmskp(ij,ka)
        sum_wmskp(ij,jl,ka)     = sum_wmskp(ij,jl,ka)     + wmskadd
        sum_wmskp_day(ij,jl,ka) = sum_wmskp_day(ij,jl,ka) + wmskadd
        sum_wmskp_6hr(ij,jl,ka) = sum_wmskp_6hr(ij,jl,ka) + wmskadd
      END DO
    END DO

    DO ka=1,kmaxm
      DO ij=1,ijmax
        snp_wmskp(ij,jl,ka) = delt*wmskp(ij,ka)
      END DO
    END DO

  END IF   !-----   IF ( akwdi == 'WMSK   ')

END SUBROUTINE monit_mask
