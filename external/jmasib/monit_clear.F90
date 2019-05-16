SUBROUTINE monit_clear(totsec)

  USE com_monit,ONLY : &
    sum_2,     &
    sum_3,     &
    sum_3eta,  &
    snp_wmsk,  &
    snp_wmskp, &
    sum_wmsk,  &
    sum_wmskp

  IMPLICIT NONE

  REAL(8),INTENT(OUT) :: totsec

  totsec   = 0.0D0

  sum_2   (:,:,:)   = 0.0D0
  sum_3   (:,:,:,:) = 0.0D0
  sum_3eta(:,:,:,:) = 0.0D0

!HY  snp_wmsk (:,:,:) = 0.0D0
!HY  snp_wmskp(:,:,:) = 0.0D0
  sum_wmsk (:,:,:) = 0.0D0
  sum_wmskp(:,:,:) = 0.0D0

END SUBROUTINE monit_clear
