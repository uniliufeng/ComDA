! forcing_init_sib.F90 - 陸面の初期値を作成
! vi: set ts=72 sw=2:

subroutine forcing_init_sib( &
  tmp_cnp_nos_phy, fr_wtr_cnp_nos_phy, fr_gla_cnp_nos_phy, &! out
  tmp_cnp_snw_phy, fr_wtr_cnp_snw_phy, fr_gla_cnp_snw_phy, &! out
  tmp_grsk_phy, fr_wtr_grs_phy, fr_gla_grs_phy, &! out
  tmp_snsk_phy, &! out
  info_snow_phy, &! out
  tmp_soil_phy, fr_wtr_soil_phy, fr_gla_soil_phy, &! out
  tmp_snsl_phy, fr_wtr_snsl_phy, fr_gla_snsl_phy, &! out
  tmp_snow_phy, wtr_snow_phy, gla_snow_phy, &! out
  rho_snow_inv_phy, age_snow_phy, &! out
  eng_snow_bucket_phy, h2o_snow_bucket_phy) ! out
!
  use forcing, only: IDIM, JDIM, FORCING_READ_ID
  use sibprm, only: IDP, ISN
  use sibcon, only: TMP_FREZ_c
  use com_jobinfo_sib0109, only: IDSTAR
  implicit none
  integer, parameter:: DOUBLE = kind(0.0d0)
!
! OUTPUT
!
  real(DOUBLE), intent(out), dimension(IDIM, JDIM):: &
    tmp_cnp_nos_phy, fr_wtr_cnp_nos_phy, fr_gla_cnp_nos_phy, &
    tmp_cnp_snw_phy, fr_wtr_cnp_snw_phy, fr_gla_cnp_snw_phy, &
    tmp_grsk_phy, fr_wtr_grs_phy, fr_gla_grs_phy, &
    tmp_snsk_phy
  integer, intent(out), dimension(IDIM, JDIM):: &
    info_snow_phy
  real(DOUBLE), intent(out), dimension(IDIM, JDIM, IDP):: &
    tmp_soil_phy, fr_wtr_soil_phy, fr_gla_soil_phy
  real(DOUBLE), intent(out), dimension(IDIM, JDIM):: &
    tmp_snsl_phy, fr_wtr_snsl_phy, fr_gla_snsl_phy
  real(DOUBLE), intent(out), dimension(IDIM, JDIM, ISN):: &
    tmp_snow_phy, wtr_snow_phy, gla_snow_phy, rho_snow_inv_phy
  real(DOUBLE), intent(out), dimension(IDIM, JDIM):: &
    age_snow_phy, eng_snow_bucket_phy, h2o_snow_bucket_phy
  !
  real(DOUBLE):: tmp_f
  real(DOUBLE):: tmp(IDIM, JDIM)
  character(len = 4):: cmark
  integer:: id(5)
  !
  ! 温度ファイルの最初の記録を読み、氷結地とそうでないところを決定する
  !
  cmark = "TEMP"
  id(:) = idstar(:)
  call forcing_read_id(cmark, tmp, id)
  !
  tmp_cnp_nos_phy(:, :) = tmp(:, :)
  tmp_f = TMP_FREZ_c + 1.0
  where (tmp >= tmp_f)
    tmp_grsk_phy(:, :) = tmp + 3
    tmp_soil_phy(:, :, 1) = tmp + 2
    tmp_soil_phy(:, :, 2) = tmp + 1
    fr_wtr_soil_phy(:, :, 1) = 0.5
    fr_wtr_soil_phy(:, :, 2) = 0.5
    fr_wtr_soil_phy(:, :, 3) = 0.5
    fr_gla_soil_phy(:, :, 1) = 0.0
    fr_gla_soil_phy(:, :, 2) = 0.0
    fr_gla_soil_phy(:, :, 3) = 0.0
  elsewhere
    tmp_grsk_phy(:, :) = tmp - 3
    tmp_soil_phy(:, :, 1) = tmp - 2
    tmp_soil_phy(:, :, 2) = tmp - 1
    fr_wtr_soil_phy(:, :, 1) = 0.0
    fr_wtr_soil_phy(:, :, 2) = 0.0
    fr_wtr_soil_phy(:, :, 3) = 0.0
    fr_gla_soil_phy(:, :, 1) = 0.5
    fr_gla_soil_phy(:, :, 2) = 0.5
    fr_gla_soil_phy(:, :, 3) = 0.5
  end where
  tmp_soil_phy(:, :, 3) = tmp
  !
  fr_wtr_cnp_nos_phy(:, :) = 0.0
  fr_wtr_grs_phy(:, :) = 0.0
  fr_gla_cnp_nos_phy(:, :) = 0.0
  fr_gla_grs_phy(:, :) = 0.0
  !
  tmp_cnp_snw_phy(:, :) = 0.0
  fr_wtr_cnp_snw_phy(:, :) = 0.0
  fr_gla_cnp_snw_phy(:, :) = 0.0
  !
  tmp_snsk_phy(:, :) = 0.0
  tmp_snsl_phy(:, :) = 0.0
  fr_wtr_snsl_phy(:, :) = 0.0
  fr_gla_snsl_phy(:, :) = 0.0
  !
  age_snow_phy(:, :) = 10 * 24 * 3600
  info_snow_phy(:, :) = -1
  tmp_snow_phy(:, :, :) = 0.0
  wtr_snow_phy(:, :, :) = 0.0
  gla_snow_phy(:, :, :) = 0.0
  rho_snow_inv_phy(:, :, :) = 0.0
  !
  eng_snow_bucket_phy(:, :) = 0.0
  h2o_snow_bucket_phy(:, :) = 0.0

end subroutine
