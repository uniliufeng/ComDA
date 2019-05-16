! forcing_main.F90 - 強制力を各ステップごとに与える
! vi: set sw=2 ts=72:

!=======================================================================
! forcing_main - 強制力を各ステップごとに与える

subroutine forcing_main(& 
  id_now_in, id_pre_in, id_now5, id_pre5, &
  delt_atm, rday, rsec, imask, glon, glat, &! IN
  u_phy, v_phy, pd_phy, ps_phy, pf_phy, tmp_phy, q_phy, &! OUT
  zmean_phy, ztemp_phy, ppli_phy, ppci_phy, &! OUT
  rvisb, rvisd, rnirb, rnird, dlwb)  ! OUT
  !
  use forcing, only: IDIM, JDIM, forcing_read_id, forcing_read_shortwave, &
    forcing_read_nearest_id, cur_file
  use date, only: date_diff_div, date_diff
  !
  implicit none
  !
  ! input
  !
  integer, intent(in):: id_now_in(5)  ! 本ステップの時刻
  integer, intent(in):: id_pre_in(5)  ! 本ステップの時刻
  integer, intent(in):: id_now5
  integer, intent(in):: id_pre5
  integer :: id_now(5)  ! 本ステップの時刻
  integer :: id_pre(5)  ! 本ステップの時刻

  real(8), intent(in):: rsec 
  real(8), intent(in):: rday 
  real(8), intent(in):: delt_atm
  !
  real(8), intent(in):: glon(IDIM, JDIM)  ! 経度（単位、度）
  real(8), intent(in):: glat(IDIM, JDIM)  ! 緯度（単位、度）
  integer, intent(in):: imask(IDIM, JDIM)     
  !
  ! output 
  !
  real(8), intent(out):: u_phy(IDIM, JDIM)  ! u
  real(8), intent(out):: v_phy(IDIM, JDIM)  ! v 
  real(8), intent(out):: pd_phy(IDIM, JDIM)  !(ps-ph)*2  
  real(8), intent(out):: ps_phy(IDIM, JDIM)  ! ハーフレベル = Ps HPa
  real(8), intent(out):: pf_phy(IDIM, JDIM)  ! フルレベル HPa 
  real(8), intent(out):: tmp_phy(IDIM, JDIM)  ! 温度 
  real(8), intent(out):: q_phy(IDIM, JDIM)  ! 比湿 kg/kg
  real(8), intent(out):: zmean_phy(IDIM, JDIM)  ! 一時間?平均天頂角
  real(8), intent(out):: ztemp_phy(IDIM, JDIM)  ! 各ステップ天頂角
  real(8), intent(out):: ppli_phy(IDIM, JDIM)  ! 大規模凝結性降水
  real(8), intent(out):: ppci_phy(IDIM, JDIM)  ! 積雲性降水
  !
  integer, parameter:: IJRAD = IDIM * JDIM
  real(8), intent(out):: rvisb(IJRAD)  ! 放射計算時可視直達
  real(8), intent(out):: rvisd(IJRAD)  ! 放射計算時可視散乱
  real(8), intent(out):: rnirb(IJRAD)  ! 放射計算時近赤直達
  real(8), intent(out):: rnird(IJRAD)  ! 放射計算時近赤散乱
  real(8), intent(out):: dlwb(IDIM, JDIM)  ! 放射計算時長波
  !
  real(8):: tprc(IDIM, JDIM)
  real(8):: cprc(IDIM, JDIM)
  real(8):: swdn(IDIM, JDIM)
  real(8):: cloudiness(IDIM, JDIM)
  real(8):: delt_tprc  ! 降水データの時間間隔
  real(8), save:: delt_tprc_last = 3600.0
  integer:: idiff_tprc(5)  ! 降水データの時間間隔
!
  id_now(:) = id_now_in(:)
  id_pre(:) = id_pre_in(:)
  id_now(5) = id_now5
  id_pre(5) = id_pre5
  !
  ! とりあえず、全部 id_now で読んでみる
  !
#ifdef DEBUG
  print "(' ',a,i4.4,4('-',i2.2))", 'forcing_main', id_now
#endif
  call forcing_read_id('MWND', u_phy, id_now)
  v_phy = 0.0
  call forcing_read_id('TEMP', tmp_phy, id_now)
  call forcing_read_id('QREF', q_phy, id_now)
  ! Pa 単位なので hPa 単位に変換
  call forcing_read_id('PRSS', ps_phy, id_now)
  ps_phy = ps_phy / 100.0_8
  ! 時間間隔はファイルから読み取り [mm/ステップ] に変換
  call forcing_read_nearest_id('CPRC', cprc, id_now, 0.0_8)
  call forcing_read_nearest_id('TPRC', tprc, id_now, 0.0_8)
  call date_diff(cur_file%last%id, cur_file%before%id, idiff_tprc)
  call date_diff_div(idiff_tprc, (/0, 0, 0, 0, 1/), delt_tprc)
  if (delt_tprc == 0.0) then
    delt_tprc = delt_tprc_last
  else
    delt_tprc_last = delt_tprc
  endif
  ppci_phy(:, :) = cprc(:, :) * delt_atm / delt_tprc
  ppli_phy(:, :) = (tprc - cprc) * delt_atm / delt_tprc

  call forcing_read_id('LWDN', dlwb, id_now)
  call forcing_read_shortwave( &
    'SWDN', 'CLD ', id_now, delt_atm, glon, glat, imask, &! IN
    swdn, rvisb, rvisd, rnirb, rnird, zmean_phy, ztemp_phy, cloudiness &! OUT
  )
  pd_phy(:, :) = 10.0_8
  pf_phy(:, :) = ps_phy(:, :) - 5.0

  !
  ! モニタ出力
  !

#ifdef MONYOS
  call forcing_monitor
#endif
  !
  ! お別れ
  !
#ifdef DEBUG
  print "(' ',a,i4.4,4('-',i2.2))", 'forcing_main ended ', id_now
#endif
  return

#ifdef MONYOS
contains

  subroutine forcing_monitor
    use prm, only: JLPHY, IJPHY
    use com_runconf, only: JCNIMNT
    use com_runconf_sib0109, only: JCN_IWL_SKIP
    use com_step, only: ICNSW
    use sib_monit, only : imonit_level , imonit_select 

    real(8), parameter:: UNITY = 1.0_8
    real(8):: precip_time = 0.0
    character(len = 7):: cmark 
    integer:: jl
    integer:: nsteps_rad
    integer:: istart

    if (abs(JCN_IWL_SKIP) == 1) then   
      ! 毎ステップ放射計算
      nsteps_rad = 1 
    else                                            ! 一時間に一度
      nsteps_rad = 3601 / DELT_ATM 
    endif

    if (precip_time == 0.0) then
      if (JCNIMNT >= 900) then
        precip_time = UNITY
      else if (JCNIMNT <= 0) then
        precip_time = UNITY
      else
        precip_time = UNITY
      endif
    endif

    if (imonit_level .ge. imonit_select) then
      do, jl = 1, JLPHY
        istart = (jl - 1) * IJPHY + 1              ! ポインタ。
        if (ICNSW == 1) then
          ! 全短波
          cmark = 'FSR' 
          call monit_add_2(cmark, swdn(istart, 1), &
            jl, DELT_ATM * NSTEPS_RAD, UNITY)     
          ! 長波
          cmark = 'FLR' 
          call monit_add_2(cmark, dlwb(istart, 1), &
            jl, DELT_ATM * NSTEPS_RAD, UNITY)     
        endif

        ! 雲量
        cmark = 'FCLD' 
        call monit_add_2(cmark, cloudiness(istart, 1), jl, DELT_ATM, UNITY)

        ! 気温
        cmark = 'FTMP' 
        call monit_add_2(cmark, tmp_phy(istart, 1), jl, DELT_ATM, UNITY)

        ! 比湿
        cmark = 'FQ' 
        call monit_add_2(cmark, q_phy(istart, 1), jl, DELT_ATM, UNITY)

        ! 地表面気圧
        cmark = 'FPS' 
        call monit_add_2(cmark, ps_phy(istart, 1), jl, DELT_ATM, UNITY)

        ! 積雲性降水
        cmark = 'FPC' 
        call monit_add_2(cmark, ppci_phy(istart, 1), jl, precip_time, UNITY)

        ! 積雲性降水
        cmark = 'FPL' 
        call monit_add_2(cmark, ppli_phy(istart, 1), jl, precip_time, UNITY)

        ! 風速
        cmark = 'FU'
        call monit_add_2(cmark, u_phy(istart, 1), jl, DELT_ATM, UNITY)

      enddo
    endif

  end subroutine
#endif

end subroutine

!=======================================================================
!
! forcing_main で出力する変数リストの登録
!

subroutine monit_regist_forcing
!
  use sib_monit, only : imonit_level , imonit_select 
!
  implicit none
  !
  character(len = 7):: cmark 
  character(len = 32):: ctitle
  character(len = 13):: cunit
  logical, save:: lfirst = .TRUE.
!
  if (lfirst) then
    write(6, *) "monit_regist_forcing 2000-08-04 toyoda"
    lfirst = .FALSE.
  endif
  !
  ! 各強制変数の登録
  !

  if (imonit_level .ge. imonit_select) then
!
  cmark = 'FSR'
  ctitle = 'shortwave radiation (atmospheric forcing)'
  cunit = 'W/m**2'
  call monit_regist_sib(cmark, ctitle, cunit)

  cmark = 'FLR'
  ctitle = 'longwave radiation (atmospheric forcing)'
  cunit = 'W/m**2'
  call monit_regist_sib(cmark, ctitle, cunit)

  cmark = 'FTMP'
  ctitle = 'temperature (atmospheric forcing)'
  cunit = 'K'
  call monit_regist_sib(cmark, ctitle, cunit)

  cmark = 'FQ'
  ctitle = 'watar vapor mixing ratio (atmospheric forcing)'
  cunit = 'kg/kg'
  call monit_regist_sib(cmark, ctitle, cunit)

  cmark = 'FPS'
  ctitle = 'surface pressure (atmospheric forcing)'
  cunit = 'hPa'
  call monit_regist_sib(cmark, ctitle, cunit)

  cmark = 'FPC'
  ctitle = 'convective precipitation (atmospheric forcing)'
  cunit = 'kg/m**2/s'
  call monit_regist_sib(cmark, ctitle, cunit)

  cmark = 'FPL'
  ctitle = 'large scale precipitation (atmospheric forcing)'
  cunit = 'kg/m**2/s'
  call monit_regist_sib(cmark, ctitle, cunit)

  cmark = 'FCLD'
  ctitle = 'cloud amount (atmospheric forcing)'
  cunit = ''
  call monit_regist_sib(cmark, ctitle, cunit)

  cmark = 'FU'
  ctitle = 'wind speed (atmospheric forcing)'
  cunit = 'm/s'
  call monit_regist_sib(cmark, ctitle, cunit)
!
  endif

end subroutine
