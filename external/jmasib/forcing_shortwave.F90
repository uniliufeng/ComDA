! forcing_shortwave.F90 - 短波放射の補間
! vi: set sw=2 ts=72:

!=======================================================================
! forcing_read_shortwave - 短波放射用内挿読み取りルーチン
!   forcing_read_id の短波放射用類似製品だったのだが仕様が爆発。
!   時刻 id 以前でもっとも遅いデータ時刻 id1 を求め、その次のデータ時刻
!   id2 との間で天頂角余弦で重み付け積分をしてから案分値を返す。
!   delt_atm, glon, glat, imask は天頂角余弦の計算に用いる。

subroutine forcing_read_shortwave( &
  cmark_swdn, cmark_cld, id, delt_atm, glon, glat, imask, &! IN
  swdn, rvisb, rvisd, rnirb, rnird, zmean_phy, ztemp_phy, cloudiness &! OUT
)
  use forcing, only: IDIM, JDIM, forcing_read_id2, forcing_read_nearest_id
  use com_step, only: icnsw
  use date
  implicit none
  character(len = 4), intent(in):: cmark_swdn
  character(len = 4), intent(in):: cmark_cld
  integer, intent(in):: id(5)
  real(8), intent(in):: delt_atm
  real(8), intent(in):: glon(IDIM, JDIM)
  real(8), intent(in):: glat(IDIM, JDIM)
  integer, intent(in):: imask(IDIM, JDIM)
  real(8), intent(out):: swdn(IDIM, JDIM)
  real(8), intent(out):: rvisb(IDIM, JDIM)
  real(8), intent(out):: rvisd(IDIM, JDIM)
  real(8), intent(out):: rnirb(IDIM, JDIM)
  real(8), intent(out):: rnird(IDIM, JDIM)
  real(8), intent(out):: zmean_phy(IDIM, JDIM)
  real(8), intent(out):: ztemp_phy(IDIM, JDIM)
  real(8), intent(out):: cloudiness(IDIM, JDIM)  ! 雲量
  real:: buf1(IDIM, JDIM)  ! 時刻の早いほうのデータ
  real:: buf2(IDIM, JDIM)  ! 時刻の遅いほうのデータ
  real(8):: weight  ! ここでは使わない
  integer:: id1(5)  ! データの早いほうの時刻
  integer:: id2(5)  ! データの遅いほうの時刻
  real(8), save:: delt_data  ! 補間すべき2つのデータの時刻の間隔 [秒]
  integer:: idiff(5)  ! 時刻の間隔
  integer:: id_base(5)  ! 時刻の起点
  integer:: idiff_span(5)  ! 時刻の間隔
  integer:: id_afterspan(5)  ! 時刻の終点
  real(8):: rday  ! 年初からの時間 [年]
  real(8):: rsec  ! 正子からの時間 [日]
  logical:: update  ! 読み取りによって buf1, buf2 が更新された場合真
  ! 天頂角余弦で重みをつけたステップごと短波入射量
  real(8), save:: sr_flux_work(IDIM, JDIM)
  ! sr_flux_work を計算した時刻。初期値はありえない日付
  integer, save:: id_sr(5) = (/0, 0, 0, 0, 0/)
  logical, save:: lfirst = .TRUE.
  integer:: nsteps_rad
! ---
  !
  ! 時刻 id 以前でもっとも遅いデータ時刻 id1 とその次のデータ時刻 id2 を
  ! 求める。id1 が新しくなったときに update が真になり、sr_flux_work が
  ! 更新される。計算時刻とデータ時刻が重なったときには更新すべきなので、
  ! id + 1[sec] の前後で検索する。データ間隔は2秒以上であることを仮定。
  !
  id_afterspan(:) = id(:)
  id_afterspan(5) = id_afterspan(5) + 1
  call forcing_read_id2(cmark_swdn, id_afterspan, &! IN
    id1, buf1, id2, buf2, weight, update) ! OUT
  if (all(id1 == id2)) update = .FALSE.
  !
  ! プログラム起動時にビンゴしてしまったら救済措置
  ! (正時 - 1sec で起動しないかぎり不要)
  !
  if (lfirst .and. .not. update) then
    call forcing_read_id2(cmark_swdn, id, &! IN
      id1, buf1, id2, buf2, weight, update) ! OUT
  endif
  if (lfirst) lfirst = .FALSE.
  if (update) then
    ! つまり idiff = id2 - id1
    call date_diff(id2, id1, idiff)
    ! つまり delt_data = idiff / 1[sec]
    call date_diff_div(idiff, (/0, 0, 0, 0, 1/), delt_data)
  endif
  if (delt_data == 0.0_8) delt_data = 3600.0_8 * 6.0_8  ! 安全策

  ! id1 から rday を求める
  ! つまり rday = (id1 - id1#yr#Jan1) / ((id1#yr + 1)#Jan1 - id1#yr#Jan1)
  id_base(:) = id1(:)
  id_base(2:5) = (/1, 1, 0, 0/)
  call date_diff(id1, id_base, idiff)
  id_afterspan(:) = id_base(:)
  id_afterspan(1) = id_afterspan(1) + 1
  call date_diff(id_afterspan, id_base, idiff_span)
  call date_diff_div(idiff, idiff_span, rday)

  ! id から rsec を求める
  ! つまり rsec = (id1 - id1#day#0am) / ((id1#day + 1)#0am - id1#day#0am)
  id_base(:) = id(:)
  call date_normalize(id_base)
  rsec = (id_base(4) * 3600.0_8 + real(id_base(5), kind=8)) / 86400.0_8
#ifdef DEBUG
  print *, 'forcing_shortwave_step', id_base, rsec
#endif

  ! 新しいデータ時刻に達したならばステップ毎平均放射を求める
  if (update) then
    call forcing_shortwave_average( &
      rday, rsec, delt_atm, delt_data, glon, glat, imask, &! IN
      buf1, &! INOUT
      sr_flux_work) ! OUT
  endif

  ! icnsw が真になったステップ毎に放射を計算する
  if (icnsw == 1) then
    call forcing_shortwave_step( &
      rday, rsec, delt_atm, sr_flux_work, glon, glat, &! IN
      nsteps_rad, swdn, zmean_phy) ! OUT
    call forcing_read_nearest_id(cmark_cld, cloudiness, id, 0.5_8)
    cloudiness(:, :) = cloudiness(:, :) / 100.0_8
    call forcing_shortwave_divide( &
      rvisb, rvisd, rnirb, rnird, &! OUT
      swdn, cloudiness, zmean_phy) ! IN
  endif
  call islscp_sunang(rday,  rsec,  glon,  glat, &! IN
    ztemp_phy) ! OUT
  ztemp_phy(:, :) = max(ztemp_phy(:, :), 0.0_8)

end subroutine

!=======================================================================
! forcing_shortwave_average - 短波入射量を天頂角で案分するための準備
!   短波入射量 rswd を時間平均余弦天頂角で割り、sr_flux_work を
!   作成する。この値によりあとで余弦天頂角 cosz をかけてフラックスを計算する。
!   積分時間は季節 rday の時刻 rsec から delt_data 秒である。

subroutine forcing_shortwave_average( &
  rday, rsec, delt_atm, delt_data, &! IN
  glon, glat, imask, &! IN
  rswd, &! INOUT
  sr_flux_work &! OUT
)
  ! cosz_min は正の値を持つことを仮定している
  use sibcon, only: cosz_min_c
  use prm, only: idim, jdim 
  use message, only: message_put
  !
  implicit none
  !
  real, intent(inout):: rswd(idim, jdim)
  real(8), intent(in):: rday 
  real(8), intent(in):: rsec
  real(8), intent(in):: delt_atm
  real(8), intent(in):: delt_data
  integer, intent(in):: imask(idim, jdim)
  real(8), intent(in):: glon(idim, jdim)
  real(8), intent(in):: glat(idim, jdim)
  !
  real(8), intent(out):: sr_flux_work(idim, jdim) 
  !
  integer:: i
  integer:: j
  integer:: ist
  integer:: nsteps_6hr 
  real(8):: dsec 
  !
  real(8):: cosz        (idim, jdim) 
  real(8):: cosz_sum    (idim, jdim) 
  !
  ! ----------------------
  ! > delt_data 秒間のステップ数 < 
  ! ----------------------
  !
  nsteps_6hr = (delt_data + 1.0) / delt_atm 
  !
  ! ------------------
  ! > データチェック <
  ! ------------------
  !   -10 以下があったらはねる。-10<rswd<0 ならばメッセージを吐いて 0 にする。
  !
  do j=1, jdim
    do i=1, idim
!     if (rswd(i, j) < -1.0) then 
      if (rswd(i, j) < -10.0) then 
        write(6, *) 'rswd forcing_shortwave_average severe warning',  &
          i,  j,  rswd(i, j), 'is modified to 0'
        rswd(i, j) = 0.0
        stop 999
      elseif (rswd(i, j) < 0.0) then 
        write(6, *) 'rswd forcing_shortwave_average warning', &
          i,  j,  rswd(i, j), 'is modified to 0'
        rswd(i, j) = 0.0
      endif
    enddo 
  enddo 
  !
  ! --------------------------------------------------
  ! > delt_data 秒間の間の各ステップでの天頂角の計算、同積算 <
  ! --------------------------------------------------
  !
  cosz_sum (:,:) = 0.d0
  !
  do ist = 1, nsteps_6hr
    dsec = rsec + delt_atm / 86400.d0 * ( ist - 0.5d0 ) 
    call islscp_sunang(rday,  dsec,  glon,  glat, &! IN
      cosz) ! OUT
    where (cosz >= cosz_min_c)
      cosz_sum = cosz_sum + cosz
    end where
  enddo
  !
  !   (3) データチェック
  !
  do j=1, jdim
    do i=1, idim
  !
  !         短波があるのに天頂角の和が負ならば、警告、短波はゼロに
      if (rswd(i, j) > 0.0 .and. cosz_sum(i, j) < cosz_min_c &
        .and. imask(i, j) /= 0) then
        call message_put("forcing_shortwave: rswd>0 cosz=0")
#       if defined(DEBUG) || defined(CHECK)
          write(6, *) i,  j,  rswd(i, j),  cosz_sum(i, j) 
#       endif
        rswd(i, j) = 0.
  !
  !         短波がないのに天頂角の和が正ならば、警告。
      elseif (rswd (i, j) <= 0.0 .and. cosz_sum(i, j) >= cosz_min_c &
        .and. imask(i, j) /= 0) then
        call message_put("forcing_shortwave: rswd=0 cosz>0")
#       if defined(DEBUG) || defined(CHECK)
          write(6, *) i,  j,  rswd(i, j),  cosz_sum(i, j) 
#       endif
      endif
    enddo
  enddo
  !
  !   (4) s/cosΘ の計算
  !       sr_flux_work * cosz がフラックスになるようにする。
  !
  !       cosz_sum があまりに小さいときは、短波を零にしてしまう。
  !       害は無かろう。多分。
  !
  !       rswd は日が沈んでいる間も含めた delt_data 秒間の平均。
  !       よって、短波総計を rswd*nsteps_6hr で求めておく。
  !       cosz/cosz_sum は、そのうちの配分比を表す。
  !
  sr_flux_work(:,:) = 0.0d0 
  where (cosz_sum >= cosz_min_c)
    sr_flux_work = rswd * nsteps_6hr / cosz_sum
  end where
  !
  return
end subroutine

!================================================================
subroutine forcing_shortwave_step( &
  rday, rsec, delt_atm, sr_flux_work, glon, glat, &! IN
  nsteps_rad, rshrt, zmean_phy) ! OUT
!
!  短波放射計算ステップの処理 
!     次の短波放射計算までの平均の天頂角・短波放射を計算する。
!     積分時間は季節 rday の時刻 rsec から 1 時間である。
!
  use prm, only: idim, jdim
    ! jcn_iwl_skip が ±1 ならば 毎ステップ放射計算
  use com_runconf_sib0109, only: jcn_iwl_skip
    ! cosz_min が正であることを仮定
  use sibcon, only: cosz_min_c
  !
  implicit none
  !
  real(8), intent(in):: rday 
  real(8), intent(in):: rsec
  real(8), intent(in):: delt_atm
  real(8), intent(in):: sr_flux_work(idim, jdim)
  real(8), intent(in):: glon(idim, jdim)
  real(8), intent(in):: glat(idim, jdim)
  real(8), intent(out):: zmean_phy(idim, jdim)
  real(8), intent(out):: rshrt(idim, jdim)
  !
  integer, intent(out):: nsteps_rad
  !
  real(8):: cosz(idim, jdim) 
  real(8):: cosz_sum(idim, jdim) 
  integer:: icosz_sum(idim, jdim) 
  integer:: ist
  real(8):: dsec
  !
  !  ステップ数
  !
  if (abs(jcn_iwl_skip) == 1) then   
    ! 毎ステップ放射計算
    nsteps_rad = 1 
  else
    ! 一時間に一度
    nsteps_rad = 3601 / delt_atm 
  endif  
  !
  cosz_sum(:, :) = 0.d0
  icosz_sum(:, :) = 0
  !
  do ist = 1, nsteps_rad
    dsec = rsec + delt_atm / 86400.d0 *( ist - 0.5d0 ) 
    !
    call islscp_sunang(rday, dsec, glon, glat, &! IN
      cosz) ! OUT
    !
    where (cosz >= cosz_min_c)
      cosz_sum = cosz_sum + cosz
      icosz_sum = icosz_sum + 1
    end where
  enddo
  !
  zmean_phy(:, :) = 0.0
  rshrt(:, :) = 0.0
  where (cosz_sum >= cosz_min_c)
    zmean_phy = cosz_sum / icosz_sum
    rshrt = sr_flux_work * cosz_sum / nsteps_rad
  end where
  !
  return
end subroutine

!=======================================================================
! forcing_shortwave_divide - 短波放射の4種分解
!   短波入射量 rshort, 雲量 cloudiness, 天頂角 zmean_phy から
!   放射4成分 rvisb, rvisd, rnirb, rnird を算出する

subroutine forcing_shortwave_divide( &
  rvisb, rvisd, rnirb, rnird, &! OUT
  rshort, cloudiness, zmean_phy &! IN
)
  use forcing, only: IDIM, JDIM
  implicit none
  !
  integer, parameter:: double = kind(0.0d0)
  !
  real(double), intent(out):: rvisb(IDIM * JDIM)
  real(double), intent(out):: rvisd(IDIM * JDIM)
  real(double), intent(out):: rnirb(IDIM * JDIM)
  real(double), intent(out):: rnird(IDIM * JDIM)
  real(double), intent(in):: rshort(IDIM * JDIM)
  real(double), intent(in):: cloudiness(IDIM * JDIM)
  real(double), intent(in):: zmean_phy(IDIM * JDIM)
  !
  real(double), parameter:: ZERO = 0.0_double
  real(double), parameter:: UNITY = 1.0_double
  real(double), parameter:: PERCENT = 0.01_double
  real(double):: dif_ratio(IDIM * JDIM)
  real(double):: vn_ratio(IDIM * JDIM)
  !
  ! Goudriaan(1977) による下向き短波放射の分解
  !
  dif_ratio(:) = 0.0683_double &
    + 0.0604_double / max(PERCENT, ZMEAN_PHY(:) - 0.0223_double)
  dif_ratio(:) = min(UNITY, max(ZERO, dif_ratio(:)))
  dif_ratio(:) = dif_ratio(:) + (1.0_double - dif_ratio(:)) * cloudiness(:)
  vn_ratio(:) = (580.0_double - cloudiness(:) * 464.0_double) &
    / ((580.0_double - cloudiness(:) * 499.0_double) &
      + (580.0_double - cloudiness(:) * 464.0_double))
  rvisb(:) = (1.0_double - dif_ratio(:)) * vn_ratio(:) * rshort(:)
  rvisd(:) = dif_ratio(:) * vn_ratio(:) * rshort(:)
  rnirb(:) = (1.0_double - dif_ratio(:)) * (1.0_double - vn_ratio(:)) * rshort(:)
  rnird(:) = dif_ratio(:) * (1.0_double - vn_ratio(:)) * rshort(:)
end subroutine

