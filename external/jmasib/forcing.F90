! forcing.F90 - ISLSCP 互換強制力入力ルーチンの共有情報
! vi: set sw=2 ts=72:
!
! 旧 islscp.F90 と islscp_file_unit.F90 の情報を両方持っている。
! forcing_open が最初の読み込み時に初期化する。

module forcing

  use prm, only: IDIM, JDIM
  implicit none

  ! 入力ファイルの数
  integer, parameter:: IDATE_LENGTH = 5
  integer, parameter:: CMARK_LENGTH = 4

  ! ある時刻の記録を読み取ったもの。バッファリングのために使う。

  type FORCING_RECORD
    ! 記録番号
    integer:: irec
    ! 年月日時秒
    integer:: id(5)   
    ! データ本体
    real, pointer:: buf(:, :)
    ! 種別名
    character(len = CMARK_LENGTH):: cmark
  end type

  ! 読み出す物理量ごとに FORCING_FILE 構造体を構成する。
  ! すべての FORCING_FILE 構造体は連結リスト構造になっていて、
  ! 読み取り時には cmark の一致するノードを探し出して利用する。

  type FORCING_FILE
    ! 次に開かれた物理量
    type(FORCING_FILE), pointer:: next
    ! 物理量名称
    character(len = CMARK_LENGTH):: cmark
    ! 装置番号
    integer:: unit
    ! バッファリング用
    type(FORCING_RECORD):: last
    type(FORCING_RECORD):: before
    ! 入出力に時刻補正を行う場合は非零に
    integer:: id_offset(IDATE_LENGTH)
    ! ファイル形式: 'MABIKI' = mabiki 出力, 'GRADS' = GrADS 形式
    character(len = 8):: form
    ! form="GRADS" の場合に使用される欄
    integer:: id_origin(IDATE_LENGTH)  ! 第一記録の時刻
    integer:: id_increment(IDATE_LENGTH)  ! 時間間隔
    integer:: tmax  ! 時刻最大値
    integer:: varlevs  ! 1時刻あたりのレベル数
    integer:: levoffset  ! cmark 変数の時刻群内でのオフセット(最初が1)
    real:: lat_origin, lat_increment  ! コントロールファイルの経緯度
    real:: lon_origin, lon_increment
    ! 周期的入出力を行うか、もし行うならば、起点と周期
    logical:: cyclic
    integer:: id_start(IDATE_LENGTH)
    integer:: id_cycle(IDATE_LENGTH)
    ! 入力ごとに統計検査を行う場合は、マスク配列
    logical, pointer:: maxmin_mask(:, :)
  end type

  ! すべての FORCING_FILE 構造体を保持する連結リスト構造の先頭
  type(FORCING_FILE), pointer, save:: FIRST_FILE

  ! cmark の一致するノードを探し出した結果の一時置き場
  type(FORCING_FILE), pointer, save:: cur_file

  !
  ! 長くなるのでここには置かないが、以下のサブルーチンがある。
  ! (こんなふうに interface だけを書いておけばモジュール手続の
  ! ように引数の型検査をコンパイラがやってくれるから、ないよりまし。)
  !

  interface

    !
    ! forcing/forcing_open.F90 所属
    !

    subroutine forcing_open(cmark, file, iostat)
      character(len = 4), intent(in):: cmark
      character(len = *), intent(in):: file
      integer, intent(out):: iostat
    end subroutine

    subroutine forcing_set_offset(year, month, day, hour, sec)
      integer, intent(in), optional:: year
      integer, intent(in), optional:: month
      integer, intent(in), optional:: day
      integer, intent(in), optional:: hour
      integer, intent(in), optional:: sec
    end subroutine

    subroutine forcing_seek(rec)
      integer, intent(in):: rec
    end subroutine

    subroutine forcing_select(cmark)
      character(len = *), intent(in):: cmark
    end subroutine

    subroutine forcing_close_files()
    end subroutine

    subroutine forcing_maxmin_region(mask)
      use prm, only: IDIM, JDIM
      logical, intent(in):: mask(IDIM, JDIM)
    end subroutine

    !
    ! forcing/forcing_read.F90 所属
    !

    subroutine forcing_read_id(cmark, data, id)
      use prm, only: IDIM, JDIM
      character(len = 4), intent(in):: cmark
      real(8), intent(out):: data(IDIM, JDIM)
      integer, intent(in):: id(5)
    end subroutine

    subroutine forcing_read_nearest_id(cmark, data, id, halfpoint)
      use prm, only: IDIM, JDIM
      character(len = 4), intent(in):: cmark
      real(8), intent(out):: data(IDIM, JDIM)
      integer, intent(in):: id(5)
      real(8), intent(in):: halfpoint
    end subroutine

    subroutine forcing_read_id2(cmark, id, &
      & id1, data1, id2, data2, weight, update)
      use prm, only: IDIM, JDIM
      character(len = 4), intent(in):: cmark
      integer, intent(in):: id(5)
      integer, intent(out):: id1(5)
      integer, intent(out):: id2(5)
      real, intent(out):: data1(IDIM, JDIM)
      real, intent(out):: data2(IDIM, JDIM)
      real(8), intent(out):: weight
      logical, intent(out):: update
    end subroutine

    subroutine forcing_read_cmark(cmark_f, data, cmark)
      use prm, only: IDIM, JDIM
      character(len = 4), intent(in):: cmark_f, cmark
      real, intent(out):: data(IDIM, JDIM)
    end subroutine

    subroutine forcing_read_real(cmark, irec, data)
      use prm, only: IDIM, JDIM
      character(len = 4), intent(in):: cmark
      integer, intent(in):: irec
      real, intent(out):: data(IDIM, JDIM)
    end subroutine

    subroutine forcing_read_int(cmark, irec, data)
      use prm, only: IDIM, JDIM
      character(len = 4), intent(in):: cmark
      integer, intent(in):: irec
      integer, intent(out):: data(IDIM, JDIM)
    end subroutine

    subroutine forcing_fetch(cmark, irec)
      character(len = 4), intent(in):: cmark
      integer, intent(in):: irec
    end subroutine

    !
    ! forcing/forcing_setup.F90 所属
    !

    subroutine forcing_setup()
    end subroutine

    !
    ! forcing/forcing_ini.F90 所属
    !

    subroutine forcing_ini(id_now, delt_atm, rsec, imask, glon, glat)
      use prm, only: IDIM, JDIM
      integer, intent(in):: id_now(5)
      real(8), intent(in):: delt_atm
      real(8), intent(in):: rsec
      integer, intent(in):: imask(IDIM, JDIM)
      real(8), intent(in):: glon(IDIM, JDIM)
      real(8), intent(in):: glat(IDIM, JDIM)
    end subroutine

    !
    ! forcing/forcing_main.F90 所属
    !
    subroutine forcing_main(& 
      id_now_in, id_pre_in, id_now5,id_pre5, &
      delt_atm, rday, rsec, imask, glon, glat, &! IN
      u_phy, v_phy, pd_phy, ps_phy, pf_phy, tmp_phy, q_phy, &! OUT
      zmean_phy, ztemp_phy, ppli_phy, ppci_phy, &! OUT
      rvisb, rvisd, rnirb, rnird, dlwb)  ! OUT
      !
      use prm, only: IDIM, JDIM
      !
      ! input
      !
      integer, intent(in):: id_now_in(5)  ! 本ステップの時刻
      integer, intent(in):: id_pre_in(5)  ! 本ステップの時刻
      integer, intent(in):: id_now5
      integer, intent(in):: id_pre5
      real(8), intent(in):: rsec 
      real(8), intent(in):: rday 
      real(8), intent(in):: delt_atm
      real(8), intent(in):: glon(IDIM, JDIM)  !  経度（単位、度）
      real(8), intent(in):: glat(IDIM, JDIM)  ! 緯度（単位、度）
      integer, intent(in):: imask(IDIM, JDIM)  !
      !
      ! output 
      !
      real(8), intent(out):: u_phy(IDIM, JDIM)  ! u
      real(8), intent(out):: v_phy(IDIM, JDIM)  ! v 
      real(8), intent(out):: pd_phy(IDIM, JDIM)  ! (ps-ph)*2  
      real(8), intent(out):: ps_phy(IDIM, JDIM)  ! ハーフレベル = pS hpA
      real(8), intent(out):: pf_phy(IDIM, JDIM)  ! フルレベル hpA 
      real(8), intent(out):: tmp_phy(IDIM, JDIM)  ! 温度 
      real(8), intent(out):: q_phy(IDIM, JDIM)  ! 比湿 KG/KG
      real(8), intent(out):: zmean_phy(IDIM, JDIM)  ! 一時間?平均天頂角
      real(8), intent(out):: ztemp_phy(IDIM, JDIM)  ! 各ステップ天頂角
      real(8), intent(out):: ppli_phy(IDIM, JDIM)  ! 大規模凝結性降水
      real(8), intent(out):: ppci_phy(IDIM, JDIM)  ! 積雲性降水
      real(8), intent(out):: rvisb(IDIM * JDIM)  ! 放射計算時可視直達
      real(8), intent(out):: rvisd(IDIM * JDIM)  ! 放射計算時可視散乱
      real(8), intent(out):: rnirb(IDIM * JDIM)  ! 放射計算時近赤直達
      real(8), intent(out):: rnird(IDIM * JDIM)  ! 放射計算時近赤散乱
      real(8), intent(out):: dlwb(IDIM * JDIM)  ! 放射計算時長波
    end subroutine

    subroutine monit_regist_forcing()
    end subroutine

    !
    ! forcing/forcing_get_geography.F90 所属
    !

    subroutine forcing_get_geography(out_glon, out_glat, out_imask)
      use prm, only: IDIM, JDIM
      integer, parameter:: DOUBLE = kind(0.0D0) ! i.e. 8 for most systems
      real(DOUBLE), intent(out):: out_glon(IDIM, JDIM)
      real(DOUBLE), intent(out):: out_glat(IDIM, JDIM)
      integer, intent(out), target:: out_imask(IDIM, JDIM)
    end subroutine

    !
    ! forcing/forcing_shortwave.F90 所属
    !

    subroutine forcing_read_shortwave( &
      cmark_swdn, cmark_cld, id, delt_atm, glon, glat, imask, &! IN
      swdn, rvisb, rvisd, rnirb, rnird, zmean_phy, ztemp_phy, cloudiness &! OUT
    )
      use prm, only: IDIM, JDIM
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
      real(8), intent(out):: cloudiness(IDIM, JDIM)
    end subroutine

  end interface

end module
