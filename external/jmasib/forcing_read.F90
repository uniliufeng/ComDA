! forcing_read.F90 - (汎用) 時系列データの1時刻読み取り
! vi: set sw=2 ts=72:

! forcing_locate_id の詳細メッセージを得るには
! #define DEBUG_LOCATE
! とせよ

! 時系列データの読み取りを行うルーチンをここに収める。
! サブルーチンはユーザが直接目にするものを上に、それから引用されるものを
! 下にならべている。

!=======================================================================
! forcing_read_id - 時系列データの指定時刻読み取り (内挿版)
!   時系列データ cmark から指定時刻 id のデータを読み取り data に格納する。
!   指定時刻のデータがなければ内挿を行った結果を data に格納する。
!   エラー (指定時刻の前後のデータが見つからない場合も含む) 時には
!   data には -999.0 が代入される

subroutine forcing_read_id(cmark, data, id)
  use forcing, only: cur_file, IDIM, JDIM
  use date
  implicit none
  character(len = 4), intent(in):: cmark
  integer, intent(in):: id(5)
  real(8), intent(out):: data(IDIM, JDIM)
  integer:: ir1, ir2
  real(8):: weight

  call forcing_locate_id(cmark, id, ir1, ir2, weight)
  if (ir1 == 0) goto 9999
  if (ir2 == 0) goto 9999
  if (ir1 == ir2) then
    data(:, :) = cur_file%last%buf(:, :)
  else
    data = weight * cur_file%last%buf + (1.0 - weight) * cur_file%before%buf
  endif


#ifdef DEBUG
  if (associated(cur_file%maxmin_mask)) then
    print "(a4,' ',i4.4,4('-',i2.2),' max=',g16.6,' min=',g16.6)", &
      cmark, id, &
      maxval(data, mask=cur_file%maxmin_mask), &
      minval(data, mask=cur_file%maxmin_mask)
  endif
#endif
  return

  9999 continue
  data = -999.0_8
end subroutine

!=======================================================================
! forcing_read_nearest_id - 時系列データの指定時刻読み取り (内挿版)
!   時系列データ cmark から指定時刻 id のデータを読み取り data に格納する。
!   指定時刻のデータがなければもっとも近い時刻のデータを data に格納する。
!   もっとも近いとは、weight が halfpoint を超えたときは後時刻をとること。
!   エラー (指定時刻の前後のデータが見つからない場合も含む) 時には
!   data には -999.0 が代入される

subroutine forcing_read_nearest_id(cmark, data, id, halfpoint)
  use forcing, only: cur_file, IDIM, JDIM
  use date
  implicit none
  character(len = 4), intent(in):: cmark
  integer, intent(in):: id(5)
  real(8), intent(out):: data(IDIM, JDIM)
  real(8), intent(in):: halfpoint
  integer:: ir1, ir2
  real(8):: weight

  call forcing_locate_id(cmark, id, ir1, ir2, weight)
  if (ir1 == 0) goto 9999
  if (ir2 == 0) goto 9999
  if (ir1 == ir2) then
    data(:, :) = cur_file%last%buf(:, :)
  else if (weight >= halfpoint) then
    data(:, :) = cur_file%last%buf(:, :)
  else
    data(:, :) = cur_file%before%buf(:, :)
  endif

#ifdef DEBUG
  if (associated(cur_file%maxmin_mask)) then
    print "(a4,' ',i4.4,4('-',i2.2),' max=',g16.6,' min=',g16.6)", &
      cmark, id, &
      maxval(data, mask=cur_file%maxmin_mask), &
      minval(data, mask=cur_file%maxmin_mask)
  endif
#endif
  return

  9999 continue
  data = -999.0_8
end subroutine


!=======================================================================
! forcing_read_id2 - 時系列データの指定時刻読み取り (非補間版)
!   時系列データ cmark から指定時刻 id の前後のデータをさがす。
!   直前時刻を id1 に、データを data1 に格納する。
!   直後時刻を id2 に、データを data2 に格納する。
!   指定時刻のデータがなければ内挿を行った結果を data に格納する。
!   時刻から算出した重みを weight に格納する。重みとは
!   (1.0 - weight) * data1 + weight * data2 が指定時刻 id のデータを
!   与えるようなもののことである。
!   エラー (指定時刻の前後のデータが見つからない場合も含む) 時には
!   id1, id2 には両者ともに 0 が、data1, data2 には -999.0 が代入される
!   

subroutine forcing_read_id2(cmark, id, id1, data1, id2, data2, weight, update)
  use forcing, only: forcing_select, cur_file, IDIM, JDIM
  use date
  implicit none
  character(len = 4), intent(in):: cmark
  integer, intent(in):: id(5)
  integer, intent(out):: id1(5)
  integer, intent(out):: id2(5)
  real, intent(out):: data1(IDIM, JDIM)
  real, intent(out):: data2(IDIM, JDIM)
  double precision, intent(out):: weight
  logical, intent(out):: update
  ! 記録番号. 1, 2 は before バッファか last バッファかに対応し、
  !  org は探索前の値である。
  integer:: ir1, ir2, ir1org, ir2org

  call forcing_select(cmark)
  if (.not. associated(cur_file)) goto 9999
  ir1org = cur_file%before%irec
  ir2org = cur_file%last%irec
  call forcing_locate_id(cmark, id, ir1, ir2, weight)
  if (ir1 == 0) goto 9999
  if (ir2 == 0) goto 9999
  if (ir1 == ir2) then
    id1 = cur_file%last%id
    id2 = cur_file%last%id
    data1(:, :) = cur_file%last%buf(:, :)
    data2(:, :) = cur_file%last%buf(:, :)
  else
    id1 = cur_file%before%id
    id2 = cur_file%last%id
    data1(:, :) = cur_file%before%buf(:, :)
    data2(:, :) = cur_file%last%buf(:, :)
  endif
  update = ((ir1 /= ir1org) .or. (ir2 /= ir2org))
  return

  9999 continue
  id1 = 0
  id2 = 0
  data1 = -999.0
  data2 = -999.0
  weight = -1
  update = .FALSE.
end subroutine

!=======================================================================
! forcing_locate_id - 時系列データのファイル位置付け
!
!   時系列データ cmark から指定時刻 id の前後のデータをさがす。
!   直前時刻を与える記録の時刻番号を ir1 に格納する。
!   直後時刻を与える記録の時刻番号を ir2 に格納する。
!   指定時刻のデータを線形内挿によって与えるための重みを weight に
!   格納する。重みとは
!    (1.0 - weight) * data[ir1] + weight * data[ir2] が
!   指定時刻 id のデータを与えるような weight のことである。

subroutine forcing_locate_id(cmark, id, ir1, ir2, weight)
  use forcing, only: forcing_select, cur_file, &
    forcing_fetch
  use date
  implicit none
  character(len = 4), intent(in):: cmark
  integer, intent(in):: id(5)
  integer, intent(out):: ir1, ir2
  double precision, intent(out):: weight
  integer:: id_tmp(5), idiff_tmp(5)
  integer:: id_before(5), id_dif_data(5), id_dif_request(5)
  integer:: cmp
  integer:: irec_cycle_tail
  logical:: found
  double precision:: recno

#ifdef DEBUG_LOCATE
  write(6, "(a,i4.4,4('-',i2.2))") "#locate_id " // cmark, id
#endif
  if (cmark /= ' ') call forcing_select(cmark)
  !
  ! 要求時刻 id_tmp を算出
  !
  id_tmp(:) = id(:) + cur_file%id_offset(:)
  if (cur_file%cyclic) then
    call date_modulo(id_tmp, cur_file%id_start, cur_file%id_cycle)
  endif

  if (cur_file%form == 'GRADS') then
    ! GrADS の場合 ... 手間はかかるが straightforward に計算できる
    call date_diff(id_tmp, cur_file%id_origin, idiff_tmp)
    call date_diff_div(idiff_tmp, cur_file%id_increment, recno)
    ir1 = floor(recno) + 1
    ir2 = ceiling(recno) + 1
    if (cur_file%cyclic) then
      id_tmp(:) = cur_file%id_origin(:) + (ir1 - 1) * cur_file%id_increment
      call date_compare(id_tmp, cur_file%id_start, cmp)
      if (cmp < 0) then
        call date_diff(id_tmp +cur_file%id_cycle, cur_file%id_origin, idiff_tmp)
        call date_diff_div(idiff_tmp, cur_file%id_increment, recno)
        ir1 = floor(recno) + 1
      endif
      id_tmp(:) = cur_file%id_origin(:) + (ir2 - 1) * cur_file%id_increment
      call date_compare(id_tmp, cur_file%id_start + cur_file%id_cycle, cmp)
      if (cmp > 0) then
        call date_diff(id_tmp -cur_file%id_cycle, cur_file%id_origin, idiff_tmp)
        call date_diff_div(idiff_tmp, cur_file%id_increment, recno)
        ir2 = ceiling(recno) + 1
      endif
    endif
  else if (cur_file%form == 'MABIKI') then
    ! Mabiki の場合 ... キャッシュ時刻を試し、ずれていたら探査をする
    ! 現キャッシュの時刻に挟まれていれば found
    call date_in_cur_file(id_tmp, found)
    if (found) goto 1000

    ! Mabiki の場合 ... キャッシュ時刻を試し、ずれていたら探査をする
    irec_cycle_tail = 0
    ! last の次 (なければ 1) の記録から順に要求時刻を探索
    do
#ifdef DEBUG_LOCATE
      write(6, *) 'locate/MABIKI/first-loop', cur_file%last%irec
#endif
      call forcing_fetch(cmark, cur_file%last%irec + 1)
      irec_cycle_tail = cur_file%before%irec
      ! ファイル終端またはエラー時: 記録#1からやりなおし
      if (cur_file%last%irec == 0) exit
      ! 要求時刻を挟んでいれば found
      call date_in_cur_file(id_tmp, found)
      if (found) goto 1000
      ! 時刻が要求時刻を超えるようなら 記録#1からやりなおし
      call date_compare(cur_file%last%id, id_tmp, cmp)
      if (cmp > 0) exit
    enddo

    ! 記録#1から要求時刻を挟むデータをさがす。
    cur_file%last%irec = 0  ! 巻き戻し
    if (.not. found) then
      do
#ifdef DEBUG_LOCATE
        write(6, *) 'locate/MABIKI/second-loop', cur_file%last%irec
#endif
        call forcing_fetch(cmark, cur_file%last%irec + 1)
        if (cur_file%last%irec == 0) exit
        call date_compare(cur_file%last%id, id_tmp, cmp)
        ! 要求時刻より遅いデータがみつかったところで探索停止
        if (cmp == 0) then
           ir1 = cur_file%last%irec
           goto 1100
        endif
        if (cmp > 0) then
          if (cur_file%cyclic .and. irec_cycle_tail > 0) then
            ir1 = irec_cycle_tail
            goto 1100
          endif
          if (cur_file%last%irec >= 2) goto 1000
        endif
        ! 要求時刻より早いデータがなければ irec_cycle_tail を使うのだが
        ! 見つかった時点でその処理はしなくてもよくする。
        if (cmp < 0) then
          irec_cycle_tail = 0
        endif
      enddo
    endif

    print *, 'forcing_locate_id(', cmark, &
      '): file has insufficient data records'
    ir1 = 0
    ir2 = 0
    weight = -1.0
    return

    ! ここにジャンプしてきたらキャッシュの時刻が要求時刻を挟んでいる
    1000 continue
    ir1 = cur_file%before%irec
    1100 continue
    ir2 = cur_file%last%irec
  endif

  call forcing_fetch(cmark, ir1)
  call forcing_fetch(cmark, ir2)
  if (ir1 == ir2) then
    weight = 1.0
  else
    if (ir1 > ir2) then
      id_before = cur_file%before%id - cur_file%id_cycle
    else
      id_before = cur_file%before%id
    endif
    call date_diff(id_tmp, id_before, id_dif_request)
    call date_diff(cur_file%last%id, id_before, id_dif_data)
    call date_diff_div(id_dif_request, id_dif_data, weight)
  endif

#ifdef DEBUG_LOCATE
  write(6, *) "#locate_id return(", ir1, ir2, weight, ")"
#endif
end subroutine

!=======================================================================
! date_in_cur_file - cur_file のキャッシュに指定時刻が見出されるか
!   指定時刻 id が cur_file のキャッシュに見出される場合真を found に。

subroutine date_in_cur_file(id, found)
  use forcing, only: cur_file
  use date
  implicit none
  integer, intent(in):: id(5)
  logical, intent(out):: found
  integer:: cmp_b, cmp_l
  if (cur_file%last%irec == 0 .or. cur_file%before%irec == 0) goto 9000
  if (.not. associated(cur_file%last%buf)) goto 9000
  if (.not. associated(cur_file%before%buf)) goto 9000
  call date_compare(id, cur_file%last%id, cmp_l)
  call date_compare(id, cur_file%before%id, cmp_b)
  found = cmp_l <= 0 .and. cmp_b >= 0
  return

  9000 continue
  found = .FALSE.
end subroutine

!=======================================================================
! forcing_read_cmark - 時系列データの指定要素読み取り
!  forcing_open によって cmark_f として開かれたファイルから
!  cmark 要素の「次の記録」を data に読み取る。
!  cmark が空の場合は cmark 検査をしない。
!  発見できなければ last%buf を空ポインタ化し data を -999.0 で埋める。

subroutine forcing_read_cmark(cmark_f, data, cmark)
  use forcing, only: cur_file, forcing_select, forcing_fetch, IDIM, JDIM
  implicit none
  character(len = 4), intent(in):: cmark_f, cmark
  real, intent(out):: data(IDIM, JDIM)
  integer:: irec
!
  call forcing_select(cmark_f)
  irec = cur_file%last%irec
  if (irec <= 0) irec = 1
  do
    call forcing_fetch(cmark_f, irec)
    if (cur_file%last%cmark == cmark) exit
    if (.not. associated(cur_file%last%buf)) then
      data = -999.0
      return
    endif
    if (cmark == ' ') exit
    irec = irec + 1
  enddo
  data = cur_file%last%buf
end subroutine

!=======================================================================
! forcing_swap - バッファの交換
!   読み取りバッファ last, before の内容を交換する。

subroutine forcing_swap
  use forcing, only: FORCING_RECORD, cur_file
  implicit none
  type(FORCING_RECORD):: tmp_buf
  tmp_buf = cur_file%last
  cur_file%last = cur_file%before
  cur_file%before = tmp_buf
end subroutine

!=======================================================================
! forcing_read_real - 実数配列の直接読み取り
!
!   forcing_open によって要素 cmark を与えるとして開かれたファイルから、
!   abs(irec) 番目の記録を読み取り実数型配列 data に格納する。
!   エラー時にはプログラムが停止する。

subroutine forcing_read_real(cmark, irec, data)
  use forcing, only: cur_file, forcing_fetch, IDIM, JDIM
  implicit none
  character(len = 4), intent(in):: cmark
  integer, intent(in):: irec
  real, intent(out):: data(IDIM, JDIM)
  call forcing_fetch(cmark, abs(irec))
  if (.not. associated(cur_file)) then
    print *, 'cmark=<', cmark, '> not associated to file.'
    stop 
  endif
  if (cur_file%last%irec == 0) stop
  data(:, :) = cur_file%last%buf(:, :)
end subroutine

!=======================================================================
! forcing_read_int - 整数配列の直接読み取り
!
!   forcing_open によって要素 cmark を与えるとして開かれたファイルから、
!   abs(irec) 番目の記録を読み取り整数型配列 data に格納する。
!   エラー時にはプログラムが停止する。

subroutine forcing_read_int(cmark, irec, data)
  use forcing, only: cur_file, forcing_fetch, IDIM, JDIM
  implicit none
  character(len = 4), intent(in):: cmark
  integer, intent(in):: irec
  integer, intent(out):: data(IDIM, JDIM)
  call forcing_fetch(cmark, -abs(irec))
  if (.not. associated(cur_file)) then
    print *, 'cmark=<', cmark, '> not associated to file.'
    stop 
  endif
  if (cur_file%last%irec == 0) stop
  data(:, :) = cur_file%last%buf(:, :)
end subroutine

!=======================================================================
! forcing_fetch - 時系列データの直接読み取り
!
!   forcing_open によって要素 cmark を与えるとして開かれたファイルから、
!   irec 番目の記録を読み取り cur_file%last に格納する。
!   irec が負の場合は符号反転し整数型読み取りを行う。
!   既存の cur_file%last は cur_file%before に保存される。
!   既存の cur_file%before は廃棄される。
!   エラーまたはファイル範囲外の場合は cur_file%last%{irec, id} を
!   ゼロクリアし、cur_file%last%buf を空ポインタ化する。
!   cmark と cur_file%last%cmark の整合性は検査しない。

subroutine forcing_fetch(cmark, irec)
  use forcing, only: cur_file, forcing_select, IDIM, JDIM
  use date
  implicit none
  character(len = 4), intent(in):: cmark
  integer, intent(in):: irec
  integer, allocatable:: ibuf(:, :)
  character(len = 12):: access
  character(len = 4):: my_cmark
  integer:: my_id(5)
  integer:: ios, recno
  ios = 0
  call forcing_select(cmark)

  if (.not. associated(cur_file)) goto 9000

  ! バッファリングで保存したものがひっかかるとうれしい
  if (cur_file%last%irec == irec) then
    return
  else if (cur_file%before%irec == irec) then
    ! before と last を交換
    call forcing_swap
    return
  endif

  ! before を廃棄し last を before に保存
  if (associated(cur_file%before%buf)) deallocate(cur_file%before%buf)
  cur_file%before = cur_file%last
  ! last の値を準備
  cur_file%last%irec = irec
  recno = abs(irec)
  allocate(cur_file%last%buf(IDIM, JDIM))

  if (cur_file%form == 'MABIKI') then
    if (irec < 0) then
      allocate(ibuf(IDIM, JDIM))
      read(unit=cur_file%unit, rec=recno, iostat=ios) &
        cur_file%last%cmark, cur_file%last%id(1: 4), ibuf
      cur_file%last%buf = real(ibuf)
      deallocate(ibuf)
    else
      read(unit=cur_file%unit, rec=recno, iostat=ios) &
        my_cmark, my_id(1: 4), cur_file%last%buf
      cur_file%last%cmark = my_cmark
      cur_file%last%id = my_id
    endif
    cur_file%last%id(5) = 0
  else if (cur_file%form == 'GRADS') then
    recno = (recno - 1) * cur_file%varlevs + cur_file%levoffset
    read(unit=cur_file%unit, rec=recno, iostat=ios) &
      cur_file%last%buf
    cur_file%last%cmark = cmark
    cur_file%last%id = cur_file%id_origin &
      + cur_file%id_increment * (cur_file%last%irec - 1)
  endif

  ! 読み取りに失敗したら
  if (ios /= 0) then
    inquire(unit=cur_file%unit, access=access)
    print *, 'read error = ', ios, cur_file%unit, recno, access
    if (ios == 443) print *, ' may caused by missing -Fport"(stduf)"'
    goto 9000
  endif

#ifdef DEBUG
  if (associated(cur_file%maxmin_mask)) then
    print "(a4,' record=',i8,' max=',g16.6,' min=',g16.6)", &
      cmark, irec, &
      maxval(cur_file%last%buf, mask=cur_file%maxmin_mask), &
      minval(cur_file%last%buf, mask=cur_file%maxmin_mask)
  endif
#endif

  call date_normalize(cur_file%last%id)
  return

  ! エラー時終了処理
  9000 continue
  cur_file%last%irec = 0
  cur_file%last%id(:) = 0
  cur_file%last%cmark = '!err'
  if (associated(cur_file%last%buf)) deallocate(cur_file%last%buf)
end subroutine
