! forcing_open.F90 - 強制力ファイルを開く (汎用)
! vi: set sw=2 ts=72:

! サブルーチンは典型的な使用順に並んでいる。

!=======================================================================
! forcing_open - 物理量 cmark を格納したファイル file を開く。 
!
subroutine forcing_open(cmark, file, iostat)
  use forcing, only: IDIM, JDIM, FORCING_FILE, first_file, cur_file
  implicit none
  character(len = 4), intent(in):: cmark
  character(len = *), intent(in):: file
  integer, intent(out):: iostat

  integer:: iunit, ios, recl_opened
  logical:: opened
  integer, save:: recl_mabiki = 0
  integer, save:: recl_grads = 0
  logical, save:: lfirst = .true.
  type(FORCING_FILE), pointer:: tmp_file

! 20050803 added
  real, pointer:: buf(:, :)
! till here

  ! ファイル情報置き場を拡張する
  if (lfirst .or. .not. associated(first_file)) then
    allocate(first_file)
    nullify(first_file%next)
    first_file%cmark = cmark
    cur_file => first_file
    lfirst = .false.
  else
    cur_file => first_file
    do
      if (cur_file%cmark == cmark) then
        print *, 'cmark=', cmark, ' already opened';
        return
      endif
      if (.not. associated(cur_file%next)) then
        allocate(cur_file%next)
        nullify(cur_file%next%next)
        cur_file%next%cmark = cmark
        cur_file => cur_file%next
        exit
      endif
      cur_file => cur_file%next
    enddo
  endif

  ! GrADS 読み取り以前にやっておきたい初期化
  cur_file%tmax = 0
  cur_file%form = ''
  cur_file%lat_origin = -89.5
  cur_file%lat_increment = 1.0
  cur_file%lon_origin = 0.5
  cur_file%lon_increment = 1.0

	
  ! mabiki, GrADS 形式のレコード長を計算
  if (recl_mabiki == 0) then
    !20050803
	!allocate(cur_file%last%buf(IDIM, JDIM))
    !inquire(iolength=recl_grads) cur_file%last%buf
    !inquire(iolength=recl_mabiki) &
    !  cmark, cur_file%last%id(1:4), cur_file%last%buf
    !deallocate(cur_file%last%buf)

	! new added 20050803
    allocate(buf(IDIM, JDIM))
    inquire(iolength=recl_grads) buf  ! recl_grads ､ﾏﾓ帛h餃
    inquire(iolength=recl_mabiki) cmark, cur_file%last%id(1:4), buf
    deallocate(buf)
	! till here

  endif

  ! 同名のファイルがすでに開かれていて recl があっていれば mabiki 形式
  inquire(file=file, opened=opened, number=iunit, recl=recl_opened)
  if (opened .and. recl_opened == recl_mabiki) then
    cur_file%form = 'MABIKI'
    cur_file%id_origin(:) = 0
    cur_file%id_increment(:) = 0
  endif

  if (cur_file%form == '') then
    ! GrADS 形式を試す
    call get_unused_unit(iunit)
    call forcing_open_grads(iunit, file, recl_grads)
  endif

  if (cur_file%form == '') then
    ! mabiki 形式を試す
    open(unit=iunit, file=file, access='DIRECT', recl=recl_mabiki, &
      & form='UNFORMATTED', action='READ', iostat=ios)
    cur_file%form = 'MABIKI'
    if (ios /= 0) cur_file%form = ''
    cur_file%id_origin(:) = 0
    cur_file%id_increment(:) = 0
  endif

  ! その他の成分の初期化
  cur_file%unit = iunit
  cur_file%last%id(:) = 0
  cur_file%last%irec = 0
  cur_file%last%cmark = ''
  nullify(cur_file%last%buf)
  cur_file%before = cur_file%last
  nullify(cur_file%maxmin_mask)
  cur_file%cyclic = .FALSE.
  cur_file%id_start(:) = 0
  cur_file%id_cycle(:) = 0
  cur_file%id_offset(:) = 0

  ! 報告
  if (cur_file%form == '') then
    write(6, *) 'forcing_open(', cmark, ', ', file, ') failed', ios
    iostat = ios
    !
    ! cur_file を null にする
    !
    if (associated(cur_file, first_file)) then
      first_file => first_file%next
      deallocate(cur_file)
    else
      tmp_file => first_file
      do
        if (.not. associated(tmp_file)) exit
        if (associated(tmp_file%next, cur_file)) then
          tmp_file%next => cur_file%next
          deallocate(cur_file)
        endif
        tmp_file => tmp_file%next
      enddo
    endif
  else
    write(6, "(a,i2)") ' forcing_open('//cmark//', '//trim(file)//') form=' &
      //cur_file%form//' unit=', cur_file%unit
    if (cur_file%cmark == 'GRADS') write(6, *) '   offset=', cur_file%levoffset
    iostat = 0
  endif
end subroutine

!=======================================================================
! forcing_set_offset - 入力の時間ずれ調整
!   現在のファイルの時間ずれを指定。指定単位だけ後のデータを読み取る。

subroutine forcing_set_offset(year, month, day, hour, sec)
  use forcing, only: cur_file
  implicit none
  integer, intent(in), optional:: year
  integer, intent(in), optional:: month
  integer, intent(in), optional:: day
  integer, intent(in), optional:: hour
  integer, intent(in), optional:: sec
  if (.not. associated(cur_file)) return
  if (present(year)) cur_file%id_offset(1) = year
  if (present(month)) cur_file%id_offset(2) = month
  if (present(day)) cur_file%id_offset(3) = day
  if (present(hour)) cur_file%id_offset(4) = hour
  if (present(sec)) cur_file%id_offset(5) = sec
end subroutine

!=======================================================================
! forcing_open_grads - GrADS 形式の読み込みを試みる
!
! コントロールファイル file を読み取り、
! 装置番号 iunit でデータファイルを開く。成功すれば cur_file に
! コントロールファイル情報を書き、 cur_file%form = 'GRADS' とする。
! 失敗すれば cur_file%form には何も書き込まれない。
! すでにファイルが開かれている場合には iunit を書き換える。
!
subroutine forcing_open_grads(iunit, file, recl_grads)
  use forcing, only: cur_file, JDIM
  implicit none
  integer, intent(inout):: iunit
  character(len = *), intent(in):: file
  integer, intent(in):: recl_grads

  integer, parameter:: LINELEN = 128
  character(len = LINELEN):: line, word, dset, dset_orig, nextword
  integer:: ios, i, levoffset, levels, iunitctl, iunit_already
  integer:: id_origin(5), id_increment(5)
  real:: lat_origin, lat_increment
  logical:: varsmode, opened
continue
  call get_unused_unit(iunitctl)
  open(unit=iunitctl, file=file, access='sequential', form='formatted', &
    & action='read', iostat=ios)
  if (ios /= 0) return

  ! コントロールファイルの各行についてループ
  varsmode = .FALSE.
  levoffset = 0
  cur_file%levoffset = 0
  do
    read(unit=iunitctl, fmt='(A)', iostat=ios) line
    if (ios /= 0) then
      close(iunitctl)
      write(6, *) 'unexpected EOF in GrADS Control file'
      return
    endif

    ! 最初の語を抜き出して大文字に変換
    call getword(line, word)
    call uppercase(word)

    if (word == 'DSET') then
      call getword(line, dset_orig)
      dset = dset_orig
      if (dset(1:1) == '^') then
        i = index(file, '/', back=.TRUE.)
        call replace_char(dset, '^', file(1:i))
      endif
    else if (word == 'TDEF') then
      call getword(line, word)
      read(unit=word, fmt='(I4)') cur_file%tmax
      call getword(line, word)
      call getword(line, word)
      call read_grads_date(word, id_origin)
      call getword(line, word)
      call read_grads_increment(word, id_increment)
    else if (word == 'XDEF') then
      call getword(line, word)
      call getword(line, word)
      call uppercase(word)
      if (word == 'LINEAR') then
        call getword(line, word)
        read(unit=word, fmt='(g4.0)') cur_file%lon_origin
        call getword(line, word)
        read(unit=word, fmt='(g4.0)') cur_file%lon_increment
      endif
    else if (word == 'YDEF') then
      call getword(line, word)
      call getword(line, word)
      call uppercase(word)
      if (word == 'LINEAR') then
        call getword(line, word)
        read(unit=word, fmt='(g4.0)') lat_origin
        call getword(line, word)
        read(unit=word, fmt='(g4.0)') lat_increment
        cur_file%lat_origin = lat_origin + (JDIM - 1) * lat_increment
        cur_file%lat_increment = - lat_increment
      endif
    else if (word == 'ZDEF' .or. word == 'TITLE' &
    .or. word == 'UNDEF' .or. word == 'OPTIONS') then
      continue
    else if (word == 'VARS') then
      varsmode = .TRUE.
    else if (word == 'ENDVARS') then
      exit
    else if (word(1:1) == '*') then
      continue
    else if (varsmode) then
      call getword(line, nextword)
      read(unit=nextword, fmt='(I4)') levels
      levoffset = levoffset + max(levels, 1)
      if (word == cur_file%cmark) cur_file%levoffset = levoffset
    else
      ! 未知の文があったら抜ける... バイナリファイルを開くとここに落ちる
      if (verify(word, ' 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ') == 0) then
        write(6, *) 'unsupported GrADS control word ', trim(word)
      else
        write(6, *) 'it seems not to be GrADS control file :', trim(file)
      endif
      close(iunit)
      return
    endif
  enddo

  close(iunitctl)
  cur_file%varlevs = levoffset
  if (cur_file%levoffset == 0) then
    write(6, *) 'cmark=', cur_file%cmark, " not found"
    return
  endif

  ! データファイルを開く
  inquire(file=dset, opened=opened, number=iunit_already)
  if (opened) then
    iunit = iunit_already
  else
    open(unit=iunit, file=dset, access='direct', recl=recl_grads, &
      & form='unformatted', action='read', iostat=ios)
    if (ios /= 0) print *, '#open error dset<', trim(dset), '> ios', ios
  endif
  if (ios /= 0) then
    ! 吉村モニタの ctl ファイルは余計なパスを出力する対策
    i = index(dset_orig, '/', back=.TRUE.)
    dset = dset_orig(i: len(dset_orig))
    i = index(file, '/', back=.TRUE.)
    call replace_char(dset, '/', file(1: i))
    print *, '#retry open dset<', trim(dset), '>'
    inquire(file=dset, opened=opened, number=iunit_already)
    if (opened) then
      iunit = iunit_already
    else
      open(unit=iunit, file=dset, access='direct', recl=recl_grads, &
        & form='unformatted', action='read', iostat=ios)
      print *, "#forcing_open failed ", ios
      if (ios /= 0) STOP
    endif
  endif

  cur_file%form = 'GRADS'
  cur_file%id_origin(:) = id_origin(:)
  cur_file%id_increment(:) = id_increment(:)
#ifdef DEBUG
  print *, cur_file%form, id_origin, id_increment
#endif
end subroutine

!=======================================================================
! forcing_seek - ファイルの位置づけ
!    cur_file の次の読み込み記録番号を rec に強制的に指定。
!    主に rewind のために使うもので、last バッファが破壊される。

subroutine forcing_seek(rec)
  use forcing, only: cur_file
  implicit none
  integer, intent(in):: rec
  ! 安全策: cur_file がまだ選択されていなければ、なにもしない。
  if (.not. associated(cur_file)) return
  ! 記録番号強制指定
  cur_file%last%irec = rec - 1
  ! バッファのクリア
  cur_file%last%id(:) = 0
  cur_file%last%cmark = ' '
  if (associated(cur_file%last%buf)) deallocate(cur_file%last%buf)
end subroutine

!=======================================================================
! forcing_select - ファイルの選択
!    forcing_open で cmark を使って開かれたファイルを cur_file に選択。
!    もし cmark == ' ' ならばなにもしない。

subroutine forcing_select(cmark)
  use forcing, only: cur_file, first_file
  implicit none
  character(len = *), intent(in):: cmark
  if (cmark == ' ') return
  ! forcing_open した順番に呼んでいくと高速。
  if (associated(cur_file)) then
    do
      if (cur_file%cmark == cmark) return
      cur_file => cur_file%next
      if (.not. associated(cur_file)) exit
    enddo
  endif
  cur_file => first_file
  do, while (associated(cur_file))
    if (cur_file%cmark == cmark) return
    cur_file => cur_file%next
  enddo
  write(6, *) "file for item ", cmark, " not opened";
  write(6, '(A)', advance='NO') "opened files are: "
  cur_file => first_file
  do, while (associated(cur_file))
    write(6, '(A)', advance='NO') "<", cur_file%cmark, ">, "
    cur_file => cur_file%next
  enddo
  write(6, *) ""
  nullify(cur_file)
end subroutine

!=======================================================================
! forcing_close_files - forcing_open で開いたファイルをすべて閉じる

subroutine forcing_close_files
  use forcing, only: forcing_file, first_file, cur_file
  implicit none
  cur_file => first_file
  do
    if (.not. associated(cur_file)) exit
    close(cur_file%unit)
    first_file => cur_file%next
    if (associated(cur_file%last%buf)) deallocate(cur_file%last%buf)
    if (associated(cur_file%before%buf)) deallocate(cur_file%before%buf)
    deallocate(cur_file)
  end do
end subroutine

!=======================================================================
! forcing_maxmin_region - すべてのファイルについて統計検査を行う範囲を指定
!                  

subroutine forcing_maxmin_region(mask)
  use forcing, only: forcing_file, first_file, cur_file, IDIM, JDIM
  implicit none
  logical, intent(in):: mask(IDIM, JDIM)
  logical, pointer:: mymask(:, :)
  allocate(mymask(IDIM, JDIM))
  mymask(:, :) = mask(:, :)
  cur_file => first_file
  do
    if (.not. associated(cur_file)) exit
    cur_file%maxmin_mask => mymask
    cur_file => cur_file%next
  end do
end subroutine
