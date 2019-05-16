! futil.F90 - 入出力関係の言語機能を補完するサブルーチン
!=======================================================================
! get_unused_unit - 未使用で利用可能な装置番号 unit を探索する。
! unified-sibs を想定して装置番号 51--69 が優先的に割り当てられる
!
subroutine get_unused_unit(unit)
  implicit none
  integer, intent(out):: unit 
  logical:: exist, opened
  do, unit = 51, 69
    inquire(unit=unit, exist=exist, opened=opened)
    if (exist .and. .not. opened) return
  enddo
  do, unit = 1, 1001, 2
    inquire(unit=unit, exist=exist, opened=opened)
    if (exist .and. .not. opened) return
  enddo
  unit = -1
end subroutine

!=======================================================================
! read_grads_date - GrADS コントロールファイルに許される時刻書式の読み取り
! string format: [hh[:mm]Z][dd]mmm[cc]yy
!   mmm は定数 mmm を参照
! idate は (/iy, imon, id, ih, isec/) からなる配列

subroutine read_grads_date(string, idate)
  implicit none
  character(len = *), intent(in):: string
  integer, intent(out):: idate(5)
  character(len = 15):: prefix, word, suffix
  character(len = 3), parameter:: mmm(12) = &
    & (/'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
    & 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/)
  integer:: i, ofs
continue
  word = string
  call uppercase(word)

  ! 月名を見出す
  idate(2) = 0
  do, i = 1, 12
    ofs = index(word, mmm(i))
    if (ofs == 0) cycle
    ! 月が見出されれば残りを解析
    idate(2) = i
    prefix = word(1: ofs-1)
    suffix = word(ofs+3: len(word))
    read(unit=suffix, fmt='(I4)') idate(1)
    if (len_trim(suffix) <= 2) then
      if (idate(1) < 50) idate(1) = idate(1) + 2000
      if (idate(1) < 100) idate(1) = idate(1) + 1900
    endif
    exit
  enddo
  if (idate(2) == 0) then
    print *, 'warning: read_grads_date error <', trim(string), '>'
    idate(:) = 0
    return
  endif

  ! Z を見出す
  word = prefix
  ofs = index(word, 'Z')
  if (ofs == 0) then
    prefix = "00:00"
  else
    prefix = word(1:ofs-1)
  endif
  suffix = word(ofs+1: len(word))
  
  if (suffix == '') then
    idate(3) = 1
  else
    read(suffix, '(I2)') idate(3)
  endif

  ! : を見出す
  word = prefix
  ofs = index(word, ':')
  if (ofs == 0) then
    prefix = word
    suffix = "00"
  else
    prefix = word(1:ofs-1)
    suffix = word(ofs+1: len(word))
  endif

  if (suffix == '') then
    idate(5) = 0
  else
    read(suffix, '(I2)') idate(5)
    idate(5) = idate(5) * 60
  endif
  if (prefix == '') then
    idate(4) = 0
  else
    read(prefix, '(I2)') idate(4)
  endif
end subroutine

!=======================================================================
! read_grads_increment - GrADS コントロールファイルの時間刻みの読み取り
! string format: xxMN|xxHR|xxDY|xxMO|xxYR
! idate は (/iy, imon, id, ih, isec/) からなる配列

subroutine read_grads_increment(string, idate)
  implicit none
  character(len = *), intent(in):: string
  integer, intent(out):: idate(5)
  character(len = 4):: word
  integer:: ofs, value
continue
  ofs = verify(string, ' -0123456789')
  word = string(1: ofs-1)
  read(word, '(I4)') value
  word = string(ofs: ofs+1)
  call uppercase(word)
  idate(:) = 0
  if (word == 'YR') then
    idate(1) = value
  else if (word == 'MO') then
    idate(2) = value
  else if (word == 'DY') then
    idate(3) = value
  else if (word == 'HR') then
    idate(4) = value
  else if (word == 'MN') then
    idate(5) = value * 60
  else
    print *, 'warning: read_grads_increment error <', trim(string), '>'
  endif
end subroutine

