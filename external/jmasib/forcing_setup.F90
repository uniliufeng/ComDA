! forcing_setup.F90 - unified-sibs の強制力ファイルを開く
! vi: set sw=2 ts=72:

!=======================================================================
! forcing_setup - unified-sibs の強制力ファイルを開く
!
subroutine forcing_setup
  use forcing, only: forcing_open, forcing_set_offset, IDIM, JDIM
  implicit none
  character(len = 11):: resolution, resolution_i
  character(len = 12):: resolution_c
  character(len = 4), dimension(9), parameter:: CMARKS = &
    (/'CLD ', 'MWND', 'TEMP', 'QREF', 'PRSS', 'TPRC', 'CPRC', 'SWDN', 'LWDN'/)
  integer:: ios, item

  resolution = '_III_JJJ_R4'
  call replace_int(resolution, 'III', IDIM)
  call replace_int(resolution, 'JJJ', JDIM)
  resolution_i = resolution
  call replace_char(resolution_i, '_R4', '_I4')
  resolution_c = resolution
  call replace_char(resolution_c, '_R4', '.ctl')

  !
  ! --- 時系列データ ---
  !

  do, item = 1, size(CMARKS)
    !
    ! GrADS ファイルに必要変数がすべてあればそれを使う 
    !
    call forcing_open(CMARKS(item), 'input/forcing' // resolution_c, ios)
    if (ios /= 0) goto 2000
  enddo
  call forcing_open('VEG ', ('input/forcing_veg' // resolution_c), ios)
  if (ios /= 0) goto 2000
  call forcing_open('LAT ', ('input/forcing_veg' // resolution_c), ios)
  call forcing_open('LON ', ('input/forcing_veg' // resolution_c), ios)

  return

  !
  ! mabiki 形式を試してみる
  !
  2000 continue
  do, item = 1, size(CMARKS)
    call forcing_open(CMARKS(item), &
      ('input/' // trim(CMARKS(item)) // resolution), ios)
    if (ios /= 0) then
      write(6, *) 'forcing file for cmark=', CMARKS(item), ' not found'
      stop 100
    endif
    if (CMARKS(item) == 'TPRC' .or. CMARKS(item) == 'CPRC') then
      call forcing_set_offset(hour = -6)
    else if (CMARKS(item) == 'LWDN') then
      call forcing_set_offset(hour = +3)
    endif
  enddo

  call forcing_open('VEG ', ('input/VEG_MAP' // resolution_i), ios)
  if (ios /= 0) then
    write(6, *) 'forcing file for cmark=VEG not found'
    stop 101
  endif
  call forcing_open('LAT ', ('input/VEG_MAP' // resolution_i), ios)
  call forcing_open('LON ', ('input/VEG_MAP' // resolution_i), ios)

end subroutine
