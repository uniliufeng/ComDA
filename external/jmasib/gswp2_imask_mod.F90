subroutine gswp2_imask_mod ( idim , jdim , ifile_log, imask )
!
implicit none
!
integer,intent(in) :: idim
integer,intent(in) :: jdim
integer,intent(in) :: ifile_log
integer :: ifile
integer,intent(inout) :: imask(idim,jdim)
!
integer :: i, j , ii
!!!integer :: icount_mask(0:100) = 0 
integer :: icount_mask(0: 30) = 0 
character(1)             :: cland    (idim,jdim) 
character(1)             :: cwork    (60) 
!
  if ( ifile_log > 0 ) then
    ifile = ifile_log
  else 
    ifile = 6
  endif
!
  do j=1,jdim
  do i=1,idim
    icount_mask(imask(i,j)) =  icount_mask(imask(i,j)) +1 
  enddo
  enddo
  write(ifile,*) 'original imask' , icount_mask(:)
    print *, 'check 1'
  write(ifile,*) '  12(Agri)->7(Grass) 13(Wet)->7 14(Water)->0 15(Ice) ->13 16(Miss)->0' 
    print *, 'check 2'
!C12 Agriculture or C3 Grassland  to  7
!D13 Persistent Wetland           to 10 
!E14 Water                        to  0 
!F15 Ice Cap and Glacier          to 13
!G16 Missing                      to  0
  icount_mask  (:) = 0 
  do j=1,jdim
  do i=1,idim
    if     ( imask(i,j) .eq. 12 ) then
      imask(i,j) = 7 
    elseif ( imask(i,j) .eq. 13 ) then
      imask(i,j) = 10 
    elseif ( imask(i,j) .eq. 14 ) then
      imask(i,j) =  0 
    elseif ( imask(i,j) .eq. 15 ) then
      imask(i,j) =  13 
    elseif ( imask(i,j) .eq. 16 ) then
      imask(i,j) =  0
    endif
    icount_mask(imask(i,j)) =  icount_mask(imask(i,j)) +1 
  enddo
  enddo
    print *, 'check 3'
  write(ifile,*) 'modified imask' , icount_mask(:)
    print *, 'check 4'
!
  do j=1,jdim
    do i=1,idim
      if     ( imask(i,j) <= 9 ) then 
        write(cland(i,j),'(I1.1)') imask(i,j) 
      elseif ( imask(i,j) == 10 ) then
        cland(i,j) = 'A'
      elseif ( imask(i,j) == 11 ) then
        cland(i,j) = 'B'
      elseif ( imask(i,j) == 12 ) then
        cland(i,j) = 'C'
      elseif ( imask(i,j) == 13 ) then ! Ice Sheet 
        cland(i,j) = 'D'
      else
        write(6,*) 'ERROR gswp2__ini ' , i , j, imask(i,j)
        stop 999
      endif
    enddo
!
    do ii=1,(idim-1)/60 +1
      cwork(:) = ' ' 
      cwork(1:min(idim-(ii-1)*60,60)) = cland((ii-1)*60+1:min(ii*60,idim),j)
      write(ifile,*) cwork(1:60)
    enddo
  enddo
!
end subroutine gswp2_imask_mod
