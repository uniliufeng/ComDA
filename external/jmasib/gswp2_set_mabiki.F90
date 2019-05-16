subroutine gswp2_set_mabiki ( idim, jdim, &
    ilon , xlon , jlat, ylat )  
!
implicit none
!
integer,intent(in)  :: idim, jdim 
integer,intent(out) :: ilon(idim)
integer,intent(out) :: jlat(jdim)
real(4),intent(out) :: XLON(IDIM)    
real(4),intent(out) :: YLAT(JDIM) 
!
integer :: i,j 
!
! Šiqî•ñ‚Ìİ’è
!   “Œ¼‚É‚Â‚¢‚Ä‚ÍA¼‚©‚ç‹l‚ß‚é
!   “ì–k‚É‚Â‚¢‚Ä‚ÍAÔ“¹‘¤‚©‚ç‹l‚ß‚éB
!
  j = jdim 
  if ( mod(j,2) /= 0 .and. j /= 1 ) then
    write(6,*) 'gswp2_set_mabiki ERROR J should be EVEN or 1, but J= ' , J 
    stop 999
  endif
!
! “Œ¼ (‚à‚Æ‚Ì I=1 ‚Í 179.5W )
!   mabiku  ilon(1)=1 xlon(1)=-179.5 
!
  do i=1,idim
    ilon (i) = (i-1) * 360 / idim + 1 
    xlon (i) = ilon(i) - 180. - 0.5 
    write(6,*) 'gswp2_set_mabiki lon ' , ilon(i) , xlon(i)  
  enddo
  write(6,*) 'longitudes'
  write(6,*) (xlon(i),i=1,idim)
!
  do i=2,idim
    if ( ilon(i) .le. ilon(i-1) ) then
      write(6,*) 'GSWP2_set_mabiki ERROR ILON ' , ilon 
      stop 999
    endif
  enddo
!
! “ì–k
!
  if (jdim .eq. 1 ) then
    jlat(1) = 90
    ylat(1) = 0.5
  else 
!          Ô“¹‚©‚ç‚Ì‡”Ô‚ª i, –k‹É‚©‚ç‚Ì‡”Ô‚ª j  
    do j=jdim/2,1,-1
      i  = (jdim/2 -j) * 180 / jdim 
      jlat (j) = 90-i
      ylat (j) = i + 0.5 
      write(6,*) 'gswp2_set_mabiki lat ' , jlat(j) , ylat(j)
    enddo
!
    do j=1,jdim/2
      jlat(jdim-j+1) = 181 - jlat(j) 
      ylat(jdim-j+1) = - ylat(j) 
      write(6,*) 'gswp2_set_mabiki lat ' , jlat(j) , ylat(j)
    enddo
!
    do j=2,jdim
      if ( jlat(j) .le. jlat(j-1) ) then
        write(6,*) 'GSWP2_set_mabiki ERROR JLAT ' , jlat
        stop 999
      endif
    enddo
  endif
  write(6,*) 'latitudes'
  write(6,*) '1 to jdim = ' , (ylat(j),j=1,jdim)
  write(6,*) 'jdim to 1 = ' , (ylat(j),j=jdim,1,-1)
!
  end subroutine gswp2_set_mabiki


