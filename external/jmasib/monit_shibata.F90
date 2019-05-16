SUBROUTINE monit_shibata                         &! Shibata Monitor
( ifu, ip, jp, sp, smax, smin, sl, xmag,   &
  ckwd, cttl, cunt, idate )

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ifu         ! Output file no.
  INTEGER,INTENT(IN) :: ip
  INTEGER,INTENT(IN) :: jp
  REAL(4),INTENT(IN) :: sp(ip,jp)   ! 2D data
  REAL(8),INTENT(IN) :: smax        ! max value
  REAL(8),INTENT(IN) :: smin        ! min value
  REAL(8),INTENT(IN) :: sl          ! base value
  REAL(8),INTENT(IN) :: xmag        ! magnitude,  sp*xmag
  CHARACTER( 7),INTENT(IN) :: ckwd  ! keyword
  CHARACTER(32),INTENT(IN) :: cttl  ! title
  CHARACTER(13),INTENT(IN) :: cunt  ! unit
  INTEGER,INTENT(IN) :: idate(5)    ! date

! REAL(8),PARAMETER :: zero = 0.0D0

  INTEGER,PARAMETER :: ig  = 100
  INTEGER,PARAMETER :: jg  = 37

  REAL(8) :: gx(ig)
  REAL(8) :: gy(jg)
  REAL(8) :: smap(ig)
  CHARACTER(2) :: smapc(ig)
  REAL(8) :: dx,dy,y,x
  INTEGER :: i,j,jy,jy1,ix,ix1

!----------     start  !!    ------------------

  dx = float(ip-1)/float(ig-1)
  dy = float(jp-1)/float(jg-1)

  gx(1) = 0.0D0
  DO i=2,ig-1
    gx(i) = dx*float(i-1)
  END DO
  gx(ig) = float(ip-1)-0.0001

  gy(1) = 0.0D0
  DO j=2,jg-1
    gy(j) = dy*float(j-1)
  END DO
  gy(jg) = float(jp-1)-0.0001

  WRITE(IFU,"(1H1)")
  WRITE(IFU,*) ckwd//' : '//cttl//' ('//cunt//') *',xmag
  WRITE(IFU,*) 'idate = ',idate
  WRITE(IFU,"(1H ,'MAX=',E10.3,' MIN=',E10.3,' SL=',E10.3)") &
    smax, smin, sl

  DO j = 1,jg
    jy = INT( gy(j) )
    y = gy(j)-jy             ! 0.0 <= y < 1.0
!    y = MAX( y,zero )
    jy = jy+1                ! 1 <= jy  < jp-1
    jy1 = jy+1               ! 2 <= jy1 < jp
!    jy1 = MIN( jy1,jp )
    DO i = 1,ig
      ix = INT( gx(i) )
      x  = gx(i) - ix        ! 0.0 <= x < 1.0
!      x = MAX( x, zero )
      ix = ix + 1            ! 1 <= ix  < ip-1
      ix1 = ix + 1           ! 2 <= ix1 < ip
!      ix1 = MIN( ix+1, ip )
      IF (     sp(ix ,jy ) < -9.0D33        &
          .OR. sp(ix1,jy ) < -9.0D33        &
          .OR. sp(ix ,jy1) < -9.0D33        &
          .OR. sp(ix1,jy1) < -9.0D33 ) THEN
        smap(i) = -9.99D33
      ELSE
        smap(i)                          &
        = xmag                           &
          *( (1.0-x)*(1.0-y)*sp(ix,jy)   &
             + x*(1.0-y)*sp(ix1,jy)      &
             + (1.0-x)*y*sp(ix,jy1)      &
             + x*y*sp(ix1,jy1)         )
      END IF 
    END DO
    CALL mntsgraph                      &
         ( ifu,ig,smap,smax,smin,sl,j,  &! In
           smapc )                       ! Out
    WRITE(ifu,"(1H ,5X,127A1)") ( smapc(i),i=1,ig )
  END DO

END SUBROUTINE monit_shibata

!+++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE mntsgraph           &
( ifu,nx,x,xmax,xmin,xsl,kkk,  &! In
  xc )                          ! Out

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ifu
  INTEGER,INTENT(IN) :: nx
  REAL(8),INTENT(IN) :: x(nx)
  REAL(8),INTENT(IN) :: xmax
  REAL(8),INTENT(IN) :: xmin
  REAL(8),INTENT(IN) :: xsl
  INTEGER,INTENT(IN) :: kkk
  CHARACTER(LEN=2),INTENT(OUT) :: xc(nx)

  CHARACTER(LEN=2),PARAMETER ::                               &
    CHA(30) = (/'1 ','  ','2 ','  ','3 ','  ','4 ','  ',  &
                '5 ','  ','6 ','  ','7 ','  ','8 ','  ',  &
                '9 ','  ','0 ','  ','* ','  ','+ ','  ',  &
                '% ','  ','@ ','  ','# ','X '/)

  CHARACTER(LEN=2),PARAMETER ::                               &
    CHB(30) = (/'  ','A ','  ','B ','  ','C ','  ','D ',  &
                '  ','E ','  ','F ','  ','G ','  ','H ',  &
                '  ','I ','  ','J ','  ','K ','  ','L ',  &
                '  ','M ','  ','N ','  ','P '/)

!  REAL(8),PARAMETER :: d30 = 30.0D0

  REAL(8),SAVE :: dx
  REAL(8),SAVE :: one_dx

  CHARACTER(LEN=2) :: ch(60)

  REAL(8) :: xg
  INTEGER :: i,ia,ib,ixg

!-----------     start   !!     ---------------

  IF ( kkk == 1 ) THEN
    dx = (xmax-xmin)/10.0D0
    IF ( dx == 0.0D0 ) dx = 1.0D0
    one_dx = 1.0D0/dx
    DO i=1,30
      ia = i+30
      ch(ia) = cha(i)
      ib = 31-i
      ch(i) = chb(ib)
    END DO
    WRITE(ifu,"(1H ,'DX=',E10.3,10X,60A1)") dx,ch
  END IF
 
  DO i=1,nx
    xg = ( x(i)-xsl )*one_dx
    IF ( xg >= 0.0D0 ) THEN
!      xg   = MIN( xg, d30 )
      ixg  = MIN( 30, INT(xg)+1 )
      xc(i)= cha(ixg)
    ELSE
!      xg   = MAX( xg, -d30 )
      ixg  = MIN( 30, INT(-xg)+1 )
      xc(i)= chb(ixg)
    END IF
  END DO

END SUBROUTINE mntsgraph
