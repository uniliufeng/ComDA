SUBROUTINE monit_grads_ctl                      &
( nf2ctl, nfzctl, nf3ctl, nfzetactl, nf3etactl, &! In
  nf6hrctl, nfdayctl, nfglobctl,                &! In
  c2name, czname, c3name, czetaname, c3etaname, &! In
  c6hrname, cdayname, cglobname,                &! In
  cosclt, idstar, idend, ktstar, ktend,         &! In
  model, resl, gw, pa, pb , glon , glat )                      ! In

!*****  make grads control files  *****

  USE prm
  USE com_monit
  USE com_runconf
  USE com_runconf_sib0109 , only : rcn_delt_cnp 
! USE com_papbfull , only : pafull, pbfull

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nf2ctl
  INTEGER,INTENT(IN) :: nfzctl
  INTEGER,INTENT(IN) :: nf3ctl
  INTEGER,INTENT(IN) :: nfzetactl
  INTEGER,INTENT(IN) :: nf3etactl
  INTEGER,INTENT(IN) :: nf6hrctl
  INTEGER,INTENT(IN) :: nfdayctl
  INTEGER,INTENT(IN) :: nfglobctl
!
  CHARACTER(*),INTENT(IN) :: c2name
  CHARACTER(*),INTENT(IN) :: czname
  CHARACTER(*),INTENT(IN) :: c3name
  CHARACTER(*),INTENT(IN) :: czetaname
  CHARACTER(*),INTENT(IN) :: c3etaname
  CHARACTER(*),INTENT(IN) :: c6hrname
  CHARACTER(*),INTENT(IN) :: cdayname
  CHARACTER(*),INTENT(IN) :: cglobname
!   
  REAL(8),INTENT(IN) :: cosclt(jmax)
  INTEGER,INTENT(IN) :: idstar(5)
  INTEGER,INTENT(IN) :: idend(5)
  INTEGER,INTENT(IN) :: ktstar
  INTEGER,INTENT(IN) :: ktend
  CHARACTER( 8),INTENT(IN) :: model
  CHARACTER(11),INTENT(IN) :: resl
  REAL(8),INTENT(IN) :: gw(jmax)
  REAL(8),INTENT(IN) :: pa(kmax+1)
  REAL(8),INTENT(IN) :: pb(kmax+1)
  real(8),intent(in) :: glon(idim) 
  real(8),intent(in) :: glat(jdim) 

  CHARACTER(3),PARAMETER ::                               &
    month(12) = (/ 'JAN','FEB','MAR','APR','MAY','JUN',   &
                   'JUL','AUG','SEP','OCT','NOV','DEC' /)


  REAL(8) :: alat(jmax)     ! latitude

  REAL(8) :: alat3(maxj_3)   ! lat. for 3D data
  REAL(8) :: alat3d(maxj_3)

  REAL(8) :: alat3e(maxj_3eta)  ! lat. for 3Deta data
  REAL(8) :: alat3ed(maxj_3eta)

! REAL(8) :: pf(kmax)

  CHARACTER(50) :: title
  CHARACTER(50) :: undef
  CHARACTER(50) :: options
  CHARACTER(50) :: xdef
  CHARACTER(50) :: xdef3
  CHARACTER(50) :: xdef3e
  CHARACTER(50) :: ydef
  CHARACTER(50) :: ydef3
  CHARACTER(50) :: ydef3e
  CHARACTER(50) :: zdef
  CHARACTER(50) :: tdef
  CHARACTER(50) :: tdef6hr
  CHARACTER(50) :: tdefday

  INTEGER :: j,j3,nj,nn,i,m

!---------     start  !!     ----------

  if ( nf2ctl.lt.0 ) then      ! for no warning pochaka
   write(6,*) 'monit_grs_ctl pa pb dummy' , pa(1) , pb(1) , &
    nf2ctl, nfzctl, nf3ctl, nfzetactl, nf3etactl, &
    nf6hrctl, nfdayctl, nfglobctl,                &
    c2name, czname, c3name, czetaname, c3etaname, &
    c6hrname, cdayname, cglobname,                &
    cosclt, idstar, idend, ktstar, ktend,         &
    model, resl, gw, pa, pb     
  endif

  DO j = 1,jmax
    alat(j) = asin( cosclt(j) )*180.0D0/3.14159265D0
  END DO

  alat3(:)  = 0.0D0
  alat3d(:) = 0.0D0

  DO j3 = 1,jmax/intj_3
    DO nj = 1,intj_3
      j = (j3-1)*intj_3 + nj
      alat3 (j3) = alat3 (j3) + alat(j)*gw(j)
      alat3d(j3) = alat3d(j3) + gw(j)
    END DO
  END DO

  alat3(:) =  alat3(:)/alat3d(:)

  alat3e(:)  = 0.0D0
  alat3ed(:) = 0.0D0

  DO j3 = 1,jmax/intj_3eta
    DO nj = 1,intj_3eta
      j = (j3-1)*intj_3eta + nj
      alat3e (j3) = alat3e (j3) + alat(j)*gw(j)
      alat3ed(j3) = alat3ed(j3) + gw(j)
    END DO
  END DO

  alat3e(:) =  alat3e(:)/alat3ed(:)

! pf(:) = pafull(:) + pbfull(:)*1000.0D0


  title = 'TITLE  '//model//' '//resl
  undef = 'UNDEF  -9.99E33'
  options = 'OPTIONS BIG_ENDIAN template'

! write( xdef(1:50),'(A,I5,A,F8.4)' )   &
!   'XDEF ',imax,' LINEAR 0 ', 360.0/imax
  write( xdef(1:50),'(A,I5,A,2F10.5)' )   &
    'XDEF ',imax,' LINEAR ', glon(1) , 360.0/imax

  write( xdef3(1:50),'(A,I5,A,F8.4,1X,F8.4)' ) &
    'XDEF ',maxi_3,' LINEAR ',                  &
    180.0/imax*(inti_3-1), 360.0/maxi_3
  write( xdef3e(1:50),'(A,I5,A,F8.4,1X,F8.4)' ) &
    'XDEF ',maxi_3eta,' LINEAR ',                  &
    180.0/imax*(inti_3eta-1), 360.0/maxi_3eta
!
! write( ydef(1:50),'(A,I5,A)' )    &
!   'YDEF ',jmax,' LEVELS' 

!! write( ydef(1:50),'(A,I5,A,1X,F8.4)' )    &
!!   'YDEF ',jmax,' LEVELS' , alat(jmax)

! write( ydef(1:50),'(A,I5,A,1X,2F8.4)' )    &
  write( ydef(1:50),'(A,I5,A)' )    &
       'YDEF ',jmax,' LEVELS' 

  write( ydef3(1:50),'(A,I5,A)' )    &
    'YDEF ',maxj_3,' LEVELS'
  write( ydef3e(1:50),'(A,I5,A)' )    &
    'YDEF ',maxj_3eta,' LEVELS'
  write( zdef(1:50),'(A,I5,A,F7.2)' ) &
    'ZDEF ',1,' LEVELS  ',p_monit(1)
!
! Šî–{“I‚É + 1 ˆ—
!
  IF ( jcnimnt >= 900 ) THEN
    WRITE( tdef(1:50),'(A,I5,A,I2,A3,I4,1X,A)' )            &
      'TDEF ',                                              &
      (idend(1)-idstar(1))*12 + idend(2)-idstar(2) + 1  ,   &
      ' LINEAR  ',idstar(3),                                &
      month(idstar(2)),idstar(1),'1MO'
  ELSE IF ( jcnimnt <= 0 ) THEN
    if ( ABS( INT((RCN_DELT_CNP+0.01)/60.)*60 - RCN_DELT_CNP ) &
         .gt. 0.1 ) then 
       WRITE( tdef(1:50),'(A,I2,A3,I4,1X,A)')                  &
         'TDEF   9999 LINEAR  ',idstar(3),                     &
         month(idstar(2)),idstar(1),'1HR'
    elseif ( idstar(3).ge.10 ) then
      WRITE( tdef(1:50),'(A,I2,A,I2,A3,I4,1X,I4,A)')          &
        'TDEF   9999 LINEAR  ',        &
        idstar(4),                     &
        ':00Z' ,                       &
        idstar(3),                     &
        month(idstar(2)),              &
        idstar(1),                            &
        INT((RCN_DELT_CNP+0.01)/60.) , 'MN'
    else 
      WRITE( tdef(1:50),'(A,I2,A,I1,A3,I4,1X,I4,A)')          &
        'TDEF   9999 LINEAR  ',        &
        idstar(4),                     &
        ':00Z' ,                       &
        idstar(3),                     &
        month(idstar(2)),              &
        idstar(1),                            &
        INT((RCN_DELT_CNP+0.01)/60.) , 'MN'
     endif
  ELSE IF ( 100*(jcnimnt/24) == (100*jcnimnt)/24 ) THEN
    WRITE( tdef(1:50),'(A,I5,A,I2,A3,I4,1X,I2,A)' )         &
      'TDEF ',(ktend-ktstar)/jcnimnt+1,' LINEAR  ',idstar(3), &
      month(idstar(2)),idstar(1),jcnimnt/24,'DY'
  ELSE
    WRITE( tdef(1:50),'(A,I5,A,I2,A3,I4,1X,I3,A)' )         &
      'TDEF ',(ktend-ktstar)/jcnimnt+1,' LINEAR  ',idstar(3), &
      month(idstar(2)),idstar(1),jcnimnt,'HR'
  END IF
  
  IF ( jcnimnt6hr >= 900 ) THEN
    WRITE( tdef6hr(1:50),'(A,I5,A,I2,A3,I4,1X,A)' )         &
      'TDEF ',                                              &
      (idend(1)-idstar(1))*12 + idend(2)-idstar(2) + 1 ,    &
      ' LINEAR  ',idstar(3),                                &
      month(idstar(2)),idstar(1),'1MO'
  ELSE IF ( jcnimnt6hr <= 0 ) THEN
    WRITE( tdef6hr(1:50),'(A,I2,A3,I4,1X,A)')                  &
      'TDEF   100 LINEAR  ',idstar(3),                         &
      month(idstar(2)),idstar(1),'1HR'
  ELSE IF ( 100*(jcnimnt6hr/24) == (100*jcnimnt6hr)/24 ) THEN
    WRITE( tdef6hr(1:50),'(A,I5,A,I2,A3,I4,1X,I2,A)' )         &
      'TDEF ',(ktend-ktstar)/jcnimnt6hr,' LINEAR  ',idstar(3), &
      month(idstar(2)),idstar(1),jcnimnt6hr/24,'DY'
  ELSE
    WRITE( tdef6hr(1:50),'(A,I5,A,I2,A3,I4,1X,I3,A)' )         &
      'TDEF ',(ktend-ktstar)/jcnimnt6hr,' LINEAR  ',idstar(3), &
      month(idstar(2)),idstar(1),jcnimnt6hr,'HR'
  END IF
  
  IF ( jcnimntday >= 900 ) THEN
    WRITE( tdefday(1:50),'(A,I5,A,I2,A3,I4,1X,A)' )         &
      'TDEF ',                                              &
      (idend(1)-idstar(1))*12 + idend(2)-idstar(2) + 1 ,    &
      ' LINEAR  ',idstar(3),                                &
      month(idstar(2)),idstar(1),'1MO'
  ELSE IF ( jcnimntday <= 0 ) THEN
    WRITE( tdefday(1:50),'(A,I2,A3,I4,1X,A)')                  &
      'TDEF   100 LINEAR  ',idstar(3),                         &
      month(idstar(2)),idstar(1),'1HR'
  ELSE IF ( 100*(jcnimntday/24) == (100*jcnimntday)/24 ) THEN
    WRITE( tdefday(1:50),'(A,I5,A,I2,A3,I4,1X,I2,A)' )         &
      'TDEF ',(ktend-ktstar)/jcnimntday,' LINEAR  ',idstar(3), &
      month(idstar(2)),idstar(1),jcnimntday/24,'DY'
  ELSE
    WRITE( tdefday(1:50),'(A,I5,A,I2,A3,I4,1X,I3,A)' )         &
      'TDEF ',(ktend-ktstar)/jcnimntday,' LINEAR  ',idstar(3), &
      month(idstar(2)),idstar(1),jcnimntday,'HR'
  END IF
!


!--------      2D DATA       -----------

  WRITE(nf2ctl,'(A)') 'DSET  ^'//TRIM(c2name)//'.dr'
  WRITE(nf2ctl,'(A)') title//' 2D DATA'
  WRITE(nf2ctl,'(A)') options
  WRITE(nf2ctl,'(A)') undef
  WRITE(nf2ctl,'(A)') xdef
  WRITE(nf2ctl,'(A)') ydef 
!!! WRITE(nf2ctl, '( (5X,5(F8.3,1X)) )' ) alat(jmax:1:-1)
  WRITE(nf2ctl, '( (5X,5(F8.3,1X)) )' ) glat(jmax:1:-1)
!! if ( jmax.gt.2 ) then
!!   WRITE(nf2ctl, '( (5X,5(F8.3,1X)) )' ) alat(jmax-1:1:-1)
!! endif
  WRITE(nf2ctl,'(A)') zdef
  WRITE(nf2ctl,'(A)') tdef
  nn = 0
  DO i = 1,mxptr_2
    IF ( ldisk_2(i) >= 1 ) nn=nn+1
  END DO
  WRITE(nf2ctl,'(A,I4)') 'VARS ',nn
  DO m=1,mxptr_2
    IF ( ldisk_2(m) >= 1 ) THEN
      WRITE(nf2ctl,'(A)')                         &
        ckwd_2(m)//'0 0 '//ctitle_2(m)//' '//cunit_2(m)//','//corder_2(m)
    END IF
  END DO
  WRITE(nf2ctl,'(A)') 'ENDVARS'

END SUBROUTINE monit_grads_ctl  
