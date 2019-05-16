SUBROUTINE monit_regist                       &
( idisk, atype, akwd, amask, atflag, aorder,  &! In
  atitle, aunit,                              &! In
  isdisk, smax, smin, sl, smag )               ! In

  USE prm
  USE com_monit
  USE com_runconf,ONLY : &!
    jcnimnt,    &!
    jcnimntday, &!
    jcnimnt6hr

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: idisk  
    ! 0 ( not write )
    ! 1 ( write )
    ! 2 ( write zonal mean only )
  CHARACTER( 4),INTENT(IN) :: atype
    ! '2D  '  ( 2-dimensional monitor )
    ! '3D  '  ( 3-dimensional monitor )
    ! '3DET'  ( 3-dimensional model level monitor )
    ! '6HR '  ( 6 hourly monitor )
    ! 'DAY '  ( daily monitor )
    ! 'GLOB'  ( global monitor )
  CHARACTER( 8)            :: akwd_8  
  INTEGER  ( 8)            :: Ikwd_8  
  EQUIVALENCE ( AKWD_8 , IKWD_8 )
  CHARACTER( 7),INTENT(IN) :: akwd  
    ! 'EVAP','ROFS'.....variable keyword
  CHARACTER( 4),INTENT(IN) :: amask
    ! 'LAND'  ( land only )
    ! 'SEA '  ( sea only )
    ! 'BOTH'  ( both land and sea )
  CHARACTER( 3),INTENT(IN) :: atflag
    ! 'SNP'  (snap shot)
    ! 'AVR'  (time average)
  CHARACTER( 2),INTENT(IN) :: aorder
    ! 'PH'...PHYSCS order
    !     ( PHYSCS/physcs3s.F, PHYSCS/igw.F, PHYSCS/sfcbnd.F )
    ! 'DY'...DYNAMICS order
    !     ( DYNAMICS/grddynam.F, MOIST-HONCHOU/gmoisth.F etc.)
  CHARACTER(32),INTENT(IN) :: atitle  ! title
  CHARACTER(13),INTENT(IN) :: aunit   ! unit

  INTEGER,INTENT(IN) :: isdisk 
    ! Shibata monitor ( 0=not write, 1=write )
  REAL(8),INTENT(IN) :: smax
    ! max value for Shibata monitor
  REAL(8),INTENT(IN) :: smin
    ! min value for Shibata monitor
  REAL(8),INTENT(IN) :: sl
    ! mid value for Shibata monitor
  REAL(8),INTENT(IN) :: smag
    ! magnitude for Shibata monitor

! INTEGER :: k,i
  INTEGER ::   i           ! 99/08/11 pochaka 

!-----------      start  !!     ------------

  IF ( amask /= 'LAND' .AND. amask /= 'SEA '        &
       .AND. amask /= 'BOTH'                 ) THEN
    WRITE(6,*) 'amask should be LAND or SEA or BOTH. ',akwd,' ',amask
    STOP 'monit_regist'
  END IF

  IF ( atflag /= 'SNP' .AND. atflag /= 'AVR' ) THEN
    WRITE(6,*) 'atflag should be SNP or AVR. ',akwd,' ',atflag
    STOP 'monit_regist'
  END IF

  IF ( aorder /= 'PH' .AND. aorder /= 'DY' ) THEN
    WRITE(6,*) 'aorder should be PH or DY. ',akwd,' ',aorder
    STOP 'monit_regist'
  END IF
!
  AKWD_8 = AKWD 

  IF ( atype == '2D  ' ) THEN
    DO i = 1,mxptr_2
!         ! 011114_1 ‚ÅŽè“–‚±‚±‚©‚ç mhosaka
!     IF ( Ikwd_8 == Ikwd_2(i) ) THEN     
      IF ( Akwd_8 == Ckwd_2(i) ) THEN
!         ! 011114_1 ‚ÅŽè“–‚±‚±‚Ü‚Å
!
        WRITE(6,*) 'The keyword was used :',atype,akwd
        STOP 'monit_regist'
      END IF
    END DO
    IF ( ( idisk >= 1 .OR. isdisk >= 1 ) .AND. jcnimnt /= 0 ) THEN
      mxptr_2 = mxptr_2 + 1
!
      IF ( mxptr_2 > limit_2 ) THEN
        WRITE(6,*) 'The number of data is too large.'
        STOP 'monit_regist'
      END IF
      ldisk_2(mxptr_2)  = idisk
      ckwd_2(mxptr_2)   = akwd
      cmask_2(mxptr_2)  = amask
      ctflag_2(mxptr_2) = atflag
      ctitle_2(mxptr_2) = atitle
      cunit_2(mxptr_2)  = aunit
      corder_2(mxptr_2) = aorder
      lsdisk_2(mxptr_2) = isdisk
      smax_2(mxptr_2)   = smax
      smin_2(mxptr_2)   = smin
      sl_2(mxptr_2)     = sl
      smag_2(mxptr_2)   = smag
      WRITE(6,*) atype,' ',akwd,' ' , atflag,' ' , atitle,'[',aunit,']'
    ELSE
      mxptr_nouse_2 = mxptr_nouse_2 + 1
      IF ( mxptr_nouse_2 > limit_2 ) THEN
        WRITE(6,*) 'The number of data is too large.'
        STOP 'monit_regist'
      END IF
      ckwd_nouse_2(mxptr_nouse_2)   = akwd
    END IF
!
  ELSE
    WRITE(6,*) 'Registration was failed ',atype,akwd,atitle
    STOP 'monit_regist'
  END IF

END SUBROUTINE monit_regist
