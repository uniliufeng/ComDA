SUBROUTINE monit_ini

  USE prm
  USE com_monit

  IMPLICIT NONE
  
!  monit_regist
!    (idisk,atype,akwd,amask,atflag,aorder,  &
!     atitle,aunit,                          &
!     isdisk, smax, smin, sl, smag )
!
!  idisk    : 0  ( not write )
!             1  ( write )
!             2  ( write zonal mean only )
!  atype    : '2D  '  ( 2-dimensional monitor )
!             '3D  '  ( 3-dimensional monitor )
!             '3DET'  ( 3-dimensional model level monitor )
!             '6HR '  ( 6 hourly monitor )
!             'DAY '  ( daily monitor )
!             'GLOB'  ( global monitor )
!  akwd     : 'EVAP','ROFS'.....variable keyword
!  amask    : 'LAND'  ( land only )
!             'SEA '  ( sea only )
!             'BOTH'  ( both land and sea )
!  atflag   : 'SNP'  (snap shot)
!             'AVR'  (time average)
!  aorder   : 'PH'...PHYSCS order
!                ( PHYSCS/physcs3s.F, PHYSCS/igw.F, PHYSCS/sfcbnd.F )
!             'DY'...DYNAMICS order
!                ( DYNAMICS/grddynam.F, MOIST-HONCHOU/gmoisth.F etc.)
!  atitle   : title
!  aunit    : unit

!  isdisk   : Shibata monitor ( 0=not write, 1=write )
!  smax     : max value for Shibata monitor
!  smin     : min value for Shibata monitor
!  sl       : mid value for Shibata monitor
!  smag     : magnitude for Shibata monitor

! REAL(8),PARAMETER :: one = 1.0D0
! REAL(8),PARAMETER :: daysec = 86400.0D0

!  REAL(8) :: smax,smin,sl,smag

! REAL(8) :: d850, d500, d200, d10

!-------------       start  !!        -------------------

!---  DUMMY DATA FOR SHIBATA MONITOR     7/16 1998 BY H.YOSHIMURA
!  smax = 100.0D0
!  smin = -100.0D0
!  sl   = 0.0D0
!  smag = ONE

!   =================================================================
!   FOR 2D DATA
!   =================================================================
!
  write(6,*) 'IO-YOSHIMURA touroku start' 
!
   CALL MONIT_REGIST_SIBMAIN
!
   CALL SET_COM_MONIT
!
  write(6,*) 'IO-YOSHIMURA touroku end' 

END SUBROUTINE monit_ini
