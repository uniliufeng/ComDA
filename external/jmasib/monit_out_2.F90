SUBROUTINE monit_out_2              &!
( imonit_2fl, imonit_2fl8 , imonit_sfl,           &! Input
  idate, model, resl, expr, cinf,   &! Input
  kt, totsec, fsec, imask, gw )      ! Input
!
  USE prm,ONLY : &
    idim,  &
    jdim,  &
    imax,  &
    jmax,  &
    jlmax, &
    jblk2, &
    jblk1
  USE com_runconf,ONLY : &
    jcnmmm 
!
  USE com_jobinfo_sib0109,ONLY : &
    iout_8byte

  USE com_monit,ONLY : &
    mxptr_2,  &
    ldisk_2,  &
    lsdisk_2, &
    ctflag_2, &
    corder_2, &
    sum_2,    &
    cmask_2,  &
    smax_2,   &
    smin_2,   &
    sl_2,     &
    smag_2,   &
    ckwd_2,   &
    ctitle_2, &
    cunit_2

    
 USE com_step,ONLY : &
   scn_delt

 use com_step_sib0109 , only :  &
   icn_SIB0109_CALC_SOIL_SNOW 

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: imonit_2fl
  INTEGER,INTENT(IN) :: imonit_2fl8
  INTEGER,INTENT(IN) :: imonit_sfl
  INTEGER,INTENT(IN) :: idate(5)
  CHARACTER(LEN= 8),INTENT(IN) :: model
  CHARACTER(LEN=11),INTENT(IN) :: resl
  CHARACTER(LEN= 4),INTENT(IN) :: expr
  CHARACTER(LEN=80),INTENT(IN) :: cinf(10)
  INTEGER,INTENT(IN) :: kt
  REAL(8),INTENT(IN) :: totsec
  REAL(8),INTENT(IN) :: fsec
  INTEGER,INTENT(IN) :: imask(idim,jdim)
  REAL(8),INTENT(IN) :: gw(jmax)

! REAL(8),PARAMETER :: crwmsk = 0.1D0

  INTEGER,SAVE :: irec = 0

  REAL(4) :: agcm (imax,jmax)
  REAL(8) :: agcm8(imax,jmax)

  REAL(8) :: one_totsec
  
  INTEGER :: i,j,k,ij,jl,lat

!-----------       start !       ------------
# ifdef DEBUG
  write(6,*) 'monit_out_2' , fsec , expr , resl , model  
# endif

! write(6,*) 'monit_out_2' , fsec , expr , resl , model , cinf , gw , kt
! if ( mod(kt,10000).eq.9999 ) write(6,*) cinf

! IF ( kt <= 0 ) RETURN     ! ppp 99/10/05 pochaka 

!-- FOR GRSMK (2-D)

  DO k=1,mxptr_2             !-----   Loop 2D  start   -----
!
    IF ( ldisk_2(k) == 0 .AND. lsdisk_2(k) == 1 ) CYCLE
!HYHY    IF ( ctflag_2(k) == 'AVR' .AND. totsec < 0.1D0 ) CYCLE
    IF ( ctflag_2(k) == 'AVR' ) THEN
      one_totsec = 1.0D0/totsec
    ELSE
      one_totsec = 1.0D0/scn_delt
!HYHY      one_totsec = 1.0D0
    END IF

    IF ( corder_2(k) == 'PH' ) THEN
!      write(6,*) 'monit_out_2 ' ,ckwd_2(k)

      DO jl = 1,jlmax
        DO j=1,jblk2
          DO i = 1,imax
            agcm8( i, (jl-1)*jblk2+j )                     &
                = sum_2( (j-1)*imax+i, jl, k )*one_totsec
            agcm( i, (jl-1)*jblk2+j ) = agcm8( i, (jl-1)*jblk2+j ) 
          END DO
        END DO
      END DO
    ELSE
      DO jl = 1,jlmax
        DO j=1,jblk1
          lat  = (jl-1)*jblk1 + j
          DO i = 1,imax
            ij = (j-1)*imax+i
            agcm8(i,     lat  )                          &
                = sum_2(ij          ,jl,k)*one_totsec
            agcm8(i,jmax-lat+1)                          &
                = sum_2(ij+imax*jblk1,jl,k)*one_totsec
            agcm(i,     lat  ) = agcm8(i,     lat  ) 
            agcm(i,jmax-lat+1) = agcm8(i,jmax-lat+1)
          END DO
        END DO
      END DO
    END IF

    DO j=1,jmax
      DO i=1,imax
        IF ( ( cmask_2(k) == 'LAND' .AND. imask(i,j) == 0 )       &
             .OR. ( cmask_2(k) == 'SEA ' .AND. imask(i,j) /= 0 ) ) THEN
          agcm8(i,j) = 0. 
          agcm(i,j) = -9.99E33
        END IF
      END DO
    END DO
!
    irec = irec + 1
    if ( iout_8byte .eq. 1 ) then
      WRITE( imonit_2fl8, REC=irec ) agcm8(1:imax,jmax:1:-1)
    endif
    WRITE( imonit_2fl , REC=irec ) agcm (1:imax,jmax:1:-1)
!   WRITE( 6,* ) agcm
!XXX    CALL monit_wrt2d(imonit_2fl,irec,imax,jmax,agcm)

    IF ( lsdisk_2(k) == 1 .AND. jcnmmm == 1 ) THEN
      CALL monit_shibata( imonit_sfl,imax,jmax,agcm,           &
                   smax_2(k),smin_2(k),sl_2(k),smag_2(k),  &
                   ckwd_2(k),ctitle_2(k),cunit_2(k),idate )
    END IF

  END DO                  !-----  Loop 2D  end  -----


END SUBROUTINE monit_out_2





