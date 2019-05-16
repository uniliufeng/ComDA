SUBROUTINE PRO4SAIL(Cab,Car,Cbrown,Cw,Cm,N,a,lai,q,tts,tto,psi,psoil,skyl,resh,resv)

!	Jean-Baptiste Féret (feret@ipgp.fr)
!	Institut de Phyysique du Globe de Paris
!	Space and Planetary Geophysics
!	October 2009

!	this model PRO4SAIL is based on a version provided by
!	Wout Verhoef 
!	NLR	
!	April/May 2003,
!	Improved and extended version of SAILH model that avoids numerical singularities
!	and works more efficiently if only few parameters change.
!
!	The module "staticvar" is meant for internal use in order to keep
!	several variables static between successive calls of the subroutine.
!	For calculation of vertical flux profiles inside the canopy, the quantities 
!	m, sb, sf, ks, and rinf are needed by the calling routine, so in this case 
!	this module must also be declared in the calling routine.

	USE ANGLE
	USE spectrum_width_P5
	USE dataparm
	USE lidf_parm
	USE output_PROSPECT
	USE outSAIL
	USE rfltrn
	USE spectral
	USE flag_util
	USE staticvar
	USE dataSpec_P5B
	IMPLICIT NONE

!!!!!!!!		INPUT / OUTPUT		!!!!!!!!
REAL(KIND=8),INTENT(in) :: N, Cab,Car,Cbrown,Cw,Cm
REAL(KIND=8),INTENT(in) :: a,lai,psoil
REAL(KIND=8),INTENT(in) :: q
REAL(KIND=8),INTENT(in) :: skyl,tts,tto,psi
REAL(KIND=8),INTENT(out):: resh(nw),resv(nw)
!!!!!!!!	END INPUT / OUTPUT		!!!!!!!!

REAL(KIND=8) :: chi_o,chi_s
REAL(KIND=8) :: ctl,sto,sts,tanto,tants,cospsi,koli,ksli,bf,bfli
REAL(KIND=8) :: alf
REAL(KIND=8) :: f1,f2,g1(nw),g2(nw)
REAL(KIND=8) :: fhot,fint,frho,ftau
REAL(KIND=8) :: J1ko(nw),J1ks(nw),J2ko(nw),J2ks(nw),Jfunc3
REAL(KIND=8) :: ttl
REAL(KIND=8) :: sigf(nw),sobli,sofli
REAL(KIND=8) :: T1(nw),T2(nw),T3(nw),Tv1(nw),Tv2(nw)
REAL(KIND=8) :: x1,x2,y1,y2,z

REAL(KIND=8) :: e1(nw),e2(nw),rinf2(nw),re(nw),denom(nw)
REAL(KIND=8) :: Ps(nw),Qs(nw),Pv(nw),Qv(nw)
REAL(KIND=8) :: dn(nw)

! ----------in main--------
! ANGLE CONVERSION
pi=ATAN(1.)*4.
rd=pi/180.

! PROSPECT output
ALLOCATE (LRT(nw,3))

ALLOCATE (sb(nw),sf(nw),vb(nw),vf(nw),w(nw))
ALLOCATE (m(nw),m2(nw))
ALLOCATE (att(nw))
ALLOCATE (sigb(nw))
ALLOCATE (rinf(nw))
ALLOCATE (rho(nw),tau(nw))
ALLOCATE (PARdiro(nw),PARdifo(nw))

ALLOCATE(tsd(nw),tdd(nw),tdo(nw),rsd(nw),rdd(nw),rso(nw),rdo(nw))
ALLOCATE(rddt(nw),rsdt(nw),rdot(nw),rsodt(nw),rsost(nw),rsot(nw),rsos(nw),rsod(nw))

! LEAF ANGLE DISTRIBUTION
na = 18
ALLOCATE (thetal(na))
ALLOCATE (lidf(na))
DO i=1,na
	thetal(i)=2.5+5*(i-1)
ENDDO
ALLOCATE (rsoil(nw,2),rsoil0(nw))
!------in main------
DO i=1,9
	flag(i)=.not.init_completed
ENDDO

IF (init_completed) THEN

	!	Detect which inputs have changed (if it is not the first time)
	delta_lai	= (lai.ne.lai_old)
	delta_hot	= (q.ne.q_old)
	delta_geom	= (tts.ne.tts_old).or.(tto.ne.tto_old).or.&
					(psi.ne.psi_old)
	delta_soil	= (psoil.ne.psoil_old)
	delta_lidf	= (a.ne.angl_old)
	delta_leaf	= (N.ne.N_old).or.(Cab.ne.Cab_old).or.(Car.ne.Car_old).or.&
					(Cbrown.ne.Cbrown_old).or.(Cw.ne.Cw_old).or.(Cm.ne.Cm_old)
	delta_skyl	= (skyl.ne.skyl_old)

	!	Raise the flags for the modules to be executed
	flag(1) = delta_geom
	flag(2) = delta_lidf
	flag(3) = delta_geom.or.delta_lidf
	flag(4) = flag(3).or.delta_leaf
	flag(5) = flag(4).or.delta_lai.or.delta_soil
	flag(6) = flag(3).or.delta_lai.or.delta_hot
	flag(7) = delta_leaf
	flag(8) = delta_skyl
	flag(9) = delta_soil
ENDIF

	init_completed=.true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	LEAF OPTICAL PROPERTIES	!!												!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!												!!
IF(flag(7)) THEN															!!
	CALL prospect_5B(N,Cab,Car,Cbrown,Cw,Cm,LRT)							!!
	rho	=	LRT(:,2)														!!
	tau	=	LRT(:,3)														!!
ENDIF																		!!
																			!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!												!!
!!	LEAF ANGLE DISTRIBUTION	!!												!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!												!!
IF(flag(2)) THEN															!!
	CALL calc_LIDF(a,lidf)													!!
ENDIF																		!!
																			!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!												!!
!!	direct / diffuse light	!!												!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!												!!
! Es = direct																!!
! Ed = diffuse																!!
IF(flag(8)) THEN															!!
	! PAR direct															!!
	PARdiro	=	(1-skyl/100.)*Es											!!
	! PAR diffus															!!
	PARdifo	=	(skyl/100.)*Ed												!!
ENDIF																		!!
																			!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!											!!
!!	Soil Reflectance Properties	!!											!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!											!!
! rsoil1 = dry soil															!!
! rsoil2 = wet soil															!!
IF(flag(9)) THEN															!!
	rsoil0=psoil*Rsoil1+(1-psoil)*Rsoil2									!!
ENDIF																		!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(flag(1)) THEN
	cts		= COS(rd*tts)
	cto		= COS(rd*tto)
	ctscto	= cts*cto
	tants	= TAN(rd*tts)
	tanto	= TAN(rd*tto)
	cospsi	= COS(rd*psi)
	dso		= SQRT(tants*tants+tanto*tanto-2.*tants*tanto*cospsi)
ENDIF

! angular distance, compensation of shadow length
IF (flag(3)) THEN
	!	Calculate geometric factors associated with extinction and scattering 
	!	Initialise sums
	ks	= 0
	ko	= 0
	bf	= 0
	sob	= 0
	sof	= 0

	!	Weighted sums over LIDF
	DO i=1,na
		ttl = thetal(i)	! leaf inclination discrete values
		ctl = COS(rd*ttl)
		!	SAIL volume scattering phase function gives interception and portions to be 
		!	multiplied by rho and tau
		CALL volscatt(tts,tto,psi,ttl,chi_s,chi_o,frho,ftau)

		!********************************************************************************
		!*                   SUITS SYSTEM COEFFICIENTS 
		!*
		!*	ks  : Extinction coefficient for direct solar flux
		!*	ko  : Extinction coefficient for direct observed flux
		!*	att : Attenuation coefficient for diffuse flux
		!*	sigb : Backscattering coefficient of the diffuse downward flux
		!*	sigf : Forwardscattering coefficient of the diffuse upward flux
		!*	sf  : Scattering coefficient of the direct solar flux for downward diffuse flux
		!*	sb  : Scattering coefficient of the direct solar flux for upward diffuse flux
		!*	vf   : Scattering coefficient of upward diffuse flux in the observed direction
		!*	vb   : Scattering coefficient of downward diffuse flux in the observed direction
		!*	w   : Bidirectional scattering coefficient
		!********************************************************************************

		!	Extinction coefficients
		ksli = chi_s/cts
		koli = chi_o/cto

		!	Area scattering coefficient fractions
		sobli	= frho*pi/ctscto
		sofli	= ftau*pi/ctscto
		bfli	= ctl*ctl
		ks	= ks+ksli*lidf(i)
		ko	= ko+koli*lidf(i)
		bf	= bf+bfli*lidf(i)
		sob	= sob+sobli*lidf(i)
		sof	= sof+sofli*lidf(i)

	ENDDO
	!	Geometric factors to be used later with rho and tau
	sdb	= 0.5*(ks+bf)
	sdf	= 0.5*(ks-bf)
	dob	= 0.5*(ko+bf)
	dof	= 0.5*(ko-bf)
	ddb	= 0.5*(1.+bf)
	ddf	= 0.5*(1.-bf)
ENDIF

IF(flag(4)) THEN
	!	Here rho and tau come in
	sigb= ddb*rho+ddf*tau
	sigf= ddf*rho+ddb*tau
	att	= 1.-sigf
	m2=(att+sigb)*(att-sigb)
	WHERE (m2.LT.0) 
		m2=0
	ENDWHERE
	m=SQRT(m2)
	sb	= sdb*rho+sdf*tau
	sf	= sdf*rho+sdb*tau
	vb	= dob*rho+dof*tau
	vf	= dof*rho+dob*tau
	w	= sob*rho+sof*tau
ENDIF

IF (flag(5)) THEN
	!	Here the LAI comes in
	!   Outputs for the case LAI = 0
	IF (lai.le.0) THEN
		tss		= 1.
		too		= 1.
		tsstoo	= 1.
		rdd		= 0.
		tdd		= 1.
		rsd		= 0.
		tsd		= 0.
		rdo		= 0.
		tdo		= 0.
		rso		= 0.
		rsos	= 0.
		rsod	= 0.

		rddt	= rsoil0
		rsdt	= rsoil0
		rdot	= rsoil0
		rsodt	= 0.
		rsost	= rsoil0
		rsot	= rsoil0
		RETURN
	ENDIF

	!	Other cases (LAI > 0)
	e1		= EXP(-m*lai)
	e2		= e1*e1
	rinf	= (att-m)/sigb
	rinf2	= rinf*rinf
	re		= rinf*e1
	denom	= 1.-rinf2*e2

	CALL Jfunc1(ks,m,lai,J1ks)
	CALL Jfunc2(ks,m,lai,J2ks)
	CALL Jfunc1(ko,m,lai,J1ko)
	CALL Jfunc2(ko,m,lai,J2ko)

	Ps = (sf+sb*rinf)*J1ks
	Qs = (sf*rinf+sb)*J2ks
	Pv = (vf+vb*rinf)*J1ko
	Qv = (vf*rinf+vb)*J2ko

	rdd	= rinf*(1.-e2)/denom
	tdd	= (1.-rinf2)*e1/denom
	tsd	= (Ps-re*Qs)/denom
	rsd	= (Qs-re*Ps)/denom
	tdo	= (Pv-re*Qv)/denom
	rdo	= (Qv-re*Pv)/denom

	tss	= EXP(-ks*lai)
	too	= EXP(-ko*lai)
	z	= Jfunc3(ks,ko,lai)
	g1	= (z-J1ks*too)/(ko+m)
	g2	= (z-J1ko*tss)/(ks+m)

	Tv1 = (vf*rinf+vb)*g1
	Tv2 = (vf+vb*rinf)*g2
	T1	= Tv1*(sf+sb*rinf)
	T2	= Tv2*(sf*rinf+sb)
	T3	= (rdo*Qs+tdo*Ps)*rinf

	!	Multiple scattering contribution to bidirectional canopy reflectance
	rsod = (T1+T2-T3)/(1.-rinf2)

ENDIF

IF (flag(6)) THEN
	!	Treatment of the hotspot-effect
	alf=1e6
	!	Apply correction 2/(K+k) suggested by F.-M. Bréon
	IF (q.gt.0.) THEN
		alf=(dso/q)*2./(ks+ko)
	ENDIF
	IF (alf.GT.200.) THEN	!inserted H. Bach 1/3/04
		alf=200.
	ENDIF
	IF (alf.eq.0.) THEN
		!	The pure hotspot - no shadow
		tsstoo = tss
		sumint = (1-tss)/(ks*lai)
	ELSE
		!	Outside the hotspot
		fhot=lai*SQRT(ko*ks)
		!	Integrate by exponential Simpson method in 20 steps
		!	the steps are arranged according to equal partitioning
		!	of the slope of the joint probability function
		x1=0.
		y1=0.
		f1=1.
		fint=(1.-EXP(-alf))*.05
		sumint=0.

		DO i=1,20
			IF (i.lt.20) THEN
				x2=-LOG(1.-i*fint)/alf
			ELSE
				x2=1.
			ENDIF
			y2=-(ko+ks)*lai*x2+fhot*(1.-EXP(-alf*x2))/alf 
			f2=EXP(y2)
			sumint=sumint+(f2-f1)*(x2-x1)/(y2-y1)
			x1=x2
			y1=y2
			f1=f2
		ENDDO
		tsstoo=f1
	ENDIF
ENDIF

!	Bidirectional reflectance
!	Single scattering contribution
	rsos = w*lai*sumint
!	Total canopy contribution
	rso=rsos+rsod

!	Interaction with the soil

dn=1.-rsoil0*rdd

rddt=rdd+tdd*rsoil0*tdd/dn
rsdt=rsd+(tsd+tss)*rsoil0*tdd/dn
rdot=rdo+tdd*rsoil0*(tdo+too)/dn

rsodt=rsod+((tss+tsd)*tdo+(tsd+tss*rsoil0*rdd)*too)*rsoil0/dn
rsost=rsos+tsstoo*rsoil0
rsot=rsost+rsodt

resh	= (rddt*PARdifo+rsdt*PARdiro)/(PARdiro+PARdifo)
resv	= (rdot*PARdifo+rsot*PARdiro)/(PARdiro+PARdifo)

	!	Before returning, save current parameters as old ones 

	N_old		= N
	Cab_old		= Cab
	Car_old		= Car
	Cbrown_old		= Cbrown
	Cw_old		= Cw
	Cm_old		= Cm
	angl_old	= a
	lai_old		= lai
	q_old		= q
	skyl_old	= skyl
	tts_old		= tts
	tto_old		= tto
	psi_old		= psi
	psoil_old	= psoil

100 FORMAT(i8,2000f15.8)

!---deallocate
DEALLOCATE(LRT,sb,sf,vb,vf)
DEALLOCATE(w,m,m2,att,sigb,rinf)
DEALLOCATE(rho, tau, PARdiro, PARdifo)
DEALLOCATE (tsd,tdd,tdo)
DEALLOCATE(rsd, rdd, rso, rdo, rddt, rsdt, rdot, rsodt, rsost)
DEALLOCATE(rsot,rsos,rsod)
DEALLOCATE (thetal,lidf)
DEALLOCATE (rsoil,rsoil0)
!---------

RETURN
END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE Jfunc1(k,l,t,Jout)

	USE dataparm
	IMPLICIT NONE

!	J1 function with avoidance of singularity problem
!	
REAL(KIND=8),INTENT(in) :: k,l(nw),t
REAL(KIND=8),INTENT(out) :: Jout(nw)
REAL(KIND=8) :: del(nw)


del=(k-l)*t
WHERE (ABS(del)>1e-3)
	Jout=(EXP(-l*t)-EXP(-k*t))/(k-l)
ELSEWHERE
	Jout=0.5*t*(EXP(-k*t)+EXP(-l*t))*(1.-del*del/12.)
END WHERE

END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE Jfunc2(k,l,t,Jout)

	USE dataparm
	IMPLICIT NONE

!	J2 function

REAL(KIND=8),INTENT(in) :: k,l(nw),t
REAL(KIND=8),INTENT(out) :: Jout(nw)
REAL(KIND=8) :: del(nw)

Jout=(1.-EXP(-(k+l)*t))/(k+l)

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


FUNCTION Jfunc3(k,l,t)

	USE dataparm
	IMPLICIT NONE

!	J2 function

REAL(KIND=8) :: k,l,t
REAL(KIND=8) :: Jfunc3

Jfunc3=(1.-EXP(-(k+l)*t))/(k+l)

END
