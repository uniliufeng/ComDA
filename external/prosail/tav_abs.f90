!     ******************************************************************
!     TAV(teta,ref) computation of the transmittivity at the leaf 
!     surface for a given incidence solid angle. teta is the incidence
!     solid angle (in radian). The average angle that works in most 
!     cases is 40deg*pi/180. ref is the refaction index.
!
!     Jacquemoud S., 1992
! 
!     ******************************************************************
!     STERN F., 1964, Transmission of isotropic radiation across an
!     interface between two dielectrics, Appl.Opt., Vol.3, 1:111-113
!     ALLEN W.A., 1973, Transmission of isotropic light across a
!     dielectric surface in two and three dimensions, J.Opt.Soc.Am.,
!     Vol.63, 6:664-666
!	  FERET et al. (2008), PROSPECT-4 and 5: Advances in the Leaf Optical
!	  Properties Model Separating Photosynthetic Pigments, Remote Sensing of
!	  Environment
!     ******************************************************************

SUBROUTINE tav_abs(teta,refr,res)

	USE spectrum_width_P5
	IMPLICIT NONE

REAL(KIND=8),INTENT(in) :: teta,refr(nl)
REAL(KIND=8),INTENT(out) :: res(nl)
REAL(KIND=8) teta0
REAL(KIND=8) refr2(nl)
REAL(KIND=8) a(nl),b(nl),b1(nl),b2(nl),k(nl)
REAL(KIND=8) ts(nl),tp(nl),tp1(nl),tp2(nl),tp3(nl),tp4(nl),tp5(nl)
REAL(KIND=8) pi

! ANGLE CONVERSION
pi	=	ATAN(1.)*4.

IF (teta.eq.0.) THEN
	res=4.*refr/((refr+1.)*(refr+1.))
	RETURN
ENDIF

refr2=	refr*refr
a	=	(refr+1.)*(refr+1.)/2.
k	=	-((refr2-1.)*(refr2-1.))/4.
teta0=	pi*teta/180.

IF (teta0.eq.pi/2.) THEN
	b1=0.
ELSE
	b1=DSQRT((SIN(teta0)*SIN(teta0)-(refr2+1.)/2.)**2+k)
ENDIF

b2	=	SIN(teta0)*SIN(teta0)-(refr2+1.)/2.
b	=	b1-b2
ts	=	(k*k/(6.*b*b*b)+k/b-b/2.)-(k*k/(6.*a*a*a)+k/a-a/2.)
tp1	=	-2.*refr2*(b-a)/((refr2+1.)*(refr2+1.))
tp2	=	-2.*refr2*(refr2+1.)*DLOG(b/a)/((refr2-1.)*(refr2-1.))
tp3	=	refr2*(1./b-1./a)/2.
tp4	=	16.*refr2*refr2*(refr2*refr2+1.)*DLOG((2.*(refr2+1.)*b-(refr2-1.)*(refr2-1.))/ &
(2.*(refr2+1.)*a-(refr2-1.)*(refr2-1.)))/((refr2+1.)*(refr2+1.)*(refr2+1.)*(refr2-1.)*(refr2-1.))
tp5	=	16.*refr2*refr2*refr2*(1./(2.*(refr2+1.)*b-((refr2-1.)*(refr2-1.)))-1./(2.*(refr2&
+1.)*a-(refr2-1.)*(refr2-1.)))/((refr2+1.)*(refr2+1.)*(refr2+1.))
tp	=	tp1+tp2+tp3+tp4+tp5
res	=	(ts+tp)/(2.*SIN(teta0)*SIN(teta0))

RETURN
END
