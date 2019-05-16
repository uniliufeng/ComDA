SUBROUTINE prospect_5B(N,Cab,Car,Cbrown,Cw,Cm,LRT)

! prospect_5B (including carotenoids and brown pigments)
! version 5b (october, 20th 2009)
! subroutines required: tav, dataSpec_P5B

! ***********************************************************************
! Jacquemoud S., Baret F. (1990), PROSPECT: a model of leaf optical
! properties spectra, Remote Sens. Environ., 34:75-91.
! Féret et al. (2008), PROSPECT-4 and 5: Advances in the Leaf Optical
! Properties Model Separating Photosynthetic Pigments, Remote Sensing of
! Environment, 112:3030-3043
! The specific absorption coefficient corresponding to brown pigment is
! provided by Frederic Baret (EMMAH, INRA Avignon, baret@avignon.inra.fr)
! and used with his autorization.
! ***********************************************************************

	USE dataSpec_P5B
	IMPLICIT NONE

REAL(KIND=8),INTENT(IN) :: N,Cab,Car,Cbrown,Cw,Cm
REAL(KIND=8),INTENT(OUT) :: LRT(nl,3)

REAL(KIND=8) :: k(nl),k1(nl),k2(nl),k3(nl),k4(nl),k40(nl)
REAL(KIND=8) :: mask1(nl),mask2(nl),mask3(nl),mask4(nl)
REAL(KIND=8) :: k0(nl),y0(nl)
REAL(KIND=8) :: x(nl),y(nl)
REAL(KIND=8) :: x1(nl),x2(nl),x3(nl),x4(nl),x5(nl),x6(nl)
REAL(KIND=8) :: x_2(nl),x_4(nl),y_2(nl),y_3(nl),y_4(nl)
REAL(KIND=8) :: r(nl),t(nl),ra(nl),ta(nl)
REAL(KIND=8) :: delta(nl),beta(nl),va(nl),vb(nl),s11(nl),s12(nl),s2(nl),s3(nl)
REAL(KIND=8) :: k_2(nl),t_2(nl),r_2(nl),vbNN(nl),vbNN1(nl),va1(nl)
REAL(KIND=8) :: eps
REAL(KIND=8) :: t1(nl),t2(nl),t12(nl),nref2(nl)
REAL(KIND=8) :: ang1,ang2

INTEGER l

k=(Cab*aCab+Car*aCar+Cbrown*aBrown+Cw*aCw+Cm*aCm)/N

! Create masks:
mask1=0;mask2=0;mask3=0;mask4=0
WHERE (k.le.0.0) 
	mask1=1
ELSEWHERE (k.gt.4.0.and.k.le.85)
	mask2=1
ELSEWHERE (k.gt.85.0)
	mask3=1
ELSEWHERE (k.gt.0.0.and.k.le.4.0)
	mask4=1
ENDWHERE

k1=mask1
k2=k*mask2
k3=k*mask3
k4=k*mask4
k40=k4+mask1+mask2+mask3

! ***********************************************************************
! reflectance and transmittance of one layer
! ***********************************************************************
! Allen W.A., Gausman H.W., Richardson A.J., Thomas J.R. (1969),
! Interaction of isotropic ligth with a compact plant leaf, J. Opt.
! Soc. Am., 59(10):1376-1379.
! ***********************************************************************

! ******************************************************************
! exponential integral - NAG - S13AAF
! ******************************************************************

! small k - argument evaluation
	x_4 = 0.5 * k4 - 1.0
	y_4 = (((((((((((((((-3.60311230482612224d-13 &
		*x_4+3.46348526554087424d-12)*x_4-2.99627399604128973d-11) &
		*x_4+2.57747807106988589d-10)*x_4-2.09330568435488303d-9) &
		*x_4+1.59501329936987818d-8)*x_4-1.13717900285428895d-7) &
		*x_4+7.55292885309152956d-7)*x_4-4.64980751480619431d-6) &
		*x_4+2.63830365675408129d-5)*x_4-1.37089870978830576d-4) &
		*x_4+6.47686503728103400d-4)*x_4-2.76060141343627983d-3) &
		*x_4+1.05306034687449505d-2)*x_4-3.57191348753631956d-2) &
		*x_4+1.07774527938978692d-1)*x_4-2.96997075145080963d-1
	y_4 = (y_4*x_4+8.64664716763387311d-1)*x_4 + 7.42047691268006429d-1
	y_4 = y_4 - LOG(k40)

! large k - asymptotic test
	x_2 = 14.5 / (k2+3.25) - 1.0
	y_2 = (((((((((((((((-1.62806570868460749d-12 &
		*x_2-8.95400579318284288d-13)*x_2-4.08352702838151578d-12) &
		*x_2-1.45132988248537498d-11)*x_2-8.35086918940757852d-11) &
		*x_2-2.13638678953766289d-10)*x_2-1.10302431467069770d-9) &
		*x_2-3.67128915633455484d-9)*x_2-1.66980544304104726d-8) &
		*x_2-6.11774386401295125d-8)*x_2-2.70306163610271497d-7) &
		*x_2-1.05565006992891261d-6)*x_2-4.72090467203711484d-6) &
		*x_2-1.95076375089955937d-5)*x_2-9.16450482931221453d-5) &
		*x_2-4.05892130452128677d-4)*x_2-2.14213055000334718d-3
	y_2 = ((y_2*x_2-1.06374875116569657d-2)*x_2-8.50699154984571871d-2)*x_2 + &
		9.23755307807784058d-1

WHERE (k2.eq.0)
	y_2=0
ELSEWHERE
	y_2 = EXP(-k2) * y_2 / k2
ENDWHERE
WHERE (k4.eq.0) 
	y_4 =0
ENDWHERE

! asymptotic range
y_3 = 0.0

y0=y_2+y_3+y_4
k0=k2+k3+k4
k0=(1.-k0)*DEXP(-k0)+k0**2*y0
k=k1+k0

ang1=90;ang2=40
CALL tav_abs(ang1,nrefrac,t1)
CALL tav_abs(ang2,nrefrac,t2)
t12	=	t1*t1
nref2=nrefrac*nrefrac

k_2=k*k
x1=1-t1
x2=t12*k_2*(nref2-t1)
x3=t12*k*nref2
x4=nref2*nref2-k_2*(nref2-t1)*(nref2-t1)
x5=t2/t1
x6=x5*(t1-1)+1-t2

r=x1+x2/x4
t=x3/x4
ra=x5*r+x6
ta=x5*t
! ***********************************************************************
! reflectance and transmittance of N layers
! ***********************************************************************
! Stokes G.G. (1862), On the intensity of the light reflected from
! or transmitted through a pile of plates, Proc. Roy. Soc. Lond.,
! 11:545-556.
! ***********************************************************************

t_2		=	t*t
r_2		=	r*r
delta=(t_2-r_2-1)*(t_2-r_2-1)-4*r_2

! correction of unstability if no absorption and N=1...
eps=1e-16
WHERE (delta .le. 0)
	delta=eps
ENDWHERE
beta=(1+r_2-t_2-SQRT(delta))/(2*r)
va=(1+r_2-t_2+SQRT(delta))/(2*r)

! numerical unstabilities for VERY important absorptions
! (k>12)
! theorically impossible... (0.1cm water or 600 µg/cm² Chl...)
WHERE (beta-r.le.0)
	vb=SQRT(beta*(va-r)/(va*eps))
ELSEWHERE
	vb=SQRT(beta*(va-r)/(va*(beta-r)))
ENDWHERE

vbNN	=	vb**(N-1)
vbNN1	=	vb**(-(N-1))
va1		=	va**(-1)

s11=ra*(va*vbNN-va1*vbNN1)
s12=(ta*t-ra*r)*(vbNN-vbNN1)
s2=ta*(va-va1)
s3=va*vbNN-va1*vbNN1-r*(vbNN-vbNN1)

WHERE ((s11+s12)*s2*s3.ne. 0)
	LRT(:,2) = (s11+s12)/s3
	LRT(:,3) = s2/s3
ELSEWHERE
	LRT(:,2) = 999.
	LRT(:,3) = 999.
ENDWHERE
FORALL(l=1:nl) LRT(l,1) = l+399

END SUBROUTINE
