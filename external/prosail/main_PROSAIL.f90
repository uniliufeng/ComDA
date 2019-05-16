!*************************************************************************
!*                                                                       *
	PROGRAM main_PROSAIL

	! 10 20 2009
	! This program allows modeling reflectance data from canopy
	! - modeling leaf optical properties with PROSPECT-5 (feret et al. 2008)
	! - modeling leaf inclination distribution function with Campbell.f
	! (Ellipsoidal distribution function caracterised by the average leaf 
	! inclination angle in degree)
	! - modeling canopy reflectance with 4SAIL (Verhoef et al., 2007)

	! Verhoef et al. (2007) Unified Optical-Thermal Four-Stream Radiative
	! Transfer Theory for Homogeneous Vegetation Canopies, IEEE TRANSACTIONS 
	! ON GEOSCIENCE AND REMOTE SENSING, VOL. 45, NO. 6, JUNE 2007
	! Féret et al. (2008), PROSPECT-4 and 5: Advances in the Leaf Optical
	! Properties Model Separating Photosynthetic Pigments, REMOTE SENSING OF 
	! ENVIRONMENT
	! The specific absorption coefficient corresponding to brown pigment is
	! provided by Frederic Baret (EMMAH, INRA Avignon, baret@avignon.inra.fr)
	! and used with his autorization.
	! the model PRO4SAIL is based on a version provided by
	!	Wout Verhoef
	!	NLR	
	!	April/May 2003

!*                                                                       *
!*************************************************************************

	USE ANGLE			! defines pi & rad conversion
	USE dataparm		! defines nb of spectra
	USE lidf_parm		! angles taken into account for LADF
	USE spectral		! spectral properties of soil and atmosphere
	USE staticvar		! static variables kept in memory for optimization
	USE flag_util		! flags for optimization
	USE output_PROSPECT	! output variables of PROSPECT
	USE outSAIL			! output variables of SAIL
	USE rfltrn			! internal variables of SAIL
	USE spectrum_width_P5	! nb opf wavelength
	IMPLICIT NONE

! LEAF BIOCHEMISTRY
REAL(KIND=8) :: N,Cab,Car,Cbrown,Cw,Cm
! CANOPY
REAL(KIND=8) :: lai,angl,psoil
REAL(KIND=8) :: skyl,hspot,ihot
REAL(KIND=8) :: tts,tto,psi
REAL(KIND=8),ALLOCATABLE,SAVE :: resh(:),resv(:)

INTEGER :: i

! SPECTRUM WIDTH
!nw	= nl

! SAIL output
! resh : hemispherical reflectance
! resv : directional reflectance
ALLOCATE (resh(nl),resv(nl))


! INITIAL PARAMETERS

Cab		=	30.		! chlorophyll content (µg.cm-2) 
Car		=	10.		! carotenoid content (µg.cm-2)
Cbrown	=	0.0		! brown pigment content (arbitrary units)
Cw		=	0.015	! EWT (cm)
Cm		=	0.009	! LMA (g.cm-2)
N		=	1.2		! structure coefficient
LAI		=	2.		! leaf area index
angl	=	50.		! average leaf angle (°)
psoil	=	1.		! soil coefficient
skyl	=	70.		! % diffuse/direct radiation
hspot	=	0.2		! hot spot
ihot	=	1.0		! 
tts		=	30.		! solar zenith angle (°)
tto		=	0.		! observer zenith angle (°)
psi		=	0.		! azimuth (°)

CALL PRO4SAIL(Cab,Car,Cbrown,Cw,Cm,N,angl,lai,hspot,tts,tto,psi,psoil,skyl,resh,resv)

OPEN (unit=11,file='Refl_CAN_P5B.txt')
	DO i=1,nl
		WRITE(11,100) i+399,resh(i),resv(i)
	ENDDO
CLOSE(11)

100 FORMAT(i8,2000f15.8)

STOP
END
