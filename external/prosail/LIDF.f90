!********************************************************************************
!*                          Campbell.f                            
!*     
!*    Computation of the leaf angle distribution function value (freq) 
!*    Ellipsoidal distribution function caracterised by the average leaf 
!*    inclination angle in degree (ala)                                     
!*    Campbell 1986                                                      
!*                                                                              
!********************************************************************************

SUBROUTINE campbell(ala,freq)

	REAL(KIND=8) :: ala,freq(18),dum
	REAL(KIND=8) :: excent,tl1(18),tl2(18),x1,x2,alpha
	REAL(KIND=8) :: x12,x22,almx1,almx2,alpx1,alpx2,sum
	REAL(KIND=8) :: tx1(18),tx2(18)
	REAL(KIND=8) :: x(18)
	INTEGER i,n
	DATA tx2 /0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85/
	DATA tx1 /5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90/

	n=18
	DO i=1,n
		x(i)=(tx2(i)+tx1(i))/2.
	ENDDO
	DO i=1,n
		tl1(i)=tx1(i)*ATAN(1.)/45.
		tl2(i)=tx2(i)*ATAN(1.)/45.
	ENDDO
	excent=EXP(-1.6184e-5*ala**3+2.1145e-3*ala**2-1.2390e-1*ala+3.2491)
	sum = 0.
	DO i=1,n
		x1  = excent/(SQRT(1.+excent**2*TAN(tl1(i))**2))
		x2  = excent/(SQRT(1.+excent**2*TAN(tl2(i))**2))
		IF (excent.eq.1.) THEN
			freq(i) = ABS(COS(tl1(i))-COS(tl2(i)))
		ELSE
			alpha  = excent/SQRT(ABS(1.-excent**2))
			alpha2 = alpha**2
			x12 = x1**2
			x22 = x2**2
			IF (excent.gt.1) THEN
			   alpx1 = SQRT(alpha2+x12)
			   alpx2 = SQRT(alpha2+x22)
			   dum   = x1*alpx1+alpha2*LOG(x1+alpx1)
			   freq(i) = ABS(dum-(x2*alpx2+alpha2*LOG(x2+alpx2)))
			ELSE
			   almx1 = SQRT(alpha2-x12)
			   almx2 = SQRT(alpha2-x22)
			   dum   = x1*almx1+alpha2*ASIN(x1/alpha)
			   freq(i) = ABS(dum-(x2*almx2+alpha2*ASIN(x2/alpha)))
			ENDIF
		 ENDIF
		 sum = sum+freq(i)
	  ENDDO
	DO i=1,n
	   freq(i)=freq(i)/sum	!*100.
	ENDDO

RETURN
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE calc_LIDF(alpha,freqvar)
	IMPLICIT NONE

REAL(KIND=8) :: opt,bemu,benu,alphadeg,alpha,freqvar(18),freq(18),pi
INTEGER(KIND=4) :: ia,icou,k
!...............................................................................
!     Call leaf angle distribution fonction                  
!...............................................................................
	
	alphadeg=alpha

	CALL campbell(alphadeg,freq)
	DO ia=1,18
		freqvar(ia)=freq(ia)
	ENDDO

RETURN
END SUBROUTINE
