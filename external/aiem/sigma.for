c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!c
c         Advanced Integral Equation Model                  c
c                                                           c
c  Program computes bistatics scattering coefficients       c
c  from 3d finitely conducting surface                      c
c  (only for single scattering)                             c
c                                                           c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!c
c                                                           c
c  Parameters:                                              c
c                                                           c
c     er: surface relative dielectric constant              c
c     kl: normalized surface correlation length             c
c     ks: normalized surface rms height                     c
c                                                           c
c     itype: select what type of surface correlation        c
c           =1  Gaussian correlation function               c
c           =2  exponential correlation function            c
c           =3  transformed exponential correlation         c
c           =4  power-law spectrum correlation function
c           =5  x-power correlated surface
c           =6  x-exponential  correlated surface
c           =7  exponential-like  correlated surface        c
c     theta: incident angle in deg                          c
c     thetas: scattering angle in deg                       c
c     phis: scattering azimuth angle in deg                 c
c     (phi=0 : incident azimuth angle in deg)               c
c                                                           c
c  -----------------------------------------------------    c
c  sigma0: array for IEM scattering coefficients            c
c            sigma(1) ---> VV polarization                  c
c            sigma(2) ---> HH polarization                  c
c            sigma(3) ---> HV polarization                  c
c            sigma(4) ---> VH polarization                  c
c                                                           c
c***********************************************************c
c        Approximations of Fresnel reflection coeff.        c
c        -------------------------------------------        c
c   Reference: IGARSS'98                                    c
c     A transition model for the reflection coefficient     c
c     in surface scattering                                 c
c                                                           c
c    irc=1: Fresnel reflection coeff. approxmated by        c
c           R(incident_angle)                               c
c    irc=2: Fresnel reflection coeff. approxmated by        c
c           R(specular_angle)                               c
c    irc=3: Fresnel reflection coeff. approxmated by        c
c           R(Transition)                                   c
c    irc=4: calculate irc=1,2,3 simultaneously              c
c                                                           c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!c
c                                                           c
c                                   Sept. 18, 2001          c
c                                                           c
c                                             T. D. Wu      c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!c 

!!-- modified by Liangxu Wang, 2010.
!!-- for mixed programming with C++.
!! theta,thetas,phis, unit: radian
!! itran is irc.
!! err,eri is the real and imag of the DIELECTRIC CONSTANT

c*******************************************************************c
c subroutine sigma calculates the scattering coefficients           c
c*******************************************************************c
      subroutine sigma(kstmp,kl,theta,thetas,phis,sigma0,itran,err,
     &eri,itype)
c	!DEC$ ATTRIBUTES DLLEXPORT::SIGMA
      implicit real*8 (a-h,k,o-z)                     !!功率普函数的类型
      integer itran,itype
      real*8 err,eri
      complex*16 er,ur,rv,rh,rvh,stem,steml,fvv,fhh,fvh,fhv
      complex*16 rvi,rhi,rvhi,rvl,rhl,rvhl
      complex*16 rvtran,rhtran,rvhtran
      complex*16 Fahh,Favv,Favh,Fahv
	complex*16 Fbhh,Fbvv,Fbvh,Fbhv
	complex*16 A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16
	complex*16 B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,B13,B14,B15,B16
	complex*16 C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16
	complex*16 D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14,D15,D16
	complex*16 Fahh1,Fahh2,Fahh3,Fahh4
	complex*16 Favv1,Favv2,Favv3,Favv4
	complex*16 Favh1,Favh2,Favh3,Favh4
	complex*16 Fahv1,Fahv2,Fahv3,Fahv4
	complex*16 Fbhh1,Fbhh2,Fbhh3,Fbhh4
	complex*16 Fbvv1,Fbvv2,Fbvv3,Fbvv4
	complex*16 Fbvh1,Fbvh2,Fbvh3,Fbvh4
	complex*16 Fbhv1,Fbhv2,Fbhv3,Fbhv4
	complex*16 kc1,kc2,kc3,kc4,kc5,kc6,kc7,kc8
	complex*16 expkc1,expkc2,expc1,expc2,expc3,expc4
	complex*16 expkcaup,expkcadn,expkcbup,expkcbdn
	complex*16 expcauau,expcadau,expcauad,expcadad
      complex*16 expcbuau,expcbdau,expcbuad,expcbdad
	complex*16 expcaubu,expcadbu,expcaubd,expcadbd
	complex*16 expcbubu,expcbdbu,expcbubd,expcbdbd
      complex*16 qq,qqs,qqt,qqts
      complex*16 kterm(4),kcterm(4),cterm(4)
	real*8 w(5000)
	dimension sigma0(4)
	common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
	common /cont2/ ks2,iterm
	common /spectra/ w
      common /eu/er,ur
!      common /surface/itype
      common /rhhvv/rh,rv,rvh
	common /tran/tranv,tranh,tranv0,tranh0,rvtran,rhtran
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real*8 pIndex, xIndex,ap,bp
	common /p/pIndex, xIndex,ap,bp
	real*8 fn,L,K
	common /cor/fn,L,K
	EXTERNAL Corr,CorrExp

	er=dcmplx(err,eri)
	ur=1.0	
c------------------------------------------------------------------c
C      The below is to set parameteres for itype= 4, 5,6,7         c
c------------------------------------------------------------------c
	if (itype .eq. 4) then  !如果是power law，则计算出ap和bp，取p＝2
           
		pIndex=2.0d0
		call gamma(pIndex-0.5d0,temp1)
		call gamma(pIndex,temp2)
		ap	=	temp1/temp2
		bp	=	compute_Bp(pIndex)
	else if (itype .eq. 5) then  !如果是x-power，取1.5
		xIndex=1.5d0              
	else if (itype .eq. 6) then
		xIndex=1.6d0
	else if (itype .eq. 7) then
		xIndex=1d0
	endif
c------------------------------------------------------------------c
C      the above is to set parameteres for itype= 4, 5,6,7         c
c------------------------------------------------------------------c 


c	write(*,*) kstmp,kl,theta,thetas,phis,itran
C==== To faciliate the compuation efficiency, some clauses are commented====C
      if (itran.eq.4) then  !!4是什么情况？为什么5指的是后向散射？
	   ks=0.00001
	else
	   ks=kstmp
	end if
C==== To faciliate the compuation efficiency, some clauses are commented====C
	torlant=1.0d-16
      cs=dcos(theta) !入射角余玄
      css=dcos(thetas) !散射角余玄
      si=sin(theta)
      sis=sin(thetas)
      csfs=cos(phis) !散射方位角余玄
      sfs=sin(phis)

	cs2=cs*cs
	css2=css*css
	si2=si*si
	sis2=sis*sis
      ks2=ks*ks
      kl2=kl*kl
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
c reflection coefficients based on the incident angle   c
c ==> R(theta)                                          c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
      stem=cdsqrt(er*ur-si2) !这里的ur是什么，前面设的是1
      rvi=(er*cs-stem)/(er*cs+stem)
      rhi=(ur*cs-stem)/(ur*cs+stem)
	rvhi=(rvi-rhi)/2.0     !论文中47页两种情况下反射系数的
	                        !开根号，相加
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
c reflection coefficients based on the specular angle   c
c ==> R(theta_specular)                                 c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
      csl=dsqrt(1.0d0+cs*css-si*sis*csfs)/1.41421356d0
	sil=dsqrt(1.0d0-csl*csl)
	steml=cdsqrt(er*ur-sil*sil)
      rvl=(er*csl-steml)/(er*csl+steml)
      rhl=(ur*csl-steml)/(ur*csl+steml)
	rvhl=(rvl-rhl)/2.0    !基于specular angle的反射系数
	                      !公式在哪???
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
c     Reflection coefficients rv, rh, rvh               c
c           for kirchhoff field coefficients            c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
      go to (101,102,103,104,105), itran
101		rh=rhi
         rv=rvi
	   rvh=rvhi
      go to 100
102      rh=rhl
         rv=rvl
	   rvh=rvhl
      go to 100
103      rh=rhl
         rv=rvl
	   rvh=rvhl
      go to 100
104      rh=rhl
         rv=rvl
	   rvh=rvhl
      go to 100
105      rh=rhtran
         rv=rvtran
	   rvh=rvhtran
100   continue
c---------------------------------------------c
c     kirchhoff field coefficients            c
c---------------------------------------------c
      fvv=2.0d0*rv*(si*sis-(1.0d0+cs*css)*csfs)/(css+cs)
      fhh=-2.0d0*rh*(si*sis-(1.0d0+cs*css)*csfs)/(css+cs)
	fvh=-2.0d0*rvh*sfs
	fhv=2.0d0*rvh*sfs
c---------------------------------------------c
c compute roughness spectrum w(n)             c
c---------------------------------------------c
      iterm=1
	tempold=0.0
c	temp=(ks2*(cs+css)**2.0)
c	do while(dabs(temp-tempold).gt.torlant)
c		tempold=temp
c	   iterm=iterm+1
c	   fiterm=float(iterm)
c        temp=tempold*(ks2*(cs+css)**2.0)/fiterm    
c	end do
	iterm=2.72*(ks2*(cs+css)**2.0)+136
c     if (itran.eq.4) iterm=1

C==== To faciliate the compuation efficiency, some clauses are commented====C
c      do 109,n=1,1000
c      w(n)=0.0
c109   continue
C==== To faciliate the compuation efficiency, some clauses are commented====C

	ql=kl*dsqrt((-si+sis*csfs)**2.0+(0.+sis*sfs)**2.0)
c-----Gaussian correlated surface ------------------------------------------c
      if( itype.eq.1 ) then
		do 110,n=1,iterm
			fn=float(n)  
			w(n)=kl2*dexp(-ql*ql/(4.0*fn))/(2*fn)
110		continue
c-----Exponential correlated surface ---------------------------------------c
      else if( itype.eq.2 ) then
		do 111,n=1,iterm
			fn=float(n)  
			w(n)=(kl/fn)**2.0*(1.0+(ql/fn)**2.0)**(-1.5)
111			continue
c-----1.5-power correlated surface -----------------------------------------c
      else if( itype.eq.3 ) then  
		do 112,n=1,iterm
			fn=float(n)
			e = 1.5 * fn - 1.0d0
			y = 1.5 * fn
			gam = ( alogam(y) ) !gamma function (1.5n)
			if ( ql .eq. 0.0d0 ) then
				w(n) = kl*kl/(3.0*fn-2.0)
			endif
			if ( dmod( fn, 2.0d0) .eq. 0.0d0) then
				m = 1.5 * fn - 1.0            !fn is even;
				bk = log(BESSK( m, ql))       !integer order, check OK
			else
				m = 1.5 * fn - 1.0  - 0.5     !fn is odd
				bk = log(BesselK( m, ql))     !fractional order
			endif
			out = kl * kl * ( ql / 2.0d0) ** e
			w(n)=  out * dexp(bk - gam ) 
112		continue
c-----Generalized power law sprectrum -------------------------------------c
	else if( itype.eq.4 ) then  
		do 113,n=1,iterm
			fn=float(n)
			fp=0.5*(1+(1.5/pIndex)**2)
			apBp2=(ap/bp)**2
			fpn2=(fn**fp)**2
			w(n)=0.5*kl2/fpn2*(pIndex-1)*apBp2*
     &			(1+apBp2*ql*ql/(4*fpn2))**(-pIndex)
113		continue
c-----X-power correlation function ----------------------------------------c
	else if( itype.eq.5 ) then  
		do 114,n=1,iterm
			fn=float(n)
			u=xIndex*fn-1.0
	
			if (u .lt. 90d0 ) then
				call CIKVA(u,dcmplx(ql,0d0),VM,Res)	
			else
				call CIKLV(u,dcmplx(ql,0d0),Res)		
			endif
			call gamma(fn*xIndex, R_gamma)      
			if (IsNAN(Res) ) then
 				w(n)=0.0d0
			else
				w(n)=kl*kl*(ql/2.0d0)**u/R_gamma*Res
			endif
114		continue
c-----X expoential correlation function ----------------------------------c
	else if( itype.eq.6 ) then
		do 115,n=1,iterm
			fn=float(n)
			L=kl
			K=ql/kl	
	call quanc8(Corr,0d0,1.0d03,0.001d0,0.001d0,result,ert,nofun,fg)
			w(n)=result/1.0d305
115		continue
c-----Expoential-like correlation function -----------c
	else if( itype.eq.7 ) then
		do 116,n=1,iterm
			fn=float(n)
			K=ql/kl
C	w(n)=W_expLike(fn,xIndex,kl,K)
			L=kl
			K=ql/kl
	call quanc8(CorrExp,0d0,1.d4,1d-3,1d-3,result,ert,nofun,fg)
			if (isnan(RESULT)) then
				w(n)=0
			else
				w(n)=RESULT/1d305
			endif
116		continue

      endif	
	if (w(1) .lt. 0 ) then
		continue
	end if

c------------------------------------------c
c    compute kirchhoff term                c
c------------------------------------------c
      sum=0.0
	temp=1.0
	temp0=(ks2*(cs+css)**2.0)
      do 3 n=1,iterm
	   fn=float(n)
         temp=temp*temp0/fn
	   sum=sum+temp*w(n)
3     continue  
      expk=dexp(-ks2*(css+cs)**2.0)*sum
      kterm(1)=0.5*expk*cdabs(fvv)**2
      kterm(2)=0.5*expk*cdabs(fhh)**2
      kterm(3)=0.5*expk*cdabs(fhv)**2
      kterm(4)=0.5*expk*cdabs(fvh)**2
c
c-------------------------------------------------------c
c  end of kirchhoff term computation                    c
c-------------------------------------------------------c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
c     Reflection coefficients rv, rh, rvh               c
c           for complementary field coefficients        c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
        go to (201,202,203,204,205), itran
201          rh=rhi
             rv=rvi
	       rvh=rvhi
        go to 200
202          rh=rhl
             rv=rvl
	       rvh=rvhl
        go to 200
203          rh=rhl
             rv=rvl
	       rvh=rvhl
        go to 200
204          rh=rhl
             rv=rvl
	       rvh=rvhl
        go to 200
205          rh=rhtran
             rv=rvtran
	       rvh=rvhtran
200   continue
c
c--------------------------------------------------------------c
       ex1=dexp(-ks2*(cs2+css2+cs*css))
       ex2=dexp(-ks2*(cs2+css2))
       qq=cs
       qqt=cdsqrt(er-si2)
       qqs=css
       qqts=cdsqrt(er-sis2)
c--------------------------------------------------------------c
c------------------------------------------c
c   compute cross term                     c
c------------------------------------------c

      Favv1=Favv(-si,0.0d0,qq,qq)
	Favv2=Favv(-sis*csfs,-sis*sfs,qqs,qqs)
      Favv3=Favv(-si,0.0d0,-qq,qq)
	Favv4=Favv(-sis*csfs,-sis*sfs,-qqs,qqs)
      Fbvv1=Fbvv(-si,0.0d0,qqt,qqt)
	Fbvv2=Fbvv(-sis*csfs,-sis*sfs,qqts,qqts)
      Fbvv3=Fbvv(-si,0.0d0,-qqt,qqt)
	Fbvv4=Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts)
      Fahh1=Fahh(-si,0.0d0,qq,qq)
	Fahh2=Fahh(-sis*csfs,-sis*sfs,qqs,qqs)
      Fahh3=Fahh(-si,0.0d0,-qq,qq)
	Fahh4=Fahh(-sis*csfs,-sis*sfs,-qqs,qqs)
      Fbhh1=Fbhh(-si,0.0d0,qqt,qqt)
	Fbhh2=Fbhh(-sis*csfs,-sis*sfs,qqts,qqts)
      Fbhh3=Fbhh(-si,0.0d0,-qqt,qqt)
	Fbhh4=Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts)
      Fahv1=Fahv(-si,0.0d0,qq,qq)
	Fahv2=Fahv(-sis*csfs,-sis*sfs,qqs,qqs)
      Fahv3=Fahv(-si,0.0d0,-qq,qq)
	Fahv4=Fahv(-sis*csfs,-sis*sfs,-qqs,qqs)
      Fbhv1=Fbhv(-si,0.0d0,qqt,qqt)
      Fbhv2=Fbhv(-sis*csfs,-sis*sfs,qqts,qqts)
      Fbhv3=Fbhv(-si,0.0d0,-qqt,qqt)
	Fbhv4=Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts)
      Favh1=Favh(-si,0.0d0,qq,qq)
	Favh2=Favh(-sis*csfs,-sis*sfs,qqs,qqs)
      Favh3=Favh(-si,0.0d0,-qq,qq)
	Favh4=Favh(-sis*csfs,-sis*sfs,-qqs,qqs)
      Fbvh1=Fbvh(-si,0.0d0,qqt,qqt)
	Fbvh2=Fbvh(-sis*csfs,-sis*sfs,qqts,qqts)
      Fbvh3=Fbvh(-si,0.0d0,-qqt,qqt)
	Fbvh4=Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts)

	kc1=expkc1(qq)
	kc2=expkc2(qqs)
	kc3=expkc1(-qq)
	kc4=expkc2(-qqs)
	kc5=expkc1(qqt)
	kc6=expkc2(qqts)
	kc7=expkc1(-qqt)
	kc8=expkc2(-qqts)


      expkcaup=dreal(dconjg(fvv)*Favv1)*kc1
     & +dreal(dconjg(fvv)*Favv2)
     & *kc2
      expkcadn=dreal(dconjg(fvv)*Favv3)*kc3
     & +dreal(dconjg(fvv)*Favv4)
     & *kc4
      expkcbup=dreal(dconjg(fvv)*Fbvv1)*kc5
     & +dreal(dconjg(fvv)*Fbvv2)
     & *kc6
      expkcbdn=dreal(dconjg(fvv)*Fbvv3)*kc7
     & +dreal(dconjg(fvv)*Fbvv4)
     & *kc8

	kcterm(1)=0.25*ex1*(expkcaup+expkcadn+expkcbup+expkcbdn)

      expkcaup=dreal(dconjg(fhh)*Fahh1)* kc1 
     & +dreal(dconjg(fhh)*Fahh2)
     & * kc2 
      expkcadn=dreal(dconjg(fhh)*Fahh3)* kc3
     & +dreal(dconjg(fhh)*Fahh4)
     & * kc4
      expkcbup=dreal(dconjg(fhh)*Fbhh1)* kc5
     & +dreal(dconjg(fhh)*Fbhh2)
     & * kc6
      expkcbdn=dreal(dconjg(fhh)*Fbhh3)*kc7
     & +dreal(dconjg(fhh)*Fbhh4)
     & *kc8
	kcterm(2)=0.25*ex1*(expkcaup+expkcadn+expkcbup+expkcbdn)

      expkcaup=dreal(dconjg(fhv)*Fahv1)*kc1
     & +dreal(dconjg(fhv)*Fahv2)
     & *kc2
      expkcadn=dreal(dconjg(fhv)*Fahv3)*kc3
     & +dreal(dconjg(fhv)*Fahv4)
     & *kc4
      expkcbup=dreal(dconjg(fhv)*Fbhv1)*kc5
     & +dreal(dconjg(fhv)*Fbhv2)
     & *kc6
      expkcbdn=dreal(dconjg(fhv)*Fbhv3)*kc7
     & +dreal(dconjg(fhv)*Fbhv4)
     & *kc8
	kcterm(3)=0.25*ex1*(expkcaup+expkcadn+expkcbup+expkcbdn)


      expkcaup=dreal(dconjg(fvh)*Favh1)*kc1
     & +dreal(dconjg(fvh)*Favh2)
     & *kc2
      expkcadn=dreal(dconjg(fvh)*Favh3)*kc3
     & +dreal(dconjg(fvh)*Favh4)
     & *kc4
      expkcbup=dreal(dconjg(fvh)*Fbvh1)*kc5
     & +dreal(dconjg(fvh)*Fbvh2)
     & *kc6
      expkcbdn=dreal(dconjg(fvh)*Fbvh3)*kc7
     & +dreal(dconjg(fvh)*Fbvh4)
     & *kc8
	kcterm(4)=0.25*ex1*(expkcaup+expkcadn+expkcbup+expkcbdn)
c
c--------------------------------------c
c end of computation of cross terms    c
c--------------------------------------c
c
c--------------------------------------c
c evaluate  complementary term         c
c--------------------------------------c
	A1=expc1(qq,qq)              
	A2=expc2(qq,qqs)             
	A3=expc3(qqs,qq)             
	A4=expc4(qqs,qqs)            
                                     
	A5=expc1(-qq,qq)             
	A6=expc2(-qq,qqs)            
	A7=expc3(-qqs,qq)            
	A8=expc4(-qqs,qqs)           
                                     
	A9=expc1(qq,-qq)             
	A10=expc2(qq,-qqs)           
	A11=expc3(qqs,-qq)           
	A12=expc4(qqs,-qqs)          
                                     
	A13=expc1(-qq,-qq)           
	A14=expc2(-qq,-qqs)          
	A15=expc3(-qqs,-qq)          
	A16=expc4(-qqs,-qqs)         
                                     


	B1=expc1(qqt,qq)          
	B2=expc2(qqt,qqs)         
	B3=expc3(qqts,qq)         
	B4=expc4(qqts,qqs)        
                                  
	B5=expc1(-qqt,qq)         
	B6=expc2(-qqt,qqs)        
	B7=expc3(-qqts,qq)        
	B8=expc4(-qqts,qqs)       
                                  
	B9=expc1(qqt,-qq)         
	B10=expc2(qqt,-qqs)       
	B11=expc3(qqts,-qq)       
	B12=expc4(qqts,-qqs)      
                                  
	B13=expc1(-qqt,-qq)       
	B14=expc2(-qqt,-qqs)      
	B15=expc3(-qqts,-qq)      
	B16=expc4(-qqts,-qqs)     


	C1=expc1(qq,qqt)                                    
	C2=expc2(qq,qqts) 
	C3=expc3(qqs,qqt) 
	C4=expc4(qqs,qqts)

	C5=expc1(-qq,qqt)                                   
	C6=expc2(-qq,qqts) 
	C7=expc3(-qqs,qqt) 
	C8=expc4(-qqs,qqts) 

	C9=expc1(qq,-qqt)                                   
	C10=expc2(qq,-qqts) 
	C11=expc3(qqs,-qqt) 
	C12=expc4(qqs,-qqts) 

	C13=expc1(-qq,-qqt)                                   
	C14=expc2(-qq,-qqts) 
	C15=expc3(-qqs,-qqt) 
	C16=expc4(-qqs,-qqts) 



	D1=expc1(qqt,qqt)        
	D2=expc2(qqt,qqts)       
	D3=expc3(qqts,qqt)       
	D4=expc4(qqts,qqts)      
                                 
	D5=expc1(-qqt,qqt)       
	D6=expc2(-qqt,qqts)      
	D7=expc3(-qqts,qqt)      
	D8=expc4(-qqts,qqts)     
                                 
	D9=expc1(qqt,-qqt)       
	D10=expc2(qqt,-qqts)     
	D11=expc3(qqts,-qqt)     
	D12=expc4(qqts,-qqts)    
                                 
	D13=expc1(-qqt,-qqt)     
	D14=expc2(-qqt,-qqts)    
	D15=expc3(-qqts,-qqt)    
	D16=expc4(-qqts,-qqts)   


      expcauau=Favv1                                                                                                                                     
     &  *dconjg(Favv1)*A1                                                                                                                      
     &  +Favv1                                                                                                                                           
     &  *dconjg(Favv2)*A2                                                                                                                     
     &  +Favv2                                                                                                                                           
     &  *dconjg(Favv1)*A3                                                                                                                    
     &  +Favv2                                                                                                                                           
     &  *dconjg(Favv2)*A4                                                                                                                    
      expcadau=Favv3                                                                                                                                     
     &  *dconjg(Favv1)*A5                                                                                                                     
     &  +Favv3                                                                                                                                           
     &  *dconjg(Favv2)*A6                                                                                                                    
     &  +Favv4                                                                                                                                           
     &  *dconjg(Favv1)*A7                                                                                                                    
     &  +Favv4                                                                                                                                           
     &  *dconjg(Favv2)*A8                                                                                                                   
      expcauad=Favv1                                                                                                                                     
     &  *dconjg(Favv3)*A9                                                                                                                     
     &  +Favv1                                                                                                                                  
     &  *dconjg(Favv4)*A10                                                                                               
     &  +Favv2                                                                                                         
     &  *dconjg(Favv3)*A11                                  
     &  +Favv2                                                         
     &  *dconjg(Favv4)*A12                                 
      expcadad=Favv3                                                   
     &  *dconjg(Favv3)*A13                                  
     &  +Favv3                                                         
     &  *dconjg(Favv4)*A14                                 
     &  +Favv4                                                         
     &  *dconjg(Favv3)*A15                                 
     &  +Favv4                                                         
     &  *dconjg(Favv4)*A16                                
      expcbuau=Fbvv1                                                   
     &  *dconjg(Favv1)*B1                                   
     &  +Fbvv1                                                         
     &  *dconjg(Favv2)*B2                                  
     &  +Fbvv2                                                         
     &  *dconjg(Favv1)*B3                                  
     &  +Fbvv2                                                         
     &  *dconjg(Favv2)*B4                                
      expcbdau=Fbvv3                                                   
     &  *dconjg(Favv1)*B5                                  
     &  +Fbvv3                                                         
     &  *dconjg(Favv2)*B6                                 
     &  +Fbvv4                                                         
     &  *dconjg(Favv1)*B7                                 
     &  +Fbvv4                                                         
     &  *dconjg(Favv2)*B8                                
      expcbuad=Fbvv1                                                   
     &  *dconjg(Favv3)*B9                                  
     &  +Fbvv1                                                         
     &  *dconjg(Favv4)*B10                                 
     &  +Fbvv2                                                         
     &  *dconjg(Favv3)*B11                                 
     &  +Fbvv2                                                         
     &  *dconjg(Favv4)*B12                                
      expcbdad=Fbvv3                                                   
     &  *dconjg(Favv3)*B13                                 
     &  +Fbvv3                                                         
     &  *dconjg(Favv4)*B14                                
     &  +Fbvv4                                                         
     &  *dconjg(Favv3)*B15                                
     &  +Fbvv4                                                         
     &  *dconjg(Favv4)*B16                               
      expcaubu=Favv1                                                   
     &  *dconjg(Fbvv1)*C1                                   
     &  +Favv1                                                         
     &  *dconjg(Fbvv2)*C2                                  
     &  +Favv2                                                         
     &  *dconjg(Fbvv1)*C3                                  
     &  +Favv2                                                         
     &  *dconjg(Fbvv2)*C4                                 
      expcadbu=Favv3                                                   
     &  *dconjg(Fbvv1)*C5                                  
     &  +Favv3                                                         
     &  *dconjg(Fbvv2)*C6                                 
     &  +Favv4                                                         
     &  *dconjg(Fbvv1)*C7                                 
     &  +Favv4                                                         
     &  *dconjg(Fbvv2)*C8                                
      expcaubd=Favv1                                                   
     &  *dconjg(Fbvv3)*C9                                  
     &  +Favv1                                                         
     &  *dconjg(Fbvv4)*C10                                 
     &  +Favv2                                                         
     &  *dconjg(Fbvv3)*C11                                 
     &  +Favv2                                                         
     &  *dconjg(Fbvv4)*C12                                
      expcadbd=Favv3                                                   
     &  *dconjg(Fbvv3)*C13                                 
     &  +Favv3                                                         
     &  *dconjg(Fbvv4)*C14                                
     &  +Favv4                                                         
     &  *dconjg(Fbvv3)*C15                                
     &  +Favv4                                                         
     &  *dconjg(Fbvv4)*C16                               
      expcbubu=Fbvv1                                                   
     &  *dconjg(Fbvv1)*D1                                  
     &  +Fbvv1                                                         
     &  *dconjg(Fbvv2)*D2                                 
     &  +Fbvv2                                                         
     &  *dconjg(Fbvv1)*D3                                 
     &  +Fbvv2                                                         
     &  *dconjg(Fbvv2)*D4                                
      expcbdbu=Fbvv3                                                   
     &  *dconjg(Fbvv1)*D5                                 
     &  +Fbvv3                                                         
     &  *dconjg(Fbvv2)*D6                                
     &  +Fbvv4                                                         
     &  *dconjg(Fbvv1)*D7                                
     &  +Fbvv4                                                         
     &  *dconjg(Fbvv2)*D8                               
      expcbubd=Fbvv1                                                   
     &  *dconjg(Fbvv3)*D9                                 
     &  +Fbvv1                                                         
     &  *dconjg(Fbvv4)*D10                                
     &  +Fbvv2                                                         
     &  *dconjg(Fbvv3)*D11                                
     &  +Fbvv2                                                         
     &  *dconjg(Fbvv4)*D12                               
      expcbdbd=Fbvv3                                                   
     &  *dconjg(Fbvv3)*D13                                
     &  +Fbvv3                                                         
     &  *dconjg(Fbvv4)*D14                               
     &  +Fbvv4                                                         
     &  *dconjg(Fbvv3)*D15                               
     &  +Fbvv4                                                         
     &  *dconjg(Fbvv4)*D16                              
      cterm(1)=0.03125*ex2*(expcauau+expcadau+expcauad+expcadad        
     &                     +expcbuau+expcbdau+expcbuad+expcbdad        
     &                     +expcaubu+expcadbu+expcaubd+expcadbd        
     &                     +expcbubu+expcbdbu+expcbubd+expcbdbd)       
c                                                                      
      expcauau=Fahh1                                                   
     &  *dconjg(Fahh1)*A1                                    
     &  +Fahh1                                                         
     &  *dconjg(Fahh2)*A2                                   
     &  +Fahh2                                                         
     &  *dconjg(Fahh1)*A3                                  
     &  +Fahh2                                                         
     &  *dconjg(Fahh2)*A4                                  
      expcadau=Fahh3                                                   
     &  *dconjg(Fahh1)*A5                                   
     &  +Fahh3                                                         
     &  *dconjg(Fahh2)*A6                                  
     &  +Fahh4                                                         
     &  *dconjg(Fahh1)*A7                                  
     &  +Fahh4                                                         
     &  *dconjg(Fahh2)*A8                                 
      expcauad=Fahh1                                                   
     &  *dconjg(Fahh3)*A9                                   
     &  +Fahh1                                                         
     &  *dconjg(Fahh4)*A10                                  
     &  +Fahh2                                                         
     &  *dconjg(Fahh3)*A11                                  
     &  +Fahh2                                                         
     &  *dconjg(Fahh4)*A12                                 
      expcadad=Fahh3                                                   
     &  *dconjg(Fahh3)*A13                                  
     &  +Fahh3                                                         
     &  *dconjg(Fahh4)*A14                                 
     &  +Fahh4                                                         
     &  *dconjg(Fahh3)*A15                                 
     &  +Fahh4                                                         
     &  *dconjg(Fahh4)*A16                                
      expcbuau=Fbhh1                                                   
     &  *dconjg(Fahh1)*B1                                   
     &  +Fbhh1                                                         
     &  *dconjg(Fahh2)*B2                                  
     &  +Fbhh2                                                         
     &  *dconjg(Fahh1)*B3                                  
     &  +Fbhh2                                                         
     &  *dconjg(Fahh2)*B4                                
      expcbdau=Fbhh3                                                   
     &  *dconjg(Fahh1)*B5                                  
     &  +Fbhh3                                                         
     &  *dconjg(Fahh2)*B6                                 
     &  +Fbhh4                                                         
     &  *dconjg(Fahh1)*B7                                 
     &  +Fbhh4                                                         
     &  *dconjg(Fahh2)*B8                                
      expcbuad=Fbhh1                                                   
     &  *dconjg(Fahh3)*B9                                  
     &  +Fbhh1                                                         
     &  *dconjg(Fahh4)*B10                                 
     &  +Fbhh2                                                         
     &  *dconjg(Fahh3)*B11                                 
     &  +Fbhh2                                                         
     &  *dconjg(Fahh4)*B12                                
      expcbdad=Fbhh3                                                   
     &  *dconjg(Fahh3)*B13                                 
     &  +Fbhh3                                                         
     &  *dconjg(Fahh4)*B14                                
     &  +Fbhh4                                                         
     &  *dconjg(Fahh3)*B15                                
     &  +Fbhh4                                                         
     &  *dconjg(Fahh4)*B16                               
      expcaubu=Fahh1                                                   
     &  *dconjg(Fbhh1)*C1                                   
     &  +Fahh1                                                         
     &  *dconjg(Fbhh2)*C2                                  
     &  +Fahh2                                                         
     &  *dconjg(Fbhh1)*C3                                  
     &  +Fahh2                                                         
     &  *dconjg(Fbhh2)*C4                                 
      expcadbu=Fahh3                                                   
     &  *dconjg(Fbhh1)*C5                                  
     &  +Fahh3                                                         
     &  *dconjg(Fbhh2)*C6                                 
     &  +Fahh4                                                         
     &  *dconjg(Fbhh1)*C7                                 
     &  +Fahh4                                                         
     &  *dconjg(Fbhh2)*C8                                
      expcaubd=Fahh1                                                   
     &  *dconjg(Fbhh3)*C9                                  
     &  +Fahh1                                                         
     &  *dconjg(Fbhh4)*C10                                 
     &  +Fahh2                                                         
     &  *dconjg(Fbhh3)*C11                                 
     &  +Fahh2                                                         
     &  *dconjg(Fbhh4)*C12                                
      expcadbd=Fahh3                                                   
     &  *dconjg(Fbhh3)*C13                                 
     &  +Fahh3                                                         
     &  *dconjg(Fbhh4)*C14                                
     &  +Fahh4                                                         
     &  *dconjg(Fbhh3)*C15                                
     &  +Fahh4                                                         
     &  *dconjg(Fbhh4)*C16                               
      expcbubu=Fbhh1                                                   
     &  *dconjg(Fbhh1)*D1                                  
     &  +Fbhh1                                                         
     &  *dconjg(Fbhh2)*D2                                 
     &  +Fbhh2                                                         
     &  *dconjg(Fbhh1)*D3                                 
     &  +Fbhh2                                                         
     &  *dconjg(Fbhh2)*D4                                
      expcbdbu=Fbhh3                                                   
     &  *dconjg(Fbhh1)*D5                                 
     &  +Fbhh3                                                         
     &  *dconjg(Fbhh2)*D6                                
     &  +Fbhh4                                                         
     &  *dconjg(Fbhh1)*D7                                
     &  +Fbhh4                                                         
     &  *dconjg(Fbhh2)*D8                               
      expcbubd=Fbhh1                                                   
     &  *dconjg(Fbhh3)*D9                                 
     &  +Fbhh1                                                         
     &  *dconjg(Fbhh4)*D10                                
     &  +Fbhh2                                                         
     &  *dconjg(Fbhh3)*D11                                
     &  +Fbhh2                                                         
     &  *dconjg(Fbhh4)*D12                               
      expcbdbd=Fbhh3                                                   
     &  *dconjg(Fbhh3)*D13                                
     &  +Fbhh3                                                         
     &  *dconjg(Fbhh4)*D14                               
     &  +Fbhh4                                                         
     &  *dconjg(Fbhh3)*D15                               
     &  +Fbhh4                                                         
     &  *dconjg(Fbhh4)*D16                              
      cterm(2)=0.03125*ex2*(expcauau+expcadau+expcauad+expcadad        
     &                     +expcbuau+expcbdau+expcbuad+expcbdad        
     &                     +expcaubu+expcadbu+expcaubd+expcadbd        
     &                     +expcbubu+expcbdbu+expcbubd+expcbdbd)       
c                                                                      
      expcauau=Fahv1                                                   
     &  *dconjg(Fahv1)*A1                                    
     &  +Fahv1                                                         
     &  *dconjg(Fahv2)*A2                                   
     &  +Fahv2                                                         
     &  *dconjg(Fahv1)*A3                                  
     &  +Fahv2                                                         
     &  *dconjg(Fahv2)*A4                                  
      expcadau=Fahv3                                                   
     &  *dconjg(Fahv1)*A5                                   
     &  +Fahv3                                                         
     &  *dconjg(Fahv2)*A6                                  
     &  +Fahv4                                                         
     &  *dconjg(Fahv1)*A7                                  
     &  +Fahv4                                                         
     &  *dconjg(Fahv2)*A8                                 
      expcauad=Fahv1                                                   
     &  *dconjg(Fahv3)*A9                                   
     &  +Fahv1                                                         
     &  *dconjg(Fahv4)*A10                                  
     &  +Fahv2                                                         
     &  *dconjg(Fahv3)*A11                                  
     &  +Fahv2                                                         
     &  *dconjg(Fahv4)*A12                                 
      expcadad=Fahv3                                                   
     &  *dconjg(Fahv3)*A13                                  
     &  +Fahv3                                                         
     &  *dconjg(Fahv4)*A14                                 
     &  +Fahv4                                                         
     &  *dconjg(Fahv3)*A15                                 
     &  +Fahv4                                                         
     &  *dconjg(Fahv4)*A16                                
      expcbuau=Fbhv1                                                   
     &  *dconjg(Fahv1)*B1                                   
     &  +Fbhv1                                                         
     &  *dconjg(Fahv2)*B2                                  
     &  +Fbhv2                                                         
     &  *dconjg(Fahv1)*B3                                  
     &  +Fbhv2                                                         
     &  *dconjg(Fahv2)*B4                                
      expcbdau=Fbhv3                                                   
     &  *dconjg(Fahv1)*B5                                  
     &  +Fbhv3                                                         
     &  *dconjg(Fahv2)*B6                                 
     &  +Fbhv4                                                         
     &  *dconjg(Fahv1)*B7                                 
     &  +Fbhv4                                                         
     &  *dconjg(Fahv2)*B8                                
      expcbuad=Fbhv1                                                   
     &  *dconjg(Fahv3)*B9                                  
     &  +Fbhv1                                                         
     &  *dconjg(Fahv4)*B10                                 
     &  +Fbhv2                                                         
     &  *dconjg(Fahv3)*B11                                 
     &  +Fbhv2                                                         
     &  *dconjg(Fahv4)*B12                                
      expcbdad=Fbhv3                                                   
     &  *dconjg(Fahv3)*B13                                 
     &  +Fbhv3                                                         
     &  *dconjg(Fahv4)*B14                                
     &  +Fbhv4                                                         
     &  *dconjg(Fahv3)*B15                                
     &  +Fbhv4                                                         
     &  *dconjg(Fahv4)*B16                               
      expcaubu=Fahv1                                                   
     &  *dconjg(Fbhv1)*C1                                   
     &  +Fahv1                                                         
     &  *dconjg(Fbhv2)*C2                                  
     &  +Fahv2                                                         
     &  *dconjg(Fbhv1)*C3                                  
     &  +Fahv2                                                         
     &  *dconjg(Fbhv2)*C4                                 
      expcadbu=Fahv3                                                   
     &  *dconjg(Fbhv1)*C5                                  
     &  +Fahv3                                                         
     &  *dconjg(Fbhv2)*C6                                 
     &  +Fahv4                                                         
     &  *dconjg(Fbhv1)*C7                                 
     &  +Fahv4                                                         
     &  *dconjg(Fbhv2)*C8                                
      expcaubd=Fahv1                                                   
     &  *dconjg(Fbhv3)*C9                                  
     &  +Fahv1                                                         
     &  *dconjg(Fbhv4)*C10                                 
     &  +Fahv2                                                         
     &  *dconjg(Fbhv3)*C11                                 
     &  +Fahv2                                                         
     &  *dconjg(Fbhv4)*C12                                
      expcadbd=Fahv3                                                   
     &  *dconjg(Fbhv3)*C13                                 
     &  +Fahv3                                                         
     &  *dconjg(Fbhv4)*C14                                
     &  +Fahv4                                                         
     &  *dconjg(Fbhv3)*C15                                
     &  +Fahv4                                                         
     &  *dconjg(Fbhv4)*C16                               
      expcbubu=Fbhv1                                                   
     &  *dconjg(Fbhv1)*D1                                  
     &  +Fbhv1                                                         
     &  *dconjg(Fbhv2)*D2                                 
     &  +Fbhv2                                                         
     &  *dconjg(Fbhv1)*D3                                 
     &  +Fbhv2                                                         
     &  *dconjg(Fbhv2)*D4                                
      expcbdbu=Fbhv3                                                   
     &  *dconjg(Fbhv1)*D5                                 
     &  +Fbhv3                                                         
     &  *dconjg(Fbhv2)*D6                                
     &  +Fbhv4                                                         
     &  *dconjg(Fbhv1)*D7                                
     &  +Fbhv4                                                         
     &  *dconjg(Fbhv2)*D8                               
      expcbubd=Fbhv1                                                   
     &  *dconjg(Fbhv3)*D9                                 
     &  +Fbhv1                                                         
     &  *dconjg(Fbhv4)*D10                                
     &  +Fbhv2                                                         
     &  *dconjg(Fbhv3)*D11                                
     &  +Fbhv2                                                         
     &  *dconjg(Fbhv4)*D12                               
      expcbdbd=Fbhv3                                                   
     &  *dconjg(Fbhv3)*D13                                
     &  +Fbhv3                                                         
     &  *dconjg(Fbhv4)*D14                               
     &  +Fbhv4                                                         
     &  *dconjg(Fbhv3)*D15                               
     &  +Fbhv4                                                         
     &  *dconjg(Fbhv4)*D16                              
      cterm(3)=0.03125*ex2*(expcauau+expcadau+expcauad+expcadad        
     &                     +expcbuau+expcbdau+expcbuad+expcbdad        
     &                     +expcaubu+expcadbu+expcaubd+expcadbd        
     &                     +expcbubu+expcbdbu+expcbubd+expcbdbd)       
c                                                                      
      expcauau=Favh1                                                   
     &  *dconjg(Favh1)*A1                                    
     &  +Favh1                                                         
     &  *dconjg(Favh2)*A2                                   
     &  +Favh2                                                         
     &  *dconjg(Favh1)*A3                                  
     &  +Favh2                                                         
     &  *dconjg(Favh2)*A4                                  
      expcadau=Favh3                                                   
     &  *dconjg(Favh1)*A5                                   
     &  +Favh3                                                         
     &  *dconjg(Favh2)*A6                                  
     &  +Favh4                                                         
     &  *dconjg(Favh1)*A7                                  
     &  +Favh4                                                         
     &  *dconjg(Favh2)*A8                                 
      expcauad=Favh1                                                   
     &  *dconjg(Favh3)*A9                                   
     &  +Favh1                                                         
     &  *dconjg(Favh4)*A10                                  
     &  +Favh2                                                         
     &  *dconjg(Favh3)*A11                                  
     &  +Favh2                                                         
     &  *dconjg(Favh4)*A12                                 
      expcadad=Favh3                                                   
     &  *dconjg(Favh3)*A13                                  
     &  +Favh3                                                         
     &  *dconjg(Favh4)*A14                                 
     &  +Favh4                                                         
     &  *dconjg(Favh3)*A15                                 
     &  +Favh4                                                         
     &  *dconjg(Favh4)*A16                                
      expcbuau=Fbvh1                                                   
     &  *dconjg(Favh1)*B1                                   
     &  +Fbvh1                                                         
     &  *dconjg(Favh2)*B2                                  
     &  +Fbvh2                                                         
     &  *dconjg(Favh1)*B3                                  
     &  +Fbvh2                                                         
     &  *dconjg(Favh2)*B4                                
      expcbdau=Fbvh3                                                   
     &  *dconjg(Favh1)*B5                                  
     &  +Fbvh3                                                         
     &  *dconjg(Favh2)*B6                                 
     &  +Fbvh4                                                         
     &  *dconjg(Favh1)*B7                                 
     &  +Fbvh4                                                         
     &  *dconjg(Favh2)*B8                                
      expcbuad=Fbvh1                                                   
     &  *dconjg(Favh3)*B9                                  
     &  +Fbvh1                                                         
     &  *dconjg(Favh4)*B10                                 
     &  +Fbvh2                                                         
     &  *dconjg(Favh3)*B11                                 
     &  +Fbvh2                                                         
     &  *dconjg(Favh4)*B12                                
      expcbdad=Fbvh3                                                   
     &  *dconjg(Favh3)*B13                                 
     &  +Fbvh3                                                         
     &  *dconjg(Favh4)*B14                                
     &  +Fbvh4                                                         
     &  *dconjg(Favh3)*B15                                
     &  +Fbvh4                                                         
     &  *dconjg(Favh4)*B16                               
      expcaubu=Favh1                                                   
     &  *dconjg(Fbvh1)*C1                                   
     &  +Favh1                                                         
     &  *dconjg(Fbvh2)*C2                                  
     &  +Favh2                                                         
     &  *dconjg(Fbvh1)*C3                                  
     &  +Favh2                                                         
     &  *dconjg(Fbvh2)*C4                                 
      expcadbu=Favh3                                                   
     &  *dconjg(Fbvh1)*C5                                  
     &  +Favh3                                                         
     &  *dconjg(Fbvh2)*C6                                 
     &  +Favh4                                                         
     &  *dconjg(Fbvh1)*C7                                 
     &  +Favh4                                                         
     &  *dconjg(Fbvh2)*C8                                
      expcaubd=Favh1                                                   
     &  *dconjg(Fbvh3)*C9                                  
     &  +Favh1                                                         
     &  *dconjg(Fbvh4)*C10                                 
     &  +Favh2                                                         
     &  *dconjg(Fbvh3)*C11                                 
     &  +Favh2                                                         
     &  *dconjg(Fbvh4)*C12                                
      expcadbd=Favh3                                                   
     &  *dconjg(Fbvh3)*C13                                 
     &  +Favh3                                                         
     &  *dconjg(Fbvh4)*C14                                
     &  +Favh4                                                         
     &  *dconjg(Fbvh3)*C15                                
     &  +Favh4                                                         
     &  *dconjg(Fbvh4)*C16                               
      expcbubu=Fbvh1                                                   
     &  *dconjg(Fbvh1)*D1                                  
     &  +Fbvh1                                                         
     &  *dconjg(Fbvh2)*D2                                 
     &  +Fbvh2                                                         
     &  *dconjg(Fbvh1)*D3                                 
     &  +Fbvh2                                                         
     &  *dconjg(Fbvh2)*D4                                
      expcbdbu=Fbvh3                                                   
     &  *dconjg(Fbvh1)*D5                                 
     &  +Fbvh3                                                         
     &  *dconjg(Fbvh2)*D6                                
     &  +Fbvh4                                                         
     &  *dconjg(Fbvh1)*D7                                
     &  +Fbvh4                                                         
     &  *dconjg(Fbvh2)*D8                               
      expcbubd=Fbvh1                                                   
     &  *dconjg(Fbvh3)*D9                                 
     &  +Fbvh1                                                         
     &  *dconjg(Fbvh4)*D10                                
     &  +Fbvh2                                                         
     &  *dconjg(Fbvh3)*D11                                
     &  +Fbvh2                                                         
     &  *dconjg(Fbvh4)*D12                               
      expcbdbd=Fbvh3                                                   
     &  *dconjg(Fbvh3)*D13                                
     &  +Fbvh3                                                         
     &  *dconjg(Fbvh4)*D14                               
     &  +Fbvh4                                                         
     &  *dconjg(Fbvh3)*D15                               
     &  +Fbvh4                                                         
     &  *dconjg(Fbvh4)*D16                              
      cterm(4)=0.03125*ex2*(expcauau+expcadau+expcauad+expcadad        
     &                     +expcbuau+expcbdau+expcbuad+expcbdad        
     &                     +expcaubu+expcadbu+expcaubd+expcadbd        
     &                     +expcbubu+expcbdbu+expcbubd+expcbdbd)       
    
c                                                                      
c---------------------------------------------------------------c
c end of computation of complementary terms                     c
c---------------------------------------------------------------c
c
      sigma0(1)=dreal(kterm(1)+kcterm(1)+cterm(1))
      sigma0(2)=dreal(kterm(2)+kcterm(2)+cterm(2))
	sigma0(3)=dreal(kterm(3)+kcterm(3)+cterm(3))
	sigma0(4)=dreal(kterm(4)+kcterm(4)+cterm(4))
c
c----------------------------------------------------------------c
c     calculate transition functions                             c
c----------------------------------------------------------------c
      if(itran.eq.3) then
         tranv=dreal(cterm(1)/(kterm(1)+kcterm(1)+cterm(1)))
         tranh=dreal(cterm(2)/(kterm(2)+kcterm(2)+cterm(2)))
	end if
      if(itran.eq.4) then
         tranv0=dreal(cterm(1)/(kterm(1)+kcterm(1)+cterm(1)))
         tranh0=dreal(cterm(2)/(kterm(2)+kcterm(2)+cterm(2)))
	   rvtran=rvi+(rvl-rvi)*(1.0-tranv/tranv0)
	   rhtran=rhi+(rhl-rhi)*(1.0-tranh/tranh0)
	   rvhtran=(rvtran-rhtran)/2.0
	end if
c----------------------------------------------------------------c
      return
      end
c----------------------------------------------------------------c
      function expkc1(q)
      implicit real*8 (a-h,k,o-z)
	complex*16 expkc1,q,temp,sum
	complex*16 temp0
	real*8 w(5000)
	common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
	common /cont2/ ks2,iterm
	common /spectra/ w
      sum=0.0
	temp=1.0
	temp0=ks2*(css-q)*(css+cs)
      do 80 n=1,iterm
	   fn=float(n)
         temp=temp*(temp0)/fn
	   sum=sum+temp*w(n)
80    continue  
      expkc1=cdexp(-ks2*(q*q-q*css+q*cs))*sum
      return
	end
c
c------------------------------------------------------------------c
      function expkc2(q)
      implicit real*8 (a-h,k,o-z)
	complex*16 expkc2,q,temp,sum
	real*8 w(5000)
	complex*16 temp0
	common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
	common /cont2/ ks2,iterm
	common /spectra/ w
      sum=0.0
	temp=1.0
	temp0=ks2*(cs+q)*(css+cs)
       do 81 n=1,iterm
	   fn=float(n)
         temp=temp*(temp0)/fn
	   sum=sum+temp*w(n)
81     continue
      expkc2=cdexp(-ks2*(q*q-q*css+q*cs))*sum
      return
	end
c
c------------------------------------------------------------------c
      function expc1(q,qp)
      implicit real*8 (a-h,k,o-z)
	complex*16 expc1,q,qp,temp,sum
	real*8 w(5000)
	complex*16 temp0
	common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
	common /cont2/ ks2,iterm
	common /spectra/ w
      sum=0.0
	temp=1.0
	temp0=ks2*(css-q)*(css-qp)
       do 82 n=1,iterm
	   fn=float(n)
         temp=temp*(temp0)/fn
	   sum=sum+temp*w(n)
82     continue
      expc1=cdexp(-ks2*(q*q+qp*qp-css*(q+qp)+cs*(q+qp)))*sum
      return
	end
c
c------------------------------------------------------------------c
      function expc2(q,qp)
      implicit real*8 (a-h,k,o-z)
	complex*16 expc2,q,qp,temp,sum
	real*8 w(5000)
	complex*16 temp0
	common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
	common /cont2/ ks2,iterm
	common /spectra/ w
      sum=0.0
	temp=1.0
	temp0=ks2*(css-q)*(cs+qp)
       do 83 n=1,iterm
	   fn=float(n)
         temp=temp*(temp0)/fn
	   sum=sum+temp*w(n)
83     continue
      expc2=cdexp(-ks2*(q*q+qp*qp-css*(q+qp)+cs*(q+qp)))*sum
      return
	end
c
c------------------------------------------------------------------c
      function expc3(q,qp)
      implicit real*8 (a-h,k,o-z)
	complex*16 expc3,q,qp,temp,sum
	real*8 w(5000)
	complex*16 temp0
	common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
	common /cont2/ ks2,iterm
	common /spectra/ w
      sum=0.0
	temp=1.0
	temp0=ks2*(cs+q)*(css-qp)
       do 80 n=1,iterm
	   fn=float(n)
         temp=temp*(temp0)/fn
	   sum=sum+temp*w(n)
80     continue
      expc3=cdexp(-ks2*(q*q+qp*qp-css*(q+qp)+cs*(q+qp)))*sum
      return
	end
c
c------------------------------------------------------------------c
      function expc4(q,qp)
      implicit real*8 (a-h,k,o-z)
	complex*16 expc4,q,qp,temp,sum
	real*8 w(5000)
	complex*16 temp0
	common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
	common /cont2/ ks2,iterm
	common /spectra/ w
      sum=0.0
	temp=1.0
	temp0=ks2*(cs+q)*(cs+qp)
       do 80 n=1,iterm
	   fn=float(n)
         temp=temp*(temp0)/fn
	   sum=sum+temp*w(n)
80     continue
      expc4=cdexp(-ks2*(q*q+qp*qp-css*(q+qp)+cs*(q+qp)))*sum
      return
	end
c
c-------------------------------------------------------------------c
      function Favv(u,v,q,qfix)
      implicit real*8 (a-h,k,o-z)
      complex*16 er,ur,Favv,rh,rv,rvh,q,qfix
      complex*16 c1,c2,c3,c4,c5,c6
	complex*16 zx,zy,zxp,zyp
	complex*16 rpv,rmv,av,bv
	common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
      common /eu/er,ur
      common /rhhvv/rh,rv,rvh
      kxu=si+u
	ksxu=sis*csfs+u
	kyv=v
	ksyv=sis*sfs+v
	if((dabs(dreal(css-q)).lt.0.00000001).or.
     &   (dabs(dreal(cs+q)).lt.0.00000001)) then
	 c1=0.0
	 c2=0.0
	 c3=0.0
	 c4=0.0
	 c5=0.0
	 c6=0.0
	else
      zx=(-ksxu)/(css-q)
      zy=(-ksyv)/(css-q)
      zxp=(kxu)/(cs+q)
      zyp=(kyv)/(cs+q)
      c1=-csfs*(-1.0-zx*zxp)+sfs*zxp*zy
      c2=-csfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp-si*v*zx*zyp)
     &   +sfs*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp)
      c3=-csfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)
     &   +sfs*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)
      c4=-css*sfs*(-si*zyp+cs*zx*zyp)-csfs*css*(-cs-si*zxp-cs*zy*zyp)
     &   +sis*(-cs*zx-si*zx*zxp-si*zy*zyp)
      c5=-css*sfs*(-v*zx+v*zxp)-csfs*css*(q+u*zxp+v*zy)
     &   +sis*(q*zx+u*zx*zxp+v*zxp*zy)
	c6=-css*sfs*(-u*zyp+q*zx*zyp)-csfs*css*(v*zyp-q*zy*zyp)
     &   +sis*(v*zx*zyp-u*zy*zyp)
	end if
c
	rpv=1.0+rv
	rmv=1.0-rv
      av=rpv/qfix
	bv=rmv/qfix
      Favv=bv*(-rpv*c1+rmv*c2+rpv*c3)+av*(rmv*c4+rpv*c5+rmv*c6)
      return
      end
c
c-------------------------------------------------------------------c
      function Fahh(u,v,q,qfix)
      implicit real*8 (a-h,k,o-z)
      complex*16 er,ur,Fahh,rh,rv,rvh,q,qfix
      complex*16 c1,c2,c3,c4,c5,c6
	complex*16 zx,zy,zxp,zyp
	complex*16 rph,rmh,ah,bh
	common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
      common /eu/er,ur
      common /rhhvv/rh,rv,rvh
      kxu=si+u
	ksxu=sis*csfs+u
	kyv=v
	ksyv=sis*sfs+v
	if((dabs(dreal(css-q)).lt.0.00000001).or.
     &   (dabs(dreal(cs+q)).lt.0.00000001)) then
	 c1=0.0
	 c2=0.0
	 c3=0.0
	 c4=0.0
	 c5=0.0
	 c6=0.0
	else
      zx=(-ksxu)/(css-q)
      zy=(-ksyv)/(css-q)
      zxp=(kxu)/(cs+q)
      zyp=(kyv)/(cs+q)
      c1=-csfs*(-1.0-zx*zxp)+sfs*zxp*zy
      c2=-csfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp-si*v*zx*zyp)
     &   +sfs*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp)
      c3=-csfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)
     &   +sfs*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)
      c4=-css*sfs*(-si*zyp+cs*zx*zyp)-csfs*css*(-cs-si*zxp-cs*zy*zyp)
     &   +sis*(-cs*zx-si*zx*zxp-si*zy*zyp)
      c5=-css*sfs*(-v*zx+v*zxp)-csfs*css*(q+u*zxp+v*zy)
     &   +sis*(q*zx+u*zx*zxp+v*zxp*zy)
	c6=-css*sfs*(-u*zyp+q*zx*zyp)-csfs*css*(v*zyp-q*zy*zyp)
     &   +sis*(v*zx*zyp-u*zy*zyp)
	end if
c
	rph=1.0+rh
	rmh=1.0-rh
      ah=rph/qfix
	bh=rmh/qfix
      Fahh=-bh*(-rph*c1+rmh*c2+rph*c3)-ah*(rmh*c4+rph*c5+rmh*c6)
      return
      end
c
c-------------------------------------------------------------------c
      function Fbvv(u,v,q,qfix)
      implicit real*8 (a-h,k,o-z)
      complex*16 er,ur,Fbvv,rh,rv,rvh,q,qfix
      complex*16 c1,c2,c3,c4,c5,c6
	complex*16 zx,zy,zxp,zyp
	complex*16 rpv,rmv,av,bv
	common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
      common /eu/er,ur
      common /rhhvv/rh,rv,rvh
      kxu=si+u
	ksxu=sis*csfs+u
	kyv=v
	ksyv=sis*sfs+v
      zx=(-ksxu)/(css-q)
      zy=(-ksyv)/(css-q)
      zxp=(kxu)/(cs+q)
      zyp=(kyv)/(cs+q)
      c1=-csfs*(-1.0-zx*zxp)+sfs*zxp*zy
      c2=-csfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp-si*v*zx*zyp)
     &   +sfs*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp)
      c3=-csfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)
     &   +sfs*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)
      c4=-css*sfs*(-si*zyp+cs*zx*zyp)-csfs*css*(-cs-si*zxp-cs*zy*zyp)
     &   +sis*(-cs*zx-si*zx*zxp-si*zy*zyp)
      c5=-css*sfs*(-v*zx+v*zxp)-csfs*css*(q+u*zxp+v*zy)
     &   +sis*(q*zx+u*zx*zxp+v*zxp*zy)
	c6=-css*sfs*(-u*zyp+q*zx*zyp)-csfs*css*(v*zyp-q*zy*zyp)
     &   +sis*(v*zx*zyp-u*zy*zyp)
	rpv=1.0+rv
	rmv=1.0-rv
      av=rpv/qfix
	bv=rmv/qfix
      Fbvv=av*(rpv*c1-rmv*c2-rpv*c3/er)-bv*(rmv*c4*er+rpv*c5+rmv*c6)
      return
      end
c
c-------------------------------------------------------------------c
      function Fbhh(u,v,q,qfix)
      implicit real*8 (a-h,k,o-z)
      complex*16 er,ur,Fbhh,rh,rv,rvh,q,qfix
      complex*16 c1,c2,c3,c4,c5,c6
	complex*16 zx,zy,zxp,zyp
	complex*16 rph,rmh,ah,bh
	common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
      common /eu/er,ur
      common /rhhvv/rh,rv,rvh
      kxu=si+u
	ksxu=sis*csfs+u
	kyv=v
	ksyv=sis*sfs+v
      zx=(-ksxu)/(css-q)
      zy=(-ksyv)/(css-q)
      zxp=(kxu)/(cs+q)
      zyp=(kyv)/(cs+q)
      c1=-csfs*(-1.0-zx*zxp)+sfs*zxp*zy
      c2=-csfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp-si*v*zx*zyp)
     &   +sfs*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp)
      c3=-csfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)
     &   +sfs*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)
      c4=-css*sfs*(-si*zyp+cs*zx*zyp)-csfs*css*(-cs-si*zxp-cs*zy*zyp)
     &   +sis*(-cs*zx-si*zx*zxp-si*zy*zyp)
      c5=-css*sfs*(-v*zx+v*zxp)-csfs*css*(q+u*zxp+v*zy)
     &   +sis*(q*zx+u*zx*zxp+v*zxp*zy)
	c6=-css*sfs*(-u*zyp+q*zx*zyp)-csfs*css*(v*zyp-q*zy*zyp)
     &   +sis*(v*zx*zyp-u*zy*zyp)
	rph=1.0+rh
	rmh=1.0-rh
      ah=rph/qfix
	bh=rmh/qfix
      Fbhh=ah*(-rph*c1*er+rmh*c2+rph*c3)+bh*(rmh*c4+rph*c5+rmh*c6/er)
      return
      end
c
c-------------------------------------------------------------------c
      function Fahv(u,v,q,qfix)
      implicit real*8 (a-h,k,o-z)
      complex*16 er,ur,Fahv,rh,rv,rvh,q,qfix
      complex*16 b1,b2,b3,b4,b5,b6
	complex*16 zx,zy,zxp,zyp
	complex*16 rp,rm,a,b
	common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
      common /eu/er,ur
      common /rhhvv/rh,rv,rvh
      kxu=si+u
	ksxu=sis*csfs+u
	kyv=v
	ksyv=sis*sfs+v
	if((dabs(dreal(css-q)).lt.0.00000001).or.
     &   (dabs(dreal(cs+q)).lt.0.00000001)) then
	 b1=0.0
	 b2=0.0
	 b3=0.0
	 b4=0.0
	 b5=0.0
	 b6=0.0
	else
      zx=(-ksxu)/(css-q)
      zy=(-ksyv)/(css-q)
      zxp=(kxu)/(cs+q)
      zyp=(kyv)/(cs+q)
      b1=-css*sfs*(-1.0-zx*zxp)-sis*zy-csfs*css*zxp*zy
	b2=-css*sfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp 
     &    -si*v*zx*zyp)+ 
     & sis*(-cs*q*zy-q*si*zxp*zy+q*si*zx*zyp-cs*u*zx*zyp-cs*v*zy*zyp)- 
     & csfs*css*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp)
      b3=-css*sfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)- 
     &   csfs*css*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)+ 
     &   sis*(-si*v*zx+cs*v*zx*zxp+si*u*zy-cs*u*zxp*zy)
	b4=-csfs*(-si*zyp+cs*zx*zyp)+sfs*(-cs-si*zxp-cs*zy*zyp)
      b5=-csfs*(-v*zx+v*zxp)+sfs*(q+u*zxp+v*zy)
	b6=-csfs*(-u*zyp+q*zx*zyp)+sfs*(v*zyp-q*zy*zyp)
	end if
c
	rp=1.0+rvh
	rm=1.0-rvh
      a=rp/qfix
	b=rm/qfix
      Fahv=b*(rp*b1-rm*b2-rp*b3)+a*(rm*b4+rp*b5+rm*b6)
      return
      end
c
c-------------------------------------------------------------------c
      function Favh(u,v,q,qfix)
      implicit real*8 (a-h,k,o-z)
      complex*16 er,ur,Favh,rh,rv,rvh,q,qfix
      complex*16 b1,b2,b3,b4,b5,b6
	complex*16 zx,zy,zxp,zyp
	complex*16 rp,rm,a,b
	common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
      common /eu/er,ur
      common /rhhvv/rh,rv,rvh
      kxu=si+u
	ksxu=sis*csfs+u
	kyv=v
	ksyv=sis*sfs+v
	if((dabs(dreal(css-q)).lt.0.00000001).or.
     &   (dabs(dreal(cs+q)).lt.0.00000001)) then
	 b1=0.0
	 b2=0.0
	 b3=0.0
	 b4=0.0
	 b5=0.0
	 b6=0.0
	else
      zx=(-ksxu)/(css-q)
      zy=(-ksyv)/(css-q)
      zxp=(kxu)/(cs+q)
      zyp=(kyv)/(cs+q)
      b1=-css*sfs*(-1.0-zx*zxp)-sis*zy-csfs*css*zxp*zy
	b2=-css*sfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp 
     &    -si*v*zx*zyp)+ 
     & sis*(-cs*q*zy-q*si*zxp*zy+q*si*zx*zyp-cs*u*zx*zyp-cs*v*zy*zyp)- 
     & csfs*css*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp)
      b3=-css*sfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)- 
     &   csfs*css*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)+ 
     &   sis*(-si*v*zx+cs*v*zx*zxp+si*u*zy-cs*u*zxp*zy)
	b4=-csfs*(-si*zyp+cs*zx*zyp)+sfs*(-cs-si*zxp-cs*zy*zyp)
      b5=-csfs*(-v*zx+v*zxp)+sfs*(q+u*zxp+v*zy)
	b6=-csfs*(-u*zyp+q*zx*zyp)+sfs*(v*zyp-q*zy*zyp)
      end if
c
	rp=1.0+rvh
	rm=1.0-rvh
      a=rp/qfix
	b=rm/qfix
      Favh=b*(rp*b4+rm*b5+rp*b6)-a*(-rm*b1+rp*b2+rm*b3)
      return
      end
c
c-------------------------------------------------------------------c
      function Fbhv(u,v,q,qfix)
      implicit real*8 (a-h,k,o-z)
      complex*16 er,ur,Fbhv,rh,rv,rvh,q,qfix
      complex*16 b1,b2,b3,b4,b5,b6
	complex*16 zx,zy,zxp,zyp
	complex*16 rp,rm,a,b
	common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
      common /eu/er,ur
      common /rhhvv/rh,rv,rvh
      kxu=si+u
	ksxu=sis*csfs+u
	kyv=v
	ksyv=sis*sfs+v
      zx=(-ksxu)/(css-q)
      zy=(-ksyv)/(css-q)
      zxp=(kxu)/(cs+q)
      zyp=(kyv)/(cs+q)
      b1=-css*sfs*(-1.0-zx*zxp)-sis*zy-csfs*css*zxp*zy
	b2=-css*sfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp 
     &    -si*v*zx*zyp)+ 
     & sis*(-cs*q*zy-q*si*zxp*zy+q*si*zx*zyp-cs*u*zx*zyp-cs*v*zy*zyp)- 
     & csfs*css*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp)
      b3=-css*sfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)- 
     &   csfs*css*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)+ 
     &   sis*(-si*v*zx+cs*v*zx*zxp+si*u*zy-cs*u*zxp*zy)
	b4=-csfs*(-si*zyp+cs*zx*zyp)+sfs*(-cs-si*zxp-cs*zy*zyp)
      b5=-csfs*(-v*zx+v*zxp)+sfs*(q+u*zxp+v*zy)
	b6=-csfs*(-u*zyp+q*zx*zyp)+sfs*(v*zyp-q*zy*zyp)
	rp=1.0+rvh
	rm=1.0-rvh
      a=rp/qfix
	b=rm/qfix
      Fbhv=a*(-rp*b1+rm*b2+rp*b3/er)-b*(rm*b4*er+rp*b5+rm*b6)
      return
      end
c
c-------------------------------------------------------------------c
      function Fbvh(u,v,q,qfix)
      implicit real*8 (a-h,k,o-z)
      complex*16 er,ur,Fbvh,rh,rv,rvh,q,qfix
      complex*16 b1,b2,b3,b4,b5,b6
	complex*16 zx,zy,zxp,zyp
	complex*16 rp,rm,a,b
	common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
      common /eu/er,ur
      common /rhhvv/rh,rv,rvh
      kxu=si+u
	ksxu=sis*csfs+u
	kyv=v
	ksyv=sis*sfs+v
      zx=(-ksxu)/(css-q)
      zy=(-ksyv)/(css-q)
      zxp=(kxu)/(cs+q)
      zyp=(kyv)/(cs+q)
      b1=-css*sfs*(-1.0-zx*zxp)-sis*zy-csfs*css*zxp*zy
	b2=-css*sfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp 
     &    -si*v*zx*zyp)+ 
     & sis*(-cs*q*zy-q*si*zxp*zy+q*si*zx*zyp-cs*u*zx*zyp-cs*v*zy*zyp)- 
     & csfs*css*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp)
      b3=-css*sfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)- 
     &   csfs*css*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)+ 
     &   sis*(-si*v*zx+cs*v*zx*zxp+si*u*zy-cs*u*zxp*zy)
	b4=-csfs*(-si*zyp+cs*zx*zyp)+sfs*(-cs-si*zxp-cs*zy*zyp)
      b5=-csfs*(-v*zx+v*zxp)+sfs*(q+u*zxp+v*zy)
	b6=-csfs*(-u*zyp+q*zx*zyp)+sfs*(v*zyp-q*zy*zyp)
	rp=1.0+rvh
	rm=1.0-rvh
      a=rp/qfix
	b=rm/qfix
      Fbvh=-a*(rp*b4+rm*b5+rp*b6/er)+b*(-rm*b1*er+rp*b2+rm*b3)
      return
      end
c
c
c-------------------------------------------------------------------c
c
      function bessj0(x)
      real*8 x,y,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     &s1,s2,s3,s4,s5,s6,bessj0
      data r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,
     &651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,
     &s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0,
     &59272.64853d0,267.8532712d0,1.d0/
      data p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     &-.2073370639d-5,.2093887211d-6/,q1,q2,q3,q4,q5/-.1562499995d-1,
     &.1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      if(abs(x).lt.8.)then
	y=x**2
      bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     &  /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
	else
	ax=abs(x)
	z=8./ax
	y=z**2
	xx=ax-.785398164
      bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y
     &  *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
	endif
      return
      end

      subroutine shadowg(ti,ts,s,shfct)
      implicit real*8 (a-h,o-z)
c       if(ts.ge.ti) then
c         arg=ts
c         else
c           arg=ti
c        endif
	arg=ti
	if(arg.eq.0.0) then
	  shfct=1.0
	   return
	endif
       u=1.0/dtan(arg)
       pi=acos(-1.0d0)
       et=u/(sqrt(2.0)*s)
       if(et.ge.20) then
	   shfct=1.0
	    return
	 endif
       f1=dsqrt(2.0/pi)*s*dexp(-et*et)/u
	f2=erfc(et)
       f=(f1-f2)/2.0
       shfct=1.0/(1.0+f)
c       shfct=(1.0-0.5*erfc(et))/(1.0+f)
       return
       end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
c subroutine calculates shadowing function     c
c (Tsang et al 1985, pp.95                     c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
      subroutine shadow(back,ti,ts,shfct)
      implicit real*8 (a-h,o-z)
      logical back
      common /al/cl,effslop
      s=effslop
	if(ti.eq.0.0.and.ts.eq.0.0) then
	  shfct=1.0
	   return
	endif
       if(ti.eq.0.0) then
	    ui=1.0d+30
	  else
	    ui=1.0/dtan(ti)
       endif
       if(ts.eq.0.0) then
	    us=1.0d+30
	   else
	    us=1.0/dtan(ts)
       endif
       pi=acos(-1.0d0)
       eti=ui/(sqrt(2.0)*s)
       ets=us/(sqrt(2.0)*s)
	      
       if(.not.back) then
	f1i=sqrt(2.0/pi)*s/ui*exp(-eti*eti)
       f1s=sqrt(2.0/pi)*s/us*exp(-ets*ets)
       f2i=erfc(eti)
       f2s=erfc(ets)
       fi=(f1i-f2i)/2.0
       fs=(f1s-f2s)/2.0
	 shfct=1.0/(1.0+fi+fs)
       else
       if(ts.ge.ti) then
	 f1s=sqrt(2.0/pi)*s/us*exp(-ets*ets)
	 f2s=erfc(ets)
	 fs=(f1s-f2s)/2.0
	  shfct=(1.0-0.5*f2s)/(1.0+fs) 
	 else
	   f1i=sqrt(2.0/pi)*s/ui*exp(-eti*eti)
	   f2i=erfc(eti)
	   fi=(f1i-f2i)/2.0
	  shfct=(1.0-0.5*f2i)/(fi+1.0) 
	endif
	endif
       return
       end

      FUNCTION erfcc(x)
      REAL*8 erfcc,x
      REAL*8 t,z
      z=dabs(x)
      t=1.d0/(1.d0+0.5d0*z)
      erfcc=t*exp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+t*
     *(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*(-1.13520398d0+t*
     *(1.48851587d0+t*(-.82215223d0+t*.17087277d0)))))))))
      if (x.lt.0.d0) erfcc=2.d0-erfcc
      return
      END

      function erfc(x)
      implicit real*8 (a-h,o-z)
      if(x.lt.0.d0)then
        erfc=1.d0+gammp(.5d0,x**2.d0)
      else
        erfc=gammq(.5d0,x**2.d0)
      endif
      return
      end

      function gammp(a,x)
      implicit real*8 (a-h,o-z)
      if(x.lt.0.d0.or.a.le.0.d0) pause
      if(x.lt.a+1.d0)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.d0-gammcf
      endif
      return
      end

      function gammq(a,x)
      implicit real*8 (a-h,o-z)
      if(x.lt.0.d0.or.a.le.0.d0)pause
      if(x.lt.a+1.d0)then
        call gser(gamser,a,x,gln)
        gammq=1.-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      end

      subroutine gser(gamser,a,x,gln)
      implicit real*8 (a-h,o-z)
      parameter (itmax=100,eps=1.d-18)
      gln=gammln(a)
      if(x.le.0.d0)then
        if(x.lt.0.d0)pause
        gamser=0.d0
        return
      endif
      ap=a
      sum=1.d0/a
      del=sum
      do 11 n=1,itmax
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(dabs(del).lt.dabs(sum)*eps)go to 1
11    continue
      pause 'a too large, itmax too small'
1     gamser=sum*dexp(-x+a*dlog(x)-gln)
      return
      end

      subroutine gcf(gammcf,a,x,gln)
      implicit real*8 (a-h,o-z)
      parameter (itmax=100,eps=1.d-18)
      gln=gammln(a)
      gold=0.0d0
      a0=1.0d0
      a1=x
      b0=0.0d0
      b1=1.0d0
      fac=1.0d0
      do 11 n=1,itmax
        an=dfloat(n)
        ana=an-a
        a0=(a1+a0*ana)*fac
        b0=(b1+b0*ana)*fac
        anf=an*fac
        a1=x*a0+anf*a1
        b1=x*b0+anf*b1
        if(a1.ne.0.d0)then
          fac=1.0d0/a1
          g=b1*fac
          if(dabs((g-gold)/g).lt.eps)go to 1
          gold=g
        endif
11    continue
      pause 'a too large, itmax too small'
1     gammcf=dexp(-x+a*dlog(x)-gln)*g
      return
      end

      function gammln(xx)
      implicit real*8 (a-h,o-z)
      dimension cof(6)
      data cof,stp/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      data half,one,fpf/0.5D0,1.0D0,5.5D0/
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*dlog(tmp)-tmp
      ser=one
      do 11 j=1,6
        x=x+one
        ser=ser+cof(j)/x
11    continue
      gammln=tmp+dlog(stp*ser)
      return
      end
c
c
c
      function BesselK(n,x)
c***************************************************
c  Modified Bessel function of order n+0.5
c
c  Input parameters:
c    n : integer part of order
c    x : real parameter
c
c***************************************************

      implicit real*8 (a - h, k, o - z)
      if(x.eq.0.0d0) then
	print*,'BesselK: Singularity encountered !'
	return
      endif
      PI = 4.0d0*atan(1.0d0)
      cons = dsqrt(PI/(2.0d0*x))
      K0 = cons*exp(-x)
      if(n.eq.0) then
	BesselK = K0
	return
      endif
      K1 = K0*(1.0d0+1.0d0/x)
      if(n.eq.1) then
	BesselK = K1
	return
      endif
      K2 = K0*(1.0d0+3.0d0/x+3.0d0/x/x)
      if(n.eq.2) then
	BesselK = K2
	return
      endif
      fn0 = K1*cons
      fn1 = -K2*cons
      do 10 i=2,n-1
	fn2 = fn0-float(2*i+1)/x*fn1
	fn0 = fn1
	fn1 = fn2
10    continue
      BesselK = dabs(fn2/cons)
      return
      end
c
c
c
      function alogam( x )
c
c evaluates natural logarithm of GAMMA(x)
c for x > 0
c
      implicit real*8 ( a - h, o - z)
      pi = acos (-1.0d0)
      a1 = dlog( 2 * pi) / 2.0d0
      a2 = 1.0 / 1680.0d0
      a3 = 1.0 / 1260.0d0
      a4 = 1.0 / 360.0d0
      a5 = 1.0 / 12.0d0
      alogam = 0.0d0
      ifault = 1
      if (x .le. 0.0d0) return
      ifault = 0
      y = x
      f = 0.0d0
      if (y .ge. 7.0d0) go to 30
      f = y  
10    y = y + 1.0d0
      if (y .ge. 7.0d0) go to 20
      f = f * y
      go to 10
20    f = -dlog(f)
30    z = 1.0d0 / ( y * y)
      alogam = f + (y - 0.5d0) * dlog(y) - y + a1
     &         + ((( -a2 * z + a3) * z - a4) * z
     &         + a5) / y
      return
      end

      FUNCTION BESSK(N,X)
      implicit real*8 (a-h,o-z)
      IF (N.ge.2) go to 1
      if(n.eq.0) bessk = bessk0(x)
      if(n.eq.1) bessk = bessk1(x)
      return
1     continue
      TOX=2.0/X
      BKM=BESSK0(X)
      BK=BESSK1(X)
      DO 11 J=1,N-1
	BKP=BKM+J*TOX*BK
	BKM=BK
	BK=BKP
11    CONTINUE
      BESSK=BK
      RETURN
      END

      FUNCTION BESSK0(X)
      implicit real*8 (a-h,o-z)
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,0.23069756D0,
     *    0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,0.2189568D-1,
     *    -0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
      IF (X.LE.2.0d0) THEN
	Y=X*X/4.0d0
	BESSK0=(-dLOG(X/2.0d0)*BESSI0(X))+(P1+Y*(P2+Y*(P3+
     *        Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
	Y=(2.0/X)
	BESSK0=(dEXP(-X)/dSQRT(X))*(Q1+Y*(Q2+Y*(Q3+
     *        Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END

      FUNCTION BESSK1(X)
      implicit real*8 (a-h,o-z)
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,0.15443144D0,-0.67278579D0,
     *    -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,-0.3655620D-1,
     *    0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/
      IF (X.LE.2.0d0) THEN
	Y=X*X/4.0d0
	BESSK1=(dLOG(X/2.0)*BESSI1(X))+(1.0/X)*(P1+Y*(P2+
     *      Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
	Y=2.0/X
	BESSK1=(dEXP(-X)/dSQRT(X))*(Q1+Y*(Q2+Y*(Q3+
     *      Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END

      FUNCTION BESSI0(X)
      implicit real*8 (a-h,o-z)
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,1.2067492D
     *0,
     *    0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     *    0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     *    0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF (dABS(X).LT.3.75) THEN
	Y=(X/3.75)**2
	BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
	AX=ABS(X)
	Y=3.75/AX
	BESSI0=(dEXP(AX)/dSQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4
     *      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END

      FUNCTION BESSI1(X)
      implicit real*8 (a-h,o-z)
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,
     *    0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,
     *    -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,
     *    -0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF (dABS(X).LT.3.75) THEN
	Y=(X/3.75)**2
	BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
	AX=dABS(X)
	Y=3.75/AX
	BESSI1=(dEXP(AX)/dSQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+
     *      Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END




c************************************************************************************c

C        SUBROUTINE CJYVA(V,Z,VM,CBJ,CDJ,CBY,CDY)
        SUBROUTINE CIKVA(V,Z,VM,Res)
C
C       ============================================================
C       Purpose: Compute the modified Bessel functions Iv(z), Kv(z)
C                and their derivatives for an arbitrary order and
C                complex argument
C       Input :  z --- Complex argument
C                v --- Real order of Iv(z) and Kv(z)
C                      ( v = n+v0, n = 0,1,2,..., 0 v0 < 1 )
C       Output:  CBI(n) --- In+v0(z)
C                CDI(n) --- In+v0'(z)
C                CBK(n) --- Kn+v0(z)
C                CDK(n) --- Kn+v0'(z)
C                VM --- Highest order computed
C       Routines called:
C            (1) GAMMA for computing the gamma function
C            (2) MSTA1 and MSTA2 for computing the starting 
C                point for backward recurrence
C       ============================================================
C
        IMPLICIT DOUBLE PRECISION (A,G,P,R,V,W)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBI(0:1000),CDI(0:1000),CBK(0:1000),CDK(0:1000)
        PI=3.141592653589793D0
        CI=(0.0D0,1.0D0)
        A0=CDABS(Z)
        Z1=Z
        Z2=Z*Z
        N=INT(V)
        V0=V-N
        PIV=PI*V0
        VT=4.0D0*V0*V0
        IF (N.EQ.0) N=1
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBI(K)=0.0D0
              CDI(K)=0.0D0
              CBK(K)=-1.0D+300
10            CDK(K)=1.0D+300
           IF (V0.EQ.0.0) THEN
              CBI(0)=(1.0D0,0.0D0)
              CDI(1)=(0.5D0,0.0D0)
           ENDIF
           VM=V
           RETURN
        ENDIF
        K0=14
        IF (A0.GE.35.0) K0=10
        IF (A0.GE.50.0) K0=8
        IF (REAL(Z).LT.0.0) Z1=-Z
        IF (A0.LT.18.0) THEN
           IF (V0.EQ.0.0) THEN
              CA1=(1.0D0,0.0D0)
           ELSE
              V0P=1.0D0+V0
              CALL GAMMA(V0P,GAP)
              CA1=(0.5D0*Z1)**V0/GAP
           ENDIF
           CI0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 15 K=1,50
              CR=0.25D0*CR*Z2/(K*(K+V0))
              CI0=CI0+CR
              IF (CDABS(CR).LT.CDABS(CI0)*1.0D-15) GO TO 20
15         CONTINUE
20         CBI0=CI0*CA1
        ELSE
           CA=CDEXP(Z1)/CDSQRT(2.0D0*PI*Z1)
           CS=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 25 K=1,K0
              CR=-0.125D0*CR*(VT-(2.0D0*K-1.0D0)**2.0)/(K*Z1)
25            CS=CS+CR
           CBI0=CA*CS
        ENDIF
        M=MSTA1(A0,200)
        IF (M.LT.N) THEN
           N=M
        ELSE
           M=MSTA2(A0,N,15)
        ENDIF
        CF2=(0.0D0,0.0D0)
        CF1=(1.0D-100,0.0D0)
        DO 30 K=M,0,-1
           CF=2.0D0*(V0+K+1.0D0)/Z1*CF1+CF2
           IF (K.LE.N) CBI(K)=CF
           CF2=CF1
30         CF1=CF
        CS=CBI0/CF
        DO 35 K=0,N
35         CBI(K)=CS*CBI(K)
        IF (A0.LE.9.0) THEN
           IF (V0.EQ.0.0) THEN
              CT=-CDLOG(0.5D0*Z1)-0.5772156649015329D0
              CS=(0.0D0,0.0D0)
              W0=0.0D0
              CR=(1.0D0,0.0D0)
              DO 40 K=1,50
                 W0=W0+1.0D0/K
                 CR=0.25D0*CR/(K*K)*Z2
                 CP=CR*(W0+CT)
                 CS=CS+CP
                 IF (K.GE.10.AND.CDABS(CP/CS).LT.1.0D-15) GO TO 45
40            CONTINUE
45            CBK0=CT+CS
           ELSE
              V0N=1.0D0-V0
              CALL GAMMA(V0N,GAN)
              CA2=1.0D0/(GAN*(0.5D0*Z1)**V0)
              CA1=(0.5D0*Z1)**V0/GAP
              CSU=CA2-CA1
              CR1=(1.0D0,0.0D0)
              CR2=(1.0D0,0.0D0)
              DO 50 K=1,50
                 CR1=0.25D0*CR1*Z2/(K*(K-V0))
                 CR2=0.25D0*CR2*Z2/(K*(K+V0))
                 CSU=CSU+CA2*CR1-CA1*CR2
                 WS=CDABS(CSU)
                 IF (K.GE.10.AND.DABS(WS-WS0)/WS.LT.1.0D-15) GO TO 55
                 WS0=WS
50            CONTINUE
55            CBK0=0.5D0*PI*CSU/DSIN(PIV)
           ENDIF
        ELSE
           CB=CDEXP(-Z1)*CDSQRT(0.5D0*PI/Z1)
           CS=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 60 K=1,K0
              CR=0.125D0*CR*(VT-(2.0D0*K-1.0D0)**2.0)/(K*Z1)
60            CS=CS+CR
           CBK0=CB*CS
        ENDIF
        CBK1=(1.0D0/Z1-CBI(1)*CBK0)/CBI(0)
        CBK(0)=CBK0
        CBK(1)=CBK1
        CG0=CBK0
        CG1=CBK1
        DO 65 K=2,N
           CGK=2.0D0*(V0+K-1.0D0)/Z1*CG1+CG0
           CBK(K)=CGK
           CG0=CG1
65         CG1=CGK
        IF (REAL(Z).LT.0.0) THEN
           DO 70 K=0,N
              CVK=CDEXP((K+V0)*PI*CI)
              IF (DIMAG(Z).LT.0.0D0) THEN
                 CBK(K)=CVK*CBK(K)+PI*CI*CBI(K)
                 CBI(K)=CBI(K)/CVK
              ELSE IF (DIMAG(Z).GT.0.0) THEN
                 CBK(K)=CBK(K)/CVK-PI*CI*CBI(K)
                 CBI(K)=CVK*CBI(K)
              ENDIF
70         CONTINUE
        ENDIF
        CDI(0)=V0/Z*CBI(0)+CBI(1)
        CDK(0)=V0/Z*CBK(0)-CBK(1)
        DO 75 K=1,N
           CDI(K)=-(K+V0)/Z*CBI(K)+CBI(K-1)
75         CDK(K)=-(K+V0)/Z*CBK(K)-CBK(K-1)
        VM=N+V0
	  Res=CBK(INT(V))
        RETURN
        END


        SUBROUTINE GAMMA(X,GA)
C
C       ==================================================
C       Purpose: Compute gamma function a(x)
C       Input :  x  --- Argument of a(x)
C                       ( x is not equal to 0,-1,-2, ...)
C       Output:  GA --- a(x)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END


        INTEGER FUNCTION MSTA1(X,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward  
C                recurrence such that the magnitude of    
C                Jn(x) at that point is about 10^(-MP)
C       Input :  x     --- Argument of Jn(x)
C                MP    --- Value of magnitude
C       Output:  MSTA1 --- Starting point   
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        N0=INT(1.1*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20             
           NN=N1-(N1-N0)/(1.0D0-F0/F1)                  
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
 10        F1=F
 20     MSTA1=NN
        RETURN
        END


        INTEGER FUNCTION MSTA2(X,N,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that all Jn(x) has MP
C                significant digits
C       Input :  x  --- Argument of Jn(x)
C                n  --- Order of Jn(x)
C                MP --- Significant digit
C       Output:  MSTA2 --- Starting point
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1*A0)+1
c	modified by zhaokaiguang
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
10         F1=F
20      MSTA2=NN+10
        RETURN
        END

        REAL*8 FUNCTION ENVJ(N,X)
        DOUBLE PRECISION X
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END





C.....................................................................C
C    To calcualte Bp using bisection methods from the nonlinar        C
C       equation(15) listed by Qi Li, et al.,2000					    C
C.....................................................................C
	FUNCTION compute_Bp(p)
      implicit real*8 (a-h,o-z)
	real*8 p, ya,yb
	dimension rt(20)    

      a=.11d0
	b=1.0d0
	eps=1.0d-14
	h=0.03D0
	nr=0
	x=a

	call Bp_Diff(x,yb,p)
10	x=x+h
	if(x.gt.b+h) goto 30
	ya=yb
	call Bp_Diff(x,yb,p)
	if(dabs(yb).lt.1.0e-14) goto 20
	if(ya*yb.gt.0.0) goto 10
	xa=x-h
	xb=x
15	x=0.5*(xa+xb)
	call Bp_Diff(x,y,p)
	if(dabs(y).lt.1.0e-14) goto 20
	if(ya*y.gt.0.0)xa=x
	if(ya*y.le.0.0)xb=x
	if(dabs((xa-xb)/x).gt.eps) goto 15
20	nr=nr+1
	rt(nr)=x
	goto 10
30	continue
	compute_Bp=rt(1)
	END 
C.....................................................................C
	subroutine Bp_Diff(x,y,p)
	implicit real*8 (a-h,k,t,o-s,u-z)
	complex*16 z

	call gamma(p-0.5d0, temp1)
	call gamma(p	, temp2)
	ap	=	temp1/temp2
	z	=	2.0d0*dcmplx(x,0)/ap
	if( p .lt. 90) then
      call	CIKVA(p-1.0d0,z,VM,Res)
	else
	call	CIKLV(p-1.0d0,z,Res)
	endif
      call	gamma(p-1.0d0, temp1)	
	y	=	(2*x/ap)**(p-1.0d0)*Res
     &		- temp1*dexp(-1.d0)*(2.0)**(p-2.0d0)
	End
C.....................................................................C
C  End calcualting Bp using bisection methods from the nonlinar       C
C       equation(15) listed by Qi Li, et al.,2000					    C
C.....................................................................C


C.....................................................................C
C  Integrand used to compute nth spectrum of X-expoential Correlation function
C.....................................................................C
	REAL*8 FUNCTION Corr (X)
	REAL*8 X
	REAL*8 pIndex, xIndex, ap,bp
	common /p/pIndex, xIndex,ap,bp
	real*8 fn,L,K,ru
	common /cor/ fn, L,K
	Corr=dexp( dlog(1.0d305)-(dabs(X)/L)**xIndex * fn) * bessj0(K*X)*X
	RETURN
	END

C.....................................................................C
C  Integrand used to compute nth spectrum of expoential-like Correlation function
C.....................................................................C

	REAL*8 FUNCTION CorrExp(X)
	REAL*8 X
	REAL*8 pIndex, xIndex, ap,bp
	common /p/pIndex, xIndex,ap,bp
	real*8 fn,L,K,ru
	common /cor/ fn, L,K


      
	ru= - ( dabs(X) / L )* (1-dexp(-X/xIndex)  )   * fn

	CorrExp=dexp(dlog(1d305) + ru)*bessj0(K*X)*X

	RETURN
	END

	




	REAL*8 FUNCTION W_expLike(fn,xIndex,L,K)
      implicit real*8 (a-h,k,o-z)
	REAL*8 fn, xIndex,K,L,Le
	Real*8 HF
	INTEGER m
	
	tmp3=0d0

	do 11, m=0, 160
	Le=xIndex*L/(fn*xIndex+m*L)
	call gamma(m+2d0, gam)
	A=(2d0+m)/2d0
	B=(3d0+m)/3d0
	C=1.0D0
	E=-K*Le
	call HYGFX(A,B,C,E,HF)
	tmp1=gam*HF*Le*Le

	
	do 12,i=0, m
	fi=float(i)
	if( m .eq. 0 )then
	tmp2=1.0d0
	else 
		if (i .eq. 0) then
		tmp2=1.0d0
		else
		tmp2=tmp2*(fn*Le/L)/fi
		endif
	endif
12	continue

	tmp3=tmp3+tmp1*tmp2

11	continue
	
	W_expLike=tmp3
	RETURN
	END










C=============Guassian Hypergeometric function =========================C

        SUBROUTINE HYGFX(A,B,C,X,HF)
C
C       ====================================================
C       Purpose: Compute hypergeometric function F(a,b,c,x)
C       Input :  a --- Parameter
C                b --- Parameter
C                c --- Parameter, c <> 0,-1,-2,...
C                x --- Argument   ( x < 1 )
C       Output:  HF --- F(a,b,c,x)
C       Routines called:
C            (1) GAMMA for computing gamma function
C            (2) PSI for computing psi function
C       ====================================================
C
        IMPLICIT Real*8 (A-H,O-Z)
        LOGICAL L0,L1,L2,L3,L4,L5
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        L0=C.EQ.INT(C).AND.C.LT.0.0
        L1=1.0D0-X.LT.1.0D-15.AND.C-A-B.LE.0.0
        L2=A.EQ.INT(A).AND.A.LT.0.0
        L3=B.EQ.INT(B).AND.B.LT.0.0
        L4=C-A.EQ.INT(C-A).AND.C-A.LE.0.0
        L5=C-B.EQ.INT(C-B).AND.C-B.LE.0.0
        IF (L0.OR.L1) THEN
           WRITE(*,*)'The hypergeometric series is divergent'
           RETURN
        ENDIF
        EPS=1.0D-15
        IF (X.GT.0.95) EPS=1.0D-8
        IF (X.EQ.0.0.OR.A.EQ.0.0.OR.B.EQ.0.0) THEN
           HF=1.0D0
           RETURN
        ELSE IF (1.0D0-X.EQ.EPS.AND.C-A-B.GT.0.0) THEN
           CALL GAMMA(C,GC)
           CALL GAMMA(C-A-B,GCAB)
           CALL GAMMA(C-A,GCA)
           CALL GAMMA(C-B,GCB)
           HF=GC*GCAB/(GCA*GCB)
           RETURN
        ELSE IF (1.0D0+X.LE.EPS.AND.DABS(C-A+B-1.0).LE.EPS) THEN
           G0=DSQRT(PI)*2.0D0**(-A)
           CALL GAMMA(C,G1)
           CALL GAMMA(1.0D0+A/2.0-B,G2)
           CALL GAMMA(0.5D0+0.5*A,G3)
           HF=G0*G1/(G2*G3)
           RETURN
        ELSE IF (L2.OR.L3) THEN
           IF (L2) NM=INT(ABS(A))
           IF (L3) NM=INT(ABS(B))
           HF=1.0D0
           R=1.0D0
           DO 10 K=1,NM
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
10            HF=HF+R
           RETURN
        ELSE IF (L4.OR.L5) THEN
           IF (L4) NM=INT(ABS(C-A))
           IF (L5) NM=INT(ABS(C-B))
           HF=1.0D0
           R=1.0D0
           DO 15 K=1,NM
              R=R*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*X
15            HF=HF+R
           HF=(1.0D0-X)**(C-A-B)*HF
           RETURN
        ENDIF
        AA=A
        BB=B
        X1=X
        IF (X.LT.0.0D0) THEN
           X=X/(X-1.0D0)
           IF (C.GT.A.AND.B.LT.A.AND.B.GT.0.0) THEN
              A=BB
              B=AA
           ENDIF
           B=C-B
        ENDIF
        IF (X.GE.0.75D0) THEN
           GM=0.0D0
           IF (DABS(C-A-B-INT(C-A-B)).LT.1.0D-15) THEN
              M=INT(C-A-B)
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(A+M,GAM)
              CALL GAMMA(B+M,GBM)
              CALL PSI(A,PA)
              CALL PSI(B,PB)
              IF (M.NE.0) GM=1.0D0
              DO 30 J=1,ABS(M)-1
30               GM=GM*J
              RM=1.0D0
              DO 35 J=1,ABS(M)
35               RM=RM*J
              F0=1.0D0
              R0=1.0D0
              R1=1.0D0
              SP0=0.D0
              SP=0.0D0
              IF (M.GE.0) THEN
                 C0=GM*GC/(GAM*GBM)
                 C1=-GC*(X-1.0D0)**M/(GA*GB*RM)
                 DO 40 K=1,M-1
                    R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(K-M))*(1.0-X)
40                  F0=F0+R0
                 DO 45 K=1,M
45                  SP0=SP0+1.0D0/(A+K-1.0)+1.0/(B+K-1.0)-1.0/K
                 F1=PA+PB+SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 55 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 50 J=1,M
50                     SM=SM+(1.0D0-A)/((J+K)*(A+J+K-1.0))+1.0/
     &                    (B+J+K-1.0)
                    RP=PA+PB+2.0D0*EL+SP+SM+DLOG(1.0D0-X)
                    R1=R1*(A+M+K-1.0D0)*(B+M+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 60
55                  HW=F1
60               HF=F0*C0+F1*C1
              ELSE IF (M.LT.0) THEN
                 M=-M
                 C0=GM*GC/(GA*GB*(1.0D0-X)**M)
                 C1=-(-1)**M*GC/(GAM*GBM*RM)
                 DO 65 K=1,M-1
                    R0=R0*(A-M+K-1.0D0)*(B-M+K-1.0)/(K*(K-M))*(1.0-X)
65                  F0=F0+R0
                 DO 70 K=1,M
70                  SP0=SP0+1.0D0/K
                 F1=PA+PB-SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 80 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 75 J=1,M
75                     SM=SM+1.0D0/(J+K)
                    RP=PA+PB+2.0D0*EL+SP-SM+DLOG(1.0D0-X)
                    R1=R1*(A+K-1.0D0)*(B+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 85
80                  HW=F1
85               HF=F0*C0+F1*C1
              ENDIF
           ELSE
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(C-A,GCA)
              CALL GAMMA(C-B,GCB)
              CALL GAMMA(C-A-B,GCAB)
              CALL GAMMA(A+B-C,GABC)
              C0=GC*GCAB/(GCA*GCB)
              C1=GC*GABC/(GA*GB)*(1.0D0-X)**(C-A-B)
              HF=0.0D0
              R0=C0
              R1=C1
              DO 90 K=1,250
                 R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(A+B-C+K))*(1.0-X)
                 R1=R1*(C-A+K-1.0D0)*(C-B+K-1.0)/(K*(C-A-B+K))
     &              *(1.0-X)
                 HF=HF+R0+R1
                 IF (DABS(HF-HW).LT.DABS(HF)*EPS) GO TO 95
90               HW=HF
95            HF=HF+C0+C1
           ENDIF
        ELSE
           A0=1.0D0
           IF (C.GT.A.AND.C.LT.2.0D0*A.AND.
     &         C.GT.B.AND.C.LT.2.0D0*B) THEN
              A0=(1.0D0-X)**(C-A-B)
              A=C-A
              B=C-B
           ENDIF
           HF=1.0D0
           R=1.0D0
           DO 100 K=1,250
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
              HF=HF+R
              IF (DABS(HF-HW).LE.DABS(HF)*EPS) GO TO 105
100           HW=HF
105        HF=A0*HF
        ENDIF
        IF (X1.LT.0.0D0) THEN
           X=X1
           C0=1.0D0/(1.0D0-X)**AA
           HF=C0*HF
        ENDIF
        A=AA
        B=BB
        IF (K.GT.120) WRITE(*,115)
115     FORMAT(1X,'Warning! You should check the accuracy')
        RETURN
        END



        SUBROUTINE PSI(X,PS)
C
C       ======================================
C       Purpose: Compute Psi function
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       ======================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        S=0.0D0
        IF (X.EQ.INT(X).AND.X.LE.0.0) THEN
           PS=1.0D+300
           RETURN
        ELSE IF (XA.EQ.INT(XA)) THEN
           N=XA
           DO 10 K=1 ,N-1
10            S=S+1.0D0/K
           PS=-EL+S
        ELSE IF (XA+.5.EQ.INT(XA+.5)) THEN
           N=XA-.5
           DO 20 K=1,N
20            S=S+1.0/(2.0D0*K-1.0D0)
           PS=-EL+2.0D0*S-1.386294361119891D0
        ELSE
           IF (XA.LT.10.0) THEN
              N=10-INT(XA)
              DO 30 K=0,N-1
30               S=S+1.0D0/(XA+K)
              XA=XA+N
           ENDIF
           X2=1.0D0/(XA*XA)
           A1=-.8333333333333D-01
           A2=.83333333333333333D-02
           A3=-.39682539682539683D-02
           A4=.41666666666666667D-02
           A5=-.75757575757575758D-02
           A6=.21092796092796093D-01
           A7=-.83333333333333333D-01
           A8=.4432598039215686D0
           PS=DLOG(XA)-.5D0/XA+X2*(((((((A8*X2+A7)*X2+
     &        A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1)
           PS=PS-S
        ENDIF
        IF (X.LT.0.0) PS=PS-PI*DCOS(PI*X)/DSIN(PI*X)-1.0D0/X
        RETURN
        END
C=============Bessel Function for large order=========================C

       SUBROUTINE CIKLV(V,Z,Res)
C
C       =====================================================
C       Purpose: Compute modified Bessel functions Iv(z) and
C                Kv(z) and their derivatives with a complex
C                argument and a large order
C       Input:   v --- Order of Iv(z) and Kv(z)
C                z --- Complex argument
C       Output:  CBIV --- Iv(z)
C                CDIV --- Iv'(z)
C                CBKV --- Kv(z)
C                CDKV --- Kv'(z)
C       Routine called:
C                CJK to compute the expansion coefficients
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CF(12),A(91)
        PI=3.141592653589793D0
        KM=12
        CALL CJK(KM,A)
        DO 30 L=1,0,-1
           V0=V-L
           CWS=CDSQRT(1.0D0+(Z/V0)*(Z/V0))
           CETA=CWS+CDLOG(Z/V0/(1.0D0+CWS))
           CT=1.0D0/CWS
           CT2=CT*CT
           DO 15 K=1,KM
              L0=K*(K+1)/2+1
              LF=L0+K
              CF(K)=A(LF)
              DO 10 I=LF-1,L0,-1
10               CF(K)=CF(K)*CT2+A(I)
15            CF(K)=CF(K)*CT**K
           VR=1.0D0/V0
           CSI=(1.0D0,0.0D0)
           DO 20 K=1,KM
20            CSI=CSI+CF(K)*VR**K
           CBIV=CDSQRT(CT/(2.0D0*PI*V0))*CDEXP(V0*CETA)*CSI
           IF (L.EQ.1) CFI=CBIV
           CSK=(1.0D0,0.0D0)
           DO 25 K=1,KM
25            CSK=CSK+(-1)**K*CF(K)*VR**K
           CBKV=CDSQRT(PI*CT/(2.0D0*V0))*CDEXP(-V0*CETA)*CSK
           IF (L.EQ.1) CFK=CBKV
30      CONTINUE
        CDIV=CFI-V/Z*CBIV
        CDKV=-CFK-V/Z*CBKV
	  Res=real(CBKV)
        RETURN
        END


        SUBROUTINE CJK(KM,A)
C
C       ========================================================
C       Purpose: Compute the expansion coefficients for the
C                asymptotic expansion of Bessel functions 
C                with large orders
C       Input :  Km   --- Maximum k
C       Output:  A(L) --- Cj(k) where j and k are related to L 
C                         by L=j+1+[k*(k+1)]/2; j,k=0,1,...,Km
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(*)
        A(1)=1.0D0
        F0=1.0D0
        G0=1.0D0
        DO 10 K=0,KM-1
           L1=(K+1)*(K+2)/2+1
           L2=(K+1)*(K+2)/2+K+2
           F=(0.5D0*K+0.125D0/(K+1))*F0
           G=-(1.5D0*K+0.625D0/(3.0*(K+1.0D0)))*G0
           A(L1)=F
           A(L2)=G
           F0=F
10         G0=G
        DO 15 K=1,KM-1
           DO 15 J=1,K
              L3=K*(K+1)/2+J+1
              L4=(K+1)*(K+2)/2+J+1
              A(L4)=(J+0.5D0*K+0.125D0/(2.0*J+K+1.0))*A(L3)
     &             -(J+0.5D0*K-1.0+0.625D0/(2.0*J+K+1.0))*A(L3-1)
15         CONTINUE
        RETURN
        END


      subroutine quanc8(fun,a,b,abserr,relerr,result,errest,nofun,flag)
c
      double precision fun, a, b, abserr, relerr, result, errest, flag
      integer nofun
c
c   estimate the integral of fun(x) from a to b
c   to a user provided tolerance.
c   an automatic adaptive routine based on
c   the 8-panel newton-cotes rule.
c
c   input ..
c
c   fun     the name of the integrand function subprogram fun(x).
c   a       the lower limit of integration.
c   b       the upper limit of integration.(b may be less than a.)
c   relerr  a relative error tolerance. (should be non-negative)
c   abserr  an absolute error tolerance. (should be non-negative)
c
c   output ..
c
c   result  an approximation to the integral hopefully satisfying the
c           least stringent of the two error tolerances.
c   errest  an estimate of the magnitude of the actual error.
c   nofun   the number of function values used in calculation of result.
c   flag    a reliability indicator.  if flag is zero, then result
c           probably satisfies the error tolerance.  if flag is
c           xxx.yyy , then  xxx = the number of intervals which have
c           not converged and  0.yyy = the fraction of the interval
c           left to do when the limit on  nofun  was approached.
c
      double precision w0,w1,w2,w3,w4,area,x0,f0,stone,step,cor11,temp
      double precision qprev,qnow,qdiff,qleft,esterr,tolerr
      double precision qright(31),f(16),x(16),fsave(8,30),xsave(8,30)
      double precision dabs,dmax1
      integer levmin,levmax,levout,nomax,nofin,lev,nim,i,j
c
c   ***   stage 1 ***   general initialization
c   set constants.
c
      levmin = 1
      levmax = 30
      levout = 6
      nomax = 5000
      nofin = nomax - 8*(levmax-levout+2**(levout+1))
c
c   trouble when nofun reaches nofin
c
      w0 =   3956.0d0 / 14175.0d0
      w1 =  23552.0d0 / 14175.0d0
      w2 =  -3712.0d0 / 14175.0d0
      w3 =  41984.0d0 / 14175.0d0
      w4 = -18160.0d0 / 14175.0d0
c
c   initialize running sums to zero.
c
      flag = 0.0d0
      result = 0.0d0
      cor11  = 0.0d0
      errest = 0.0d0
      area   = 0.0d0
      nofun = 0
      if (a .eq. b) return
c
c   ***   stage 2 ***   initialization for first interval
c
      lev = 0
      nim = 1
      x0 = a
      x(16) = b
      qprev  = 0.0d0
      f0 = fun(x0)
      stone = (b - a) / 16.0d0
      x(8)  =  (x0  + x(16)) / 2.0d0
      x(4)  =  (x0  + x(8))  / 2.0d0
      x(12) =  (x(8)  + x(16)) / 2.0d0
      x(2)  =  (x0  + x(4))  / 2.0d0
      x(6)  =  (x(4)  + x(8))  / 2.0d0
      x(10) =  (x(8)  + x(12)) / 2.0d0
      x(14) =  (x(12) + x(16)) / 2.0d0
      do 25 j = 2, 16, 2
         f(j) = fun(x(j))
   25 continue
      nofun = 9
c
c   ***   stage 3 ***   central calculation
c   requires qprev,x0,x2,x4,...,x16,f0,f2,f4,...,f16.
c   calculates x1,x3,...x15, f1,f3,...f15,qleft,qright,qnow,qdiff,area.
c
   30 x(1) = (x0 + x(2)) / 2.0d0
      f(1) = fun(x(1))
      do 35 j = 3, 15, 2
         x(j) = (x(j-1) + x(j+1)) / 2.0d0
         f(j) = fun(x(j))
   35 continue
      nofun = nofun + 8
      step = (x(16) - x0) / 16.0d0
      qleft  =  (w0*(f0 + f(8))  + w1*(f(1)+f(7))  + w2*(f(2)+f(6))
     1  + w3*(f(3)+f(5))  +  w4*f(4)) * step
      qright(lev+1)=(w0*(f(8)+f(16))+w1*(f(9)+f(15))+w2*(f(10)+f(14))
     1  + w3*(f(11)+f(13)) + w4*f(12)) * step
      qnow = qleft + qright(lev+1)
      qdiff = qnow - qprev
      area = area + qdiff
c
c   ***   stage 4 *** interval convergence test
c
      esterr = dabs(qdiff) / 1023.0d0
      tolerr = dmax1(abserr,relerr*dabs(area)) * (step/stone)
      if (lev .lt. levmin) go to 50
      if (lev .ge. levmax) go to 62
      if (nofun .gt. nofin) go to 60
      if (esterr .le. tolerr) go to 70
c
c   ***   stage 5   ***   no convergence
c   locate next interval.
c
   50 nim = 2*nim
      lev = lev+1
c
c   store right hand elements for future use.
c
      do 52 i = 1, 8
         fsave(i,lev) = f(i+8)
         xsave(i,lev) = x(i+8)
   52 continue
c
c   assemble left hand elements for immediate use.
c
      qprev = qleft
      do 55 i = 1, 8
         j = -i
         f(2*j+18) = f(j+9)
         x(2*j+18) = x(j+9)
   55 continue
      go to 30
c
c   ***   stage 6   ***   trouble section
c   number of function values is about to exceed limit.
c
   60 nofin = 2*nofin
      levmax = levout
      flag = flag + (b - x0) / (b - a)
      go to 70
c
c   current level is levmax.
c
   62 flag = flag + 1.0d0
c
c   ***   stage 7   ***   interval converged
c   add contributions into running sums.
c
   70 result = result + qnow
      errest = errest + esterr
      cor11  = cor11  + qdiff / 1023.0d0
c
c   locate next interval.
c
   72 if (nim .eq. 2*(nim/2)) go to 75
      nim = nim/2
      lev = lev-1
      go to 72
   75 nim = nim + 1
      if (lev .le. 0) go to 80
c
c   assemble elements required for the next interval.
c
      qprev = qright(lev)
      x0 = x(16)
      f0 = f(16)
      do 78 i = 1, 8
         f(2*i) = fsave(i,lev)
         x(2*i) = xsave(i,lev)
   78 continue
      go to 30
c
c   ***   stage 8   ***   finalize and return
c
   80 result = result + cor11
c
c   make sure errest not less than roundoff level.
c
      if (errest .eq. 0.0d0) return
   82 temp = dabs(result) + errest
      if (temp .ne. dabs(result)) return
      errest = 2.0d0*errest
      go to 82
      end

      subroutine quagen(z,wt,n)     
      integer n
      real*8 z(n),wt(n)
      real*8 p(513),c1(513),c2(513)
      real*8 pi,pi4,pdir,xnow,const,epsilon
      logical found
      data p(1),epsilon,nmax,ip /1.d0,1.d-10,10,512/

      if(n.le.0.or.n.gt.ip) then
	 write(*,*) 'input order out of range'
       stop
      endif
c
c     calculate the coefficients for the legendre poly. recursive formula
c   
      do 10 i=2,n
      c1(i)=(2*i-1)/dfloat(i)
      c2(i)=(i-1)/dfloat(i)
 10   continue

c
c     initial constants
c
      n1=n+1
      pi=acos(-1.d0)
      pi4=pi/4
      const=1.d0/(n+0.5)
c
c     determine the number of roots(nr) needed to be calculated
c
      n2=n/2
      if(n2*2.eq.n) then
       ncal=n2
      else
	ncal=n2+1
      endif
c 
c     main loop begins here
c   
      do i=1,ncal
       k=n-i+1
       ncount=0
       xinc=1.d0
c
c     use newton's method and a good initial guess to locate the root
c
       xnow=cos((i*pi-pi4)*const)
       found=.false.
      do while (.not.found)
	 ncount=ncount+1
	  p(2)=xnow
c
c     the following loop calculate p_n(x) using recursive formula
c
       do 20 j=2,n
	p(j+1)=c1(j)*xnow*p(j)-c2(j)*p(j-1) 
 20    continue

c
c     the derivate of p_n(x) can be calculated from p_n(x) and p_n-1(x)
c
       pdir=n*(p(n)-xnow*p(n1))/(1.d0-xnow*xnow)
	if(abs(xinc).le.epsilon.or.ncount.gt.nmax) then
	  found=.true.
	 z(k)=xnow
	 z(i)=-xnow
       wt(k)=2.d0/(1.d0-z(k)*z(k))/(pdir*pdir)
       wt(i)=wt(k)
      endif
      xinc=-p(n1)/pdir
      xnow=xnow+xinc
      enddo
      enddo
      return
      end
