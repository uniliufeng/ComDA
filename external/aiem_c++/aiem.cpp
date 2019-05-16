
#include "aiem.h"

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
//        Advanced Integral Equation Model                  //
//                                                          //
//			Program computes emissivity                     //
//		from 3d finitely conducting surface                 //
//			( only for single scattering )                  //
//                                                          //
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
//                                                          //
// Parameters:                                              //
//                                                          //
//    er: surface relative dielectric constant              //
//    kl: normalized surface correlation length             //
//    ks: normalized surface rms height                     //
//                                                          //
//    itype: select what type of surface correlation        //
//          =1  Gaussian correlation function               //
//          =2  exponential correlation function            //
//          =3  transformed exponential correlation         //
//                                                          //
//    theta: incident angle in deg                          //
//    (phi=0 : incident azimuth angle in deg)               //
//                                                          //
//                                                          //
//**********************************************************//
//       Approximations of Fresnel reflection coeff.        //
//       -------------------------------------------        //
//                                                          //
//   irc=1: Fresnel reflection coeff. approxmated by        //
//          R(incident_angle)                               //
//   irc=2: Fresnel reflection coeff. approxmated by        //
//          R(specular_angle)                               //
//                                                          //
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
//                                                          //
//                                  Sept. 18, 2001          //
//                                                          //
//                                            T. D. Wu      //
//           rewrite by Jin Rui	in may 2005					//
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

double pi=acos(-1.0);

void AIEM::emissivity()
{
    int number,irc,npp,nss,is,j;
    double ww11[513],zz11[513],ww22[513],zz22[513];
    double wt[513],zt[513];
    double sigma0[5];
    double freq,cl,sig,err,eri,tmp1,kl,ks,cs,si2,theta,thetas,phis,aa11,aa22,bb11,bb22,
           as11,bs11,as22,bs22,reflvv,reflhh,cohh,covv,hh,vv,hv,vh,emisvv,emishh;
    complex<double> stem,rvi,rhi;
    FILE *in,*out;
    FILE *temp;
    double pi=acos(-1.0);

    //npp=128;
    //nss=128;
    npp=16;         //npp=32
    nss=16;         //nss=32
    ur=1.0;
//----------------------------------------------------------//
//            open files                                    //
//----------------------------------------------------------//

    if((in=fopen("0P10.in","r"))==NULL)
        cout<<"0P10.in opend failed!"<<endl;
    else
        cout<<"0P10.in opend succesful!"<<endl;

    if((out=fopen("AIEM0P10_simple.txt","w"))==NULL)
        cout<<"AIEM0P10.txt opend failed!"<<endl;
    else
        cout<<"AIEM0P10.txt opend succesful!"<<endl;



//----------------------------------------------------------//
//         input   Parameters                               //
//----------------------------------------------------------//

    fscanf(in,"%d",&number);                //number of computation
    cout<<number<<endl;

    while(!feof(in))
    {
        //f1

        fscanf(in,"%lf %lf %lf %lf %lf %d %lf",&freq,&cl,&sig,&err,&eri,&itype,&tmp1);

//      freq:	frequency in GHz
//      cl:		correlation length in wavelength
//      sig:	rms height in wavelength
//		err:	real part of permittivity
//      eri:	imaginary part of permittivity
//      itype=1     Gaussian correlated surface
//      itype=2     Exponential correlated surface
//      itype=3     Transformed exponential correlated surface
//		tmp1:       temp variable??

        kl=2.0*pi*cl;           //kl=wavenumber*correlation length=2*pi*correlation length/wavelenght
        ks=2.0*pi*sig;          //ks=wavenumber*rms height=2*pi*rms height/wavelenght

        er=complex<double>(err, eri);            //permittivity(dielectric constant) of object
        ur=complex<double>(1.0,0.0);             //permeability of object

        theta=50.0;		// incident angle
        phis=0.0;		// scattering azimuth angle

//----------------------------------------------------------//
//       Approximations of Fresnel reflection coeff.        //
//----------------------------------------------------------//
        irc=1;      // approxmated by R(theta)
//      irc=2       // approxmated by R(theta_specular)


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//reflection coefficients based on the incident angle    //
//==> R(theta)                                           //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
        cs=cos(theta*pi/180.0);
        si2=pow(sin(theta*pi/180.0),2.0);
        stem=sqrt(er*ur-si2);
        rvi=(er*cs-stem)/(er*cs+stem);
        rhi=(ur*cs-stem)/(ur*cs+stem);

//------------------------------------------//
//             main program                 //
//------------------------------------------//

//----------------------------------------- //
//   set up integration limits              //
//----------------------------------------- //
        aa11=-3.14;
        bb11=3.14;
        aa22=0.0;
        bb22=1.57;
//----------------------------------------//
// generates zeros and absissa            //
//----------------------------------------//

        quagen(zt,wt,npp);                //子程序改变zt和wt值
        as11=(bb11-aa11)/2.0;             //3.14
        bs11=(bb11+aa11)/2.0;             //0
        for(int i=1; i<=npp; i++)
        {
            ww11[i]=wt[i]*as11;
            zz11[i]=zt[i]*as11+bs11;
        }

        quagen(zt,wt,nss);
        as22=(bb22-aa22)/2.0;              //0.785
        bs22=(bb22+aa22)/2.0;              //0.785
        for(int i=1; i<=nss; i++)
        {
            ww22[i]=wt[i]*as22;
            zz22[i]=zt[i]*as22+bs22;
        }

        reflvv=0.0;           //incoherent term ???
        reflhh=0.0;           //incoherent term ???

//-------------------------------------------------------------//
//         coherent component                                  //
//-------------------------------------------------------------//

        //reference: liou et. al. IEEE, 2001,39(1):129-135
        cohh=(pow(abs(rhi),2.0))*exp(-pow(2.0*cos(theta*pi/180.0)*ks,2.0));   //specular coherent term corrected by roughness factor
        covv=(pow(abs(rvi),2.0))*exp(-pow(2.0*cos(theta*pi/180.0)*ks,2.0));   //specular coherent term corrected by roughness factor

        for(int j=1; j<=npp; j++)                           //npp=128 散射方位角离散化为128份
        {
            //f2
            phis=zz11[j]*180.0/pi;                    //phis: -pi~pi
            if(phis==180.0)
                phis=179.99;

            for(is=1; is<=nss; is++)
            {
                //f3
                thetas=zz22[is]*180.0/pi;				//specular(scatter) angle??
                if(thetas==90.0)                      //thetas: 0~pi/2
                    thetas=89.999;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//         initialize array                   //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

                hh=0.0;                       //bistatic scatter coef.
                vv=0.0;
                hv=0.0;
                vh=0.0;

                if(thetas==theta)
                    thetas=thetas+0.001;      //not compute the specular direction


//---------------------------------------------------------------//
                sigma(ks,kl,theta,thetas,phis,sigma0,irc);       // change sigma0
//---------------------------------------------------------------//

                if(sigma0[1]>0.0) vv=sigma0[1];
                if(sigma0[2]>0.0) hh=sigma0[2];
                if(sigma0[3]>0.0) hv=sigma0[3];
                if(sigma0[4]>0.0) vh=sigma0[4];

                //双基散射系数在半空间的球面积分  Gauss求积法
                reflhh=reflhh+ww11[j]*ww22[is]*(hh+vh)*sin(thetas*pi/180.0)/(4.0*pi*cos(theta*pi/180.0));   //ww11??  ww22??
                reflvv=reflvv+ww11[j]*ww22[is]*(vv+hv)*sin(thetas*pi/180.0)/(4.0*pi*cos(theta*pi/180.0));

            }//f3
            //cout<<j<<endl;
        }//f2

        emisvv=1.0-(reflvv+covv);
        emishh=1.0-(reflhh+cohh);

        fprintf(out,"%10.4f	%10.4f",emisvv,emishh);
        cout<<emisvv<<" "<<emishh<<endl;

    }//f1
}

//********************************************************************//
//  subroutine sigma calculates the bistatic scattering coefficients  //
//********************************************************************//

void AIEM::sigma(double ks,double kl,double theta,double thetas,double phis,double sigma0[],int irc)
{
    //sigma

    complex<double> stem,steml,fvv,fhh,fvh,fhv;
    complex<double> rvi,rhi,rvhi,rvl,rhl,rvhl;
    complex<double> Fahh,Favv,Favh,Fahv;
    complex<double> Fbhh,Fbvv,Fbvh,Fbhv;
    complex<double> expkc1,expkc2,expc1,expc2,expc3,expc4;
    complex<double> expkcaup,expkcadn,expkcbup,expkcbdn;
    complex<double> expcauau,expcadau,expcauad,expcadad;
    complex<double> expcbuau,expcbdau,expcbuad,expcbdad;
    complex<double> expcaubu,expcadbu,expcaubd,expcadbd;
    complex<double> expcbubu,expcbdbu,expcbubd,expcbdbd;
    complex<double> qq,qqs,qqt,qqts;
    complex<double> kterm[5],kcterm[5],cterm[5];
    double torlant,kl2,csl,sil,tempold,temp,fiterm,fn,ql,e,y,gam,m,bk,out,sum,expk,ex1,ex2,term1,term2;

    double pi=acos(-1.0);
    torlant=1.0e-16;

    cs=cos(theta*pi/180.0);
    css=cos(thetas*pi/180.0);
    si=sin(theta*pi/180.0);
    sis=sin(thetas*pi/180.0);
    csfs=cos(phis*pi/180.0);
    sfs=sin(phis*pi/180.0);
    cs2=cs*cs;
    css2=css*css;
    si2=si*si;
    sis2=sis*sis;

    ks2=ks*ks;
    kl2=kl*kl;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//reflection coefficients based on the incident angle    //
//==> R(theta)                                           //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    stem=sqrt(er*ur-si2);
    rvi=(er*cs-stem)/(er*cs+stem);                                 //Frensel equation
    rhi=(ur*cs-stem)/(ur*cs+stem);                                 //Frensel equation
    rvhi=(rvi-rhi)/2.0;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//reflection coefficients based on the specular angle    //
//==> R(theta_specular)                                  //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    csl=sqrt(1.0+cs*css-si*sis*csfs)/sqrt(2.0);          //计算入射矢量和散射矢量的夹角的余旋。因为两个矢量方向相反，所以在点乘（cs*css-si*sis*csfs）之后还要＋1
    sil=sqrt(1.0-csl*csl);
    steml=sqrt(er*ur-sil*sil);
    rvl=(er*csl-steml)/(er*csl+steml);                   //specular angle==scatter angle和incident angle的夹角
    rhl=(ur*csl-steml)/(ur*csl+steml);
    rvhl=(rvl-rhl)/2.0;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//    Reflection coefficients rv, rh, rvh                //
//          for kirchhoff field coefficients             //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    if(irc==1)
    {
        rh=rhi;
        rv=rvi;
        rvh=rvhi;
    }
    if(irc==2)
    {
        rh=rhl;
        rv=rvl;
        rvh=rvhl;
    }


//--------------------------------------------//
//    kirchhoff field coefficients            //
//--------------------------------------------//

    fvv=2.0*rv*(si*sis-(1.0+cs*css)*csfs)/(css+cs);                //2A.4 in <<microwave remote sensing>>, ulaby et. al
    fhh=-2.0*rh*(si*sis-(1.0+cs*css)*csfs)/(css+cs);               //2A.5
    fvh=-2.0*rvh*sfs;                                              //2A.6
    fhv=2.0*rvh*sfs;                                               //2A.7

//--------------------------------------------//
//    compute roughness spectrum w(n)         //
//--------------------------------------------//
    iterm=1;
    tempold=0.0;
    temp=(ks2*pow((cs+css),2.0));

    while(fabs(temp-tempold)>torlant)     //certain the value of iterm??
    {
        tempold=temp;
        iterm=iterm+1;
        fiterm=float(iterm);
        temp=tempold*(ks2*pow((cs+css),2.0))/fiterm;
    }

    //printf("the iterm is %d\n",iterm);

    for(int n=1; n<=1001; n++)
    {
        w[n]=0.0;
    }

    for(int n=1; n<=iterm; n++)
    {
        //f5
        fn=float(n);
        ql=kl*sqrt(pow((-si+sis*csfs),2.0)+pow((0.0+sis*sfs),2.0));

        if((itype>4)||(itype<1))
            itype=1;

        if(itype==1)          //Gaussian correlated surface
            w[n]=kl2*exp(-ql*ql/(4.0*fn))/(2*fn);                //2B.4  in ulaby et. al

        if(itype==2)         //exponential correlated surface
        {
            term1=pow((kl/fn),2.0);
            term2=pow((1.0+pow(ql/fn,2.0)),1.5);
            w[n]=pow(kl/fn,2.0)/term2;                          //2B.14	 in ulaby et. al
        }

        if(itype==3)         //1.5-power correlated surface
        {
            //itype=3

            e=1.5*fn-1.0;
            y=1.5*fn;
            gam=alogam(y);          //!gamma function (1.5n)

            if(ql==0.0)
                w[n]=kl*kl/(3.0*fn-2.0);
            else
            {
                //else
                if(fmod(fn,2.0)==0.0)
                {
                    m=1.5*fn-1.0;              //fn is even;
                    bk=log(BESSK(m,ql));       //int order, check OK
                }
                else
                {
                    m=1.5*fn-1.0-0.5;            //fn is odd
                    bk=log(BesselK(m,ql));       //!fractional order
                }

                out=kl*kl*pow((ql/2.0),e);
                w[n]=out*exp(bk-gam);                 //2B.24???  in ulaby et. al
            }//else

        }//itype=3

    }//f5


//-----------------------------------------//
//   compute kirchhoff term                //
//-----------------------------------------//
    sum=0.0;
    temp=1.0;

    for(int n=1; n<=iterm; n++)
    {
        fn=float(n);
        temp=temp*(ks2*pow((cs+css),2.0))/fn;
        sum=sum+temp*w[n];
    }

    expk=exp(-ks2*pow((css+cs),2.0))*sum;
    kterm[1]=0.5*expk*pow(abs(fvv),2);                       //equation5.3 in Fung PP233
    kterm[2]=0.5*expk*pow(abs(fhh),2);
    kterm[3]=0.5*expk*pow(abs(fhv),2);
    kterm[4]=0.5*expk*pow(abs(fvh),2);

//-----------------------------------------//
// end of kirchhoff term computation       //
//-----------------------------------------//

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//    Reflection coefficients rv, rh, rvh                //
//          for complementary field coefficients         //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    if(irc==1)
    {
        rh=rhi;
        rv=rvi;
        rvh=rvhi;
    }

    if(irc==2)
    {
        rh=rhl;
        rv=rvl;
        rvh=rvhl;
    }
//-------------------------------------------------------------//
    ex1=exp(-ks2*(cs2+css2+cs*css));
    ex2=exp(-ks2*(cs2+css2));
    qq=cs;
    qqt=sqrt(er-si2);
    qqs=css;
    qqts=sqrt(er-sis2);
//-------------------------------------------------------------//

    //-----------------------------------------//
    //  compute cross term                     //
    //-----------------------------------------//

    //kcterm(1)	  vv

    expkcaup=(conj(fvv)*AIEM::Favv(-si,0.0,qq,qq)).real()*AIEM::expkc1(qq)
             +(conj(fvv)*AIEM::Favv(-sis*csfs,-sis*sfs,qqs,qqs)).real()*AIEM::expkc2(qqs);

    expkcadn=(conj(fvv)*AIEM::Favv(-si,0.0,-qq,qq)).real()*AIEM::expkc1(-qq)
             +(conj(fvv)*AIEM::Favv(-sis*csfs,-sis*sfs,-qqs,qqs)).real()*AIEM::expkc2(-qqs);

    expkcbup=(conj(fvv)*AIEM::Fbvv(-si,0.0,qqt,qqt)).real()*AIEM::expkc1(qqt)
             +(conj(fvv)*AIEM::Fbvv(-sis*csfs,-sis*sfs,qqts,qqts)).real()*AIEM::expkc2(qqts);

    expkcbdn=(conj(fvv)*AIEM::Fbvv(-si,0.0,-qqt,qqt)).real()*AIEM::expkc1(-qqt)
             +(conj(fvv)*AIEM::Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts)).real()*AIEM::expkc2(-qqts);

    kcterm[1]=0.25*ex1*(expkcaup+expkcadn+expkcbup+expkcbdn);

    //kcterm(2)  hh
    expkcaup=(conj(fhh)*AIEM::Fahh(-si,0.0,qq,qq)).real()*AIEM::expkc1(qq)
             +(conj(fhh)*AIEM::Fahh(-sis*csfs,-sis*sfs,qqs,qqs)).real()*AIEM::expkc2(qqs);

    expkcadn=(conj(fhh)*AIEM::Fahh(-si,0.0,-qq,qq)).real()*AIEM::expkc1(-qq)
             +(conj(fhh)*AIEM::Fahh(-sis*csfs,-sis*sfs,-qqs,qqs)).real()*AIEM::expkc2(-qqs);

    expkcbup=(conj(fhh)*AIEM::Fbhh(-si,0.0,qqt,qqt)).real()*AIEM::expkc1(qqt)
             +(conj(fhh)*AIEM::Fbhh(-sis*csfs,-sis*sfs,qqts,qqts)).real()*AIEM::expkc2(qqts);

    expkcbdn=(conj(fhh)*AIEM::Fbhh(-si,0.0,-qqt,qqt)).real()*AIEM::expkc1(-qqt)
             +(conj(fhh)*AIEM::Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts)).real()*AIEM::expkc2(-qqts);

    kcterm[2]=0.25*ex1*(expkcaup+expkcadn+expkcbup+expkcbdn);

    //kcterm(3)  hv
    expkcaup=(conj(fhv)*AIEM::Fahv(-si,0.0,qq,qq)).real()*AIEM::expkc1(qq)
             +(conj(fhv)*AIEM::Fahv(-sis*csfs,-sis*sfs,qqs,qqs)).real()*AIEM::expkc2(qqs);

    expkcadn=(conj(fhv)*AIEM::Fahv(-si,0.0,-qq,qq)).real()*AIEM::expkc1(-qq)
             +(conj(fhv)*AIEM::Fahv(-sis*csfs,-sis*sfs,-qqs,qqs)).real()*AIEM::expkc2(-qqs);

    expkcbup=(conj(fhv)*AIEM::Fbhv(-si,0.0,qqt,qqt)).real()*AIEM::expkc1(qqt)
             +(conj(fhv)*AIEM::Fbhv(-sis*csfs,-sis*sfs,qqts,qqts)).real()*AIEM::expkc2(qqts);

    expkcbdn=(conj(fhv)*AIEM::Fbhv(-si,0.0,-qqt,qqt)).real()*AIEM::expkc1(-qqt)
             +(conj(fhv)*AIEM::Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts)).real()*AIEM::expkc2(-qqts);

    kcterm[3]=0.25*ex1*(expkcaup+expkcadn+expkcbup+expkcbdn);

    //kcterm(4)  vh
    expkcaup=(conj(fvh)*AIEM::Favh(-si,0.0,qq,qq)).real()*AIEM::expkc1(qq)
             +(conj(fvh)*AIEM::Favh(-sis*csfs,-sis*sfs,qqs,qqs)).real()*AIEM::expkc2(qqs);

    expkcadn=(conj(fvh)*AIEM::Favh(-si,0.0,-qq,qq)).real()*AIEM::expkc1(-qq)
             +(conj(fvh)*AIEM::Favh(-sis*csfs,-sis*sfs,-qqs,qqs)).real()*AIEM::expkc2(-qqs);

    expkcbup=(conj(fvh)*AIEM::Fbvh(-si,0.0,qqt,qqt)).real()*AIEM::expkc1(qqt)
             +(conj(fvh)*AIEM::Fbvh(-sis*csfs,-sis*sfs,qqts,qqts)).real()*AIEM::expkc2(qqts);

    expkcbdn=(conj(fvh)*AIEM::Fbvh(-si,0.0,-qqt,qqt)).real()*AIEM::expkc1(-qqt)
             +(conj(fvh)*AIEM::Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts)).real()*AIEM::expkc2(-qqts);

    kcterm[4]=0.25*ex1*(expkcaup+expkcadn+expkcbup+expkcbdn);

//------------------------------------------------------------------//
//                  end of computation of cross terms               //
//------------------------------------------------------------------//

//------------------------------------------------------------------//
//                     evaluate  complementary term                 //
//------------------------------------------------------------------//

    //cterm(1) vv
    expcauau=AIEM::Favv(-si,0.0,qq,qq)
             *conj(AIEM::Favv(-si,0.0,qq,qq))*AIEM::expc1(qq,qq)
             +AIEM::Favv(-si,0.0,qq,qq)
             *conj(AIEM::Favv(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc2(qq,qqs)
             +AIEM::Favv(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Favv(-si,0.0,qq,qq))*AIEM::expc3(qqs,qq)
             +AIEM::Favv(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Favv(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc4(qqs,qqs);

    expcadau=AIEM::Favv(-si,0.0,-qq,qq)
             *conj(AIEM::Favv(-si,0.0,qq,qq))*AIEM::expc1(-qq,qq)
             +AIEM::Favv(-si,0.0,-qq,qq)
             *conj(AIEM::Favv(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc2(-qq,qqs)
             +AIEM::Favv(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Favv(-si,0.0,qq,qq))*AIEM::expc3(-qqs,qq)
             +AIEM::Favv(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Favv(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc4(-qqs,qqs);

    expcauad=AIEM::Favv(-si,0.0,qq,qq)
             *conj(AIEM::Favv(-si,0.0,-qq,qq))*AIEM::expc1(qq,-qq)
             +AIEM::Favv(-si,0.0,qq,qq)
             *conj(AIEM::Favv(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc2(qq,-qqs)
             +AIEM::Favv(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Favv(-si,0.0,-qq,qq))*AIEM::expc3(qqs,-qq)
             +AIEM::Favv(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Favv(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc4(qqs,-qqs);

    expcadad=AIEM::Favv(-si,0.0,-qq,qq)
             *conj(AIEM::Favv(-si,0.0,-qq,qq))*AIEM::expc1(-qq,-qq)
             +AIEM::Favv(-si,0.0,-qq,qq)
             *conj(AIEM::Favv(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc2(-qq,-qqs)
             +AIEM::Favv(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Favv(-si,0.0,-qq,qq))*AIEM::expc3(-qqs,-qq)
             +AIEM::Favv(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Favv(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc4(-qqs,-qqs);

    expcbuau=AIEM::Fbvv(-si,0.0,qqt,qqt)
             *conj(AIEM::Favv(-si,0.0,qq,qq))*AIEM::expc1(qqt,qq)
             +AIEM::Fbvv(-si,0.0,qqt,qqt)
             *conj(AIEM::Favv(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc2(qqt,qqs)
             +AIEM::Fbvv(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Favv(-si,0.0,qq,qq))*AIEM::expc3(qqts,qq)
             +AIEM::Fbvv(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Favv(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc4(qqts,qqs);

    expcbdau=AIEM::Fbvv(-si,0.0,-qqt,qqt)
             *conj(AIEM::Favv(-si,0.0,qq,qq))*AIEM::expc1(-qqt,qq)
             +AIEM::Fbvv(-si,0.0,-qqt,qqt)
             *conj(AIEM::Favv(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc2(-qqt,qqs)
             +AIEM::Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Favv(-si,0.0,qq,qq))*AIEM::expc3(-qqts,qq)
             +AIEM::Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Favv(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc4(-qqts,qqs);

    expcbuad=AIEM::Fbvv(-si,0.0,qqt,qqt)
             *conj(AIEM::Favv(-si,0.0,-qq,qq))*AIEM::expc1(qqt,-qq)
             +AIEM::Fbvv(-si,0.0,qqt,qqt)
             *conj(AIEM::Favv(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc2(qqt,-qqs)
             +AIEM::Fbvv(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Favv(-si,0.0,-qq,qq))*AIEM::expc3(qqts,-qq)
             +AIEM::Fbvv(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Favv(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc4(qqts,-qqs);

    expcbdad=AIEM::Fbvv(-si,0.0,-qqt,qqt)
             *conj(AIEM::Favv(-si,0.0,-qq,qq))*AIEM::expc1(-qqt,-qq)
             +AIEM::Fbvv(-si,0.0,-qqt,qqt)
             *conj(AIEM::Favv(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc2(-qqt,-qqs)
             +AIEM::Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Favv(-si,0.0,-qq,qq))*AIEM::expc3(-qqts,-qq)
             +AIEM::Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Favv(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc4(-qqts,-qqs);

    expcaubu=AIEM::Favv(-si,0.0,qq,qq)
             *conj(AIEM::Fbvv(-si,0.0,qqt,qqt))*AIEM::expc1(qq,qqt)
             +AIEM::Favv(-si,0.0,qq,qq)
             *conj(AIEM::Fbvv(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc2(qq,qqts)
             +AIEM::Favv(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fbvv(-si,0.0,qqt,qqt))*AIEM::expc3(qqs,qqt)
             +AIEM::Favv(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fbvv(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc4(qqs,qqts);

    expcadbu=AIEM::Favv(-si,0.0,-qq,qq)
             *conj(AIEM::Fbvv(-si,0.0,qqt,qqt))*AIEM::expc1(-qq,qqt)
             +AIEM::Favv(-si,0.0,-qq,qq)
             *conj(AIEM::Fbvv(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc2(-qq,qqts)
             +AIEM::Favv(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fbvv(-si,0.0,qqt,qqt))*AIEM::expc3(-qqs,qqt)
             +AIEM::Favv(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fbvv(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc4(-qqs,qqts);

    expcaubd=AIEM::Favv(-si,0.0,qq,qq)
             *conj(AIEM::Fbvv(-si,0.0,-qqt,qqt))*AIEM::expc1(qq,-qqt)
             +AIEM::Favv(-si,0.0,qq,qq)
             *conj(AIEM::Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc2(qq,-qqts)
             +AIEM::Favv(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fbvv(-si,0.0,-qqt,qqt))*AIEM::expc3(qqs,-qqt)
             +AIEM::Favv(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc4(qqs,-qqts);

    expcadbd=AIEM::Favv(-si,0.0,-qq,qq)
             *conj(AIEM::Fbvv(-si,0.0,-qqt,qqt))*AIEM::expc1(-qq,-qqt)
             +AIEM::Favv(-si,0.0,-qq,qq)
             *conj(AIEM::Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc2(-qq,-qqts)
             +AIEM::Favv(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fbvv(-si,0.0,-qqt,qqt))*AIEM::expc3(-qqs,-qqt)
             +AIEM::Favv(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc4(-qqs,-qqts);

    expcbubu=AIEM::Fbvv(-si,0.0,qqt,qqt)
             *conj(AIEM::Fbvv(-si,0.0,qqt,qqt))*AIEM::expc1(qqt,qqt)
             +AIEM::Fbvv(-si,0.0,qqt,qqt)
             *conj(AIEM::Fbvv(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc2(qqt,qqts)
             +AIEM::Fbvv(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fbvv(-si,0.0,qqt,qqt))*AIEM::expc3(qqts,qqt)
             +AIEM::Fbvv(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fbvv(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc4(qqts,qqts);

    expcbdbu=AIEM::Fbvv(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fbvv(-si,0.0,qqt,qqt))*AIEM::expc1(-qqt,qqt)
             +AIEM::Fbvv(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fbvv(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc2(-qqt,qqts)
             +AIEM::Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fbvv(-si,0.0,qqt,qqt))*AIEM::expc3(-qqts,qqt)
             +AIEM::Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fbvv(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc4(-qqts,qqts);

    expcbubd=AIEM::Fbvv(-si,0.0,qqt,qqt)
             *conj(AIEM::Fbvv(-si,0.0,-qqt,qqt))*AIEM::expc1(qqt,-qqt)
             +AIEM::Fbvv(-si,0.0,qqt,qqt)
             *conj(AIEM::Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc2(qqt,-qqts)
             +AIEM::Fbvv(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fbvv(-si,0.0,-qqt,qqt))*AIEM::expc3(qqts,-qqt)
             +AIEM::Fbvv(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc4(qqts,-qqts);

    expcbdbd=AIEM::Fbvv(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fbvv(-si,0.0,-qqt,qqt))*AIEM::expc1(-qqt,-qqt)
             +AIEM::Fbvv(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc2(-qqt,-qqts)
             +AIEM::Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fbvv(-si,0.0,-qqt,qqt))*AIEM::expc3(-qqts,-qqt)
             +AIEM::Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc4(-qqts,-qqts);

    cterm[1]=0.03125*ex2*(expcauau+expcadau+expcauad+expcadad
                          +expcbuau+expcbdau+expcbuad+expcbdad
                          +expcaubu+expcadbu+expcaubd+expcadbd
                          +expcbubu+expcbdbu+expcbubd+expcbdbd);     //16 terms


    //cterm(2)
    expcauau=AIEM::Fahh(-si,0.0,qq,qq)
             *conj(AIEM::Fahh(-si,0.0,qq,qq))*AIEM::expc1(qq,qq)
             +AIEM::Fahh(-si,0.0,qq,qq)
             *conj(AIEM::Fahh(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc2(qq,qqs)
             +AIEM::Fahh(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fahh(-si,0.0,qq,qq))*AIEM::expc3(qqs,qq)
             +AIEM::Fahh(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fahh(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc4(qqs,qqs);

    expcadau=AIEM::Fahh(-si,0.0,-qq,qq)
             *conj(AIEM::Fahh(-si,0.0,qq,qq))*AIEM::expc1(-qq,qq)
             +AIEM::Fahh(-si,0.0,-qq,qq)
             *conj(AIEM::Fahh(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc2(-qq,qqs)
             +AIEM::Fahh(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fahh(-si,0.0,qq,qq))*AIEM::expc3(-qqs,qq)
             +AIEM::Fahh(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fahh(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc4(-qqs,qqs);

    expcauad=AIEM::Fahh(-si,0.0,qq,qq)
             *conj(AIEM::Fahh(-si,0.0,-qq,qq))*AIEM::expc1(qq,-qq)
             +AIEM::Fahh(-si,0.0,qq,qq)
             *conj(AIEM::Fahh(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc2(qq,-qqs)
             +AIEM::Fahh(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fahh(-si,0.0,-qq,qq))*AIEM::expc3(qqs,-qq)
             +AIEM::Fahh(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fahh(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc4(qqs,-qqs);

    expcadad=AIEM::Fahh(-si,0.0,-qq,qq)
             *conj(AIEM::Fahh(-si,0.0,-qq,qq))*AIEM::expc1(-qq,-qq)
             +AIEM::Fahh(-si,0.0,-qq,qq)
             *conj(AIEM::Fahh(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc2(-qq,-qqs)
             +AIEM::Fahh(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fahh(-si,0.0,-qq,qq))*AIEM::expc3(-qqs,-qq)
             +AIEM::Fahh(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fahh(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc4(-qqs,-qqs);

    expcbuau=AIEM::Fbhh(-si,0.0,qqt,qqt)
             *conj(AIEM::Fahh(-si,0.0,qq,qq))*AIEM::expc1(qqt,qq)
             +AIEM::Fbhh(-si,0.0,qqt,qqt)
             *conj(AIEM::Fahh(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc2(qqt,qqs)
             +AIEM::Fbhh(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fahh(-si,0.0,qq,qq))*AIEM::expc3(qqts,qq)
             +AIEM::Fbhh(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fahh(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc4(qqts,qqs);

    expcbdau=AIEM::Fbhh(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fahh(-si,0.0,qq,qq))*AIEM::expc1(-qqt,qq)
             +AIEM::Fbhh(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fahh(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc2(-qqt,qqs)
             +AIEM::Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fahh(-si,0.0,qq,qq))*AIEM::expc3(-qqts,qq)
             +AIEM::Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fahh(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc4(-qqts,qqs);

    expcbuad=AIEM::Fbhh(-si,0.0,qqt,qqt)
             *conj(AIEM::Fahh(-si,0.0,-qq,qq))*AIEM::expc1(qqt,-qq)
             +AIEM::Fbhh(-si,0.0,qqt,qqt)
             *conj(AIEM::Fahh(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc2(qqt,-qqs)
             +AIEM::Fbhh(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fahh(-si,0.0,-qq,qq))*AIEM::expc3(qqts,-qq)
             +AIEM::Fbhh(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fahh(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc4(qqts,-qqs);

    expcbdad=AIEM::Fbhh(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fahh(-si,0.0,-qq,qq))*AIEM::expc1(-qqt,-qq)
             +AIEM::Fbhh(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fahh(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc2(-qqt,-qqs)
             +AIEM::Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fahh(-si,0.0,-qq,qq))*AIEM::expc3(-qqts,-qq)
             +AIEM::Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fahh(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc4(-qqts,-qqs);

    expcaubu=AIEM::Fahh(-si,0.0,qq,qq)
             *conj(AIEM::Fbhh(-si,0.0,qqt,qqt))*AIEM::expc1(qq,qqt)
             +AIEM::Fahh(-si,0.0,qq,qq)
             *conj(AIEM::Fbhh(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc2(qq,qqts)
             +AIEM::Fahh(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fbhh(-si,0.0,qqt,qqt))*AIEM::expc3(qqs,qqt)
             +AIEM::Fahh(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fbhh(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc4(qqs,qqts);

    expcadbu=AIEM::Fahh(-si,0.0,-qq,qq)
             *conj(AIEM::Fbhh(-si,0.0,qqt,qqt))*AIEM::expc1(-qq,qqt)
             +AIEM::Fahh(-si,0.0,-qq,qq)
             *conj(AIEM::Fbhh(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc2(-qq,qqts)
             +AIEM::Fahh(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fbhh(-si,0.0,qqt,qqt))*AIEM::expc3(-qqs,qqt)
             +AIEM::Fahh(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fbhh(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc4(-qqs,qqts);

    expcaubd=AIEM::Fahh(-si,0.0,qq,qq)
             *conj(AIEM::Fbhh(-si,0.0,-qqt,qqt))*AIEM::expc1(qq,-qqt)
             +AIEM::Fahh(-si,0.0,qq,qq)
             *conj(AIEM::Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc2(qq,-qqts)
             +AIEM::Fahh(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fbhh(-si,0.0,-qqt,qqt))*AIEM::expc3(qqs,-qqt)
             +AIEM::Fahh(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc4(qqs,-qqts);

    expcadbd=AIEM::Fahh(-si,0.0,-qq,qq)
             *conj(AIEM::Fbhh(-si,0.0,-qqt,qqt))*AIEM::expc1(-qq,-qqt)
             +AIEM::Fahh(-si,0.0,-qq,qq)
             *conj(AIEM::Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc2(-qq,-qqts)
             +AIEM::Fahh(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fbhh(-si,0.0,-qqt,qqt))*AIEM::expc3(-qqs,-qqt)
             +AIEM::Fahh(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc4(-qqs,-qqts);

    expcbubu=AIEM::Fbhh(-si,0.0,qqt,qqt)
             *conj(AIEM::Fbhh(-si,0.0,qqt,qqt))*AIEM::expc1(qqt,qqt)
             +AIEM::Fbhh(-si,0.0,qqt,qqt)
             *conj(AIEM::Fbhh(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc2(qqt,qqts)
             +AIEM::Fbhh(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fbhh(-si,0.0,qqt,qqt))*AIEM::expc3(qqts,qqt)
             +AIEM::Fbhh(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fbhh(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc4(qqts,qqts);

    expcbdbu=AIEM::Fbhh(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fbhh(-si,0.0,qqt,qqt))*AIEM::expc1(-qqt,qqt)
             +AIEM::Fbhh(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fbhh(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc2(-qqt,qqts)
             +AIEM::Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fbhh(-si,0.0,qqt,qqt))*AIEM::expc3(-qqts,qqt)
             +AIEM::Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fbhh(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc4(-qqts,qqts);

    expcbubd=AIEM::Fbhh(-si,0.0,qqt,qqt)
             *conj(AIEM::Fbhh(-si,0.0,-qqt,qqt))*AIEM::expc1(qqt,-qqt)
             +AIEM::Fbhh(-si,0.0,qqt,qqt)
             *conj(AIEM::Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc2(qqt,-qqts)
             +AIEM::Fbhh(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fbhh(-si,0.0,-qqt,qqt))*AIEM::expc3(qqts,-qqt)
             +AIEM::Fbhh(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc4(qqts,-qqts);

    expcbdbd=AIEM::Fbhh(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fbhh(-si,0.0,-qqt,qqt))*AIEM::expc1(-qqt,-qqt)
             +AIEM::Fbhh(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc2(-qqt,-qqts)
             +AIEM::Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fbhh(-si,0.0,-qqt,qqt))*AIEM::expc3(-qqts,-qqt)
             +AIEM::Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc4(-qqts,-qqts);

    cterm[2]=0.03125*ex2*(expcauau+expcadau+expcauad+expcadad
                          +expcbuau+expcbdau+expcbuad+expcbdad
                          +expcaubu+expcadbu+expcaubd+expcadbd
                          +expcbubu+expcbdbu+expcbubd+expcbdbd);

    //cterm(3)
    expcauau=AIEM::Fahv(-si,0.0,qq,qq)
             *conj(AIEM::Fahv(-si,0.0,qq,qq))*AIEM::expc1(qq,qq)
             +AIEM::Fahv(-si,0.0,qq,qq)
             *conj(AIEM::Fahv(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc2(qq,qqs)
             +AIEM::Fahv(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fahv(-si,0.0,qq,qq))*AIEM::expc3(qqs,qq)
             +AIEM::Fahv(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fahv(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc4(qqs,qqs);

    expcadau=AIEM::Fahv(-si,0.0,-qq,qq)
             *conj(AIEM::Fahv(-si,0.0,qq,qq))*AIEM::expc1(-qq,qq)
             +AIEM::Fahv(-si,0.0,-qq,qq)
             *conj(AIEM::Fahv(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc2(-qq,qqs)
             +AIEM::Fahv(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fahv(-si,0.0,qq,qq))*AIEM::expc3(-qqs,qq)
             +AIEM::Fahv(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fahv(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc4(-qqs,qqs);

    expcauad=AIEM::Fahv(-si,0.0,qq,qq)
             *conj(AIEM::Fahv(-si,0.0,-qq,qq))*AIEM::expc1(qq,-qq)
             +AIEM::Fahv(-si,0.0,qq,qq)
             *conj(AIEM::Fahv(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc2(qq,-qqs)
             +AIEM::Fahv(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fahv(-si,0.0,-qq,qq))*AIEM::expc3(qqs,-qq)
             +AIEM::Fahv(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fahv(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc4(qqs,-qqs);

    expcadad=AIEM::Fahv(-si,0.0,-qq,qq)
             *conj(AIEM::Fahv(-si,0.0,-qq,qq))*AIEM::expc1(-qq,-qq)
             +AIEM::Fahv(-si,0.0,-qq,qq)
             *conj(AIEM::Fahv(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc2(-qq,-qqs)
             +AIEM::Fahv(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fahv(-si,0.0,-qq,qq))*AIEM::expc3(-qqs,-qq)
             +AIEM::Fahv(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fahv(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc4(-qqs,-qqs);

    expcbuau=AIEM::Fbhv(-si,0.0,qqt,qqt)
             *conj(AIEM::Fahv(-si,0.0,qq,qq))*AIEM::expc1(qqt,qq)
             +AIEM::Fbhv(-si,0.0,qqt,qqt)
             *conj(AIEM::Fahv(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc2(qqt,qqs)
             +AIEM::Fbhv(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fahv(-si,0.0,qq,qq))*AIEM::expc3(qqts,qq)
             +AIEM::Fbhv(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fahv(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc4(qqts,qqs);

    expcbdau=AIEM::Fbhv(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fahv(-si,0.0,qq,qq))*AIEM::expc1(-qqt,qq)
             +AIEM::Fbhv(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fahv(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc2(-qqt,qqs)
             +AIEM::Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fahv(-si,0.0,qq,qq))*AIEM::expc3(-qqts,qq)
             +AIEM::Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fahv(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc4(-qqts,qqs);

    expcbuad=AIEM::Fbhv(-si,0.0,qqt,qqt)
             *conj(AIEM::Fahv(-si,0.0,-qq,qq))*AIEM::expc1(qqt,-qq)
             +AIEM::Fbhv(-si,0.0,qqt,qqt)
             *conj(AIEM::Fahv(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc2(qqt,-qqs)
             +AIEM::Fbhv(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fahv(-si,0.0,-qq,qq))*AIEM::expc3(qqts,-qq)
             +AIEM::Fbhv(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fahv(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc4(qqts,-qqs);

    expcbdad=AIEM::Fbhv(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fahv(-si,0.0,-qq,qq))*AIEM::expc1(-qqt,-qq)
             +AIEM::Fbhv(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fahv(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc2(-qqt,-qqs)
             +AIEM::Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fahv(-si,0.0,-qq,qq))*AIEM::expc3(-qqts,-qq)
             +AIEM::Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fahv(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc4(-qqts,-qqs);

    expcaubu=AIEM::Fahv(-si,0.0,qq,qq)
             *conj(AIEM::Fbhv(-si,0.0,qqt,qqt))*AIEM::expc1(qq,qqt)
             +AIEM::Fahv(-si,0.0,qq,qq)
             *conj(AIEM::Fbhv(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc2(qq,qqts)
             +AIEM::Fahv(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fbhv(-si,0.0,qqt,qqt))*AIEM::expc3(qqs,qqt)
             +AIEM::Fahv(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fbhv(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc4(qqs,qqts);

    expcadbu=AIEM::Fahv(-si,0.0,-qq,qq)
             *conj(AIEM::Fbhv(-si,0.0,qqt,qqt))*AIEM::expc1(-qq,qqt)
             +AIEM::Fahv(-si,0.0,-qq,qq)
             *conj(AIEM::Fbhv(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc2(-qq,qqts)
             +AIEM::Fahv(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fbhv(-si,0.0,qqt,qqt))*AIEM::expc3(-qqs,qqt)
             +AIEM::Fahv(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fbhv(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc4(-qqs,qqts);

    expcaubd=AIEM::Fahv(-si,0.0,qq,qq)
             *conj(AIEM::Fbhv(-si,0.0,-qqt,qqt))*AIEM::expc1(qq,-qqt)
             +AIEM::Fahv(-si,0.0,qq,qq)
             *conj(AIEM::Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc2(qq,-qqts)
             +AIEM::Fahv(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fbhv(-si,0.0,-qqt,qqt))*AIEM::expc3(qqs,-qqt)
             +AIEM::Fahv(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc4(qqs,-qqts);

    expcadbd=AIEM::Fahv(-si,0.0,-qq,qq)
             *conj(AIEM::Fbhv(-si,0.0,-qqt,qqt))*AIEM::expc1(-qq,-qqt)
             +AIEM::Fahv(-si,0.0,-qq,qq)
             *conj(AIEM::Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc2(-qq,-qqts)
             +AIEM::Fahv(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fbhv(-si,0.0,-qqt,qqt))*AIEM::expc3(-qqs,-qqt)
             +AIEM::Fahv(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc4(-qqs,-qqts);

    expcbubu=AIEM::Fbhv(-si,0.0,qqt,qqt)
             *conj(AIEM::Fbhv(-si,0.0,qqt,qqt))*AIEM::expc1(qqt,qqt)
             +AIEM::Fbhv(-si,0.0,qqt,qqt)
             *conj(AIEM::Fbhv(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc2(qqt,qqts)
             +AIEM::Fbhv(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fbhv(-si,0.0,qqt,qqt))*AIEM::expc3(qqts,qqt)
             +AIEM::Fbhv(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fbhv(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc4(qqts,qqts);

    expcbdbu=AIEM::Fbhv(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fbhv(-si,0.0,qqt,qqt))*AIEM::expc1(-qqt,qqt)
             +AIEM::Fbhv(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fbhv(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc2(-qqt,qqts)
             +AIEM::Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fbhv(-si,0.0,qqt,qqt))*AIEM::expc3(-qqts,qqt)
             +AIEM::Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fbhv(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc4(-qqts,qqts);

    expcbubd=AIEM::Fbhv(-si,0.0,qqt,qqt)
             *conj(AIEM::Fbhv(-si,0.0,-qqt,qqt))*AIEM::expc1(qqt,-qqt)
             +AIEM::Fbhv(-si,0.0,qqt,qqt)
             *conj(AIEM::Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc2(qqt,-qqts)
             +AIEM::Fbhv(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fbhv(-si,0.0,-qqt,qqt))*AIEM::expc3(qqts,-qqt)
             +AIEM::Fbhv(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc4(qqts,-qqts);

    expcbdbd=AIEM::Fbhv(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fbhv(-si,0.0,-qqt,qqt))*AIEM::expc1(-qqt,-qqt)
             +AIEM::Fbhv(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc2(-qqt,-qqts)
             +AIEM::Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fbhv(-si,0.0,-qqt,qqt))*AIEM::expc3(-qqts,-qqt)
             +AIEM::Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc4(-qqts,-qqts);

    cterm[3]=0.03125*ex2*(expcauau+expcadau+expcauad+expcadad
                          +expcbuau+expcbdau+expcbuad+expcbdad
                          +expcaubu+expcadbu+expcaubd+expcadbd
                          +expcbubu+expcbdbu+expcbubd+expcbdbd);

    //cterm(4)
    expcauau=AIEM::Favh(-si,0.0,qq,qq)
             *conj(AIEM::Favh(-si,0.0,qq,qq))*AIEM::expc1(qq,qq)
             +AIEM::Favh(-si,0.0,qq,qq)
             *conj(AIEM::Favh(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc2(qq,qqs)
             +AIEM::Favh(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Favh(-si,0.0,qq,qq))*AIEM::expc3(qqs,qq)
             +AIEM::Favh(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Favh(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc4(qqs,qqs);

    expcadau=AIEM::Favh(-si,0.0,-qq,qq)
             *conj(AIEM::Favh(-si,0.0,qq,qq))*AIEM::expc1(-qq,qq)
             +AIEM::Favh(-si,0.0,-qq,qq)
             *conj(AIEM::Favh(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc2(-qq,qqs)
             +AIEM::Favh(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Favh(-si,0.0,qq,qq))*AIEM::expc3(-qqs,qq)
             +AIEM::Favh(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Favh(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc4(-qqs,qqs);

    expcauad=AIEM::Favh(-si,0.0,qq,qq)
             *conj(AIEM::Favh(-si,0.0,-qq,qq))*AIEM::expc1(qq,-qq)
             +AIEM::Favh(-si,0.0,qq,qq)
             *conj(AIEM::Favh(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc2(qq,-qqs)
             +AIEM::Favh(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Favh(-si,0.0,-qq,qq))*AIEM::expc3(qqs,-qq)
             +AIEM::Favh(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Favh(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc4(qqs,-qqs);

    expcadad=AIEM::Favh(-si,0.0,-qq,qq)
             *conj(AIEM::Favh(-si,0.0,-qq,qq))*AIEM::expc1(-qq,-qq)
             +AIEM::Favh(-si,0.0,-qq,qq)
             *conj(AIEM::Favh(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc2(-qq,-qqs)
             +AIEM::Favh(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Favh(-si,0.0,-qq,qq))*AIEM::expc3(-qqs,-qq)
             +AIEM::Favh(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Favh(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc4(-qqs,-qqs);

    expcbuau=AIEM::Fbvh(-si,0.0,qqt,qqt)
             *conj(AIEM::Favh(-si,0.0,qq,qq))*AIEM::expc1(qqt,qq)
             +AIEM::Fbvh(-si,0.0,qqt,qqt)
             *conj(AIEM::Favh(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc2(qqt,qqs)
             +AIEM::Fbvh(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Favh(-si,0.0,qq,qq))*AIEM::expc3(qqts,qq)
             +AIEM::Fbvh(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Favh(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc4(qqts,qqs);

    expcbdau=AIEM::Fbvh(-si,0.0,-qqt,qqt)
             *conj(AIEM::Favh(-si,0.0,qq,qq))*AIEM::expc1(-qqt,qq)
             +AIEM::Fbvh(-si,0.0,-qqt,qqt)
             *conj(AIEM::Favh(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc2(-qqt,qqs)
             +AIEM::Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Favh(-si,0.0,qq,qq))*AIEM::expc3(-qqts,qq)
             +AIEM::Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Favh(-sis*csfs,-sis*sfs,qqs,qqs))*AIEM::expc4(-qqts,qqs);

    expcbuad=AIEM::Fbvh(-si,0.0,qqt,qqt)
             *conj(AIEM::Favh(-si,0.0,-qq,qq))*AIEM::expc1(qqt,-qq)
             +AIEM::Fbvh(-si,0.0,qqt,qqt)
             *conj(AIEM::Favh(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc2(qqt,-qqs)
             +AIEM::Fbvh(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Favh(-si,0.0,-qq,qq))*AIEM::expc3(qqts,-qq)
             +AIEM::Fbvh(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Favh(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc4(qqts,-qqs);

    expcbdad=AIEM::Fbvh(-si,0.0,-qqt,qqt)
             *conj(AIEM::Favh(-si,0.0,-qq,qq))*AIEM::expc1(-qqt,-qq)
             +AIEM::Fbvh(-si,0.0,-qqt,qqt)
             *conj(AIEM::Favh(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc2(-qqt,-qqs)
             +AIEM::Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Favh(-si,0.0,-qq,qq))*AIEM::expc3(-qqts,-qq)
             +AIEM::Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Favh(-sis*csfs,-sis*sfs,-qqs,qqs))*AIEM::expc4(-qqts,-qqs);

    expcaubu=AIEM::Favh(-si,0.0,qq,qq)
             *conj(AIEM::Fbvh(-si,0.0,qqt,qqt))*AIEM::expc1(qq,qqt)
             +AIEM::Favh(-si,0.0,qq,qq)
             *conj(AIEM::Fbvh(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc2(qq,qqts)
             +AIEM::Favh(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fbvh(-si,0.0,qqt,qqt))*AIEM::expc3(qqs,qqt)
             +AIEM::Favh(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fbvh(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc4(qqs,qqts);

    expcadbu=AIEM::Favh(-si,0.0,-qq,qq)
             *conj(AIEM::Fbvh(-si,0.0,qqt,qqt))*AIEM::expc1(-qq,qqt)
             +AIEM::Favh(-si,0.0,-qq,qq)
             *conj(AIEM::Fbvh(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc2(-qq,qqts)
             +AIEM::Favh(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fbvh(-si,0.0,qqt,qqt))*AIEM::expc3(-qqs,qqt)
             +AIEM::Favh(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fbvh(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc4(-qqs,qqts);

    expcaubd=AIEM::Favh(-si,0.0,qq,qq)
             *conj(AIEM::Fbvh(-si,0.0,-qqt,qqt))*AIEM::expc1(qq,-qqt)
             +AIEM::Favh(-si,0.0,qq,qq)
             *conj(AIEM::Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc2(qq,-qqts)
             +AIEM::Favh(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fbvh(-si,0.0,-qqt,qqt))*AIEM::expc3(qqs,-qqt)
             +AIEM::Favh(-sis*csfs,-sis*sfs,qqs,qqs)
             *conj(AIEM::Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc4(qqs,-qqts);

    expcadbd=AIEM::Favh(-si,0.0,-qq,qq)
             *conj(AIEM::Fbvh(-si,0.0,-qqt,qqt))*AIEM::expc1(-qq,-qqt)
             +AIEM::Favh(-si,0.0,-qq,qq)
             *conj(AIEM::Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc2(-qq,-qqts)
             +AIEM::Favh(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fbvh(-si,0.0,-qqt,qqt))*AIEM::expc3(-qqs,-qqt)
             +AIEM::Favh(-sis*csfs,-sis*sfs,-qqs,qqs)
             *conj(AIEM::Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc4(-qqs,-qqts);

    expcbubu=AIEM::Fbvh(-si,0.0,qqt,qqt)
             *conj(AIEM::Fbvh(-si,0.0,qqt,qqt))*AIEM::expc1(qqt,qqt)
             +AIEM::Fbvh(-si,0.0,qqt,qqt)
             *conj(AIEM::Fbvh(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc2(qqt,qqts)
             +AIEM::Fbvh(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fbvh(-si,0.0,qqt,qqt))*AIEM::expc3(qqts,qqt)
             +AIEM::Fbvh(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fbvh(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc4(qqts,qqts);

    expcbdbu=AIEM::Fbvh(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fbvh(-si,0.0,qqt,qqt))*AIEM::expc1(-qqt,qqt)
             +AIEM::Fbvh(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fbvh(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc2(-qqt,qqts)
             +AIEM::Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fbvh(-si,0.0,qqt,qqt))*AIEM::expc3(-qqts,qqt)
             +AIEM::Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fbvh(-sis*csfs,-sis*sfs,qqts,qqts))*AIEM::expc4(-qqts,qqts);

    expcbubd=AIEM::Fbvh(-si,0.0,qqt,qqt)
             *conj(AIEM::Fbvh(-si,0.0,-qqt,qqt))*AIEM::expc1(qqt,-qqt)
             +AIEM::Fbvh(-si,0.0,qqt,qqt)
             *conj(AIEM::Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc2(qqt,-qqts)
             +AIEM::Fbvh(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fbvh(-si,0.0,-qqt,qqt))*AIEM::expc3(qqts,-qqt)
             +AIEM::Fbvh(-sis*csfs,-sis*sfs,qqts,qqts)
             *conj(AIEM::Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc4(qqts,-qqts);

    expcbdbd=AIEM::Fbvh(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fbvh(-si,0.0,-qqt,qqt))*AIEM::expc1(-qqt,-qqt)
             +AIEM::Fbvh(-si,0.0,-qqt,qqt)
             *conj(AIEM::Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc2(-qqt,-qqts)
             +AIEM::Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fbvh(-si,0.0,-qqt,qqt))*AIEM::expc3(-qqts,-qqt)
             +AIEM::Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts)
             *conj(AIEM::Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts))*AIEM::expc4(-qqts,-qqts);

    cterm[4]=0.03125*ex2*(expcauau+expcadau+expcauad+expcadad
                          +expcbuau+expcbdau+expcbuad+expcbdad
                          +expcaubu+expcadbu+expcaubd+expcadbd
                          +expcbubu+expcbdbu+expcbubd+expcbdbd);      //   1/32=0.03125

//---------------------------------------------//
//end of computation of complementary terms    //
//---------------------------------------------//

    sigma0[1]=(kterm[1]+kcterm[1]+cterm[1]).real();     //kirchhoff term + Cross term + Complementary terms
    sigma0[2]=(kterm[2]+kcterm[2]+cterm[2]).real();
    sigma0[3]=(kterm[3]+kcterm[3]+cterm[3]).real();
    sigma0[4]=(kterm[4]+kcterm[4]+cterm[4]).real();


//----------------------------------------------//
}//sigma

//----------------------------------------------//

complex<double> AIEM::expkc1(const complex<double> q)    //used in cross term
{

    complex<double> rexpkc1,temp,sum;
    int n;
    double fn;

    sum=0.0;
    temp=1.0;

    for(n=1; n<=iterm; n++)
    {
        fn=float(n);
        temp=temp*(ks2*(css-q)*(css+cs))/fn;
        sum=sum+temp*w[n];
    }
    rexpkc1=exp(-ks2*(q*q-q*css+q*cs))*sum;
    return rexpkc1;

}

//-----------------------------------------------------------------//


complex<double> AIEM::expkc2(const complex<double> q)       //used in cross term
{

    complex<double> rexpkc2,temp,sum;
    int n;
    double fn;

    sum=0.0;
    temp=1.0;

    for(n=1; n<=iterm; n++)
    {
        fn=float(n);
        temp=temp*(ks2*(cs+q)*(css+cs))/fn;
        sum=sum+temp*w[n];
    }
    rexpkc2=exp(-ks2*(q*q-q*css+q*cs))*sum;
    return rexpkc2;
}

//-----------------------------------------------------------------//

complex<double> AIEM::expc1(const complex<double> q,const complex<double> qp)
{

    complex<double> rexpc1,temp,sum;
    int n;
    double fn;

    sum=0.0;
    temp=1.0;

    for(n=1; n<=iterm; n++)
    {
        fn=float(n);
        temp=temp*(ks2*(css-q)*(css-qp))/fn;
        sum=sum+temp*w[n];
    }
    rexpc1=exp(-ks2*(q*q+qp*qp-css*(q+qp)+cs*(q+qp)))*sum;
    return rexpc1;
}

//-----------------------------------------------------------------//

complex<double> AIEM::expc2(const complex<double> q,const complex<double> qp)
{
    complex<double> rexpc2,temp,sum;
    int n;
    double fn;

    sum=complex<double>(0.0, 0.0);
    temp=complex<double>(1.0,0.0);


    for(n=1; n<=iterm; n++)
    {
        fn=float(n);
        temp=temp*(ks2*(css-q)*(cs+qp))/fn;
        sum=sum+temp*w[n];
    }

    rexpc2=exp(-ks2*(q*q+qp*qp-css*(q+qp)+cs*(q+qp)))*sum;
    return rexpc2;
}

//-----------------------------------------------------------------//

complex<double> AIEM::expc3(const complex<double> q,const complex<double> qp)
{
    complex<double> rexpc3,temp,sum;
    int n;
    double fn;

    sum=0.0;
    temp=1.0;

    for(n=1; n<=iterm; n++)
    {
        fn=float(n);
        temp=temp*(ks2*(cs+q)*(css-qp))/fn;
        sum=sum+temp*w[n];
    }
    rexpc3=exp(-ks2*(q*q+qp*qp-css*(q+qp)+cs*(q+qp)))*sum;
    return rexpc3;
}

//-----------------------------------------------------------------//

complex<double> AIEM::expc4(const complex<double> q,const complex<double> qp)
{

    complex<double> rexpc4,temp,sum;
    int n;
    double fn;

    sum=0.0;
    temp=1.0;

    for(n=1; n<=iterm; n++)
    {
        fn=float(n);
        temp=temp*(ks2*(cs+q)*(cs+qp))/fn;
        sum=sum+temp*w[n];
    }
    rexpc4=exp(-ks2*(q*q+qp*qp-css*(q+qp)+cs*(q+qp)))*sum;
    return rexpc4;
}

//------------------------------------------------------------------//

complex<double> AIEM::Favv(const double u,const double v,const complex<double> q,const complex<double> qfix)
{

    //reference: A5
    complex<double> rFavv;
    complex<double> c1,c2,c3,c4,c5,c6;
    complex<double> zx,zy,zxp,zyp;
    complex<double> rpv,rmv,av,bv;
    double kxu,ksxu,kyv,ksyv;

    kxu=si+u;
    ksxu=sis*csfs+u;
    kyv=v;
    ksyv=sis*sfs+v;

    if((fabs((css-q).real())<0.000001)||(fabs((cs+q).real())<0.000001))
    {
        c1=0.0;
        c2=0.0;
        c3=0.0;
        c4=0.0;
        c5=0.0;
        c6=0.0;
    }
    else
    {
        zx=(-ksxu)/(css-q);
        zy=(-ksyv)/(css-q);
        zxp=(kxu)/(cs+q);
        zyp=(kyv)/(cs+q);
        c1=-csfs*(-1.0-zx*zxp)+sfs*zxp*zy;
        c2=-csfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp-si*v*zx*zyp)
           +sfs*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp);
        c3=-csfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)
           +sfs*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy);
        c4=-css*sfs*(-si*zyp+cs*zx*zyp)-csfs*css*(-cs-si*zxp-cs*zy*zyp)
           +sis*(-cs*zx-si*zx*zxp-si*zy*zyp);
        c5=-css*sfs*(-v*zx+v*zxp)-csfs*css*(q+u*zxp+v*zy)
           +sis*(q*zx+u*zx*zxp+v*zxp*zy);
        c6=-css*sfs*(-u*zyp+q*zx*zyp)-csfs*css*(v*zyp-q*zy*zyp)
           +sis*(v*zx*zyp-u*zy*zyp);
    }

    rpv=1.0+rv;
    rmv=1.0-rv;
    av=rpv/qfix;
    bv=rmv/qfix;
    rFavv=bv*(-rpv*c1+rmv*c2+rpv*c3)+av*(rmv*c4+rpv*c5+rmv*c6);
    return rFavv;
}

//------------------------------------------------------------------//

complex<double> AIEM::Fahh(const double u,const double v,const complex<double> q,const complex<double> qfix)
{

    //reference:A7
    complex<double> rFahh;
    complex<double> c1,c2,c3,c4,c5,c6;
    complex<double> zx,zy,zxp,zyp;
    complex<double> rph,rmh,ah,bh;
    double kxu,ksxu,kyv,ksyv;

    kxu=si+u;
    ksxu=sis*csfs+u;
    kyv=v;
    ksyv=sis*sfs+v;

    if((fabs((css-q).real())<0.000001)||(fabs((cs+q).real())<0.000001))
    {
        c1=0.0;
        c2=0.0;
        c3=0.0;
        c4=0.0;
        c5=0.0;
        c6=0.0;
    }
    else
    {
        zx=(-ksxu)/(css-q);
        zy=(-ksyv)/(css-q);
        zxp=(kxu)/(cs+q);
        zyp=(kyv)/(cs+q);
        c1=-csfs*(-1.0-zx*zxp)+sfs*zxp*zy;
        c2=-csfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp-si*v*zx*zyp)
           +sfs*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp);
        c3=-csfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)
           +sfs*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy);
        c4=-css*sfs*(-si*zyp+cs*zx*zyp)-csfs*css*(-cs-si*zxp-cs*zy*zyp)
           +sis*(-cs*zx-si*zx*zxp-si*zy*zyp);
        c5=-css*sfs*(-v*zx+v*zxp)-csfs*css*(q+u*zxp+v*zy)
           +sis*(q*zx+u*zx*zxp+v*zxp*zy);
        c6=-css*sfs*(-u*zyp+q*zx*zyp)-csfs*css*(v*zyp-q*zy*zyp)
           +sis*(v*zx*zyp-u*zy*zyp);
    }


    rph=1.0+rh;
    rmh=1.0-rh;
    ah=rph/qfix;
    bh=rmh/qfix;
    rFahh=-bh*(-rph*c1+rmh*c2+rph*c3)-ah*(rmh*c4+rph*c5+rmh*c6);
    return rFahh;

}


//------------------------------------------------------------------//

complex<double> AIEM::Fbvv(const double u,const double v,const complex<double> q,const complex<double> qfix)
{

    //reference A6
    complex<double> rFbvv;
    complex<double> c1,c2,c3,c4,c5,c6;
    complex<double> zx,zy,zxp,zyp;
    complex<double> rpv,rmv,av,bv;
    double kxu,ksxu,kyv,ksyv;

    kxu=si+u;
    ksxu=sis*csfs+u;
    kyv=v;
    ksyv=sis*sfs+v;
    zx=(-ksxu)/(css-q);
    zy=(-ksyv)/(css-q);
    zxp=(kxu)/(cs+q);
    zyp=(kyv)/(cs+q);

    c1=-csfs*(-1.0-zx*zxp)+sfs*zxp*zy;
    c2=-csfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp-si*v*zx*zyp)
       +sfs*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp);
    c3=-csfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)
       +sfs*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy);
    c4=-css*sfs*(-si*zyp+cs*zx*zyp)-csfs*css*(-cs-si*zxp-cs*zy*zyp)
       +sis*(-cs*zx-si*zx*zxp-si*zy*zyp);
    c5=-css*sfs*(-v*zx+v*zxp)-csfs*css*(q+u*zxp+v*zy)
       +sis*(q*zx+u*zx*zxp+v*zxp*zy);
    c6=-css*sfs*(-u*zyp+q*zx*zyp)-csfs*css*(v*zyp-q*zy*zyp)
       +sis*(v*zx*zyp-u*zy*zyp);
    rpv=1.0+rv;
    rmv=1.0-rv;
    av=rpv/qfix;
    bv=rmv/qfix;
    rFbvv=av*(rpv*c1-rmv*c2-rpv*c3/er)-bv*(rmv*c4*er+rpv*c5+rmv*c6);
    return rFbvv;
}


//------------------------------------------------------------------//

complex<double> AIEM::Fbhh(const double u,const double v,const complex<double> q,const complex<double> qfix)
{

    //reference: A8
    complex<double> rFbhh;
    complex<double> c1,c2,c3,c4,c5,c6;
    complex<double> zx,zy,zxp,zyp;
    complex<double> rph,rmh,ah,bh;
    double kxu,ksxu,kyv,ksyv;

    kxu=si+u;
    ksxu=sis*csfs+u;
    kyv=v;
    ksyv=sis*sfs+v;
    zx=(-ksxu)/(css-q);
    zy=(-ksyv)/(css-q);
    zxp=(kxu)/(cs+q);
    zyp=(kyv)/(cs+q);
    c1=-csfs*(-1.0-zx*zxp)+sfs*zxp*zy;
    c2=-csfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp-si*v*zx*zyp)
       +sfs*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp);
    c3=-csfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)
       +sfs*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy);
    c4=-css*sfs*(-si*zyp+cs*zx*zyp)-csfs*css*(-cs-si*zxp-cs*zy*zyp)
       +sis*(-cs*zx-si*zx*zxp-si*zy*zyp);
    c5=-css*sfs*(-v*zx+v*zxp)-csfs*css*(q+u*zxp+v*zy)
       +sis*(q*zx+u*zx*zxp+v*zxp*zy);
    c6=-css*sfs*(-u*zyp+q*zx*zyp)-csfs*css*(v*zyp-q*zy*zyp)
       +sis*(v*zx*zyp-u*zy*zyp);
    rph=1.0+rh;
    rmh=1.0-rh;
    ah=rph/qfix;
    bh=rmh/qfix;
    rFbhh=ah*(-rph*c1*er+rmh*c2+rph*c3)+bh*(rmh*c4+rph*c5+rmh*c6/er);
    return rFbhh;
}

//------------------------------------------------------------------//

complex<double> AIEM::Fahv(const double u,const double v,const complex<double> q,const complex<double> qfix)
{

    complex<double> rFahv;
    complex<double> b1,b2,b3,b4,b5,b6;
    complex<double> zx,zy,zxp,zyp;
    complex<double> rp,rm,a,b;
    double kxu,ksxu,kyv,ksyv;

    kxu=si+u;
    ksxu=sis*csfs+u;
    kyv=v;
    ksyv=sis*sfs+v;


    if((fabs((css-q).real())<0.000001)||(fabs((cs+q).real())<0.000001))
    {
        b1=0.0;
        b2=0.0;
        b3=0.0;
        b4=0.0;
        b5=0.0;
        b6=0.0;
    }
    else
    {
        zx=(-ksxu)/(css-q);
        zy=(-ksyv)/(css-q);
        zxp=(kxu)/(cs+q);
        zyp=(kyv)/(cs+q);
        b1=-css*sfs*(-1.0-zx*zxp)-sis*zy-csfs*css*zxp*zy;
        b2=-css*sfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp-si*v*zx*zyp)+
           sis*(-cs*q*zy-q*si*zxp*zy+q*si*zx*zyp-cs*u*zx*zyp-cs*v*zy*zyp)-
           csfs*css*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp);
        b3=-css*sfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)-
           csfs*css*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)+
           sis*(-si*v*zx+cs*v*zx*zxp+si*u*zy-cs*u*zxp*zy);
        b4=-csfs*(-si*zyp+cs*zx*zyp)+sfs*(-cs-si*zxp-cs*zy*zyp);
        b5=-csfs*(-v*zx+v*zxp)+sfs*(q+u*zxp+v*zy);
        b6=-csfs*(-u*zyp+q*zx*zyp)+sfs*(v*zyp-q*zy*zyp);
    }

    rp=1.0+rvh;
    rm=1.0-rvh;
    a=rp/qfix;
    b=rm/qfix;
    rFahv=b*(rp*b1-rm*b2-rp*b3)+a*(rm*b4+rp*b5+rm*b6);
    return rFahv;
}

//------------------------------------------------------------------//

complex<double> AIEM::Favh(const double u,const double v,const complex<double> q,const complex<double> qfix)
{

    complex<double> rFavh;
    complex<double> b1,b2,b3,b4,b5,b6;
    complex<double> zx,zy,zxp,zyp;
    complex<double> rp,rm,a,b;

    double kxu,ksxu,kyv,ksyv;

    kxu=si+u;
    ksxu=sis*csfs+u;
    kyv=v;
    ksyv=sis*sfs+v;

    if((fabs((css-q).real())<0.000001)||(fabs((cs+q).real())<0.000001))
    {
        b1=0.0;
        b2=0.0;
        b3=0.0;
        b4=0.0;
        b5=0.0;
        b6=0.0;
    }
    else
    {
        zx=(-ksxu)/(css-q);
        zy=(-ksyv)/(css-q);
        zxp=(kxu)/(cs+q);
        zyp=(kyv)/(cs+q);
        b1=-css*sfs*(-1.0-zx*zxp)-sis*zy-csfs*css*zxp*zy;
        b2=-css*sfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp
                     -si*v*zx*zyp)+sis*(-cs*q*zy-q*si*zxp*zy+q*si*zx*zyp-cs*u*zx*zyp-cs*v*zy*zyp)
           -csfs*css*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp);
        b3=-css*sfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)-
           csfs*css*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)+
           sis*(-si*v*zx+cs*v*zx*zxp+si*u*zy-cs*u*zxp*zy);
        b4=-csfs*(-si*zyp+cs*zx*zyp)+sfs*(-cs-si*zxp-cs*zy*zyp);
        b5=-csfs*(-v*zx+v*zxp)+sfs*(q+u*zxp+v*zy);
        b6=-csfs*(-u*zyp+q*zx*zyp)+sfs*(v*zyp-q*zy*zyp);
    }

    rp=1.0+rvh;
    rm=1.0-rvh;
    a=rp/qfix;
    b=rm/qfix;
    rFavh=b*(rp*b4+rm*b5+rp*b6)-a*(-rm*b1+rp*b2+rm*b3);
    return rFavh;
}


//------------------------------------------------------------------//

complex<double> AIEM::Fbhv(const double u,const double v,const complex<double> q,const complex<double> qfix)
{

    complex<double> rFbhv;
    complex<double> b1,b2,b3,b4,b5,b6;
    complex<double> zx,zy,zxp,zyp;
    complex<double> rp,rm,a,b;
    double kxu,ksxu,kyv,ksyv;

    kxu=si+u;
    ksxu=sis*csfs+u;
    kyv=v;
    ksyv=sis*sfs+v;
    zx=(-ksxu)/(css-q);
    zy=(-ksyv)/(css-q);
    zxp=(kxu)/(cs+q);
    zyp=(kyv)/(cs+q);
    b1=-css*sfs*(-1.0-zx*zxp)-sis*zy-csfs*css*zxp*zy;
    b2=-css*sfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp
                 -si*v*zx*zyp)+sis*(-cs*q*zy-q*si*zxp*zy+q*si*zx*zyp-cs*u*zx*zyp-cs*v*zy*zyp)
       -csfs*css*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp);
    b3=-css*sfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)-
       csfs*css*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)+
       sis*(-si*v*zx+cs*v*zx*zxp+si*u*zy-cs*u*zxp*zy);
    b4=-csfs*(-si*zyp+cs*zx*zyp)+sfs*(-cs-si*zxp-cs*zy*zyp);
    b5=-csfs*(-v*zx+v*zxp)+sfs*(q+u*zxp+v*zy);
    b6=-csfs*(-u*zyp+q*zx*zyp)+sfs*(v*zyp-q*zy*zyp);
    rp=1.0+rvh;
    rm=1.0-rvh;
    a=rp/qfix;
    b=rm/qfix;
    rFbhv=a*(-rp*b1+rm*b2+rp*b3/er)-b*(rm*b4*er+rp*b5+rm*b6);
    return rFbhv;
}

//------------------------------------------------------------------//
complex<double> AIEM::Fbvh(const double u,const double v,const complex<double> q,const complex<double> qfix)
{
    complex<double> rFbvh;
    complex<double> b1,b2,b3,b4,b5,b6;
    complex<double> zx,zy,zxp,zyp;
    complex<double> rp,rm,a,b;

    double kxu,ksxu,kyv,ksyv;

    kxu=si+u;
    ksxu=sis*csfs+u;
    kyv=v;
    ksyv=sis*sfs+v;
    zx=(-ksxu)/(css-q);
    zy=(-ksyv)/(css-q);
    zxp=(kxu)/(cs+q);
    zyp=(kyv)/(cs+q);
    b1=-css*sfs*(-1.0-zx*zxp)-sis*zy-csfs*css*zxp*zy;
    b2=-css*sfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp-si*v*zx*zyp)
       +sis*(-cs*q*zy-q*si*zxp*zy+q*si*zx*zyp-cs*u*zx*zyp-cs*v*zy*zyp)
       -csfs*css*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp);
    b3=-css*sfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)
       -csfs*css*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)
       +sis*(-si*v*zx+cs*v*zx*zxp+si*u*zy-cs*u*zxp*zy);
    b4=-csfs*(-si*zyp+cs*zx*zyp)+sfs*(-cs-si*zxp-cs*zy*zyp);
    b5=-csfs*(-v*zx+v*zxp)+sfs*(q+u*zxp+v*zy);
    b6=-csfs*(-u*zyp+q*zx*zyp)+sfs*(v*zyp-q*zy*zyp);
    rp=1.0+rvh;
    rm=1.0-rvh;
    a=rp/qfix;
    b=rm/qfix;
    rFbvh=-a*(rp*b4+rm*b5+rp*b6/er)+b*(-rm*b1*er+rp*b2+rm*b3);
    return rFbvh;
}

//------------------------------------------------------------------//

double AIEM::bessj0(const double x)
{
    double y,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6,rbessj0;
    double ax,z,xx;
    r1=57568490574.0;
    r2=-13362590354.0;
    r3=651619640.70;
    r4=-11214424.180;
    r5=77392.330170;
    r6=-184.90524560;
    s1=57568490411.0;
    s2=1029532985.0;
    s3=9494680.7180;
    s4=59272.648530;
    s5=267.85327120;
    s6=1.0;
    p1=1.0;
    p2=-0.1098628627e-2;
    p3=0.2734510407e-4;
    p4=-0.2073370639e-5;
    p5=0.2093887211e-6;
    q1=-0.1562499995e-1;
    q2=0.1430488765e-3;
    q3=-0.6911147651e-5;
    q4=0.7621095161e-6;
    q5=-0.934945152e-7;

    if(fabs(x)<8.0)
    {
        y=x*x;
        rbessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))));
    }
    else
    {
        ax=fabs(x);
        z=8.0/ax;
        y=z*z;
        xx=ax-0.785398164;
        rbessj0=sqrt(0.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))));
    }
    return rbessj0;
}



//	    void AIEM::shadowg(const double ti,const double ts,const double s,const double shfct)   //没有调用??
//		{

//      if(ts>=ti) then
//        arg=ts
//        else
//        arg=ti
//       endif

//		    double arg,u,pi,et,f1,f2,f;

//			arg=ti;

//			if(arg==0.0)
//				shfct=1.0;
//			u=1.0/tan(arg);
//			pi=acos(-1.0);
//			et=u/(sqrt(2.0)*s);

//			if(et>=20)
//				shfct=1.0;
//		    f1=sqrt(2.0/pi)*s*exp(-et*et)/u;
//			f2=erfc(et);
//          f=(f1-f2)/2.0;
//			shfct=1.0/(1.0+f);
//      shfct=(1.0-0.5*erfc(et))/(1.0+f)
//       return
//	  }



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//   subroutine calculates shadowing function   //
//        (Tsang et al 1985, pp.95)             //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

//		void AIEM::shadow(const bool back,const double ti,const double ts, double shfct)
//		{
//			double s,ui,us,pi,eti,ets,f1i,f1s,f2i,f2s,fi,fs;
//			s=effslop;

//			if(ti==0.0&&ts==0.0)
//				shfct=1.0;

//			if(ti==0.0)
//				ui=1.0+30;
//			else
//				ui=1.0/tan(ti);

//			if(ts==0.0)
//				us=1.0+30;
//			else
//				us=1.0/tan(ts);

//			pi=acos(-1.0);
//			eti=ui/(sqrt(2.0)*s);
//			ets=us/(sqrt(2.0)*s);

//			if(!back)
//			{
//				f1i=sqrt(2.0/pi)*s/ui*exp(-eti*eti);
//				f1s=sqrt(2.0/pi)*s/us*exp(-ets*ets);
//				f2i=erfc(eti);
//				f2s=erfc(ets);
//				fi=(f1i-f2i)/2.0;
//				fs=(f1s-f2s)/2.0;
//				shfct=1.0/(1.0+fi+fs);
//			}
//			else
//				if(ts>=ti)
//				{
//					f1s=sqrt(2.0/pi)*s/us*exp(-ets*ets);
//					f2s=erfc(ets);
//					fs=(f1s-f2s)/2.0;
//					shfct=(1.0-0.5*f2s)/(1.0+fs);
//				}
//				else
//				{
//					f1i=sqrt(2.0/pi)*s/ui*exp(-eti*eti);
//	                f2i=erfc(eti);
//					fi=(f1i-f2i)/2.0;
//					shfct=(1.0-0.5*f2i)/(fi+1.0);
//				}

//			}



//		double AIEM::erfcc(const double x)
//		{
//			double rerfcc,t,z;
//			z=fabs(x);
//			t=1.0/(1.0+0.5*z);
//			rerfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+
//					t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*
//					(1.48851587+t*(-.82215223+t*0.17087277)))))))));
//
//			if (x<0.0)
//				rerfcc=2.0-rerfcc;
//			return rerfcc;
//		}
//
//		double AIEM::erfc(const double x)
//		{
//			double rerfc;
//			if(x<0.0)
//				rerfc=1.0+gammp(0.5,x*x);
//			else
//				rerfc=gammq(0.5,x*x);
//			return rerfc;
//		}

//		double AIEM::gammp(const double a,const double x)
//		{
//			//if(x<0.0||a<=0.0) pause;    //???????

//			double rgammp;
//
//			if(x<(a+1.0))
//			{
//				gser(gamser,a,x,gln);
//				rgammp=gamser;
//			}
//			else
//			{
//				gcf(gammcf,a,x,gln);
//				rgammp=1.0-gammcf;
//			}
//			return rgammp;
//		}


//		double AIEM::gammq(const double a,const double x)
//		{

//if(x<0.0||a<=0.0) pause;              //??????

//			double rgammq;

//			if(x<(a+1.0))
//			{
//				gser(gamser,a,x,gln);
//				rgammq=1.0-gamser;
//			}
//			else
//			{
//				gcf(gammcf,a,x,gln);
//				rgammq=gammcf;
//			 }
//			return rgammq;

//		}


//		void AIEM::gser(double gamser,const double a,const double x,double gln)
//		{
//			const int itmax=100;
//			const double eps=1.0e-18;
//			double ap,sum,del;
//			int n;

//			gln=gammln[a];

//			if(x<=0.0)
//			{
//if(x<0.00);
//				gamser=0.0;
//			}

//			ap=a;
//			sum=1.0/a;
//			del=sum;

//			for(n=1;n<=itmax;n++)
//			{
//				ap=ap+1.0;
//				del=del*x/ap;
//				sum=sum+del;
//				if(fabs(del)<fabs(sum)*eps)
//					break;
//			}

//			if(n==itmax+1)
//				cout<<"a too large, itmax too small"<<endl;

//			gamser=sum*exp(-x+a*log(x)-gln);
//          }


//	    void AIEM::gcf(double gammcf,const int a,const double x,double gln)
//		{

//			const int itmax=100;
//			const eps=1.0e-18;
//			double gold,a0,a1,b0,b1,fac,an,ana,anf,g;
//			int n;


//			gln=gammln[a];
//			gold=0.0;
//			a0=1.0;
//			a1=x;
//			b0=0.0;
//			b1=1.0;
//			fac=1.0;

//			for(n=1;n<=itmax;n++)
//			{
//				an=float(n);
//				ana=an-a;
//				a0=(a1+a0*ana)*fac;
//				b0=(b1+b0*ana)*fac;
//				anf=an*fac;
//				a1=x*a0+anf*a1;
//				b1=x*b0+anf*b1;

//				if(a1!=0.0)
//				{
//					fac=1.0/a1;
//					g=b1*fac;
//					if(fabs((g-gold)/g)<eps)
//						break;
//					else
//						gold=g;
//				}
//			}
//				if(n==itmax+1)
//					cout<<"a too large, itmax too small"<<endl;
//				gammcf=exp(-x+a*log(x)-gln)*g;

//		}


//	    double AIEM::gammln(const double xx)
//		{
//			double x,tmp,ser,rgammln;
//			double cof[6]={76.18009173,-86.50532033,24.01409822,-1.231739516,0.120858003e-2,-0.536382e-5};
//			double stp=2.50662827465;
//			double half=0.5,one=1.0,fpf=5.5;
//			x=xx-one;
//			tmp=x+fpf;
//			tmp=(x+half)*log(tmp)-tmp;
//			ser=one;
//			for(int j=1;j<=6;j++)
//			{
//				x=x+one;
//				ser=ser+cof[j]/x;
//			}
//			rgammln=tmp+log(stp*ser);
//			return rgammln;                  //gammln返回的仅是一个double数值，但是？？
//		}




double AIEM::BesselK(const int n,const double x)
{
//**************************************************
// Modified Bessel function of order n+0.5
// Input parameters:
//   n : int part of order
//   x : real parameter
//**************************************************

    double cons,K0,K1,K2,fn0,fn1,fn2,rBesselK;
    double pi=acos(-1.0);

    if(x==0.0)
        cout<<"BesselK: Singularity encountered !"<<endl;

    cons = sqrt(pi/(2.0*x));
    K0 = cons*exp(-x);
    if(n==0)
    {
        rBesselK = K0;
        return rBesselK;
    }

    K1 = K0*(1.0+1.0/x);

    if(n==1)
    {
        rBesselK = K1;
        return rBesselK;
    }

    K2 = K0*(1.0+3.0/x+3.0/x/x);

    if(n==2)
    {
        rBesselK = K2;
        return rBesselK;
    }
    fn0 = K1*cons;
    fn1 = -K2*cons;

    for(int i=2; i<n; i++)
    {
        fn2 = fn0-float(2*i+1)/x*fn1;
        fn0 = fn1;
        fn1 = fn2;
    }
    rBesselK = fabs(fn2/cons);
    return rBesselK;
}


double AIEM::alogam(const double x)                //????
{
    //evaluates natural logarithm of GAMMA(x)
    //for x > 0

    double a1,a2,a3,a4,a5,y,f,z,ralogam;
    int ifault;
    double pi=acos(-1.0);


    a1=log(2*pi)/2.0;
    a2=1.0/1680.0;
    a3=1.0/1260.0;
    a4=1.0/360.0;
    a5=1.0/12.0;
    ralogam=0.0;
    ifault=1;

    if(x<=0.0)
        return ralogam;
    else
    {

        ifault=0;
        y=x;
        f=0.0;

        if(y>=7.0)
        {
            z=1.0/(y*y);
            ralogam=f+(y-0.5)*log(y)-y+a1+((((-a2)*z+a3)*z-a4)*z+a5)/y;
        }
        else
        {
            f = y;
            y = y + 1.0;
            while(y<7.0)
            {
                f=f*y;
                y=y+1.0;
            }
            f = -log(f);
            z = 1.0 / ( y * y);
            ralogam = f+(y-0.5)*log(y)-y+a1+((((-a2)*z+a3)*z-a4)*z+a5)/y;

        }
    }
    return ralogam;
}



double AIEM::BESSK(const double N, const double X)
{
    double TOX,BKM,BK,BKP,rBESSK;
    if(N>=2)
    {
        TOX=2.0/X;
        BKM=BESSK0(X);
        BK=BESSK1(X);

        for(int J=1; J<N; J++)
        {
            BKP=BKM+J*TOX*BK;
            BKM=BK;
            BK=BKP;
        }
        rBESSK=BK;
        return rBESSK;
    }
    else
    {
        if(N==0) rBESSK = AIEM::BESSK0(X);
        if(N==1) rBESSK = AIEM::BESSK1(X);
        return rBESSK;
    }

}


double AIEM::BESSK0(const double X)
{
    double P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Y,rBESSK0;

    P1=-0.57721566;
    P2=0.42278420;
    P3=0.23069756;
    P4=0.3488590e-1;
    P5=0.262698e-2;
    P6=0.10750e-3;
    P7=0.74e-5;
    Q1=1.25331414;
    Q2=-0.7832358e-1;
    Q3=0.2189568e-1;
    Q4=-0.1062446e-1;
    Q5=0.587872e-2;
    Q6=-0.251540e-2;
    Q7=0.53208e-3;

    if (X<=2.0)
    {
        Y=X*X/4.0;
        rBESSK0=(-log(X/2.00)*BESSI0(X))+(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))));
    }
    else
    {
        Y=(2.0/X);
        rBESSK0=(exp(-X)/sqrt(X))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))));
    }
    return rBESSK0;
}


double AIEM::BESSK1(const double X)
{
    double P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Y,rBESSK1;
    P1=1.0;
    P2=0.15443144;
    P3=-0.67278579;
    P4=-0.18156897;
    P5=-0.1919402e-1;
    P6=-0.110404e-2;
    P7=-0.4686e-4;
    Q1=1.25331414;
    Q2=0.23498619;
    Q3=-0.3655620e-1;
    Q4=0.1504268e-1;
    Q5=-0.780353e-2;
    Q6=0.325614e-2;
    Q7=-0.68245e-3;

    if (X<=2.0)
    {
        Y=X*X/4.0;
        rBESSK1=(log(X/2.0)*BESSI1(X))+(1.0/X)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))));
    }
    else
    {
        Y=2.0/X;
        rBESSK1=(exp(-X)/sqrt(X))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))));
    }
    return rBESSK1;
}



double AIEM::BESSI0(const double X)
{
    double P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Y,rBESSI0,AX;
    P1=1.0;
    P2=3.5156229;
    P3=3.0899424;
    P4=1.2067492;
    P5=0.2659732;
    P6=0.360768e-1;
    P7=0.45813e-2;
    Q1=0.39894228;
    Q2=0.1328592e-1;
    Q3=0.225319e-2;
    Q4=-0.157565e-2;
    Q5=0.916281e-2;
    Q6=-0.2057706e-1;
    Q7=0.2635537e-1;
    Q8=-0.1647633e-1;
    Q9=0.392377e-2;

    if (fabs(X)<3.75)
    {
        Y=pow((X/3.75),2);
        rBESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))));
    }
    else
    {
        AX=fabs(X);
        Y=3.75/AX;
        rBESSI0=(exp(AX)/sqrt(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))));
    }
    return rBESSI0;
}



double AIEM::BESSI1(const double X)
{
    double P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Y,AX,rBESSI1;
    P1=0.5;
    P2=0.87890594;
    P3=0.51498869;
    P4=0.15084934;
    P5=0.02658733;
    P6=0.301532e-2;
    P7=0.32411e-3;
    Q1=0.39894228;
    Q2=-0.3988024e-1;
    Q3=-0.362018e-2;
    Q4=0.163801e-2;
    Q5=-0.1031555e-1;
    Q6=0.2282967e-1;
    Q7=-0.2895312e-1;
    Q8=0.1787654e-1;
    Q9=-0.420059e-2;

    if (fabs(X)<3.75)
    {
        Y=pow((X/3.75),2);
        rBESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))));
    }
    else
    {
        AX=fabs(X);
        Y=3.75/AX;
        rBESSI1=(exp(AX)/sqrt(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))));
    }
    return rBESSI1;
}



void AIEM::quagen(double z[],double wt[],const int n)  //press et. al: P112   Gauss-Legendre
{
    //quagen

    double p[513],c1[513],c2[513];
    double pi4,pdir,xnow,const1,epsilon,xinc,n2;
    int nmax,ip,n1,ncal,k,ncount,i,j;
    bool found;
    double pi=acos(-1.0);

    p[1]=1.0;
    epsilon=1.0e-10;
    nmax=10;
    ip=512;

    if(n<=0||n>ip)
        cout<<"input order out of range"<<endl;

    // calculate the coefficients for the legendre poly. recursive formula

    for(i=2; i<=n; i++)
    {
        c1[i]=(2*i-1)/float(i);    //1.5 1.667 1.75...
        c2[i]=(i-1)/float(i);      //0.5 0.667 0.75...
    }


    // initial constants

    n1=n+1;
    pi4=pi/4;
    const1=1.0/(n+0.5);

    // determine the number of roots(nr) needed to be calculated

    n2=n/2;
    if(n2*2==n)
        ncal=n2;
    else
        ncal=n2+1;

    // main loop begins here
    for(i=1; i<=ncal; i++)
    {
        //f1
        k=n-i+1;
        ncount=0;
        xinc=1.0;

        // use newton's method and a good initial guess to locate the root
        xnow=cos((i*pi-pi4)*const1);
        found=false;

        while(!found)
        {
            //while
            ncount=ncount+1;
            p[2]=xnow;

            for(j=2; j<=n; j++) //the following loop calculate p_n(x) using recursive formula
                p[j+1]=c1[j]*xnow*p[j]-c2[j]*p[j-1];

            // the derivate of p_n(x) can be calculated from p_n(x) and p_n-1(x)

            pdir=n*(p[n]-xnow*p[n1])/(1.0-xnow*xnow);

            if(fabs(xinc)<=epsilon||ncount>nmax)
            {
                found=true;
                z[k]=xnow;
                z[i]=-xnow;
                wt[k]=2.0/(1.0-z[k]*z[k])/(pdir*pdir);
                wt[i]=wt[k];
            }

            xinc=-p[n1]/pdir;
            xnow=xnow+xinc;
        }//while
    }//f1
}//quegen
