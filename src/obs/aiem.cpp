#include "aiem.h"

using namespace ldas;
using namespace std;

AIEM::AIEM(const AIEM& other)
{
    //copy ctor
    soil=other.soilparameter();
    m_correlation=other.corelation();
    fresnel_type=other.fresnel();
    m_mode=other.mode();
    npp=128;
    nss=128;
}

AIEM& AIEM::operator=(const AIEM& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    soil=rhs.soilparameter();
    m_correlation=rhs.corelation();
    fresnel_type=rhs.fresnel();
    m_mode=rhs.mode();
    npp=128;
    nss=128;
    return *this;
}

AIEM::AIEM(const SoilParameter& sp):soil(sp),m_mode(AIEM_PASSIVE),m_correlation(AIEM_GAUSSIAN),fresnel_type(AIEM_INCIDENT_ANGLE),npp(128),nss(128)
{

}

void AIEM::gauss_integral(const int m_npp, const int m_nss)
{
    npp					= m_npp;
    nss					= m_nss;
}

void AIEM::run()
{
    switch(m_mode)
    {
    case AIEM_ACTIVE:
        active();
        break;
    case AIEM_PASSIVE:
        passive();
        break;
    default:
        break;
    }
}
void AIEM::passive()
{
    // local variables
    double sigma0[4];
    double *ww11, *zz11, *ww22, *zz22;
    double *wt, *zt;
    double theta=soil.theta*PI_NUMBER/180;
    complex<double> ur = complex<double>(1.0, 0.0);

    ww11	= new double[npp];
    zz11	= new double[npp];
    ww22	= new double[npp];
    zz22	= new double[npp];
    wt		= new double[npp>nss?npp:nss];
    zt		= new double[npp>nss?npp:nss];

    double k	= 2.0*PI_NUMBER*soil.frequency*1.0e9/LIGHT_SPEED;	// wavenumber (m^-1)
    double kl	= soil.cl *0.01*k;					// cl *0.01 (cm->m)
    double ks	= soil.rms*0.01*k;					// rms*0.01 (cm->m)

    // Fresnel reflection coefficients based on the illumination angle
    double cs		= std::cos(theta);
    double si2		= std::sin(theta) * std::sin(theta);
    complex<double> stem	= sqrt(er*ur-si2);
    Rf.h	= abs( (ur*cs-stem)/(ur*cs+stem) );
    Rf.v	= abs( (er*cs-stem)/(er*cs+stem) );
    Rf.h	= Rf.h * Rf.h;
    Rf.v	= Rf.v * Rf.v;
    Ef.h	= 1.0 - Rf.h;
    Ef.v	= 1.0 - Rf.v;

    // coherent component
    Rco.h=Rf.h*std::exp(-pow((2.*std::cos(theta)*ks),2.0));
    Rco.v=Rf.v*std::exp(-pow((2.*std::cos(theta)*ks),2.0));

    // non-coherent component

    // set up integration limits
    double aa11=-3.14;
    double bb11=3.14;
    double aa22=0.0;
    double bb22=1.57;

    // generates zeros and abscissa
    quagen_(zt, wt, &npp);
    double as11=(bb11-aa11)/2.0;
    double bs11=(bb11+aa11)/2.0;
    for (int m=0; m<npp; m++)
    {
        ww11[m]=wt[m]*as11;
        zz11[m]=zt[m]*as11+bs11;
    }
    quagen_(zt, wt, &nss);
    double as22=(bb22-aa22)/2.0;
    double bs22=(bb22+aa22)/2.0;
    for (int m=0; m<nss; m++)
    {
        ww22[m]=wt[m]*as22;
        zz22[m]=zt[m]*as22+bs22;
    }

    Rnc.h = 0.0;
    Rnc.v = 0.0;

    int ft=int(fresnel_type);
    int ct=int(m_correlation);

    for (int i=0; i<npp; i++)
    {
        double phis=zz11[i];
        if(phis==PI_NUMBER) phis-=0.00001;
        for (int j=0; j<nss; j++)
        {
            double thetas=zz22[j];
            if(thetas==PI_NUMBER/2.0) thetas -= 0.00001;
            // initialize array
            double hh=0.0;
            double vv=0.0;
            double hv=0.0;
            double vh=0.0;

            if(thetas==theta) thetas += 0.00001;

            double real = er.real();
            double imag = er.imag();

            sigma_(&ks, &kl, &theta, &thetas, &phis, sigma0, &ft, &real, &imag, &ct);

            if(sigma0[0]>0.0) vv=sigma0[0];
            if(sigma0[1]>0.0) hh=sigma0[1];
            if(sigma0[2]>0.0) hv=sigma0[2];
            if(sigma0[3]>0.0) vh=sigma0[3];
            Rnc.h += ww11[i]*ww22[j]*(hh+vh)*std::sin(thetas)/(4.0*PI_NUMBER*std::cos(theta));
            Rnc.v += ww11[i]*ww22[j]*(vv+hv)*std::sin(thetas)/(4.0*PI_NUMBER*std::cos(theta));
        }
    }

    // reflectivity & emissivity
    R.h	= Rco.h + Rnc.h;
    R.v	= Rco.v + Rnc.v;
    E.h = 1.0 - R.h;
    E.v	= 1.0 - R.v;
    delete ww11,ww22,zz11,zz22,wt,zt;
}
void AIEM::active()
{
    double slope=soil.rms/soil.cl;
    double zs = soil.rms*soil.rms/soil.cl;
    if (m_correlation==AIEM_GAUSSIAN)
        slope/=1.414;
    double theta=soil.theta*PI_NUMBER/180;

    double ks=PI_NUMBER*soil.frequency*soil.rms/15.0;
    double kl=PI_NUMBER*soil.frequency*soil.cl/15.0;

    double co=std::cos(theta);
    double si=std::sin(theta);
    double si2=si*si;

    complex<double> rh=(sqrt(er)-1.0)/(sqrt(er)+1.0);
    //r0是水平、垂直极化的法向入射时的
    complex<double> r0=pow(abs(rh),2);
    //fresnel反射系数，见任论文46页,再转化成dB
    complex<double> r0db=10.0*log10(r0);

    complex<double> siem=co+sqrt(er-si2);
    complex<double> ssem=co-sqrt(er-si2);
    //cs=cos(theta)
    //这里是计算vv的极化幅度
    complex<double> avv1=(er-1.0)*(si2-er*(1.0+si2));
    complex<double> avv2=pow((er*co+sqrt(er-si2)),2);
    complex<double> avv=pow(abs(avv1/avv2),2);
    complex<double> avvr0=sqrt(avv*r0);

    rh=ssem/siem;
    //ah是水平极化fresnel反射系数
    complex<double> ah=pow(abs(rh),2);
    ssem=(er-1.0)*(si2-er*(1.0+si2));
    siem=er*co+sqrt(er-si2);
    complex<double> ssem1=er*co-sqrt(er-si2);
    complex<double> rv=ssem1/siem;
    //rvr,见任论文47页，垂直极化方式的fresnel反射系数，上面的ah是水平极化的
    complex<double> rvr=pow(abs(rv),2);
//这一部分是计算fresnel反射系数
//这里是计算rv+rh
    complex<double> rvrh=ah+rvr;


    int istart,iend;
    switch(fresnel_type)
    {
    case AIEM_TRANSITION:
        istart=3;
        iend=5;
        break;
    case AIEM_SIMULTANEOUS:
        istart=1;
        iend=5;
        break;
    default:
        istart=int(fresnel_type);
        iend=int(fresnel_type);
        break;
    }
    //为什么如此取值？
    //转换为弧度单位
    double phis=179.01*PI_NUMBER/180;
    double thetas=theta+0.00001;
    double hh,vv,sigma0[4];
    int ct=int(m_correlation);

    //std::cout<<ks<<" "<<kl<<" "<<theta<<" "<<thetas<<" "<<phis<<" "<<er.real()<<" "<<er.imag()<<" "<<ct<<std::endl;

    for (int i=istart; i<=iend; i++)
    {
        sigma_(&ks,&kl,&theta,&thetas,&phis,sigma0,&i,&er.real(),&er.imag(),&ct);
        //std::cout<<10*std::log10(sigma0[0])<<" "<<10*std::log10(sigma0[1])<<std::endl;
    }
    //could this value less then 0?
    //it should be reviewed here!
    if(sigma0[0]>0.0)
        E.v=10*std::log10(sigma0[0]);
    if(sigma0[1]>0.0)
        E.h=10*std::log10(sigma0[1]);
    else
        E.h=100.0;
    //std::cout<<hh<<" "<<vv<<std::endl;
}
