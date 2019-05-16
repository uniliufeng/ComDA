#include <complex>
#include <cmath>
#include <iostream>
#include <cstdio>
using namespace std;

class AIEM{
public:
		void emissivity();
		void sigma(double ks,double kl,double theta,double thetas,double phis,double sigma0[],int irc);
		complex<double> expkc1(const complex<double> q);
		complex<double> expkc2(const complex<double> q);
		complex<double> expc1(const complex<double> q,const complex<double> qp);
		complex<double> expc2(const complex<double> q,const complex<double> qp);
		complex<double> expc3(const complex<double> q,const complex<double> qp);
		complex<double> expc4(const complex<double> q,const complex<double> qp);
		complex<double> Favv(const double u,const double v,const complex<double> q,const complex<double> qfix);
		complex<double> Fahh(const double u,const double v,const complex<double> q,const complex<double> qfix);
		complex<double> Fbvv(const double u,const double v,const complex<double> q,const complex<double> qfix);
		complex<double> Fbhh(const double u,const double v,const complex<double> q,const complex<double> qfix);
		complex<double> Fahv(const double u,const double v,const complex<double> q,const complex<double> qfix);
	   	complex<double> Favh(const double u,const double v,const complex<double> q,const complex<double> qfix);
		complex<double> Fbhv(const double u,const double v,const complex<double> q,const complex<double> qfix);
		complex<double> Fbvh(const double u,const double v,const complex<double> q,const complex<double> qfix);
		double bessj0(const double x);
		void shadow(const bool back,const double ti,const double ts, double shfct);
		double erfcc(const double x);
		double erfc(const double x);
		double gammp(const double a,const double x);
		double gammq(const double a,const double x);
		void gser(double gamser,const double a,const double x,double gln);
		void gcf(double gammcf,const int a,const double x,double gln);
		double gammln(const double xx);
		double BesselK(const int n,const double x);
		double alogam(const double x);
		double BESSK(const double N, const double X);
		double BESSK0(const double X);
		double BESSK1(const double X);
		double BESSI0(const double X);
		double BESSI1(const double X);
		void quagen(double z[],double wt[],const int n);

		complex<double> er,ur,rh,rv,rvh;
		int itype,iterm;
		double si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2,ks2,cl,effslop,gamser,gln,gammcf; 
		double w[1001];
		double pi;
		
				
};