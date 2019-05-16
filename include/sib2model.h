#ifndef __LDAS_SIB2MODEL_H
#define __LDAS_SIB2MODEL_H
#include <iostream>
using namespace std;
#include <fstream>
#include <cmath>
//#include "modeloperator.hpp"
#include "exception.h"
#include "model.h"
#include "constant.h"

namespace ldas
{
class SiB2Config:public ModelConfig
{
public:
    SiB2Config():ModelConfig()
    {
        //time_step=3600; // 1 hour
        //output_span=1;
        //restart_span=240;//10 days
        forcing_path="data2";
        vegetation_path="data1";
        output_path="./";
    }
    SiB2Config& operator=(const SiB2Config& other);
};
/** time constant variable */
class SiB2Parameter
{
    //todo
};
/** time variant variable */
class SiB2Variable
{
    //todo
};
class SiB2Restart
{
    //
};
class SiB2Model:public BaseModel
{
public:
    /** Default constructor */
    SiB2Model();
    /** Default destructor */
    virtual ~SiB2Model();
    /** Copy constructor
     *  \param other Object to copy from
     */
    SiB2Model(const SiB2Model& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    SiB2Model& operator=(const SiB2Model& other);
    virtual void run();
    virtual void config(const SiB2Config& cfg);
    virtual SiB2Config config() const;
    //virtual void output();
    //virtual void run();
protected:
    void initParam(const char* filename);
    SiB2Config m_config;
private:
    // prognostic variables
    double tc, tg, td, qa, capac[2], snoww[2], www[3];

    // physical constants
    static const double g=9.81;
    //static const double pie=3.14159265;
    double timcon, cpair, rhoair, psy, hlat, vkc,
           snomel, stefan, tf, clai, cw, snofac, asnow, rcp, kappa, epsfac;
    double po2m, pco2m;

    // vegetation : static, dynamic, derived parameters
    int ivtype, istype;									// gridij
    double z2, z1, vcover, chil, tran[2][2], ref[2][2], rootd, ph1, ph2, phc,
           effcon, gradm, binter, respcp, atheta, btheta,
           trda, trdm, trop, tpsa, tpsb, tpsc, slti, hlti, shti, hhti;  //vstate
    double zlt, green, fparc;							// vdyijt
    double z0d, dd, cc1, cc2, vmax0, gmudmu;			// vderiv

    // soils : space-varying, type-dependent parameters
    double sodep, soref[2];								// soilij
    double bee, phsat, poros, satco, slope, zdepth[3];	// soils
    // input atmospheric and site data
    double em, tm, um, zm, psur, ppc, ppl, radn[3][2],
           sunang, swdown, rnetm, cloud, bps;				// atmos

    // site parameters specific to 1-d model operation
    double corb1, corb2, ha, g1, g2, g3, ztz0, zwind, zmet, ht; // caerod
    double zlong, zlat, salb[2][2], rab[2][3][2];			// site
    double dtt;
    int itrunk, ilw, niter, iter, ispare;				// steps
    double time, year, month, day, hour;				// govern

    // variables returned to g.c.m. ( et(kg), h(w m-2), runoff(m),
    //                                lw(w m-2), drag(kg m-1 s-2))
    double etmass, hflux, roff, zlwup, drag;			// donor

    // variables calculated from above and ambient conditions
    double z0, d, rbc, rdc;								// rause
    double ra, rb, rd;									// aerorx
    double tgs, ta, ea, etc, etgs, getc, getgs ,u2, ustar;	// grads
    double albedo[2][2][2], radfac[2][2][2], radt[2], thermk, exrain, tgeff;
    // radabs
    double rst, rstfac[4], rsoil, cog1, cog2, hr, fc, fg;	// surfrs
    double satcap[2], wc, wg, canex, areas;				// hydrol
    double ccx, cg, csoil;								// stores
    double dtc, dtg, dtd, dth, dqm;						// delts
    double assimn, respc, respg, pco2i,gsh2o;			// carbio

    // heat fluxes : c-canopy, g-ground, t-trans, e-evap  in j m-2
    double ec, eg, hc, hg, chf, shf, ect, eci, egi, egs,
           ecmass, egmass, heaten;							// flux

    // snow variables
    double tsnow, rsnow;								// snow

    double rngdtg, rngdtc, rncdtg, rncdtc;
    double hgdtg,  hgdtc,  hgdth,  hcdtg,  hcdtc,  hcdth;
    double egdtg,  egdtc,  egdqm,  ecdtg,  ecdtc,  ecdqm;
    double deadtc, deadtg, demdqm, deadqm, aag, aac, aam, bbg, bbc, bbm;
    // flxdif
private:
    ///init the SiB2 varable arrays
    void init();
	/// initialization of physical constants of the SiB2
    void const2();

    ///Reading of meteorological data. preparation of forcing variables  data required to run sib
    void driver(ifstream& iu, const int isnow, int& nymd);
    ///solar zenith angle computation; downcoming radiation at bottom
    void radc2();

    /// energy and water balance check
    void balan(int iplace, const int nymd);

    ///calculation of  interception and drainage of rainfall
    void inter2();
    ///alteration of aerodynamic transfer properties in case of snow   accumulation.
    void snow1();
    ///temperature change due to addition of precipitation
    void adjust(double& ts, const double spechc, const double capacp,
                const double snowwp, const int iveg);
    ///marginal situation: snow exists in patches at temperature tf   with remaining area at temperature tg > tf.
    void patchs(const double p0);

    ///calculation of albedos via two stream approximation (direct and diffuse )   and partition of radiant energy
    void rada2();
    ///calculation of downward longwave.
    void longrn(double tranc1[2], double tranc2[2], double tranc3[2]);

    ///calculation of canopy and ground temperature increments   over time step, fluxes derived.
    void begtem();

    ///calculation of ea, ta, ra, rb, rd and soil moisture stress for the beginning   of the time step
    void endtem(const int ipbl);
    /// calculation of ustar, u2, ra and drag using Paulson's method.
    void rasite();
    /// calculation of Paulson psi-function for unstable condition
    void unstab(const double uest, const double a, const double b,
                const double argz, const double heat, double& psione, double& psitwo);
    ///calculation of Paulson psi-function for stable condition
    void stab(const double uest, const double a, const double b,
              const double argz, const double heat, double& psione, double& psitwo);
    /// calculation of ra for heat between z2 and zmet
    double rafcal(const double zl, const double uest, const double heat);
    ///the newton raphson iterative routine
    void newton(double& a1, double& y, const double finc, int& nox,
                const int nonpos, int& iwolk, const int l);

    ///calculation of rb and rd as functions of u2 and temperatures
    void rbrd();

    ///calculation of canopy photosynthetic rate using the integrated model  relating assimilation and stomatal conductance.
    void phosib();
    ///arranges successive pco2/error pairs in order of increasing pco2.
    void sortin(double eyy[6], double pco2y[6], const double range,
                const double gammas, const int ic);
    ///calculation equivalent to steps in figure 4 of SE-92A   c4 calculation based on CO-92.
    void cycalc(const double fparkk, const double vm, const double gradm,
                const double bintc, const double atheta, const double btheta,
                const double gah2o, const double gbh2o, const double gog1,
                const double gog2, const double wc, const double h2oi,
                const double h2om, const double h2osl, const double par,
                const double pco2m, const double psur, const double gammas,
                const double respc, const double respg, const double rrkk,
                const double omss, const double c3, const double c4,
                const double pco2i, double& eyy, double& gsh2o, double& assimn,
                double& h2os, double& h2oa);

    /// partial derivatives of radiative and sensible heat fluxes
    void delrn();

    ///calculation of partial derivatives of canopy and ground sensible heat  fluxes with respect to tc, tgs, and theta-m.
    void delhf();
    ///calculation of partial derivatives of canopy and ground latent heat fluxes  with respect to tc, tgs, theta-m, and qm.
    void delef();
    ///solve for time changes of pbl and sib variables using a semi-implicit scheme dtc, dtg, dth, dqm
    void sibslv();
    /// solve a linear system by gaussian elimination.
    void gauss(const double a[4][5], const int n, const int np1, double x[4]);
    ///calculation of temperature tendencies assuming no interaction with the pbl
    void dtcdtg();

    /// updating of all prognostic variables.
    void updat2();
    ///snowmelt / refreeze calculation
    void snow2();
    /// calculation of interflow, infiltration excess and loss to   groundwater .
    void run2();

    ///output of results to files
    void outer (ofstream& out, ofstream& out1, ofstream& out2, ofstream& out3,
                ofstream& out4, const int nymd );

    double amin1(double a, double b);
    double amax1(double a, double b);
    int min0(int a, int b);
    int max0(int a, int b);
    double sign(double a, double b);
    double e(double x);
    double ge(double x);
};
}
#endif // __LDAS_SIB2MODEL_H
