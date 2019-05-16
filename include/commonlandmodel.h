/**
\file Common Land Model
*/
#ifndef __LDAS_COMMONLANDMODEL_H
#define __LDAS_COMMONLANDMODEL_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
//#include "modeloperator.hpp"
#include "constant.h"
#include "model.h"
#include "landgrid.h"

namespace ldas
{
namespace colm
{
///常量定义
const int MAXSOILL     = 10;                      // number of soil layers
const int MAXSNL       = 5;                       // max number of snow layers
const int NFCON        = 9*MAXSOILL+29;           // number of time constant variables
const int NFTUNE       = 14;                      // number of clm tunable constants
const int NFVAR        = 5*(MAXSOILL+MAXSNL)+51;  // number of time varying variables
const int NFORC        =18;                       // number of forcing variables
const int NFLDV        =92;                       // number of output fluxes
const int NFLAI        =4;                        // number of leaf time varying variables
const int MAXPATCH     =25;                       // number of clm grid points
const int NLANDCATEG   =25;                       // number of land cover categories
const int NSOILCATEG   =17;                       // number of soil texture categories
const int MAXLONPOINT  =10;                       // max number of longitude points on model grid
const int MAXLATPOINT  =10;                       // max number of latitude points on model grid
const int MAXLAKEL     =10;                       // number of lake layers
const int ASSVAR       = 48;                      // assimilation variable number
const int NFORC_IN     = 9;                       // number of forcing variable from file
const int OUTFLUX      = 22;                       // flux variable number need to output
const int LAINUM = 12;
const int VARNUM = 48;
const int FLUXNUM = OUTFLUX;
///model parameter
const double denice = 917.0;      // density of ice [kg/m3]
const double denh2o = 1000.0;     // density of liquid water [kg/m3]
const double cpliq  = 4188.0;     // Specific heat of water [J/kg-K]
const double cpice  = 2117.27;    // Specific heat of ice [J/kg-K]
const double cpair  = 1004.64;    // specific heat of dry air [J/kg/K]
const double hfus   = 0.3336e6;   // latent heat of fusion for ice [J/kg]
const double hvap   = 2.5104e6;   // latent heat of evap for water [J/kg]
const double hsub   = 2.8440e6;   // latent heat of sublimation [J/kg]
const double tkair  = 0.023;      // thermal conductivity of air [W/m/k]
const double tkice  = 2.290;      // thermal conductivity of ice [W/m/k]
const double tkwat  = 0.6;        // thermal conductivity of water [W/m/k]
const double tfrz   = 273.16;     // freezing temperature [K]
const double rgas   = 287.04;     // gas constant for dry air [J/kg/K]
const double roverg = 4.71047e4;  // rw/g = (8.3144/0.018)/(9.80616)*1000. mm/K
const double rwat   = 461.296;    // gas constant for water vapor [J/(kg K)]
const double grav   = 9.80616;    // gravity constant [m/s2]
const double vonkar = 0.4;        // von Karman constant [-]
const double stefnc = 5.67e-8;    // Stefan-Boltzmann constant  [W/m2/K4]
const double MAXSNOWDP = 1.0;     // max snow depth
#define  EcoDynamics
}
using namespace ldas::colm;
using namespace std;
class CoLMConfig:public ModelConfig
{
public:
    CoLMConfig();
    virtual ~CoLMConfig() {}
    /** Read the config text file
    * \param The config text filename
    */
    virtual void open(const char* fn);
    std::string constant_path,variable_path,mask_path;
    std::string ssmi_path,amsr_path;
    std::string clay_path,sand_path;
    int row_num,col_num,sub_row_num,sub_col_num,sub_start_row_num,sub_start_col_num;
    int ensemble_size;
    double lat_span,lon_span,lat_corner,lon_corner;
};
class CoLMParameter:public ModelParameter
{
public:
    CoLMParameter();
    CoLMParameter(const int patches);
    virtual ~CoLMParameter();
    virtual void open(const char* fn);
    virtual void open(const char* fn,const DATA_TYPE dt, const int patches=0);
    virtual void save(const char* fn);
    virtual void patches(const unsigned int p);

    double** fvar;
    float* lai;
    Time time;
private:
    bool m_use_lai;
};
///time constant parameter
class CoLMConstant:public ModelParameter
{
public:
    CoLMConstant();
    CoLMConstant(const int patches);
    virtual ~CoLMConstant();
    virtual void open(const char* fn);
    virtual void save(const char* fn);
    virtual void patches(const unsigned int p);

    double ftune[NFTUNE];
    double** fcon;
};
class CoLMForcing
{
public:
    CoLMForcing();
    CoLMForcing(const int sum);
    virtual ~CoLMForcing();
    virtual void open(const char* fn);
    virtual void patches(const unsigned int p);

    float** forc;
private:
    int numpatch;
};


/** \brief Command Land Model, C++ version
*
* The fortran version's Original author : Yongjiu Dai
*/
class CommonLandModel:public BaseModel
{
public:
    /** Default constructor */
    CommonLandModel();
    /** Default destructor */
    virtual ~CommonLandModel();
    /** Copy constructor
     *  \param other Object to copy from
     */
    CommonLandModel(const CommonLandModel& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    CommonLandModel& operator=(const CommonLandModel& other);

    /** Constructor with a config text file
    * \param The config filename
    * \return A reference to this
    */
    CommonLandModel(const std::string config_filename);

    void init();

    virtual void run();
    virtual void run(const int ipatch);

    /** \brief CLM MODEL DRIVER
    * 单点单步运行: single point and single step
    * \param idate 输入变量，指定模型日期，格式为{年，日，秒}，其中日是指julian day，即是年中的天数，格林威治时间（UTC0）
    * \param deltim time step，即单步的具体时间，单位为秒
    * \param ftune tunable constants，可调参数？
    * \param fcon time constant variables (model restatrt required)
    * \param forc forcing variables
    * \param fvar time varying variables (clm restatrt required)，地表变量，输出
    * \param flux flux variables, 通量计算结果
    * \param dolai true if time for time-varying vegetation paramter
    * \param doalb true if time for surface albedo calculation
    * \param dosst true if time for update sst/ice/snow
    */
    virtual void run (double ftune[NFTUNE], double fcon[NFCON], double forc[NFORC], double fvar[NFVAR],double flux[OUTFLUX]);
    /** \brief model output, which is num*row*col length
    * \param fn the filename of the the output data
    * \param num the second dimension length of the outgrid array
    * \param outgrid a two dimension array which is outgrid[total_patches][num]
    */
    virtual void output(const char* fn, const int num, double** outgrid);
    virtual void config(const char* fn);
    int snowlayer() const;
    //virtual void config(const CoLMConfig& c);
    //virtual void config() const;
    void radiation(int jday,int msec,double lon,double lat, double swdown, double &sols, double &soll, double &solsd, double &solld);
protected:
    CoLMConfig m_config;
    CoLMParameter parameter;
    CoLMConstant constant;
    LandGrid grid;
    CoLMForcing forcing;
    double** flux_out; //output flux
    double** fvar_out; //output variable
private:
    //bool m_create_output;
    int snl;
private:
    void albland(int    itypwat,   double albsol,   double chil,      double ref[2][2], double tran[2][2],
                 double fveg,      double green,    double lai,       double sai,       double coszen,
                 double wt,        double fsno,     double scv,       double sag,       double ssw,
                 double tg,        double alb[2][2],double albg[2][2],double albv[2][2],double ssun[2][2],
                 double ssha[2][2],double &thermk,   double &extkb,     double &extkd);
    void albocean (double oro,double scv,double coszrs,double alb[2][2]);
    double sign(double a, double b);
    void combo (double &dz,   double &wliq, double &wice,double &t,double dz2,
                double wliq2,double wice2,double t2 );
    void dewfraction (double sigf,double lai,double sai,double dewmx,double ldew,
                      double &fwet,double &fdry);

    void eroot(int lb,int    nl_soil,           double trsmx0,             double porsl[MAXSOILL],     double bsw[MAXSOILL],        double phi0[MAXSOILL],
               double rootfr[MAXSOILL],  double dz[MAXSNL+MAXSOILL],double tss[MAXSNL+MAXSOILL],double wliq[MAXSNL+MAXSOILL],double rootr[MAXSOILL],
               double &etrc,              double &rstfac);
    void groundfluxes (double zlnd,  double zsno,  double hu,      double ht,     double hq,
                       double us,    double vs,    double tm,      double qm,     double rhoair,
                       double psrf,  double ur,    double thm,     double th,     double thv,
                       double tg,    double qg,    double dqgdT,   double htvp,   double fsno,
                       double sigf,  double &cgrnd,double &cgrndl, double &cgrnds,double &taux,
                       double &tauy, double &fsena,double &fevpa,  double &fseng, double &fevpg,
                       double &tref, double &qref, double &z0ma,   double &zol,   double &rib,
                       double &ustar,double &qstar,double &tstar,  double &f10m,  double &fm,
                       double &fh,   double &fq);
    void groundtem(int    itypwat,              int     lb,                    int     nl_soil,              double  dtime,                  double  capr,
                   double  cnfac,               double  csol[MAXSOILL],        double  porsl[MAXSOILL],      double  dkmg[MAXSOILL],         double  dkdry[MAXSOILL],
                   double  dksatu[MAXSOILL],    double  sigf,                  double  dz[MAXSNL+MAXSOILL],  double  z[MAXSNL+MAXSOILL],     double  zi[MAXSNL+MAXSOILL+1],
                   double  tss[MAXSNL+MAXSOILL],double  wice[MAXSNL+MAXSOILL], double  wliq[MAXSNL+MAXSOILL], double  &scv,                    double  &snowdp,
                   double  frl,                 double  dlrad,                 double  sabg,                 double  fseng,                  double  fevpg,
                   double  cgrnd,               double  htvp,                  double  emg,                  int     imelt[MAXSNL+MAXSOILL], double & sm,
                   double  &xmf,                 double  fact[MAXSNL+MAXSOILL]);
    void hCapacity (int    itypwat,        int lb,                       int   nl_soil1,                 double csol[MAXSOILL],
                    double porsl[MAXSOILL],double wice[MAXSNL+MAXSOILL], double wliq[MAXSNL+MAXSOILL],   double scv,
                    double dz[MAXSNL+MAXSOILL],           double cv[MAXSNL+MAXSOILL]);
    void hConductivity (int   itypwat,                int    lb,                     int    nl_soil1,              double dkmg[MAXSOILL],       double dkdry[MAXSOILL],
                        double dksatu[MAXSOILL],      double porsl[MAXSOILL],        double dz[MAXSNL+MAXSOILL],   double z[MAXSNL+MAXSOILL],   double zi[MAXSNL+MAXSOILL+1],
                        double tss[MAXSNL+MAXSOILL],  double wice[MAXSNL+MAXSOILL],  double wliq[MAXSNL+MAXSOILL], double tk[MAXSNL+MAXSOILL]);
    void lai_empirical(int    ivt, int    nl_soil,double rootfr[MAXSOILL],double t[MAXSOILL+MAXSNL],double &lai,
                       double &sai, double &fveg,   double &green);
    void lake (int    nl_lake ,              int    itypwat ,                double dlat,                 double dtime                 , double zlak[MAXSNL+MAXSOILL],
               double dzlak[MAXSNL+MAXSOILL],double zilak[MAXSNL+MAXSOILL+1],double hu  ,                   double ht,                     double hq ,
               double us   ,                 double vs   ,                   double tm  ,                 double qm ,                    double prc,
               double prl ,                  double rhoair  ,                double psrf,                 double sabg  ,                 double frl,
               double &tg ,                   double tlak[MAXSNL+MAXSOILL],   double wliq[MAXSNL+MAXSOILL],  double wice[MAXSNL+MAXSOILL]  ,double &scv,
               double &snowdp,                double &trad,                    double &tref,                 double &qref ,                  double &taux,
               double &tauy ,                 double &fsena,                   double &fevpa,                double &lfevpa,                 double &fseng,
               double &fevpg,                 double &olrg ,                   double &fgrnd,                double &tcrit ,                 double &emis ,
               double &z0ma ,                 double &zol ,                    double &rib ,                 double &ustar ,                 double &qstar,
               double &tstar,                 double &u10m ,                   double &v10m,                 double &f10m  ,                 double &fm ,
               double &fh   ,                 double &fq);
    void leafinterception (double dtime,double dewmx,double chil, double prc,
                           double prl,  double tm,   double scv,  double sigf,
                           double lai,  double sai,
                           double &ldew,double &pg);
    void leaftemone(double dtime , double csoilc ,double dewmx  ,double htvp   ,double lai    ,
                    double sai   , double displa ,double sqrtdi ,double z0m    ,double effcon ,
                    double vmax25 ,double slti   ,double hlti   ,double shti   ,double hhti   ,
                    double trda   ,double trdm   ,double trop   ,double gradm  ,double binter ,
                    double extkn  ,double extkb  ,double extkd  ,double hu     ,double ht 	,
                    double hq     ,double us     ,double vs     ,double thm    ,double th      ,
                    double thv    ,double qm     ,double psrf   ,double rhoair ,double par    ,
                    double sabv   ,double frl    ,double thermk ,double rstfac ,double po2m   ,
                    double pco2m  ,double sigf   ,double etrc   ,double tg     ,double qg     ,
                    double dqgdT  ,double emg    ,double & tl     ,double & ldew   ,double & taux   ,
                    double & tauy   ,double & fseng  ,double & fevpg  ,double & cgrnd  ,double & cgrndl ,
                    double & cgrnds ,double & tref   ,double & qref   ,double & rst    ,double & assim   ,
                    double & respc  ,double & fsenl  ,double & fevpl  ,double & etr    ,double & dlrad  ,
                    double &  ulrad  ,double & z0ma   ,double & zol    ,double & rib    ,double & ustar  ,
                    double & qstar  ,double & tstar  ,double & f10m   ,double & fm     ,double & fh     ,double & fq);
    void leaftemtwo (double dtime ,    double csoilc  ,    double dewmx  ,    double htvp   ,    double lai    ,
                     double sai   ,    double displa   ,   double sqrtdi  ,   double z0m    ,    double effcon ,
                     double vmax25 ,   double slti    ,    double hlti     ,  double shti    ,   double hhti   ,
                     double trda   ,   double trdm   ,     double trop    ,   double gradm    ,  double binter  ,
                     double extkn  ,   double extkb  ,     double extkd  ,    double hu      ,   double ht      ,
                     double hq      ,  double us     ,     double vs     ,    double thm    ,    double th      ,
                     double thv      , double qm      ,    double psrf   ,    double rhoair ,    double parsun ,
                     double parsha  ,  double sabvsun  ,   double sabvsha ,   double frl    ,    double fsun   ,
                     double thermk ,   double rstfac  ,    double po2m     ,  double pco2m   ,   double sigf   ,
                     double etrc   ,   double tg     ,     double qg      ,   double dqgdT    ,  double emg     ,
                     double & tlsun  ,   double & tlsha  ,     double & ldew   ,    double & taux    ,   double & tauy     ,
                     double & fseng   ,  double & fevpg  ,     double & cgrnd  ,    double & cgrndl ,    double & cgrnds  ,
                     double & tref     , double & qref    ,    double & rst    ,    double & assim  ,    double & respc  ,
                     double & fsenl   ,  double & fevpl    ,   double & etr     ,   double & dlrad  ,    double & ulrad  ,
                     double & z0ma   ,   double & zol     ,    double & rib      ,  double & ustar   ,   double & qstar  ,
                     double & tstar  ,   double & f10m   ,     double & fm      ,   double & fh       ,  double & fq  );
    void CLMMAIN (double dtime, bool doalb, bool dolai, bool dosst,
                  int nl_soil, int maxsnl, double dlon, double dlat,
                  int itypwat, int ivt, double &oro,

                  // soil information
                  double albsol, double csol[MAXSOILL], double porsl[MAXSOILL], double phi0[MAXSOILL], double bsw[MAXSOILL],
                  double dkmg[MAXSOILL], double dksatu[MAXSOILL], double dkdry[MAXSOILL], double hksati[MAXSOILL],

                  // vegetation information
                  double z0m, double displa, double sqrtdi,
                  double effcon, double vmax25, double slti,
                  double hlti, double shti, double hhti,
                  double trda, double trdm, double trop,
                  double gradm, double binter, double extkn,
                  double chil, double ref[2][2], double tran[2][2], double rootfr[MAXSOILL],

                  // atmospheric forcing
                  double frl, double sols, double soll, double solsd, double solld,
                  double pco2m, double po2m, double us, double vs, double tm, double qm,
                  double prc, double prl, double psrf, double rhoair,
                  double hu, double ht, double hq,

                  // land surface variables required for restart
                  int year, int jday, int msec,
                  double z[MAXSNL+MAXSOILL], double dz[MAXSNL+MAXSOILL],
                  double tss[MAXSNL+MAXSOILL], double wliq[MAXSNL+MAXSOILL], double wice[MAXSNL+MAXSOILL],
                  double &tg, double &tlsun, double &tlsha,
                  double &ldew, double &sag, double &scv, double &snowdp,
                  double &fveg, double &fsno, double &sigf,
                  double &green, double &lai, double &sai,
                  double &coszen, double albg[2][2], double albv[2][2],
                  double alb[2][2], double ssun[2][2], double ssha[2][2], double &thermk,
                  double &extkb, double &extkd,

                  // fluxes
                  double &taux,  double &tauy,
                  double &fsena, double &fevpa, double &lfevpa,
                  double &fsenl, double &fevpl, double &etr,
                  double &fseng, double &fevpg, double &olrg,
                  double &fgrnd, double &trad, double &tref, double &qref,
                  double &rsur, double &rnof, double &rst,
                  double &assim, double &respc,
                  double &parsun, double &parsha, double &sabvsun,
                  double &sabvsha, double &sabg, double &sabvg, double &xerr, double &zerr,

                  // TUNABLE modle constants
                  double &zlnd,   double &zsno,   double &csoilc, double &dewmx,  double &wtfact,
                  double &capr,   double &cnfac,  double &ssi,    double &wimp,   double &pondmx,
                  double &smpmax, double &smpmin, double &trsmx0, double &tcrit,

#ifndef EcoDynamics
                  // time-varying vegetation parameters from read-in file
                  double &lai_r, double &sai_r, double &green_r, double &fveg_r,
#endif

                  // additional variables required by coupling with WRF model
                  double &emis, double &z0ma, double &zol,
                  double &rib, double &ustar, double &qstar, double &tstar,
                  double &u10m, double &v10m, double &f10m, double &fm, double &fh, double &fq,
                  int    &snl);
    void meltf (int   lb,                    int    nl_soil,               double dtime,                  double fact[MAXSNL+MAXSOILL], double brr[MAXSNL+MAXSOILL] ,
                double hs,                   double dhsdT,                  double tssbef[MAXSNL+MAXSOILL],double tss[MAXSNL+MAXSOILL],  double wliq[MAXSNL+MAXSOILL],
                double wice[MAXSNL+MAXSOILL],int    imelt[MAXSNL+MAXSOILL], double &scv,                    double &snowdp, double &sm, double &xmf);
    void moninobuk(double hu,   double ht,   double hq,     double displa, double z0m,
                   double z0h,  double z0q,  double obu,    double um,     double &ustar,
                   double &temp1,double &temp2,double &temp12m,double &temp22m,double &f10m,
                   double &fm,   double &fh,   double &fq);
    void moninobukini(double ur, double th,  double thm,  double thv,double dth,
                      double dqh,double dthv,double zldis,double z0m,double &um,double &obu);
    void netsolar (int    itypwat  ,   double sigf,      double albg[2][2],  double albv[2][2],double alb[2][2],
                   double ssun[2][2],  double ssha[2][2],double sols,        double soll,      double solsd,
                   double solld,                                                                //in
                   double &parsun,    double &parsha,      double &sabvsun,   double &sabvsha,      //out
                   double &sabg,        double &sabvg );

    void newsnow (int    itypwat,                int    maxsnl,  double dtime,  double tm,double tg,
                  double pg,                     double tcrit,                                                   //in
                  double zi[MAXSNL+MAXSOILL+1],  double z[MAXSNL+MAXSOILL],     double dz[MAXSNL+MAXSOILL],
                  double tss[MAXSNL+MAXSOILL],   double wliq[MAXSNL+MAXSOILL],  double wice[MAXSNL+MAXSOILL],
                  double fiold[MAXSNL+MAXSOILL], int    &snl,                     double &sag,                     //inout
                  double &scv,			         double &snowdp,
                  double &pg_rain,               double &pg_snow);                                                                 //out
    double orb_coszen(double calday,double lon,double lat);
    double  psi(int k,double zeta);
    void qsadv(double T,double p,double &es,double &esdT,double &qs,double &qsdT);


    void seafluxes (double   oro,    double   hu,    double   ht,    double   hq,    double  us,
                    double   vs,     double   tm,    double   qm,    double   rhoair,double  psrf,
                    double   tssea,  double   &taux,  double  & tauy,  double  & fsena, double   &fevpa,
                    double   &fseng,  double   &fevpg, double   &tref,  double   &qref,  double   &z0ma,
                    double   &zol,    double   &rib,   double  &ustar,  double   &qstar, double   &tstar,
                    double   &u10m,   double   &v10m,  double   &f10m,  double   &fm,    double   &fh,
                    double   &fq,     double   &cgrndl,double   &cgrnds);
    void snowage(double dtime,double tg,double scv,double scvold,double &sag );
    void snowcompaction(int    lb,                   double dtime,               int    imelt[MAXSNL+MAXSOILL], double fiold[MAXSNL+MAXSOILL], double tss[MAXSNL+MAXSOILL],
                        double wliq[MAXSNL+MAXSOILL], double wice[MAXSNL+MAXSOILL],double dz[MAXSNL+MAXSOILL] );
    void snowfraction (double fveg,double z0m,double snowdp,double &wt,double &sigf,double &fsno);
    void snowlayerscombine(int lb,                      int     &snl,                  double z[MAXSNL+MAXSOILL],  double  dz[MAXSNL+MAXSOILL],double zi[MAXSNL+MAXSOILL+1],
                           double wliq[MAXSNL+MAXSOILL],double  wice[MAXSNL+MAXSOILL],double tss[MAXSNL+MAXSOILL],double &scv,
                           double &snowdp );

    void snowlayersdivide (int lb,                      int    &snl,                  double z[MAXSNL+MAXSOILL],  double dz[MAXSNL+MAXSOILL],double zi[MAXSNL+MAXSOILL+1],
                           double wliq[MAXSNL+MAXSOILL],double wice[MAXSNL+MAXSOILL],double tss[MAXSNL+MAXSOILL]);
    void snowwater (int   lb, double dtime, double ssi, double wimp,double pg_rain,
                    double qseva, double qsdew, double qsubl, double qfros, double  dz[MAXSNL+MAXSOILL],
                    double wice[MAXSNL+MAXSOILL],  double  wliq[MAXSNL+MAXSOILL], double &qout_snowb);

    void SOCEAN(bool    dosst,     double dtime,    double  &oro,   double  hu,    double  ht,
                double  hq,       double     us,    double  vs,     double  tm,    double  qm,
                double  rhoair,   double  psrf,     double  sabg,   double frl,   double  &tssea,
                double  tssub[4],  double &scv,     double & taux,   double & tauy,  double & fsena,
                double & fevpa,    double & lfevpa,   double & fseng,  double & fevpg, double & tref,
                double & qref,     double & z0ma,     double & zol,    double & rib,   double & ustar,
                double & qstar,    double & tstar,    double & u10m,   double & v10m,  double & f10m,
                double & fm,       double & fh,       double & fq,     double & emis,  double & olrg);
    void soilwater(int    nl_soil,              double dtime,               double wimp,             double smpmin,            double porsl[MAXSOILL],
                   double phi0[MAXSOILL],       double bsw[MAXSOILL],       double hksati[MAXSOILL], double z[MAXSNL+MAXSOILL],double dz[MAXSNL+MAXSOILL],
                   double zi[MAXSNL+MAXSOILL+1],double tss[MAXSNL+MAXSOILL],double vol_liq[MAXSOILL],double vol_ice[MAXSOILL], double eff_porosity[MAXSOILL],
                   double qinfl,                double etr,                 double rootr[MAXSOILL],  double dwat[MAXSOILL],    double hk[MAXSOILL],
                   double dhkdw[MAXSOILL] );
    void sortin( double eyy[6], double pco2y[6], double range,double gammas,int ic );
    void srftsb(int isrfty,double dtime,double fnt,double dfntdt,double snowh,double tsbsf[4]);
    void stomata(double   vmax25, double    effcon, double    slti, double    hlti, double    shti,
                 double   hhti,   double    trda,   double    trdm, double    trop, double    gradm,
                 double   binter, double    tm,     double    psrf, double    po2m, double    pco2m,
                 double    pco2a, double    ea,     double    ei,   double    tlef, double    par,
                 double    rb,    double    ra,     double    rstfac,double    cint[3],double &   assim,
                 double &   respc, double &   rst);

    void subsurfacerunoff (int    nl_soil,               double dtime,       double pondmx,         double dzmm[MAXSNL+MAXSOILL],double wliq[MAXSNL+MAXSOILL],
                           double eff_porosity[MAXSOILL],double hk[MAXSOILL],double dhkdw[MAXSOILL],double dwat[MAXSOILL],double &rsubst);

    void surfacerunoff(int    nl_soil,          double wtfact,           double wimp,                  double bsw[MAXSOILL],      double porsl[MAXSOILL],
                       double phi0[MAXSOILL],   double hksati[MAXSOILL], double z[MAXSNL+MAXSOILL],    double dz[MAXSNL+MAXSOILL],double zi[MAXSNL+MAXSOILL+1],
                       double vol_liq[MAXSOILL],double vol_ice[MAXSOILL],double eff_porosity[MAXSOILL],double gwat,               double &rsur);
    void THERMAL(int   itypwat ,          int lb      ,                    int   nl_soil,    		   double dtime  ,          	 double trsmx0 ,
                 double zlnd  ,           double zsno    ,			     double csoilc ,			   double dewmx  ,  			 double capr   ,
                 double cnfac   ,		   double csol[MAXSOILL]   ,         double porsl[MAXSOILL] ,      double phi0[MAXSOILL],		 double bsw[MAXSOILL] ,
                 double dkmg[MAXSOILL]   ,double dkdry[MAXSOILL]  ,         double dksatu[MAXSOILL],	   double lai    ,  			 double sai    ,
                 double z0m     ,         double displa  ,                 double sqrtdi ,      		   double rootfr[MAXSOILL] ,	 double effcon ,
                 double vmax25  ,         double slti    ,                 double hlti   ,  			   double shti   ,  			 double hhti   ,
                 double trda    ,         double trdm    ,                 double trop   ,      		   double gradm  ,		    	 double binter ,
                 double extkn   ,         double hu      ,                 double ht     ,   		   double hq     ,  			 double us     ,
                 double vs      ,         double tm      ,                 double qm     ,  			   double rhoair ,		     double psrf   ,
                 double pco2m   ,         double po2m    ,                 double coszen ,  			   double parsun ,   			 double parsha ,
                 double sabvsun ,         double sabvsha ,                 double sabg   ,  			   double frl    ,			     double extkb  ,
                 double extkd   ,         double thermk  ,                 double fsno   ,      		   double sigf   ,			     double dz[15] ,
                 double z[MAXSNL+MAXSOILL],double zi[MAXSNL+MAXSOILL+1] ,    double &tlsun  ,          	   double &tlsha  ,			     double tss[MAXSNL+MAXSOILL],
                 double wice[MAXSNL+MAXSOILL],double wliq[MAXSNL+MAXSOILL],  double &ldew   ,              double &scv    ,		    	 double &snowdp ,
                 int  imelt[MAXSNL+MAXSOILL], double &taux    ,              double &tauy   ,              double &fsena  ,			     double &fevpa  ,
                 double &lfevpa  ,         double &fsenl   ,                 double &fevpl  ,              double &etr    ,  			 double &fseng  ,
                 double &fevpg   ,         double &olrg    ,                 double &fgrnd  ,              double rootr[MAXSOILL] ,      double &qseva  ,
                 double &qsdew   ,         double &qsubl   ,                 double &qfros  ,              double &sm     ,			     double &tref   ,
                 double &qref    ,         double &trad    ,                 double &rst    ,              double &assim  ,		    	 double &respc  ,
                 double &errore  ,         double &emis    ,                 double &z0ma   ,              double &zol    ,			     double &rib    ,
                 double &ustar   ,         double &qstar   ,                 double &tstar  ,              double &u10m   ,			     double &v10m   ,
                 double &f10m    ,         double &fm      ,                 double &fh     ,              double &fq);
    void ticktime (double dtime, int idate[3]);
    void tridia (int n, double a[], double b[], double c[], double r[],double u[]);

    void twostream (double   chil,   double ref[2][2],  double tran[2][2], double green,       double lai,
                    double   sai,    double coszen,     double albg[2][2], double albv[2][2],  double tranc[2][2],
                    double   &thermk, double &extkb,      double &extkd,      double ssun[2][2],  double ssha[2][2] );
    void vec2xy(int  lat_points,          int    lon_points,           int  numpatch,                   int ixy_patch[MAXPATCH],
                int  jxy_patch[MAXPATCH], double wtxy_patch[MAXPATCH], int itypwat[MAXPATCH],           int  nfcon,
                int  nforc,               int    nfldv,                double  fcon[MAXPATCH][NFCON],   double forcxy[MAXLONPOINT][MAXLATPOINT][NFORC],
                double fldv[MAXPATCH][NFLDV],                          double fldxy_r[MAXLONPOINT][MAXLATPOINT][NFLDV]);

    void WATER (int    itypwat ,           int    lb,                    int    nl_soil ,             double dtime  ,               double z[MAXSNL+MAXSOILL],
                double dz[MAXSNL+MAXSOILL],double zi[MAXSNL+MAXSOILL+1], double bsw[MAXSOILL] ,       double porsl[MAXSOILL]  ,     double phi0[MAXSOILL] ,
                double hksati[MAXSOILL] ,  double rootr[MAXSOILL]   ,    double tss[MAXSNL+MAXSOILL] ,double wliq[MAXSNL+MAXSOILL] ,double  wice[MAXSNL+MAXSOILL] ,
                double pg_rain,            double sm      ,              double etr    ,              double qseva  ,               double  qsdew   ,
                double qsubl  ,            double qfros   ,              double &rsur   ,              double &rnof   ,               double   wtfact  ,
                double pondmx ,            double ssi     ,              double wimp   ,              double smpmin  );
};
}
#endif // __LDAS_COMMONLANDMODEL_H
