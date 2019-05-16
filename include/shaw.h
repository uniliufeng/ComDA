#ifndef __LDAS_SHAW_H
#define __LDAS_SHAW_H

#include <iostream>
#include <fstream>
#include <sstream>
#include "model.h"
#include "time.hpp"
#include "exception.h"
#include "constant.h"

namespace ldas
{
extern "C"
{
    void goshaw_(int* JULIAN, int* HOUR, int* YEAR, int* NHRPDT, float* WDT, float* DT, int* INITAL, int* NC,
                 int* NSP, int* NR, int* NS, float* TOLER, int LEVEL[6], int* MZCINP, int* INPH2O, int* MWATRXT,
                 int LVLOUT[15], int* IVLCBC, int* ITMPBC, int* NPLANT, float PLTHGT[8], float PLTWGT[8],
                 float PLTLAI[8], float ROOTDP[8], float DCHAR[8], float TCCRIT[8], float RSTOM0[8], float RSTEXP[8],
                 float PLEAF0[8], float RLEAF0[8], float RROOT0[8], float PCANDT[8], float CANALB[8], float* CANMA,
                 float* CANMB, float* WCMAX, int LANGLE[8], int ITYPE[8], float ZC[11], float TCDT[11],
                 float VAPCDT[11], float WCANDT[11], float ZSP[100], float DZSP[100], float RHOSP[100],
                 float TSPDT[100], float DLWDT[100], int ICESPT[100], float WLAG[11], float* STORE, float* SNOWEX,
                 float* SNOTMP, float ZR[10], float RHOR[10], float TRDT[10], float VAPRDT[10], float GMCDT[10],
                 float* GMCMAX, float* RLOAD, float* ZRTHIK, float* COVER, float* ALBRES, float* RESCOF, float* DIRRES,
                 float ZS[50], float TSDT[50], float VLCDT[50], float VICDT[50], float MATDT[50], float CONCDT[50][10],
                 int ICESDT[50], float SALTDT[50][10], float* ALBDRY, float* ALBEXP, float DGRADE[10],
                 float SLTDIF[10], float ASALT[50], float DISPER[50], float* ZMSRF, float* ZHSRF, float* ZERSRF,
                 float* ZMSP, float* ZHSP, float* HEIGHT, float* POND, float* PONDMX, float* ALATUD, float* SLOPE,
                 float* ASPECT, float* HRNOON, float* CLOUDS, float* DECLIN, float* HAFDAY, float* SUN,
                 float* TMP, float* WIN, float* HUM, float* PRE, float* SNO,float* SOIT, float* VLC, float* SOILX );
    void input_(int* NC, int* NSP, int* NR, int* NS, float* TOLER, int LEVEL[6], int* MTSTEP, int* MZCINP,
                int* MPLTGRO, int* INPH2O, int* MWATRXT, int LVLOUT[15], int* IVLCBC, int* ITMPBC, int* NPLANT,
                float PLTHGT[8], float PLTWGT[8], float PLTLAI[8], float ROOTDP[8], float DCHAR[8], float TCCRIT[8],
                float RSTOM0[8], float RSTEXP[8], float PLEAF0[8], float RLEAF0[8], float RROOT0[8], float CANALB[8],
                float* CANMA, float* CANMB, float* WCMAX, int LANGLE[8], int ITYPE[8], float ZC[11], float WCANDT[11],
                float ZSP[100], float DZSP[100], float RHOSP[100], float TSPDT[100], float DLWDT[100], int ICESPT[100],
                float* SNOTMP, float GMCDT[10], float* GMCMAX, float* RLOAD, float* ZRTHIK, float* COVER, float* ALBRES,
                float* RESCOF, float ZS[50], float TSDT[50], float VLCDT[50], float VICDT[50], float MATDT[50],
                float CONCDT[50][10], int ICESDT[50], float SALTDT[50][10], float* ALBDRY, float* ALBEXP,
                float DGRADE[10], float SLTDIF[10], float ASALT[50], float DISPER[50], float* ZMSRF, float* ZHSRF,
                float* ZERSRF, float* ZMSP, float* ZHSP, float* HEIGHT, float* PONDMX, int* JSTART, int* YRSTAR,
                int* HRSTAR, int* JEND, int* YREND, int* NHRPDT, float* WDT, float* ALATUD, float* SLOPE, float* ASPECT,
                float* HRNOON);


    void dayinp_(int* JULIAN, int* YEAR, int* MAXJUL, int* HOUR, int* NHRPDT, int* MTSTEP, int* MPLTGRO,
                 int* MWATRXT, int* NS, float* ALATUD, float* HRNOON, float SUNHOR[24], float TMPDAY[24],
                 float WINDAY[24], float HUMDAY[24], float PRECIP[24], float SNODEN[24], float SOITMP[24],
                 float VLCDAY[24], float SOILXT[24][50], int* NPLANT, float PLTHGT[8], float DCHAR[8],
                 float PLTWGT[8], float PLTLAI[8], float ROOTDP[8]);

    void cloudy_(float* CLOUDS, float* ALATUD, float* DECLIN, float* HAFDAY, float SUNHOR[24],
                 int* JULIAN, int* NHRPDT);

}
namespace shaw
{
#define LF 335000.0
#define LV 2500000.0
#define	LS 2835000.0
#define	G 9.81
#define	UGAS 8.3143
const float RHOL= 1000.0;
#define	RHOM 2650.0
const float	RHOI= 920.0;
#define	RHOA 1.25
#define	RHOOM 1300.0
#define	CL 4200.0
#define	CI 2100.0
#define	CM 900.0
#define	COM 1920.0
#define	CA 1006.0
#define	CV 1860.0
#define	CR 1900.0
#define	VONKRM 0.4
#define	VDIFF 0.0000212
#define	P0 101300.0
#define	TKL 0.57
#define	TKI 2.2
#define	TKA 0.025
#define	TKR 0.05
#define	TKSA 8.8
#define	TKSI 2.92
#define	TKCL 2.92
#define	TKOM 0.25
}
class SHAW : public BaseModel
{
public:
    /** Default constructor */
    SHAW();
    /** Default destructor */
    virtual ~SHAW();
    /** Copy constructor
     *  \param other Object to copy from
     */
    SHAW(const SHAW& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    SHAW& operator=(const SHAW& other);
    virtual void run();
    virtual void step();
    virtual void output();
    void input();
    void dayinput();
    void cloudy();
protected:
    /*parameter
    config
    forcing*/
private:
    // model start and end time : year julian day and hour
    int NHRPDT;                                // number of hours per time steps
    int MAXJUL;                                // maximum number of day in current year
    int JSTART, YRSTAR, HRSTAR, YREND, JEND;   // start and end time of model
    int JULIAN;                                // current julian day of model
    int HOUR;                                  // current hour of model
    int YEAR;                                  // current year of model

    // state variable of canopy
    float ZC[11];                         // distance of node i from top of canopy, ZC(1) must 0.0 m
    float TCDT[11];                       // canopy temperature.
    float VAPCDT[11];                     // canpy vapor pressure
    float WCANDT[10];                     // initial moisture content for dead plant material
    int NC;                               // number of node in canopy


    // state variables of snow
    float ZSP[100];                       // depth of node in snow layer (m) ??
    float DZSP[100];                      // thickness of snow layer (m)
    float RHOSP[100];                     // density of snow layer
    float TSPDT[100];                     // snow temperature (Celsius)
    float DLWDT[100];                     // depth of liquid water in snow layer (m)
    float WLAG[11];                       // snow lag ???
    int ICESPT[100];                      // flag indicating ice in snow (0 or 1)
    int NSP;                              // number of node in snow


    // state variables of residue
    float ZR[10];                         // depth of node in residue
    float RHOR[10];                       // density of residue
    float TRDT[10];                       // residue temperature
    float VAPRDT[10];                     // residue vapor pressure
    float GMCDT[10];                      // initial gravimetric water content of residue
    int NR;                               // number of node in residue


    // state variables of soil
    float ZS[50];                         // depth of node in soil
    float TSDT[50];                       // soil temperature (Celsius)
    float MATDT[50];                      // soil Matrix potential
    float VLCDT[50];                      // soil moisture (liquid+ice)
    float VICDT[50];                      // soil ice content
    int ICESDT[50];                       // flag indicating ice in soil (0 or 1)
    int NS;                               // number of node in soil


    // state variables of solute
    float CONCDT[50][10];                 // concentration of solute of soil in layer
    float SALTDT[50][10];                 // of moles of solute per kg of soil in layer
    float CONCDTT[10][50];                // REMARK: CONCDT and SALTDT transfer and receive value between C++ and Fortran
    float SALTDTT[10][50];                // CONCDTT is the transpose of CONCDT, we should visit the CONCDTT in C++
    float DGRADE[10];                     // degradation of solute
    float SLTDIF[10];                     // diffusion coefficient for solute at 0 Celsius
    float ASALT[50];                      // molecular diffusion parameter for solute in soil layer
    float DISPER[50];                     // parameter for hydrodynamic dispersion coefficient (m)
    int NSALT;                            // number of solute

    //state variables of plant
    float PLTHGT[8];                      // height of plant species j (m)
    float PLTWGT[8];                      // dry biomass of plant species j (kg/m3)
    float PLTLAI[8];                      // leaf area index of plant species j
    float ROOTDP[8];                      // effective rooting depth of plant species j (m)
    float DCHAR[8];                       // characteristic dimension of leaves of plant species j (cm)
    float TCCRIT[8];                      // temperature above which plant j will transpire
    float RSTOM0[8];                      // stomatal resistance of plant species j
    float RSTEXP[8];                      // empirical exponent relating actual stomatal resistance
    float PLEAF0[8];                      // critical leaf water potential for plant species j
    float RLEAF0[8];                      // resistence of leaves for plant species j
    float RROOT0[8];                      // resistance of roots for plant species j
    float PCANDT[8];                      // precipitation interception for plant species j
    float CANALB[8];                      // albedo of plant species j
    int NPLANT;                           // number of plant species


    // soil sink
    float SOILXT[24][50];

    // forcing data and lower boundary condition
    float SUNHOR[24];                     // total solar radiation measured on a H surface (W/m2)
    float TMPDAY[24];                     // air temperature (Celsius)
    float WINDAY[24];                     // windspeed
    float HUMDAY[24];                     // relative humidity
    float PRECIP[24];                     // precipitation
    float SNODEN[24];                     // density of newly fallen snow (g/cm3)
    float SOITMP[24];                     // temperature of soil lowest boundary
    float VLCDAY[24];                     // moisture of soil lowest boundary
    float SUN;
    float TMP;
    float WIN;
    float HUM;
    float PRE;
    float SNO;
    float SOIT;
    float VLC;


    // other constant parameter
    int LEVEL[6];                                          // debug control
    int LVLOUT[15];                                        // output control
    int	LANGLE[8];                                         // parameter specifying leaf angle for plant species j
    int ITYPE[8];                                          // Parameter specifying plant type for plant species j
    float CANMA, CANMB;                                    // parameter of water potential of dead plant material
    float SNOTMP;                                          // maximum temperature at which precipitation is snow
    float RLOAD, ZRTHIK, COVER, ALBRES, RESCOF;            // parameter of residue
    float ALBDRY, ALBEXP;                                  // parameter of soil albedo
    float ZMSP, ZHSP, ZMSRF, ZHSRF, ZERSRF, HEIGHT;        // parameter of boundary layer
    float PONDMX;                                          // maximum depth of ponding for rainfall or snowmelt (m)
    float ALATUD, SLOPE, ASPECT, HRNOON, DECLIN, HAFDAY;   // parameter related to orientation
    float CLOUDS;                                          // fraction of cloud cover
    float DT, WDT;                                         //
    float TOLER;                                           // error tolerance for convergence criteria


    // other varing parameter
    float WCMAX;                                           //
    float GMCMAX;                                          //
    float STORE;                                           //
    float SNOWEX;                                          //
    float DIRRES;                                          //
    float POND;                                            //

    // flag
    int MTSTEP, MZCINP, MPLTGRO, INPH2O, MWATRXT, IVLCBC, ITMPBC;
    int INITAL;

    // others
    float SOILX;

    // soil parameters
    float SAND[50];
    float CLAY[50];
    float RHOB[50];

};
class SHAWConfig: public ModelConfig
{
public:
    SHAWConfig() {}
    virtual ~SHAWConfig() {}
    virtual void open(const char* fn);

    int mtstep,inph2o,mwatrxt;
    std::string site_path,weather_path,moisture_path,temperature_path,soilsink_path;
    int output_frequency,comparison_flag,temperature_frequency,watercontent_frequency;
    int matricpotential_frequency,energy_frequency,waterbalance_frequency,rootwater_frequency;
    int frostdepth_frequency,salt_frequency,soilsolution_frequency,updatescreen_frequency;

    std::string soilsolution_path,salt_path,frostdepth_path,rootwater_path,waterflow_path;
    std::string waterbalance_path,energyflux_path,waterpotential_path,watercontent_path;
    std::string simulated_temperature_path,profile_path,outputgeneral_path;

private:
};
class SHAWParameter
{
public:

};
}
#endif // __LDAS_SHAW_H
