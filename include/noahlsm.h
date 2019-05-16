#ifndef __LDAS_NOAHLSM_H
#define __LDAS_NOAHLSM_H

#include "model.h"
#include "time.hpp"
#include "exception.h"
#include "constant.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>

namespace ldas
{
extern "C"
{
    void __module_sf_noahlsm_MOD_sflx(float* FFROZP,int* ICE,int* ISURBAN,float* DT,float* ZLVL,int* NSOIL,float* SLDPTH,
               int* LOCAL,
               char* LLANDUSE, char* LSOIL,
               float* LWDN,float* SOLDN,float* SOLNET,float* SFCPRS,float* PRCP,float* SFCTMP,float* Q2,float* SFCSPD,
               float* COSZ,float* PRCPRAIN,float* SOLARDIRECT,
               float* TH2,float* Q2SAT,float* DQSDT2,
               int* VEGTYP,int* SOILTYP,int* SLOPETYP,float* SHDFAC,float* SHDMIN,float* SHDMAX,
               float* ALB, float* SNOALB,float* TBOT, float* Z0BRD, float* Z0, float* EMISSI, float* EMBRD,
               float* CMC,float* T1,float STC[],float SMC[],float SH2O[],float* SNOWH,float* SNEQV,float* ALBEDO,float* CH,float* CM,
               /** ----------------------------------------------------------------------
                OUTPUTS, DIAGNOSTICS, PARAMETERS BELOW GENERALLY NOT NECESSARY WHEN
                COUPLED WITH E.G. A NWP MODEL (SUCH AS THE NOAA/NWS/NCEP MESOSCALE ETA
                MODEL).  OTHER APPLICATIONS MAY REQUIRE DIFFERENT OUTPUT VARIABLES.
                ----------------------------------------------------------------------
               */
               float* ETA,float* SHEAT, float* ETA_KINEMATIC,float* FDOWN,
               float* EC,float* EDIR,float ET[],float* ETT,float* ESNOW,float* DRIP,float* DEW,
               float* BETA,float* ETP,float* SSOIL,
               float* FLX1,float* FLX2,float* FLX3,
               float* SNOMLT,float* SNCOVR,
               float* RUNOFF1,float* RUNOFF2,float* RUNOFF3,
               float* RC,float* PC,float* RSMIN,float* XLAI,float* RCS,float* RCT,float* RCQ,float* RCSOIL,
               float* SOILW,float* SOILM,float* Q1,float SMAV[],
               int* RDLAI2D,int* USEMONALB,
               float* SNOTIME1,
               float* RIBB,
               float* SMCWLT,float* SMCDRY,float* SMCREF,float* SMCMAX,int* NROOT);
    /*void __module_ascii_io_MOD_open_forcing_file(int* iunit, char* forcing_filename, char infotext[], int* nsoil, char startdate[], char enddate[], bool* loop_for_a_while,
                            float* latitude, float* longitude,
                            int* forcing_timestep, int* noahlsm_timestep, int* ice, float* t1, float stc[], float smc[], float sh2o[], float sldpth[], float* cmc, float* snowh, float* sneqv, float* tbot,
                            int* vegtyp, int* soiltyp, int* slopetyp, float* snoalb, float* zlvl, float* zlvl_wind, float albedo_monthly[], float shdfac_monthly[],
                            float z0brd_monthly[], float lai_monthly[], bool* use_urban_module, int* isurban, float* shdmin, float* shdmax, bool* usemonalb, bool* rdlai2d, char* llanduse,
                            int* sfcdif_option, int* iz0tlnd,int a_len,int b_len);
    void __module_ascii_io_MOD_read_forcing_text(int* iunit, char* nowdate, int* forcing_timestep, float* wspd, float* u, float* v,
                            float* sfctmp, float* spechumd, float* sfcprs, float* swrad, float* lwrad, float* pcprate, int* ierr);*/
void __module_sf_noahlsm_MOD_sfcdif_off(float* zlm,float* z0,float* thz0,float* thlm,float* sfcspd,float* czil,float* akms,float* akhs,
                    int* vegtyp, int* isurban, int* iz0tlnd);
                            extern int __module_sf_noahlsm_MOD_lucats;
							extern int __module_sf_noahlsm_MOD_slcats;
							extern int __module_sf_noahlsm_MOD_slpcats;
							extern int __module_sf_noahlsm_MOD_bare,__module_sf_noahlsm_MOD_natural;
							extern int __module_sf_noahlsm_MOD_nrotbl[50];
							extern float __module_sf_noahlsm_MOD_shdtbl[50],__module_sf_noahlsm_MOD_rstbl[50],__module_sf_noahlsm_MOD_rgltbl[50],__module_sf_noahlsm_MOD_hstbl[50],__module_sf_noahlsm_MOD_snuptbl[50],
          __module_sf_noahlsm_MOD_maxalb[50],__module_sf_noahlsm_MOD_laimintbl[50],__module_sf_noahlsm_MOD_laimaxtbl[50],__module_sf_noahlsm_MOD_emissmintbl[50],
          __module_sf_noahlsm_MOD_emissmaxtbl[50],__module_sf_noahlsm_MOD_albedomintbl[50],__module_sf_noahlsm_MOD_albedomaxtbl[50],__module_sf_noahlsm_MOD_z0mintbl[50],__module_sf_noahlsm_MOD_z0maxtbl[50];
		  extern float __module_sf_noahlsm_MOD_bb[30],__module_sf_noahlsm_MOD_drysmc[30],__module_sf_noahlsm_MOD_f11[30],__module_sf_noahlsm_MOD_maxsmc[30],__module_sf_noahlsm_MOD_refsmc[30],__module_sf_noahlsm_MOD_satpsi[30],__module_sf_noahlsm_MOD_satdk[30],__module_sf_noahlsm_MOD_satdw[30],__module_sf_noahlsm_MOD_wltsmc[30],__module_sf_noahlsm_MOD_qtz[30];
		  extern float __module_sf_noahlsm_MOD_slope_data[30];
    extern float __module_sf_noahlsm_MOD_topt_data,__module_sf_noahlsm_MOD_cmcmax_data,__module_sf_noahlsm_MOD_cfactr_data,__module_sf_noahlsm_MOD_rsmax_data;
    extern float __module_sf_noahlsm_MOD_sbeta_data,__module_sf_noahlsm_MOD_fxexp_data,__module_sf_noahlsm_MOD_csoil_data,__module_sf_noahlsm_MOD_salp_data,__module_sf_noahlsm_MOD_refdk_data,__module_sf_noahlsm_MOD_refkdt_data,__module_sf_noahlsm_MOD_frzk_data,__module_sf_noahlsm_MOD_zbot_data,__module_sf_noahlsm_MOD_czil_data,__module_sf_noahlsm_MOD_smlow_data,__module_sf_noahlsm_MOD_smhigh_data,__module_sf_noahlsm_MOD_lvcoef_data;

}
namespace noah
{
const float STBOLT=5.67051E-8;
const int kztm=10001;
const float badval=-1.E36;
const float g = 9.81;
#if ( NMM_CORE == 1 )
const float r_d          = 287.04;
const float cp           = 1004.6;
#else
const float  r_d          = 287.;
const float  cp           = 7.*r_d/2.;
#endif
const float r_v          = 461.6;
const float rvovrd       = r_v/r_d;
const float p608=rvovrd-1.;
const float p1000mb=100000;
const float p0=p1000mb;
const float rcp= r_d/cp;
}
class Noahlsm : public BaseModel
{
public:
    /** Default constructor */
    Noahlsm();
    /** Default destructor */
    virtual ~Noahlsm();
    /** Copy constructor
     *  \param other Object to copy from
     */
    Noahlsm(const Noahlsm& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    Noahlsm& operator=(const Noahlsm& other);
    virtual void config(const std::string& fn);
    virtual void run();
    void init();

    /** calculate temperature
    * \param t1 skin temperature [k], input
    * \param sfctmp air temperature [k] at level ZLVL, input
    * \param sfcprs atmospheric pressure [pa] at level ZLVL, input
    * \param zlvl height [m] whre atmospheric fields are available, input
    * \param q2, specific huminity [kg/kg] at level ZLVL, input
    * \param th2, potential temperature (considering the reference pressure to be at the surface, output
    * \param t1v, virtual skin temperature, output
    * \param th2v, virtual potential temperature at ZLVL, output
    * \param rho, densith, output
    */
    void caltmp(float t1, float sfctmp, float sfcprs, float zlvl, float q2,
                float &th2, float &t1v, float &th2v, float &rho);

    /** calculate humidity
    * \param sfctmp air temperature [k] at level ZLVL, input
    * \param sfcprs atmospheric pressure [pa] at level ZLVL, input
    * \param q2sat saturated specific humidity, output
    * \param dqsdt2 slope of saturated specific huminity curve, output
    */
    void calhum(float sfctmp, float sfcprs, float &q2sat, float &dqsdt2);

    /** Given a set of 12 values, taken to be valid on the fifteenth of each month (Jan through Dec)
    * \param a12, array of 12 month
    * \param t Time variable
    * \return a value valid for the day given in <nowdate>, as an interpolation from the 12 monthly values.
    */
    float month_d(const float a12[12], const Time& t) const;
    void myjsfcinit();
    void sfcdif_myj(float zsl, float z0, float z0base, float sfcprs, float tz0, float tlow, float qz0,
                    float qlow, float sfcspd, float czil, float& rib, float& akms, float& akhs, int vegtyp, int isurban, int iz0tlnd );
    void modi_buf(char *buf);
    void  soil_veg_gen_parm();
    Time time_convert(const char startdate[12]) const;
    void sfcdif_off(float zlm,float z0,float thz0,float thlm,float sfcspd,float czil,float& akms,float& akhs,
                    int vegtyp, int isurban, int iz0tlnd);
    std::string convert_time(const Time& t) const;
    /** open and read config and forcing*/
    void open_forcing_file();
    void read_forcing_text();
protected:
    char forcing_filename[256];
    char nowdate[12]; //The date of each time step, ( YYYYMMDDHHmm ), updated in each step of the main loop TIMELOOP
    int ktime; //A counter for the timesteps in the main loop TIMELOOP.

    //
    //Useful data attributes describing the data in the initial/forcing conditions file
    //
    char infotext[4096];  // Character string returned by subroutine OPEN_FORCING_FILE, giving some possibly useful information for the user.
    float  latitude;        /// Latitude of the point ( Degrees North )
    float  longitude;        /// Longitude of the point ( Degrees East )
    bool   loop_for_a_while; /// Whether to loop the same year ad infinitum
    char startdate[12];        /// Starting date of the data ( YYYYMMDDHHmm )
    char enddate[12];          /// Ending date of the data ( YYYYMMDDHHmm )
    int forcing_timestep; /// The time interval ( seconds ) of the data in the forcing file
    int noahlsm_timestep; /// The timestep ( seconds ) to use when integrating the Noah LSM
    float albedo_monthly[12];   /// Monthly values of background (i.e., snow-free) albedo ( Fraction [0.0-1.0] )
    float shdfac_monthly[12];   /// Monthly values for green vegetation fraction ( Fraction [0.0-1.0] )
    float z0brd_monthly[12];    /// Monthly values for background (i.e., snow-free) roughness length ( m )
    float lai_monthly[12];      /// Monthly values for Leaf Area Index ( dimensionless )

    ///
    /// Various arguments to subroutine SFLX:
    ///

    float ffrozp;     /// Fraction of precip which is frozen (0.0 - 1.0).
    int ice;        /// Flag for sea-ice (1) or land (0).
    int isurban;    /// Vegetation category for urban land class.
    float dt;         /// Time step (seconds).
    float zlvl;       /// Height at which atmospheric forcing variables are taken to be valid (m)
    float zlvl_wind;  /// Height at which the wind forcing variable is taken to be valid (m)
    float* sldpth; /// Thicknesses of each soil level
    int nsoil;      /// Number of soil levels.
    bool local;      /// Not used in SFLX
    char llanduse[256];  /// Land-use dataset.  Valid values are :
    ///                               /// "USGS" (USGS 24/27 category dataset) and
    ///                               /// "MODIFIED_IGBP_MODIS_NOAH" (MODIS 20-category dataset)
    char lsoil[256];     /// Soil-category dateset.  Only "STAS" (STATSGO dataset) supported.
    float   lwdn;       /// Downward longwave radiation flux at surface (W m-2) [Forcing]
    float   soldn;      /// Downward shortwave radiation flux at surface (W m-2) [Forcing]
    float   solnet;     /// Net downward shortwave radiation flux at the surface (W m-2)
    float   sfcprs;     /// Surface atmospheric pressure (Pa) [Forcing]
    float   prcp;       /// Precipitation rate (kg m-2 s-1) [Forcing]
    float   sfctmp;     /// Air temperature (K) [Forcing]
    float   q2;         /// Surface specific humidity (kg kg-1) [Forcing]
    float   sfcspd;     /// Surface wind speed (m s-1) [Forcing]
    float   sfcu;       /// West-to-east component of the surface wind (m s-1)
    float   sfcv;       /// South-to-north component of the surface wind (m s-1)
    float   cosz;       /// Unused if we're not using urban canopy model.
    float   prcprain;   /// Unused.
    float   solardirect;/// Unused.
    float   th2;        /// Potential temperature at level ZLVL (K)
    float   t1v;        /// Virtual skin temperature (K).  Used in SFCDIF_off for computing CM and CH, but not passed to SFLX
    float   th2v;       /// Virtual potential temperature at level ZLVL (K).  Used in SFCDIF_off
    ///                     /// for computing CM and CH, but not passed to SFLX
    float   rho;        /// Air density (dummy value output from CALTMP, not passed to SFLX).
    float   q2sat;      /// Saturated specific humidity (kg kg-1)
    float   dqsdt2;     /// Slope of the Saturated specific humidity curve W.R.T. Temperature.
    int vegtyp;     /// Vegetation category.
    int soiltyp;    /// Soil category.
    int slopetyp;   /// Slope category.
    float   shdfac;     /// Shade factor (0.0-1.0).
    float   shdmin;     /// Minimum shade factor (0.0-1.0).
    float   shdmax;     /// Maximum shade factor (0.0-1.0).
    float   alb;        /// Background snow-free albedo (0.0-1.0).
    float   snoalb;     /// Maximum snow albedo over deep snow (0.0-1.0)
    float   tbot;       /// Deep-soil time-invariant temperature (K).  Representing sort of a mean annual air temperature.
    float   z0brd;      /// Background Z0 value (m).
    float   z0;         /// Roughness length (m)
    float   emissi ;    /// Surface emissivity (0.0 - 1.0).  This includes the snow-cover effect.
    float   embrd;      /// Background value (i.e., not including snow-cover effect) of surface emissivity (0.0 - 1.0)
    float   cmc ;       /// Canopy moisture content (kg m-2)
    float   t1;         /// Skin temperature (K)
    float* stc;  /// Soil temperature (K)
    float* smc;  /// Total soil moisture content (m3 m-3)
    float* sh2o; /// Liquid soil moisture content (m3 m-3)
    float* et;   /// Plant transpiration from each soil level.
    float* smav; /// Soil Moisture Availability at each level, fraction between
    /// SMCWLT (SMAV=0.0) and SMCMAX (SMAV=1.0)
    float   snowh;      /// Physical snow depth.
    float   sneqv;      /// Water equivalent of accumulated snow depth (m).
    float   albedo;     /// Surface albedo including possible snow-cover effect.  This is set in SFLX,
    ///                     /// overriding any value given; it should perhaps be INTENT(OUT) from SFLX.
    float   ch ;        /// Exchange coefficient for head and moisture (m s-1).  An initial value is needed for SFCDIF_off.
    float   cm ;        /// Exchange coefficient for momentum (m s-1).  An initial value is needed for SFCDIF_off.
    float   eta ;       /// Latent heat flux (evapotranspiration) ( W m{-2} )
    float   sheat ;     /// Sensible heat flux ( W m{-2} )
    float   etakin ;    /// Latent heat flux (evapotranspiration) ( kg m{-2} s{-1} )
    float   fdown;      /// Radiation forcing at the surface ( W m{-2} )
    float   ec  ;       /// Latent heat flux component: canopy water evaporation ( W m{-2} )
    float   edir ;      /// Latent heat flux component: direct soil evaporation ( W m{-2} )
    float   ett  ;      /// Latent heat flux component: total plant transpiration ( W m{-2} )
    float   esnow ;     /// Latent heat flux component: sublimation from (or deposition to) snowpack ( W m{-2} )
    float   drip  ;     /// Precipitation or dew falling through canopy, in excess of canopy holding capacity ( m )
    float   dew  ;      /// Dewfall (or frostfall for T<273.15) ( m )
    float   beta ;      /// Ratio of actual to potential evapotranspiration ( Fraction [0.0-1.0] )
    float   etp  ;      /// Potential evapotranspiration ( W m{-2} )
    float   ssoil ;     /// Soil heat flux ( W m{-2} )
    float   flx1 ;      /// Latent heat flux from precipitation accumulating as snow ( W m{-2} )
    float   flx2  ;     /// Latent heat flux from freezing rain converting to ice ( W m{-2} )
    float   flx3  ;     /// Latent heat flux from melting snow ( W m{-2} )
    float   snomlt ;    /// Snow melt water ( m )
    float   sncovr ;   /// Fractional snow cover ( Fraction [0.0-1.0] )
    float   runoff1 ;   /// Surface runoff, not infiltrating the soil ( m s{-1} )
    float   runoff2 ;   /// Subsurface runoff, drainage out the bottom of the last soil layer ( m s{-1} )
    float   runoff3 ;   /// Internal soil layer runoff ( m s{-1} )
    float   rc    ;     /// Canopy resistance ( s m{-1} )
    float   pc    ;     /// Plant coefficient, where PC * ETP = ETA ( Fraction [0.0-1.0] )
    float   rsmin  ;    /// Minimum canopy resistance ( s m{-1} )
    float   xlai  ;     /// Leaf area index ( dimensionless )
    float   rcs   ;     /// Incoming solar RC factor ( dimensionless )
    float   rct   ;     /// Air temperature RC factor ( dimensionless )
    float   rcq   ;     /// Atmospheric water vapor deficit RC factor ( dimensionless )
    float   rcsoil;     /// Soil moisture RC factor ( dimensionless )
    float   soilw ;     /// Available soil moisture in the root zone ( Fraction [SMCWLT-SMCMAX] )
    float   soilm ;     /// Total soil column moisture content, frozen and unfrozen ( m )
    float   q1    ;     /// Effective mixing ratio at the surface ( kg kg{-1} )
    bool rdlai2d ;   /// If RDLAI2D == .TRUE., then the XLAI value that we pass to SFLX will be used.
    ///                     /// If RDLAI2d == .FALSE., then XLAI will be computed within SFLX, from table
    ///                     /// minimum and maximum values in VEGPARM.TBL, and the current Green Vegetation Fraction.
    bool usemonalb;  /// If USEMONALB == .TRUE., then the ALB value passed to SFLX will be used as the background
    ///                     /// snow-free albedo term.  If USEMONALB == .FALSE., then ALB will be computed within SFLX
    ///                     /// from minimum and maximum values in VEGPARM.TBL, and the current Green Vegetation Fraction.
    float   snotime1;   /// Age of the snow on the ground.
    float   ribb  ;     /// Bulk Richardson number used to limit the dew/frost.
    float   smcwlt ;    /// Wilting point ( m{3} m{-3} )
    float   smcdry ;    /// Dry soil moisture threshold where direct evaporation from the top layer ends ( m{3} m{-3} )
    float   smcref ;    /// Soil moisture threshold where transpiration begins to stress ( m{3} m{-3} )
    float   smcmax ;    /// Porosity, i.e., saturated value of soil moisture ( m{3} m{-3} )
    int nroot ;     /// Number of root layers ( count )

    int z0tlnd ;   /// Option to turn on (IZ0TLND=1) or off (IZ0TLND=0) the vegetation-category-dependent
    /// calculation of the Zilitinkivich coefficient CZIL in the SFCDIF subroutines.

    int sfcdif_option ;/// Option to use previous (SFCDIF_OPTION=0) or updated (SFCDIF_OPTION=1) version of
    /// SFCDIF subroutine.

    ///
    /// Some diagnostics computed from the output of subroutine SFLX
    ///
    float qfx ;      /// Evapotranspiration ( W m{-2} )  the sum of 1) direct evaporation
    ///                 /// from soil; 2) evaporation from canopy; 3) total plant transpiration;
    ///                 /// and 4) evaporation from snowpack.  Mostly, this should be the
    ///                 /// same as ETA

    float res ;      /// Residual of the surface energy balance equation ( W m{-2} )
    float fup ;      /// Upward longwave radiation flux from the surface ( W m{-2} )
    float f ;        /// Incoming shortwave and longwave radiation flux  ( W m{-2} )

    ///
    /// Miscellaneous declarations
    ///
    int           ierr   ;          /// Error flag returned from read routines.
    int iunit ;      /// Fortran unit number for reading initial/forcing conditions file.
    bool           use_urban_module; /// Flag, set to TRUE in the initial/forcing conditions file, if the
    ///                                      /// user wants to use the urban canopy model.  Since this code does not
    ///                                      /// include the urban canopy model, a TRUE value of this flag will simply
    ///                                      /// stop the execution.
    //float, external    month_d ;         /// External function (follows this main program):  given an array (dimension 12)
    ///                                      /// representing monthly values for some parameter, return a value for
    ///                                      /// a specified date.
    float               czil  ;          /// Zilitinkevich constant, read from GENPARM.TBL and used to compute surface
    ///                                      /// exchange coefficients
    float               longwave ;       /// Longwave radiation as read from the forcing data, which is immediately
    ///                                      /// adjusted (by the emissivity factor) to set variable LWDN.

private:
    float dzeta2;
    float ztmax2;
    float* psih2;
    float* psim2;
    std::string lutype,sltype;

    int iz0tlnd;
    Time endtime;
    std::ifstream forcing;
private:
    float psphs(const float yy) const;
    float psphu(const float xx) const;
    float pspms(const float yy) const;
    float pspmu(const float xx) const;
    float pslhs(const float zz) const;
    float pslhu(const float zz) const;
    float pslms(const float zz) const;
    float pslmu(const float zz) const;

};
}
#endif // __LDAS_NOAHLSM_H
