#ifndef __LDAS_VIC_H
#define __LDAS_VIC_H

#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>

#include "model.h"
#include "constant.h"

#include "vicNl.h"

namespace ldas
{
namespace vic
{

/****************************
 **                         **
 **  STRUCTURE DEFINITIONS  **
 **                         **
 ****************************/
typedef struct
{
    int ndays;             /* number of days of data in input file */
    int indewpt;           /* input dewpoint temperature flag (0=NO, 1=YES) */
    int outhum;            /* output humidity flag            (0=VPD, 1=VP) */
    int inyear;            /* input year flag                 (0=NO, 1=YES) */
} control_struct;

typedef struct
{
    double base_elev;      /* base elevation, meters */
    double base_isoh;      /* base annual precip isohyet, cm */
    double site_lat;       /* site latitude, dec. degrees (- for south) */
    double site_elev;      /* site elevation, meters */
    double site_slp;       /* site slope, degrees */
    double site_asp;       /* site aspect, degrees */
    double site_isoh;      /* site annual precip isohyet, cm */
    double site_ehoriz;    /* site east horizon, degrees */
    double site_whoriz;    /* site west horizon, degrees */
    double tmax_lr;        /* maximum temperature lapse rate, deg C/1000m */
    double tmin_lr;        /* minimum temperature lapse rate, deg C/1000m */
} parameter_struct;

typedef struct
{
    int *year;             /* array of year values */
    int *yday;             /* array of yearday values */
    double *tmax;          /* array of base maximum temperature values */
    double *tmin;          /* array of base minimum temperature values */
    double *prcp;          /* array of base daily precipitation values */
    double *tdew;          /* array of base dewpoint temperature values */
    double *s_tmax;        /* array of site tmax values */
    double *s_tmin;        /* array of site tmin values */
    double *s_tday;        /* array of site daylight temperature values */
    double *s_prcp;        /* array of site prcp values */
    double *s_hum;         /* array of site humidity values (VPD or VP, Pa) */
    double *s_srad;        /* array of site shortwave radiation values */
    double *s_dayl;        /* array of site daylength values */
    /* start vic_change */
    double *s_tskc;	 /* array of cloudiness values */
    /* end vic_change */
} data_struct;

}
class Vic:public BaseModel
{
public:
    /** Default constructor */
    Vic();
    /** Default destructor */
    virtual ~Vic();
    /** Copy constructor
     *  \param other Object to copy from
     */
    Vic(const Vic& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    Vic& operator=(const Vic& other);
    virtual void run();
    ///结果输出
    void output();
    void init();
    void config(const char*);
protected:
    /** Variable Declarations **/
    char                     NEWCELL;
    char                     LASTREC;
    char                     MODEL_DONE;
    char                    *init_STILL_STORM;
    char                     ErrStr[MAXSTRING];
    int                      rec, i, j;
    int                      veg;
    int                      dist;
    int                      band;
    int                      Ndist;
    int                      Nveg_type;
    int                      cellnum;
    int                      index;
    int                     *init_DRY_TIME;
    int                      RUN_MODEL;
    int                      Ncells;
    int                      cell_cnt;
    int                      startrec;
    int                      ErrorFlag;
    float                    mu;
    double                   storage;
    double                   veg_fract;
    double                   band_fract;
    double                   Clake;
    dmy_struct              *dmy;
    atmos_data_struct       *atmos;
    veg_con_struct          *veg_con;
    soil_con_struct          soil_con;
    dist_prcp_struct         prcp; /* stores information about distributed  precipitation */
    filenames_struct         filenames;
    filep_struct             filep;
    lake_con_struct          lake_con;
    out_data_file_struct     *out_data_files;
    out_data_struct          *out_data;
    save_data_struct         save_data;
private:
    //from global.h
    char *version;
    char *optstring;
    char   ref_veg_over[4];
    double ref_veg_rarc[4];
    double ref_veg_rmin[4];
    double ref_veg_lai[4];
    double ref_veg_albedo[4];
    double ref_veg_rough[4];
    double ref_veg_displ[4];
    double ref_veg_wind_h[4];
    double ref_veg_RGL[4];
    double ref_veg_rad_atten[4];
    double ref_veg_wind_atten[4];
    double ref_veg_trunk_ratio[4];
    char ref_veg_ref_crop[6];
    //======
    int              NF, NR;
#if QUICK_FS
    double   temps[8];
#endif

#if LINK_DEBUG
    debug_struct debug;
#endif // LINK_DEBUG


    Error_struct Error;

public:
    param_set_struct param_set;
    global_param_struct global_param;
    veg_lib_struct *veg_lib;
    option_struct options;
protected:
    double advected_sensible_heat(double, double, double, double, double);
    void alloc_atmos(int, atmos_data_struct **);
    double arno_evap(layer_data_struct *, layer_data_struct *, double, double,
                     double, double, double, double, double, double, double,
                     double, double, double, double, double,
#if SPATIAL_FROST
                     double, double *);
#else
                     double);
#endif // SPATIAL_FROST
    unsigned char average_moisture_for_storm(double *, double *, double, double);

    int   CalcAerodynamic(char, double, double, double, double, double,
                          double *, double *, double *, double *, double *);
    void   calc_cloud_cover_fraction(atmos_data_struct *, dmy_struct *, int,
                                     int, int, double *);
    void   calc_energy_balance_error(int, double, double, double, double, double);
#if OUTPUT_FORCE_STATS
    void   calc_forcing_stats(int, atmos_data_struct *);
#endif // OUTPUT_FORCE_STATS
    void   calc_longwave(double *, double, double, double);
    void   calc_netlongwave(double *, double, double, double);
    double calc_netshort(double, int, double, double *);
    double calc_rainonly(double,double,double,double,double);
    double calc_rc(double,double,float,double,double,double,double,char);
    void   calc_root_fractions(veg_con_struct *, soil_con_struct *);
    double calc_snow_coverage(int *, double, double, double, double, double,
                              double, double, double *, double *, double *,
                              double *, double *);
    double calc_snow_ground_flux(int, int, int, int, double, double, double,
                                 double, double, double *, double *, double *,
                                 double *, energy_bal_struct *,
                                 snow_data_struct *, layer_data_struct *,
                                 layer_data_struct *, soil_con_struct *, char *);
#if QUICK_FS
    int    calc_soil_thermal_fluxes(int, double *, double *, char *, int *, double *, double *,
                                    double *, double *, double *,double *,
                                    double *, double *, double *,
                                    double *, double *, double *, double ***, int, int, int, int);
#else
    int    calc_soil_thermal_fluxes(int, double *, double *, char *, int *, double *, double *,
                                    double *, double *, double *,double *,
                                    double *, double *, double *,
                                    double *, double *, double *,
#if EXCESS_ICE
                                    double *, double *,
#endif // EXCESS_ICE
                                    int, int, int, int);
#endif // QUICK_FS
    double CalcSnowPackEnergyBalance(double Tsurf, ...);
    double CalcBlowingSnow(double, double, int, double, double, double, double,
                           double, double, double, double, double, float,
                           float, double, int, int, float, double, double, double *);
    double calc_atmos_energy_bal(double, double, double, double, double, double,
                                 double, double, double, double, double, double,
                                 double, double, double, double,
                                 double *, double *, double *, double *,
                                 double *, double *, double *, double *, char *, int *);
    double calc_surf_energy_bal(double, double, double, double, double, double,
                                double, double, double, double, double, double,
                                double, double, double, double, double, double,
                                double, double, double, double, double, double,
                                double, double, double,
                                double *, double *, double *, double *, double *,
                                double *, double *, double *, double *, double *,
                                float *, int, int,
                                int, int, int, int, int, int, int, int, int, int,
                                atmos_data_struct *, dmy_struct *,
                                energy_bal_struct *, layer_data_struct *,
                                layer_data_struct *,
                                snow_data_struct *, soil_con_struct *,
                                veg_var_struct *, veg_var_struct *, int);
    double calc_trans(double, double);
    double calc_veg_displacement(double);
    double calc_veg_height(double);
    double calc_veg_roughness(double);
    double calc_water_balance_error(int, double, double, double);
    double canopy_evap(layer_data_struct *, layer_data_struct *,
                       veg_var_struct *, veg_var_struct *, char, int, int,
                       double, double *, double, double, double, double,
                       double, double, double, double, double, double,
                       double *, double *, double *, double *,
#if SPATIAL_FROST
                       double *, float *);
#else
                       float *);
#endif
    void   check_files(filep_struct *, filenames_struct *);
    FILE  *check_state_file(char *, dmy_struct *, global_param_struct *, int, int,
                            int *);
    void   close_files(filep_struct *, out_data_file_struct *, filenames_struct *);
    filenames_struct cmd_proc(int argc, char *argv[]);
    void   collect_eb_terms(energy_bal_struct, snow_data_struct, cell_data_struct,
                            int *, int *, int *, int *, int *, double, double, double,
                            int, int, int, int, double *, double *,
#if SPATIAL_FROST
                            double *, double,
#endif
                            out_data_struct *);
    void   collect_wb_terms(cell_data_struct, veg_var_struct, snow_data_struct, lake_var_struct,
                            double, double, double, double, int, int, int, double *, out_data_struct *);
    void   compress_files(char string[]);
    void   compute_dz(double *, double *, int, double);
    void   correct_precip(double *, double, double, double, double);
    void   compute_pot_evap(int, dmy_struct *, int, int, double, double , double, double, double, double **, double *);
    void   compute_soil_layer_thermal_properties(layer_data_struct *, double *,
            double *, double *,
            double *,
#if SPATIAL_FROST
            double *,
#endif
            int);
    void   compute_treeline(atmos_data_struct *, dmy_struct *, double, double *, char *);
    out_data_struct *create_output_list();

    void   display_current_settings(int, filenames_struct *, global_param_struct *);
    int    dist_prec(atmos_data_struct *,dist_prcp_struct *,soil_con_struct *,
                     veg_con_struct *, lake_con_struct *,
                     dmy_struct *,global_param_struct *,
                     filep_struct *, out_data_file_struct *,
                     out_data_struct *, save_data_struct *,
                     int, int, char, char, char *, int *);
#if QUICK_FS
    int  distribute_node_moisture_properties(double *, double *, double *, double *,
            double *, double *, double *, double ***,
            double *, double *, double *,
            double *, double *, int, int, char);
#else
#if EXCESS_ICE
    int  distribute_node_moisture_properties(double *, double *, double *, double *,
            double *, double *, double *, double *,
            double *, double *, double *,
            double *, double *, double *,
            double *, double *, int, int, char);
#else
    int  distribute_node_moisture_properties(double *, double *, double *,
            double *, double *, double *,
            double *, double *, double *,
            double *, double *, double *,
            double *, double *, int, int, char);
#endif
#endif
    void   distribute_soil_property(double *,double,double,
                                    double **l_param,
                                    int, int, double *, double *);

    double error_calc_atmos_energy_bal(double Tcanopy, ...);
    double error_calc_atmos_moist_bal(double , ...);
    double error_calc_canopy_energy_bal(double Tsurf, ...);
    double error_calc_snow_ground_flux(double Tsurf, ...);
    double error_calc_surf_energy_bal(double Tsurf, ...);
    double ErrorSnowPackEnergyBalance(double Tsurf, ...);
    double error_print_atmos_energy_bal(double, va_list);
    double error_print_atmos_moist_bal(double, va_list);
    double error_print_canopy_energy_bal(double, va_list);
    double error_print_snow_ground_flux(double, va_list);
    double ErrorPrintSnowPackEnergyBalance(double, va_list);
    double error_print_solve_T_profile(double, va_list);
    double error_print_surf_energy_bal(double, va_list);
    double error_solve_T_profile(double Tsurf, ...);
    double estimate_dew_point(double, double, double, double, double);
#if QUICK_FS
    int estimate_layer_ice_content(layer_data_struct *, double *, double *,
                                   double *, double ***, double *,
                                   double *, double ***,
#if SPATIAL_FROST
                                   double *, double,
#endif // SPATIAL_FROST
                                   double *, double *, double *, double *,
                                   int, int, char);
#else
    int estimate_layer_ice_content(layer_data_struct *, double *, double *,
                                   double *, double *, double *, double *,
                                   double *, double *, double *,
#if SPATIAL_FROST
                                   double *, double,
#endif // SPATIAL_FROST
#if EXCESS_ICE
                                   double *, double *,
#endif // EXCESS_ICE
                                   double *, double *, double *, double *,
                                   int, int, char);
#endif
    double estimate_T1(double, double, double, double, double, double, double,
                       double, double, double, double);
    double exp_interp(double,double,double,double,double);

    double f(double, double, double, double, double, double, double, double,
             double, double, int, double *, double, double, double, double *,
             double *, double *, double *, double *, double *);
    void   fda_heat_eqn(double *, double *, int, int, ...);
    void   fdjac3(double *, double *, double *, double *, double *,
                  void (ldas::Vic::*vecfunc)(double *, double *, int, int, ...),
                  int);
    void   find_0_degree_fronts(energy_bal_struct *, double *, double *, int);
    layer_data_struct find_average_layer(layer_data_struct *, layer_data_struct *,
                                         double, double);
    void   find_sublayer_temperatures(layer_data_struct *, double *, double *,
                                      double *, double, double, int, int);
    int    finish_frozen_soil_calcs(energy_bal_struct *, layer_data_struct *,
                                    layer_data_struct *, layer_data_struct *,
                                    soil_con_struct *, int, int, double,
                                    double *, double *, double *, double *);
    void   free_atmos(int nrecs, atmos_data_struct **atmos);
    void   free_dist_prcp(dist_prcp_struct *, int);
    void   free_dmy(dmy_struct **dmy);
    void   free_vegcon(veg_con_struct **);
    void   free_veglib(veg_lib_struct **);
    void   free_out_data_files(out_data_file_struct **);
    void   free_out_data(out_data_struct **);
    int    full_energy(char, int, int, atmos_data_struct *, dist_prcp_struct *,
                       dmy_struct *, global_param_struct *, lake_con_struct *,
                       soil_con_struct *, veg_con_struct *);
    double func_aero_resist(double,double,double,double,double);
    static double func_atmos_energy_bal(Vic*,double, va_list);
    double func_atmos_moist_bal(double, va_list);
    static double func_canopy_energy_bal(Vic*, double, va_list);
    double func_snow_ground_flux(double, va_list);
    static double func_surf_energy_bal(Vic*,double, va_list);
    double get_thresh(double Tair, double SurfaceLiquidWater, double Zo_salt, int flag);
    double get_avg_temp(double, double, double *, double *, int);
    double get_dist(double, double, double, double);
    void   get_force_type(char *, int, int *);
    global_param_struct get_global_param(filenames_struct *, FILE *);
    void   get_next_time_step(int *, int *, int *, int *, int *, int);

    double hermint(double, int, double *, double *, double *, double *, double *);
    void   hermite(int, double *, double *, double *, double *, double *);
    void   HourlyT(int, int, int *, double *, int *, double *, double *);

    void   init_output_list(out_data_struct *, int, char *, int, float);
    void   initialize_atmos(atmos_data_struct *, dmy_struct *, FILE **, double,
                            double, double, double, double, double, double,
                            double, double *,
#if OUTPUT_FORCE
                            char *, out_data_file_struct *, out_data_struct *);
#else
                            char *);
#endif
    void   initialize_global();
    int   initialize_model_state(dist_prcp_struct *, dmy_struct,
                                 global_param_struct *, filep_struct,
                                 int, int, int, int,
                                 double, soil_con_struct *,
                                 veg_con_struct *, lake_con_struct,
                                 char **, int **, save_data_struct *);
    int    initialize_new_storm(cell_data_struct ***, veg_var_struct ***,
                                int, int, int, double, double);
    void   initialize_snow(snow_data_struct **, int, int);
    void   initialize_soil(cell_data_struct **, soil_con_struct *, veg_con_struct *, int);
    void   initialize_veg( veg_var_struct **, veg_con_struct *,
                           global_param_struct *, int);

    void   latent_heat_from_snow(double, double, double, double, double,
                                 double, double, double *, double *,
                                 double *, double *, double *);
    double linear_interp(double,double,double,double,double);

    cell_data_struct **make_cell_data(int, int);
    dist_prcp_struct make_dist_prcp(int);
    dmy_struct *make_dmy(global_param_struct *);
    energy_bal_struct **make_energy_bal(int);
    void make_in_and_outfiles(filep_struct *, filenames_struct *,
                              soil_con_struct *, out_data_file_struct *);
    out_data_struct *make_out_data(int);
    snow_data_struct **make_snow_data(int);
    veg_var_struct **make_veg_var(int);
    void   MassRelease(double *,double *,double *,double *);
#if EXCESS_ICE
    double maximum_unfrozen_water(double, double, double, double, double, double);
#else
    double maximum_unfrozen_water(double, double, double, double);
#endif
#if QUICK_FS
    double maximum_unfrozen_water_quick(double, double, double **);
#endif
    double modify_Ksat(double);
    void mtclim42_wrapper(int, int, double, double, double, double,
                          global_param_struct *, dmy_struct *, double *,
                          double *, double *, double *, double *, double *);

    double new_snow_density(double);
    int    newt_raph(void (ldas::Vic::*vecfunc)(double *, double *, int, int, ...),
                     double *, int);
    void   nrerror(char *);

    void   open_debug();
    FILE  *open_file(char string[], char type[]);
    FILE  *open_state_file(global_param_struct *, filenames_struct, int, int);

    void parse_output_info(filenames_struct *, FILE *, out_data_file_struct **, out_data_struct *);
    double penman(double, double, double, double, double, double, double);
    void   prepare_full_energy(int, int, int, dist_prcp_struct *,
                               soil_con_struct *, double *, double *);
    double priestley(double, double);
    int    put_data(dist_prcp_struct *, atmos_data_struct *,
                    soil_con_struct *, veg_con_struct *,
                    lake_con_struct *, out_data_file_struct *,
                    out_data_struct *, save_data_struct *,
                    dmy_struct *, int);

    double read_arcinfo_value(char *, double, double);
    int    read_arcinfo_info(char *, double **, double **, int **);
    void   read_atmos_data(FILE *, global_param_struct, int, int, double **);
    double **read_forcing_data(FILE **, global_param_struct);
    void   read_initial_model_state(FILE *, dist_prcp_struct *,
                                    global_param_struct *, int, int, int,
                                    soil_con_struct *, int, char *,
                                    int *, lake_con_struct);
    void   read_snowband(FILE *, soil_con_struct *);
    void   read_snowmodel(atmos_data_struct *, FILE *, int, int, int, int);
    soil_con_struct read_soilparam(FILE *, int);
    soil_con_struct read_soilparam_arc(FILE *, char *, int *, int *, int);
    veg_lib_struct *read_veglib(FILE *, int *);
    veg_con_struct *read_vegparam(FILE *, int, int);
    void rescale_lake_fluxes(double oldfrac,
                             double newfrac,
                             lake_var_struct *lake);
    int    redistribute_during_storm(cell_data_struct ***, veg_var_struct ***,
                                     int, int, int, double, double, double,
                                     double *);
    void   redistribute_moisture(layer_data_struct *, double *, double *,
                                 double *, double *, double *, int);
    unsigned char redistribute_moisture_for_storm(double *, double *, double,
            double, double);
    double root_brent(double, double, char *, double (*Function)(Vic*,double, va_list), ...);
    int    runoff(cell_data_struct *, cell_data_struct *,
                  energy_bal_struct *, soil_con_struct *, double *,
#if EXCESS_ICE
                  int,
#endif
#if SPATIAL_FROST
                  double *,
#endif
                  double, int, int, int, int, int);

    void set_max_min_hour(double *, int, int *, int *);
    void set_node_parameters(double *, double *, double *, double *, double *, double *,
                             double *, double *, double *, double *, double *,
                             double *, double *,
#if QUICK_FS
                             double ***,
#endif
#if EXCESS_ICE
                             double *, double *, double *, double *,
#endif
                             int, int, char);
    out_data_file_struct *set_output_defaults(out_data_struct *);
    int set_output_var(out_data_file_struct *, int, int, out_data_struct *, char *, int, char *, int, float);
    double snow_albedo(double, double, double, double, double, double, int, char);
    double snow_density(snow_data_struct *, double, double, double, double, double);
    int    snow_intercept(double, double, double, double, double, double, double,
                          double, double, double, double, double, double, double,
                          double, double,
                          double *, double *, double *, double *, double *,
                          double *, double *, double *, double *, double *,
                          double *, double *, double *, double *, double *,
                          double *, char *, int *, double *, double *, double *,
                          double *, double *, double *, float *, int, int, int,
                          int, int, int, int, layer_data_struct *,
                          layer_data_struct *, soil_con_struct *,
                          veg_var_struct *, veg_var_struct *);
    int    snow_melt(double, double, double, double, double *, double, double *, double,
                     double, double, double, double, double, double, double,
                     double, double, double, double, double, double,
                     double *, double *, double *, double *, double *, double *,
                     double *, double *, double *, double *, double *, double *,
                     int, int, int, int, snow_data_struct *, soil_con_struct *);
    static double SnowPackEnergyBalance(Vic*, double, va_list);
    double soil_conductivity(double, double, double, double, double);
    void   soil_thermal_calc(soil_con_struct *, layer_data_struct *,
                             energy_bal_struct, double *, double *, double *,
                             int, int);
    void advect_soil_veg_storage(double lakefrac,
                                 double max_newfraction,
                                 double newfraction,
                                 double *delta_moist,
                                 soil_con_struct  *soil_con,
                                 veg_con_struct   *veg_con,
                                 cell_data_struct *cell,
                                 veg_var_struct   *veg_var);
    static double soil_thermal_eqn(Vic*,double, va_list);
    double solve_snow(char                 overstory,
                      double               BareAlbedo,
                      double               LongUnderOut, // LW from understory
                      double               MIN_RAIN_TEMP,
                      double               MAX_SNOW_TEMP,
                      double               Tcanopy, // canopy air temperature
                      double               Tgrnd, // soil surface temperature
                      double               air_temp, // air temperature
                      double               density,
                      double               dp,
                      double               ice0,
                      double               longwave,
                      double               moist,
                      double               mu,
                      double               prec,
                      double               pressure,
                      double               shortwave,
                      double               snow_grnd_flux,
                      double               vp,
                      double               vpd,
                      double               wind_h,
                      double              *AlbedoUnder,
                      double              *Evap,
                      double              *Le,
                      double              *LongUnderIn, // surface incomgin LW
                      double              *NetLongSnow, // net LW at snow surface
                      double              *NetShortGrnd, // net SW reaching ground
                      double              *NetShortSnow, // net SW at snow surface
                      double              *ShortUnderIn, // surfave incoming SW
                      double              *Torg_snow,
                      double              *aero_resist,
                      double              *aero_resist_used,
                      double              *coverage, // best guess snow coverage
                      double              *delta_coverage, // cover fract change
                      double              *delta_snow_heat, // change in pack heat
                      double              *displacement,
                      double              *gauge_correction,
                      double              *melt_energy,
                      double              *out_prec,
                      double              *out_rain,
                      double              *out_snow,
                      double              *ppt,
                      double              *rainfall,
                      double              *ref_height,
                      double              *roughness,
                      double              *snow_inflow,
                      double              *snowfall,
                      double              *surf_atten,
                      double              *wind,
                      float               *root,
                      int                  INCLUDE_SNOW,
                      int                  Nnodes,
                      int                  Nveg,
                      int                  band,
                      int                  hour,
                      int                  iveg,
                      int                  day_in_year,
                      int                  dt,
                      int                  month,
                      int                  day,
                      int                  year,
                      int                  rec,
                      int                  veg_class,
                      int                 *UnderStory,
                      energy_bal_struct   *energy,
                      layer_data_struct   *layer_dry,
                      layer_data_struct   *layer_wet,
                      snow_data_struct    *snow,
                      soil_con_struct     *soil_con,
                      veg_var_struct      *veg_var_dry,
                      veg_var_struct      *veg_var_wet);
    double solve_atmos_energy_bal(double Tcanopy, ...);
    double solve_atmos_moist_bal(double , ...);
    double solve_canopy_energy_bal(double Tfoliage, ...);
    double solve_snow_ground_flux(double Tsurf, ...);
    double solve_surf_energy_bal(double Tsurf, ...);
#if QUICK_FS
    int    solve_T_profile(double *, double *, char *, int *, double *, double *,double *,
                           double *, double, double *, double *, double *,
                           double *, double *, double *, double *, double, double *, double ***,
                           int, int *, int, int, int, int);
#else
    int    solve_T_profile(double *, double *, char *, int *, double *, double *,double *,
                           double *, double, double *, double *, double *,
                           double *, double *, double *, double *, double, double *,
#if EXCESS_ICE
                           double *, double *,
#endif
                           int, int *, int, int, int, int);

#endif
    int   solve_T_profile_implicit(double *, double *, double *, double *, double *,
                                   double *, double, double *, double *, double *,
#if EXCESS_ICE
                                   double *, double *,
#endif
                                   double *, double *, double *, double *, double, int, int *,
                                   int, int, int, int,
                                   double *, double *, double *, double *);
    double StabilityCorrection(double, double, double, double, double, double);
    void   store_moisture_for_debug(int,int,double *,cell_data_struct ***,
                                    veg_var_struct ***,snow_data_struct **,
                                    soil_con_struct *);
    int    surface_fluxes(char, double, double, double, double,
#if EXCESS_ICE
                          int, double *, double *,
#endif
                          double, double, double *, double *, double **,
                          double *, double *, double *, double *,
                          double *, double *, double *, double *, double *,
                          float *, int, int, int, int, int,
                          int, int, int, int, atmos_data_struct *, dmy_struct *,
                          energy_bal_struct *, global_param_struct *,
                          cell_data_struct *, cell_data_struct *,
                          snow_data_struct *, soil_con_struct *,
                          veg_var_struct *, veg_var_struct *, float, float, float);
    double svp(double);
    double svp_slope(double);

    void transpiration(layer_data_struct *, int, int, double, double, double,
                       double, double, double, double, double, double, double,
                       double *, double *, double *, double *, double *, double *,
#if SPATIAL_FROST
                       double *,
#endif
                       float *);
    void tridag(double *,double *,double *,double *,double *,int);
    void tridiag(double *, double *, double *, double *, unsigned);
    int update_thermal_nodes(dist_prcp_struct *,
                             int, int, int, soil_con_struct *, veg_con_struct *);
    void usage(char *);

    void   vicerror(char *);
    double volumetric_heat_capacity(double,double,double);

    void write_atmosdata(atmos_data_struct *, int);
    void write_data(out_data_file_struct *, out_data_struct *, dmy_struct *, int);
    void write_debug(atmos_data_struct *, soil_con_struct *, cell_data_struct *,
                     energy_bal_struct *, snow_data_struct *, veg_var_struct *,
                     dmy_struct *, global_param_struct *,
                     double, double, int, int, int, int, int, char);
    void write_dist_prcp(dist_prcp_struct *);
#if OUTPUT_FORCE
    void write_forcing_file(atmos_data_struct *, int, out_data_file_struct *, out_data_struct *);
#endif
    void write_forcing_file(atmos_data_struct *atmos,
                            int                nrecs,
                            out_data_file_struct *out_data_files,
                            out_data_struct   *out_data);
    void write_header(out_data_file_struct *, out_data_struct *, dmy_struct *, global_param_struct);
    void write_layer(layer_data_struct *, int, int,
#if SPATIAL_FROST
                     double *,
#endif
                     double *);
    void write_model_state(dist_prcp_struct *, global_param_struct *, int,
                           int, filep_struct *, soil_con_struct *, char *,
                           int *, lake_con_struct);
    void write_snow_data(snow_data_struct, int, int);
    void write_soilparam(soil_con_struct *);
    void write_vegparam(veg_con_struct *);
    void write_vegvar(veg_var_struct *, int);

    void zero_output_list(out_data_struct *);


    double qromb(double (*sub_with_height)(double z,double es,  double Wind, double AirDens, double ZO,
                                           double EactAir,double F, double hsalt, double phi_r,
                                           double ushear, double Zrh), double es, double Wind, double AirDens, double ZO,
                 double EactAir, double F, double hsalt, double phi_r, double ushear, double Zrh,
                 double a, double b);
    static double (*funcd)(double z,double es,  double Wind, double AirDens, double ZO,
                           double EactAir,double F, double hsalt, double phi_r,
                           double ushear, double Zrh);
    static double sub_with_height(double z,double es,  double Wind, double AirDens, double ZO,
                                  double EactAir,double F, double hsalt, double phi_r,
                                  double ushear, double Zrh);
    static double transport_with_height(double z,double es,  double Wind, double AirDens, double ZO,
                                        double EactAir,double F, double hsalt, double phi_r,
                                        double ushear, double Zrh);
    double rtnewt(double x1, double x2, double xacc, double Ur, double Zr);
    void get_shear(double x, double *f, double *df, double Ur, double Zr);
    double get_prob(double Tair, double Age, double SurfaceLiquidWater, double U10);
    void shear_stress(double U10, double ZO,double *ushear, double *Zo_salt, double utshear);
    double CalcSubFlux(double EactAir, double es, double Zrh, double AirDens, double utshear,
                       double ushear, double fe, double Tsnow, double Tair, double U10,
                       double Zo_salt, double F, double *Transport);
    void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
    static double trapzd(double (*funcd)(double z,double es,  double Wind, double AirDens, double ZO,
                                         double EactAir,double F, double hsalt, double phi_r,
                                         double ushear, double Zrh), double es, double Wind, double AirDens, double ZO,
                         double EactAir, double F, double hsalt, double phi_r, double ushear,
                         double Zrh, double a, double b, int n);
    static double IceEnergyBalance(Vic*,double TSurf, va_list ap);
    double get_mean(double *, int, double);
    double get_stdev(double *, int, double, double);
    double get_sum(double *, int, double);
    double get_min(double *, int, double);
    double get_max(double *, int, double);
    double ErrorPrintIcePackEnergyBalance(double TSurf, va_list ap);
    double ErrorIcePackEnergyBalance(double Tsurf, ...);
    double CalcIcePackEnergyBalance(double Tsurf, ...);
    int ice_melt(double            z2,
                 double            aero_resist,
                 double            *aero_resist_used,
                 double            Le,
                 snow_data_struct *snow,
                 lake_var_struct  *lake,
                 int               delta_t,
                 double            displacement,
                 double            Z0,
                 double            surf_atten,
                 double            rainfall,
                 double            snowfall,
                 double            wind,
                 double            Tcutoff,
                 double            air_temp,
                 double            net_short,
                 double            longwave,
                 double            density,
                 double            pressure,
                 double            vpd,
                 double            vp,
                 double           *melt,
                 double           *save_advection,
                 double           *save_deltaCC,
                 double           *save_SnowFlux,
                 double           *save_latent,
                 double           *save_sensible,
                 double           *save_Qnet,
                 double           *save_refreeze_energy,
                 double           *save_LWnet,
                 double            fracprv);
    int ice_depth(lake_con_struct lake_con, double volume, double ice_water_eq, double *hice);
    int get_depth_from_sarea(lake_con_struct lake_con, double sarea, double *depth);
    int get_depth(lake_con_struct lake_con, double volume, double *depth);
    int get_volume(lake_con_struct lake_con, double depth, double *volume);
    int get_sarea(lake_con_struct lake_con, double depth, double *sarea);
    int initialize_lake (lake_var_struct   *lake,
                         lake_con_struct   lake_con,
                         soil_con_struct  *soil_con,
                         double            airtemp);

    /*** Subroutine prototypes ***/

    double adjflux(double, double, double ,double, double, double, double,
                   double, double, double, double *, double *);
    void advect_snow_storage(double, double, double, snow_data_struct *);
    void alblake(double, double, double *, double *, float *, float *, double, double,
                 int, int *, double, double, char *, int);
    double calc_density(double);
    void colavg (double *, double *, double *, float, double *, int, double, double);
    float dragcoeff(float, double, double);
    void eddy (int, double, double * , double *, double *, double, int, double, double);
    void energycalc(double *, double *, int, double, double,double *, double *, double *);
    void iceform (double *,double *,double ,double,double *,int, int, double, double, double *, double *, double *, double *, double *, double);
    void icerad(double,double ,double,double *, double *,double *);
    int lakeice(double *, double, double, double, double, int,
                double, double, double *, double, double, int, dmy_struct, double *, double *, double, double);
    void latsens(double,double, double, double, double, double, double, double,
                 double *, double *, double);
    float lkdrag(float, double, double, double, double);
    lake_con_struct read_lakeparam(FILE *, soil_con_struct, veg_con_struct *);
    void rescale_soil_veg_fluxes(double, double, cell_data_struct *, veg_var_struct *);
    void rescale_snow_energy_fluxes(double, double, snow_data_struct *, energy_bal_struct *);
    void rhoinit(double *, double);
    int solve_lake(double, double, double, double, double, double, double, double,
                   double, double, lake_var_struct *, lake_con_struct,
                   soil_con_struct, int, int, double, dmy_struct, double);
    double specheat (double);
    void temp_area(double, double, double, double *, double *, double *, double *, int, double *, int, double, double, double*, double *, double *);
    void tracer_mixer(double *, int *, int, double*, int, double, double, double *);
    void tridia(int, double *, double *, double *, double *, double *);
    int water_balance (lake_var_struct *, lake_con_struct, int, dist_prcp_struct *, int, int, int, double, soil_con_struct, veg_con_struct,
#if EXCESS_ICE
                       int, double,
#endif
                       double, double);
    int  water_energy_balance(int, double*, double*, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double *, double *, double *, double*, double *, double *, double *, double, double *, double *, double *, double *, double *, double);
    int water_under_ice(int, double,  double, double *, double *, double, int, double, double, double, double *, double *, double *, double *, int, double, double, double, double *);
    /********************************
     **                             **
     **    FUNCTION PROTOTYPES      **
     **                             **
     ********************************/
    int calc_tair(const ldas::vic::control_struct *ctrl, const ldas::vic::parameter_struct *p,
                  ldas::vic::data_struct *data);
    int calc_prcp(const ldas::vic::control_struct *ctrl, const ldas::vic::parameter_struct *p,
                  ldas::vic::data_struct *data);
    /* start vic_change */
    int calc_srad_humidity(const ldas::vic::control_struct *ctrl, const ldas::vic::parameter_struct *p,
                           ldas::vic::data_struct *data, double *tiny_radfract);
    /* end vic_change */
    /* start vic_change */
    int calc_srad_humidity_iterative(const ldas::vic::control_struct *ctrl,
                                     const ldas::vic::parameter_struct *p, ldas::vic::data_struct *data,
                                     double *hourly_radfract);
    /* end vic_change */
    int data_alloc(const ldas::vic::control_struct *ctrl, ldas::vic::data_struct *data);
    int data_free(const ldas::vic::control_struct *ctrl, ldas::vic::data_struct *data);
    double calc_pet(double rad, double ta, double pa, double dayl);
    double atm_pres(double elev);
    int pulled_boxcar(double *input,double *output,int n,int w,int w_flag);
    void mtclim42_init(int have_dewpt, double elevation, double annual_prcp,
                       double lat, global_param_struct *vic_global, dmy_struct *dmy,
                       double *prec, double *tmax, double *tmin, double *hourlyrad,
                       double *tiny_radfract, ldas::vic::control_struct *ctrl,
                       ldas::vic::parameter_struct *p, ldas::vic::data_struct *mtclim42_data);

    void mtclim42_to_vic(int have_dewpt, int have_shortwave, double hour_offset,
                         global_param_struct *vic_global, dmy_struct *dmy,
                         double *tiny_radfract, ldas::vic::control_struct *ctrl,
                         ldas::vic::data_struct *mtclim42_data, double *tskc, double *vp,
                         double *hourlyrad);
    void ttrim( char *c );
};
}
#endif // __LDAS_VIC_H
