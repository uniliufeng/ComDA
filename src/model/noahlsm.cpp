#include "noahlsm.h"

using namespace ldas;
using namespace ldas::noah;
using namespace std;

Noahlsm::Noahlsm()
{
    //ctor
    init();
}

Noahlsm::~Noahlsm()
{
    //dtor
    if (nsoil>0)
    {
        delete et;
        delete smav;
    }
    /*
    if (__module_sf_noahlsm_MOD_lucats>0)
    {
        delete __module_sf_noahlsm_MOD_shdtbl;
        delete __module_sf_noahlsm_MOD_nrotbl;
        delete __module_sf_noahlsm_MOD_rstbl;
        delete __module_sf_noahlsm_MOD_rgltbl;
        delete __module_sf_noahlsm_MOD_hstbl;
        delete __module_sf_noahlsm_MOD_snuptbl;
        delete __module_sf_noahlsm_MOD_maxalb;
        delete __module_sf_noahlsm_MOD_laimintbl;
        delete __module_sf_noahlsm_MOD_laimaxtbl;
        delete __module_sf_noahlsm_MOD_emissmintbl;
        delete __module_sf_noahlsm_MOD_emissmaxtbl;
        delete __module_sf_noahlsm_MOD_albedomintbl;
        delete __module_sf_noahlsm_MOD_albedomaxtbl;
        delete __module_sf_noahlsm_MOD_z0mintbl;
        delete __module_sf_noahlsm_MOD_z0maxtbl;
    }*/
    /*
    if (__module_sf_noahlsm_MOD_slcats>0)
    {
        delete __module_sf_noahlsm_MOD_bb;
        delete __module_sf_noahlsm_MOD_drysmc;
        delete __module_sf_noahlsm_MOD_f11;
        delete __module_sf_noahlsm_MOD_maxsmc;
        delete __module_sf_noahlsm_MOD_refsmc;
        delete __module_sf_noahlsm_MOD_satpsi;
        delete __module_sf_noahlsm_MOD_satdk;
        delete __module_sf_noahlsm_MOD_satdw;
        delete __module_sf_noahlsm_MOD_wltsmc;
        delete __module_sf_noahlsm_MOD_qtz;
    }*/
    /*
    if (__module_sf_noahlsm_MOD_slpcats>0)
        delete __module_sf_noahlsm_MOD_slope_data;*/
    if (nsoil>0)
    {
        delete sldpth;
        delete smc;
        delete stc;
        delete sh2o;
    }
    delete psih2;
    delete psim2;
}

Noahlsm::Noahlsm(const Noahlsm& other)
{
    //copy ctor
}

Noahlsm& Noahlsm::operator=(const Noahlsm& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}
void Noahlsm::config(const std::string& fn)
{
    strcpy(forcing_filename,fn.c_str());
    //std::cout<<forcing_filename<<std::endl;
    open_forcing_file();
    /*__module_ascii_io_MOD_open_forcing_file(&iunit, forcing_filename, infotext, &nsoil, startdate, enddate, &loop_for_a_while,
                                            &latitude, &longitude,
                                            &forcing_timestep, &noahlsm_timestep, &ice, &t1, stc, smc, sh2o, sldpth, &cmc, &snowh, &sneqv, &tbot,
                                            &vegtyp, &soiltyp, &slopetyp, &snoalb, &zlvl, &zlvl_wind, albedo_monthly, shdfac_monthly,
                                            z0brd_monthly, lai_monthly, &use_urban_module, &isurban, &shdmin, &shdmax, &usemonalb, &rdlai2d, llanduse,
                                            &sfcdif_option, &iz0tlnd,fn.length(),0);
                                            */
    //fortran版本中有错误,sfcdif_option不变为1
    sfcdif_option=0;
    dt = noahlsm_timestep;
    //std::cout<<"here?"<<std::endl;
    m_time=time_convert(startdate);
    endtime=time_convert(enddate);
    std::cout<<m_time;
    std::cout<<endtime;

    et=new float[nsoil];
    smav=new float[nsoil];
    for(int i=0; i<nsoil; i++)
    {
        et[i]=badval;
        smav[i] = badval;
    }


    ///
    /// Set up some input variables for SFLX.
    ///

    ///
    /// LLANDUSE:
    ///Currently only the USGS vegetation dataset as used in WRF is supported.
    ///

    strcpy(llanduse,"USGS");

    ///
    /// LSOIL:
    /// Currently, only the STATSGO soil dataset as used in WRF is supported.
    ///

    strcpy(lsoil,"STAS");

    ///
    /// Read our lookup tables and parameter tables:
    ///VEGPARM.TBL, SOILPARM.TBL, GENPARM.TBL
    ///

    soil_veg_gen_parm();

    ///
    /// COSZ is unused if we're not using the urban canopy model.  If we implement the
    /// urban canopy model for this simple point driver, we will need to compute a COSZ
    /// somewhere.
    ///

    cosz = badval;

    ///
    /// PRCPRAIN is unused.
    ///

    prcprain = badval;

    ///
    /// SOLARDIRECT is unused.
    ///

    solardirect = badval;

    ///
    /// Set EMISSI for our first time step.  Just a guess, but it's only for the
    /// first step.  Later time steps get EMISSI from what was set in the prior
    /// time step by SFLX.
    ///

    emissi = 0.96;

    ///
    /// For the initial ALBEDO value used in computing SOLNET, just use our
    /// snow-free value.  Subsequent timesteps use the value computed in the
    /// previous call to SFLX:
    ///

    albedo = month_d(albedo_monthly, m_time);

    ///
    /// For the initial value of Z0 (used in SFCDIF_off to compute CH and CM),
    /// just use a snow-free background value.  Subsequent timesteps use this
    /// value as computed in the previous call to SFLX:
    ///

    z0 = month_d(z0brd_monthly, m_time);

    ///
    /// Z0BRD is computed within SFLX.  But we need an initial value, so the call to
    /// SFCDIF_MYJ can do its thing.  Subsequent timesteps will recycle the Z0BRD
    /// value as returned from SFLX in the previous timestep.
    ///

    if ( sfcdif_option == 1 )
        z0brd  = z0;
    else
        z0brd = badval;


    ///
    /// CZIL is needed for the SFCDIF_OFF step.  This comes from CZIL_DATA, as read
    /// from the GENPARM.TBL file, which is how REDPRM ultimately gets it as well:
    ///

    czil = __module_sf_noahlsm_MOD_czil_data;

    ///
    ///  CM and CH, computed in subroutine SFCDIF_OFF, need initial values.  Values are
    ///  subsequently updated for each time step.  So, just take a guess at reasonable
    ///  initial values:
    ///

    ch = 1.e-4;
    cm = 1.e-4;

    if ( sfcdif_option == 1 )
        myjsfcinit();

}
void Noahlsm::init()
{
    iunit = 10;
    snotime1  = 0.0;
    ribb      = 0.0;
    float badval=-1.E36;
    sheat = badval;
    etakin = badval;
    fdown = badval;
    ec = badval;
    edir = badval;
    ett = badval;
    esnow = badval;
    drip = badval;
    dew = badval;
    beta=badval;
    t1 = badval;
    snowh = badval;
    sneqv = badval;
    etp = badval;
    ssoil = badval;
    flx1 = badval;
    flx2 = badval;
    flx3 = badval;
    snomlt = badval;
    sncovr = badval;
    runoff1 = badval;
    runoff2 = badval;
    runoff3 = badval;
    rc = badval;
    pc = badval;
    rcs = badval;
    rct = badval;
    rcq = badval;
    rcsoil = badval;
    soilw = badval;
    soilm = badval;
    q1 = badval;
    smcwlt = badval;
    smcdry = badval;
    smcref = badval;
    smcmax = badval;
    rsmin = badval;
    nroot = -999999;
    __module_sf_noahlsm_MOD_lucats=0;
    __module_sf_noahlsm_MOD_slcats=0;
    __module_sf_noahlsm_MOD_slpcats=0;
    nsoil=0;
    psih2=new float[kztm];
    psim2=new float[kztm];
	local=false;
	embrd=0;
}
void Noahlsm::run()
{
    ktime = 0;
    forcing.open(forcing_filename);
    while(m_time<endtime)
    {
        ++ktime;
        if ((loop_for_a_while) && (m_time == endtime))
        {
            std::cout<<m_time;
            //output_close();

            m_time = time_convert(startdate);
            ktime=1;
            /*__module_ascii_io_MOD_read_forcing_text(&iunit, startdate, &forcing_timestep,
                                                    &sfcspd, &sfcu, &sfcv, &sfctmp, &q2, &sfcprs, &soldn, &longwave, &prcp, &ierr);
            if (ierr == 0) throw Exception("Noahlsm","run()","Wrong input for looping a year." );*/
            //rewind(iunit);
            read_forcing_text();
        }
        strcpy(nowdate,convert_time(m_time).c_str());
        /*__module_ascii_io_MOD_read_forcing_text(&iunit, nowdate , &forcing_timestep,
                                                &sfcspd, &sfcu, &sfcv, &sfctmp, &q2, &sfcprs, &soldn, &longwave, &prcp, &ierr);*/
        read_forcing_text();
        //if (ierr != 0)
          //  throw Exception("Noahlsm","run()",":  FORCING DATA READ PROBLEM" );
        if ( (prcp > 0) && (sfctmp < 273.15) )
            ffrozp = 1.0;
        else
            ffrozp = 0.0;
        /**
        ! At each time step, using the forcing fields (and T1, the skin temperature, which
        ! gets updated by SFLX), we need to compute a few additional thermodynamic variables.
        ! Ultimately, TH2, Q2SAT and DQSDT2 get passed to SFLX;
        !             T1V and TH2V get used in SFCDIF_off but are not used by SFLX.
        !             RHO is not used in SFLX, but is used in URBAN.  It is a dummy variable
        !             as far as the non-urban code is concerned.
        */
        caltmp(t1, sfctmp, sfcprs, zlvl, q2, th2, t1v, th2v, rho);
        calhum(sfctmp, sfcprs, q2sat, dqsdt2);
        if (usemonalb)
            alb    = month_d(albedo_monthly, m_time);
        else
            alb = badval;
        if (rdlai2d)
            xlai = month_d(lai_monthly, m_time);
        else
            xlai = badval;
        shdfac = month_d(shdfac_monthly, m_time);
        if (q1 == badval)
            q1 = q2;
        /**
        ! SFCDIF_OFF computes mixing lengths for momentum and heat, CM and CH.
        ! Z0 is needed for SFCDIF_OFF.  We use the Z0 as computed in the previous
        ! timestep of SFLX, but on the first time step, we need a value of Z0.  This
        ! is set above from our background value.  The initial value may not be quite
        ! what we want, but for that one timestep, it should be OK.  Additionally,
        ! CH and CM need some values for the initial timestep.  These values are
        ! set above.
        */
        iz0tlnd=1;//why?
        if ( sfcdif_option == 0 )

            //sfcdif_off ( zlvl_wind , z0 , t1v , th2v , sfcspd , czil , cm , ch ,
              //           vegtyp , isurban , iz0tlnd );// out:  cm, ch
__module_sf_noahlsm_MOD_sfcdif_off ( &zlvl_wind , &z0 , &t1v , &th2v , &sfcspd , &czil , &cm , &ch ,
                         &vegtyp , &isurban , &iz0tlnd );// out:  cm, ch

        else if ( sfcdif_option == 1 )

            sfcdif_myj ( zlvl_wind , z0 , z0brd , sfcprs , t1 , sfctmp , q1 ,
                         q2 , sfcspd , czil , ribb , cm , ch , vegtyp , isurban , iz0tlnd );

        solnet = soldn * (1.0-albedo);
        lwdn = longwave * emissi;

        /**
             !  Call the Noah LSM routine for a single time step.
             !

             !
             ! Input:
             !
             !    FFROZP      -- Fraction of total precipitation which is frozen ( Fraction [0.0-1.0] )
             !    ICE         -- Land point (ICE==0) or sea-ice point (ICE==1) ( Integer flag 0 or 1 )
             !    ISURBAN     -- The vegetation category for Urban points
             !    DT          -- Time step ( seconds )
             !    ZLVL        -- Height of atmospheric forcing variables ( m AGL )
             !    NSOIL       -- Number of soil layers ( count )
             !    SLDPTH      -- Thickness of each soil layer ( m )
             !    LOCAL       -- Logical flag, .TRUE. to use table values for ALB, SHDFAC, and Z0BRD
             !                   .FALSE. to use values for ALB, SHDFAC, and Z0BRD as set in this driver routine
             !    LLANDUSE    -- Land-use dataset we're using.  "USGS" is the only dataset supported
             !    LSOIL       -- Soil dataset we're using.  "STAS" (for STATSGO) is the only dataset supported
             !    LWDN        -- Longwave downward radiation flux ( W m{-2} )
             !    SOLDN       -- Shortwave downward radiation flux ( W m{-2} )
             !    SOLNET      -- Shortwave net radiation flux ( W m{-2} )
             !    SFCPRS      -- Atmospheric pressure at height ZLVL m AGL ( Pa )
             !    PRCP        -- Precipitation rate ( kg m{-2} s{-1} )
             !    SFCTMP      -- Air temperature at height ZLVL m AGL ( K )
             !    Q2          -- Atmospheric mixing ratio at height ZLVL m AGL ( kg kg{-1} )
             !    SFCSPD      -- Wind speed at height ZLVL m AGL ( m s{-1} )
             !    COSZ        -- Cosine of the Solar Zenith Angle (unused in SFLX)
             !    PRCPRAIN    -- Liquid precipitation rate ( kg m{-2} s{-1} ) (unused)
             !    SOLARDIRECT -- Direct component of downward solar radiation ( W m{-2} ) (unused)
             !    TH2         -- Air potential temperature at height ZLVL m AGL ( K )
             !    Q2SAT       -- Saturation specific humidity at height ZLVL m AGL ( kg kg{-1} )
             !    DQSDT2      -- Slope of the Saturation specific humidity curve at temperature SFCTMP ( kg kg{-1} K{-1} )
             !    VEGTYP      -- Vegetation category ( index )
             !    SOILTYP     -- Soil category ( index )
             !    SLOPETYP    -- Slope category ( index )
             !    SHDFAC      -- Areal fractional coverage of green vegetation ( fraction [0.0-1.0] ).
             !                   SHDFAC will be set by REDPRM if (LOCAL == .TRUE.)
             !    SHDMIN      -- Minimum areal fractional coverage of green vegetation ( fraction [0.0-1.0] )
             !    SHDMAX      -- Maximum areal fractional coverage of green vegetation ( fraction [0.0-1.0] )
             !    ALB         -- Surface background snow-free albedo (fraction [0.0-1.0]).  ALB will
             !                   be set by REDPRM if (LOCAL == .TRUE.).
             !    SNOALB      -- Maximum deep-snow albedo. ( fraction [0.0-1.0] )
             !    TBOT        -- Constant deep-soil temperature ( K )
             !    Z0BRD       -- Background (i.e., without snow-cover effects) roughness length ( M )
             !    Z0          -- Roughness length, including snow-cover effects ( M )
             !    EMBRD       -- Background emissivity (i.e., not including snow-cover effects) ( fraction [0.0-1.0] )
             !
             ! Updated:
             !
             !    EMISSI      -- Emissivity ( fraction )
             !    CMC         -- Canopy moisture content ( kg m{-2} )
             !    T1          -- Skin temperature ( K )
             !    STC         -- Soil temperature at NSOIL levels ( K )
             !    SMC         -- Volumetric soil moisture content at NSOIL levels ( m{3} m{-3} )
             !    SH2O        -- Liquid portion of the volumetric soil moisture content at NSOIL levels ( m{3} m{-3} )
             !    SNOWH       -- Snow depth ( m )
             !    SNEQV       -- Water equivalent snow depth ( m )
             !    ALBEDO      -- Surface albedo, including any snow-cover effects ( Fraction [0.0-1.0] )
             !    CH          -- Surface exchange coefficient for heat and moisture ( m s{-1} )
             !    CM          -- Surface exchange coefficient for momentum, unused in this code ( m s{-1} )
             !    ETA         -- Latent heat flux (evapotranspiration) ( W m{-2} )
             !    SHEAT       -- Sensible heat flux ( W m{-2} )
             !    ETAKIN      -- Latent heat flux (evapotranspiration) ( kg m{-2} s{-1} )
             !    FDOWN       -- Radiation forcing at the surface ( W m{-2} )
             !    EC          -- Latent heat flux component: canopy water evaporation ( W m{-2} )
             !    EDIR        -- Latent heat flux component: direct soil evaporation ( W m{-2} )
             !    ET          -- Latent heat flux component: plant transpiration from each of NSOIL levels ( W m{-2} )
             !    ETT         -- Latent heat flux component: total plant transpiration ( W m{-2} )
             !    ESNOW       -- Latent heat flux component: sublimation from (or deposition to) snowpack ( W m{-2} )
             !    DRIP        -- Precipitation or dew falling through canopy, in excess of canopy holding capacity ( m )
             !    DEW         -- Dewfall (or frostfall for T<273.15) ( m )
             !    BETA        -- Ratio of actual to potential evapotranspiration ( Fraction [0.0-1.0] )
             !    ETP         -- Potential evapotranspiration ( W m{-2} )
             !    SSOIL       -- Soil heat flux ( W m{-2} )
             !    FLX1        -- Latent heat flux from precipitation accumulating as snow ( W m{-2} )
             !    FLX2        -- Latent heat flux from freezing rain converting to ice ( W m{-2} )
             !    FLX3        -- Latent heat flux from melting snow ( W m{-2} )
             !    SNOMLT      -- Snow melt water ( m )
             !    SNCOVR      -- Fractional snow cover ( Fraction [0.0-1.0] )
             !    RUNOFF1     -- Surface runoff, not infiltrating the soil ( m s{-1} )
             !    RUNOFF2     -- Subsurface runoff, drainage out the bottom of the last soil layer ( m s{-1} )
             !    RUNOFF3     -- Internal soil layer runoff ( m s{-1} )
             !    RC          -- Canopy resistance ( s m{-1} )
             !    PC          -- Plant coefficient, where PC * ETP = ETA ( Fraction [0.0-1.0] )
             !    RSMIN       -- Minimum canopy resistance ( s m{-1} )
             !    XLAI        -- Leaf area index ( dimensionless )
             !    RCS         -- Incoming solar RC factor ( dimensionless )
             !    RCT         -- Air temperature RC factor ( dimensionless )
             !    RCQ         -- Atmospheric water vapor deficit RC factor ( dimensionless )
             !    RCSOIL      -- Soil moisture RC factor ( dimensionless )
             !    SOILW       -- Available soil moisture in the root zone ( Fraction [SMCWLT-SMCMAX] )
             !    SOILM       -- Total soil column moisture content, frozen and unfrozen ( m )
             !    Q1          -- Effective mixing ratio at the surface ( kg kg{-1} )
             !    SMAV        -- Soil Moisture Availability at each level, fraction between SMCWLT (SMAV=0.0) and SMCMAX (SMAV=1.0)
             !    SMCWLT      -- Wilting point ( m{3} m{-3} )
             !    SMCDRY      -- Dry soil moisture threshold where direct evaporation from the top layer ends ( m{3} m{-3} )
             !    SMCREF      -- Soil moisture threshold where transpiration begins to stress ( m{3} m{-3} )
             !    SMCMAX      -- Porosity, i.e., saturated value of soil moisture ( m{3} m{-3} )
             !    NROOT       -- Number of root layers ( count )
             */
		/*std::cout<<"ffrozp:"<<ffrozp<<std::endl;
		std::cout<<"ice:"<<ice<<std::endl;
		std::cout<<"isurban:"<<isurban<<std::endl;
		std::cout<<"dt:"<<dt<<std::endl;
		std::cout<<"zlvl:"<<zlvl<<std::endl;
		std::cout<<"nsoil:"<<nsoil<<std::endl;
		//std::cout<<"sldpth:"<<sldpth<<std::endl;
		std::cout<<"local:"<<local<<std::endl;
		std::cout<<"llanduse:"<<llanduse<<std::endl;
		std::cout<<"lsoil:"<<lsoil<<std::endl;
		std::cout<<"lwdn:"<<lwdn<<std::endl;
		std::cout<<"soldn:"<<soldn<<std::endl;
		std::cout<<"solnet:"<<solnet<<std::endl;
		std::cout<<"sfcprs:"<<sfcprs<<std::endl;
		std::cout<<"prcp:"<<prcp<<std::endl;
		std::cout<<"sfctmp:"<<sfctmp<<std::endl;
		std::cout<<"q2:"<<q2<<std::endl;
		std::cout<<"sfcspd:"<<sfcspd<<std::endl;
		std::cout<<"cosz:"<<cosz<<std::endl;
		std::cout<<"prcprain:"<<prcprain<<std::endl;
		std::cout<<"solardirect:"<<solardirect<<std::endl;
		std::cout<<"th2:"<<th2<<std::endl;
		std::cout<<"q2sat:"<<q2sat<<std::endl;
		std::cout<<"dqsdt2:"<<dqsdt2<<std::endl;
		std::cout<<"vegtyp:"<<vegtyp<<std::endl;
		std::cout<<"soiltyp:"<<soiltyp<<std::endl;
		std::cout<<"slopetyp:"<<slopetyp<<std::endl;
		std::cout<<"shdfac:"<<shdfac<<std::endl;
		std::cout<<"shdmin:"<<shdmin<<std::endl;
		std::cout<<"shdmax:"<<shdmax<<std::endl;
		std::cout<<"alb:"<<alb<<std::endl;
		std::cout<<"snoalb:"<<snoalb<<std::endl;
		std::cout<<"tbot:"<<tbot<<std::endl;
		std::cout<<"z0brd:"<<z0brd<<std::endl;
		std::cout<<"z0:"<<z0<<std::endl;
		std::cout<<"emissi:"<<emissi<<std::endl;
		std::cout<<"embrd:"<<embrd<<std::endl;*/
		int local1,rdlai2d1,usemonalb1;
		local1=local?1:0;
		rdlai2d1=rdlai2d?1:0;
		usemonalb1=usemonalb?1:0;

        __module_sf_noahlsm_MOD_sflx(&ffrozp, &ice, &isurban, &dt, &zlvl, &nsoil, sldpth,
                                     &local1,
                                     llanduse, lsoil,
                                     &lwdn, &soldn, &solnet, &sfcprs, &prcp, &sfctmp, &q2, &sfcspd,
                                     &cosz, &prcprain, &solardirect,
                                     &th2, &q2sat, &dqsdt2,
                                     &vegtyp, &soiltyp, &slopetyp, &shdfac, &shdmin, &shdmax,
                                     &alb, &snoalb, &tbot, &z0brd, &z0, &emissi, &embrd,
                                     &cmc, &t1, stc, smc, sh2o, &snowh, &sneqv, &albedo, &ch, &cm,
                                     &eta, &sheat, &etakin, &fdown,
                                     &ec, &edir, et, &ett, &esnow, &drip, &dew,
                                     &beta, &etp, &ssoil,
                                     &flx1, &flx2, &flx3,
                                     &snomlt, &sncovr,
                                     &runoff1, &runoff2, &runoff3,
                                     &rc, &pc, &rsmin, &xlai, &rcs ,&rct, &rcq, &rcsoil,
                                     &soilw, &soilm, &q1, smav, &rdlai2d1, &usemonalb1, &snotime1,
                                     &ribb,
                                     &smcwlt, &smcdry, &smcref, &smcmax, &nroot);
        qfx = edir + ec + ett + esnow;
        f = solnet + lwdn;
        fup = emissi * STBOLT * std::pow(t1,4);
        res = f - sheat + ssoil - eta - fup - flx1 - flx2 - flx3;
        std::cout<<"res:"<<res<<std::endl;
        m_time.addSeconds(noahlsm_timestep);
    }
    forcing.close();
}


void Noahlsm::caltmp(float t1, // skin temperature [k]
                     float sfctmp, // air temperature [k] at level ZLVL
                     float sfcprs, // atmospheric pressure [pa] at level ZLVL
                     float zlvl, // height [m] whre atmospheric fields are available
                     float q2, // specific huminity [kg/kg] at level ZLVL
                     float &th2, // potential temperature (considering the reference pressure to be at the surface
                     float &t1v, // virtual skin temperature
                     float &th2v, // virtual potential temperature at ZLVL
                     float &rho) // densith
{
    float t2v;
    float RD = 287.04;
    th2 = sfctmp + (0.0098 * zlvl); // potential temperature
    t1v = t1 * (1.0 + 0.61 * q2);
    th2v = th2 * (1.0 + 0.61 * q2);
    t2v = sfctmp * (1.0 + 0.61 * q2);
    rho = sfcprs / (RD * t2v);
}

void Noahlsm::calhum(float sfctmp,
                     float sfcprs,
                     float &q2sat,  // saturated specific humidity
                     float &dqsdt2) // slope of saturated specific huminity curve
{
    float A2 = 17.67;
    float A3 = 273.15;
    float A4 = 29.65;
    float ELWV = 2.501E6;
    float E0 = 611.0;
    float RV = 461.0;
    float EPSILON = 0.622;
    float A23M4 = A2 * (A3 - A4);

    float es;
    es = E0 * std::exp(ELWV / RV * (1. / A3 - 1. / sfctmp));
    q2sat = EPSILON * es / (sfcprs - (1 - EPSILON) * es);
    dqsdt2 = q2sat * A23M4 / std::pow((sfctmp - A4), 2);
}

float Noahlsm::month_d(const float a12[12], const Time& t) const
{
    int nowm = t.getMonth();
    int nowd = t.getDay();
    int ndays[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    int prevm, postm;
    float factor;
    if(t.is_leap_year()==true)
    {
        ndays[2] = 29;
    }
    if(nowd==15)
    {
        return(a12[nowm-1]);
    }
    else if(nowd<15)
    {
        postm = nowm;
        prevm = nowm -1;
        if(prevm==0)
        {
            prevm=12;
        }
        factor = (1.0*ndays[prevm]- 15.0 + 1.0*nowd)/(1.0 * ndays[prevm]);
    }
    else if(nowd>15)
    {
        prevm = nowm;
        postm = nowm + 1;
        if(postm==13)
        {
            postm = 1;
        }
        factor = (1.0*nowd - 15.0)/(1.0*ndays[prevm]);
    }
    return(a12[prevm-1]*(1.0-factor) + a12[postm-1]*factor);
}

void Noahlsm::myjsfcinit()
{
    int k;
    float x,zeta1,zeta2,zrng1,zrng2,ztmax1,dzeta1;
    const float pihf=3.1415926/2.;
    const float eps=1.e-6;
    const float ztmin1 = -5.0;
    const float ztmin2=-5.0;
    /**----------------------------------------------------------------------
    !
    !***  compute surface layer integral functions
    !
    !----------------------------------------------------------------------
    */
    ztmax1=1.0;
    ztmax2=1.0;

    zrng1=ztmax1-ztmin1;
    zrng2=ztmax2-ztmin2;

    dzeta1=zrng1/(kztm-1);
    dzeta2=zrng2/(kztm-1);

    /**----------------------------------------------------------------------
    !***  function definition loop
    !----------------------------------------------------------------------
    */
    zeta1=ztmin1;
    zeta2=ztmin2;

    for(k=1; k<=kztm; k++)
    {
        /**----------------------------------------------------------------------
        !***  unstable range
        !----------------------------------------------------------------------
        */
        if(zeta2<0.)
        {
            /**----------------------------------------------------------------------
            !***  paulson 1970 functions
            !----------------------------------------------------------------------
            */
            x=std::sqrt(std::sqrt(1.-16.*zeta2));
            psim2[k-1]=-2.*log((x+1.)/2.)-log((x*x+1.)/2.)+2.*atan(x)-pihf;
            psih2[k-1]=-2.*log((x*x+1.)/2.);
        }
        /**--------------------------------------------------------------------
        !***  stable range
        !----------------------------------------------------------------------
        */
        else
        {
            /**--------------------------------------------------------------------
            !***  holtslag and de bruin 1988
            !----------------------------------------------------------------------
            */
            psim2[k-1]=0.7*zeta2+0.75*zeta2*(6.-0.35*zeta2)*std::exp(-0.35*zeta2);
            psih2[k-1]=0.7*zeta2+0.75*zeta2*(6.-0.35*zeta2)*std::exp(-0.35*zeta2);
        }
        if(k==kztm)
        {
            ztmax1=zeta1;
            ztmax2=zeta2;
        }
        zeta1=zeta1+dzeta1;
        zeta2=zeta2+dzeta2;
    }
    ztmax1=ztmax1-eps;
    ztmax2=ztmax2-eps;
}


void Noahlsm::sfcdif_myj(float zsl, float z0, float z0base, float sfcprs, float tz0, float tlow, float qz0,
                         float qlow, float sfcspd, float czil, float& rib, float& akms, float& akhs, int vegtyp, int isurban, int iz0tlnd )
{
    const int itrmx = 5;

    const float excml   = 0.0001;
    const float excms   = 0.0001;
    const float vkarman = 0.4;
    const float ztfc    = 1.0;
    const float elocp   = 2.72e6 / cp;
    const float epsu2   = 1.e-6;
    const float epsust  = 1.e-9;
    const float sqvisc  = 258.2;
    const float ric     = 0.505;
    const float epszt   = 1.e-28;

    const int kztm2 = kztm-2;

    const float wwst  = 1.2;
    const float wwst2 = wwst*wwst;
    const float ztmin2 = -5.0;

    /**   ****************************************************************
    !     *                                                              *
    !     *                       SURFACE LAYER                          *
    !     *                                                              *
    !     ****************************************************************
    */
    float thlow; /// potential temperature (k) at level zsl above ground.
    float thz0; ///potential temperature (k) at level z0 above ground.
    float thelow;
    float cwmlow; ///cloud water (assumed to be zero for hrldas)
    float ustar;
    float rlmo;
    int itr,k;
    float czil_local;
    float zilfc;

    float a,b,btgh,btgx,cxchl,cxchs,czetmax,dthv,du2,elfc,
          pshz,pshzl,psmz,psmzl,
          rdzt,
          rlogt,rlogu,rz,rzst,rzsu,simh,simm,tem,thm,
          ustark,wstar2,
          zetalt,zetalu,
          zetat,zetau,zq,zslt,zslu,zt,zu,zzil;
    /**----------------------------------------------------------------------
    ! Compute potential temperatures from input <Za>-level temperature TLOW
    ! and input Skin-level (or Z0-level) temperature.

        THLOW = TLOW * (P0/SFCPRS) ** RCP
        THZ0  = TZ0  * (P0/SFCPRS) ** RCP  ! Should actually use a pressure at Z0?

        ! For THELOW, WRF uses (1-QC*ELOCP/T)*TH.  If QC is assumed to be zero,
        ! Then THELOW is just TH ?
        THELOW = THLOW
    !----------------------------------------------------------------------
    */
    thlow = tlow * std::pow(p0/sfcprs,rcp);
    thz0  = tz0  * std::pow(p0/sfcprs,rcp);
	thelow = thlow;
    cxchl=excml/zsl;

    btgx=g/thlow;
    elfc=vkarman*btgx;

//mbkwm    if(pblh>1000.)then
//mbkwm       btgh=btgx*pblh
//mbkwm    else
    btgh=btgx*1000.;
//mbkwm    endif

    thm=(thelow+thz0)*0.5;
    tem=(tlow+tz0)*0.5;

    a=thm*p608;
    b=(elocp/tem-1.-p608)*thm;
    cwmlow = 0.0; //! mb/kwm assume no cloud water

    dthv=((thelow-thz0)*((qlow+qz0+cwmlow)*(0.5*p608)+1.)
          +(qlow-qz0+cwmlow)*a+cwmlow*b);

//!    du2=max((ulow-uz0)**2+(vlow-vz0)**2,epsu2)
    du2 = max(sfcspd*sfcspd,epsu2);
    rib=btgx*dthv*zsl/du2;

    zu=z0;
    zt=zu*ztfc;
    zslu=zsl+zu;
    rzsu=zslu/zu;
    rlogu=std::log(rzsu);
    zslt=zsl+zu; //! u,v and t are at the same level

    if ( (iz0tlnd==0) || (vegtyp == isurban) )
    {
        //! just use the incoming czil value.
        zilfc=-czil*vkarman*sqvisc;
    }
    else
    {
        //! modify czil according to chen & zhang, 2009
        //! czil = 10 ** -0.40 h, ( where h = 10*zo )
        czil_local = std::pow(10.0,  ( -0.40 * ( z0 / 0.07 ) ));
        zilfc=-czil_local*vkarman*sqvisc;
    }

    /**----------------------------------------------------------------------
    !
    !  rib modification to zilfc term
    !
    !      czetmax = 20.
    */
    czetmax = 10.;
//! stable
    if(dthv>0.)
        if (rib<ric)
            zzil=zilfc*(1.0+(rib/ric)*(rib/ric)*czetmax);
        else
            zzil=zilfc*(1.0+czetmax);
//! unstable
    else
        zzil=zilfc;

    //! first guess of ustar based on input (previous timestep) akhs and akms?
    if (btgh * akhs * dthv != 0.0)
        wstar2 = wwst2* pow(abs (btgh * akhs * dthv), 2./3.);
    else
        wstar2 = 0.0;

    ustar = max (std::sqrt (akms * std::sqrt (du2+ wstar2)),epsust);

    //land_point_iteration: do itr=1,itrmx
    for(itr=1; itr<=itrmx; itr++)
    {
        /**----------------------------------------------------------------------
        !***  zilitinkevitch fix for zt
        !----------------------------------------------------------------------
        */
        zt=max(exp(zzil*sqrt(ustar*z0base))*z0base,epszt);
        rzst=zslt/zt;
        rlogt=log(rzst);

        /**--------------------------------------------------------------------
        !***  1./monin-obukhov length-scale
        !----------------------------------------------------------------------
        */
        rlmo=elfc*akhs*dthv/pow(ustar,3);
        zetalu=zslu*rlmo;
        zetalt=zslt*rlmo;
        zetau=zu*rlmo;
        zetat=zt*rlmo;

        zetalu=min(max(zetalu,ztmin2),ztmax2);
        zetalt=min(max(zetalt,ztmin2),ztmax2);
        zetau=min(max(zetau,ztmin2/rzsu),ztmax2/rzsu);
        zetat=min(max(zetat,ztmin2/rzst),ztmax2/rzst);

        /**--------------------------------------------------------------------
        !***  land functions
        !----------------------------------------------------------------------
        */
        rz=(zetau-ztmin2)/dzeta2;
        k=int(rz);
        rdzt=rz-k;
        k=min(k,kztm2);
        k=max(k,0);
        psmz=(psim2[k+2-1]-psim2[k+1-1])*rdzt+psim2[k+1-1];

        rz=(zetalu-ztmin2)/dzeta2;
        k=int(rz);
        rdzt=rz-k;
        k=min(k,kztm2);
        k=max(k,0);
        psmzl=(psim2[k+2-1]-psim2[k+1-1])*rdzt+psim2[k+1-1];

        simm=psmzl-psmz+rlogu;

        rz=(zetat-ztmin2)/dzeta2;
        k=int(rz);
        rdzt=rz-k;
        k=min(k,kztm2);
        k=max(k,0);
        pshz=(psih2[k+2-1]-psih2[k+1-1])*rdzt+psih2[k+1-1];

        rz=(zetalt-ztmin2)/dzeta2;
        k=int(rz);
        rdzt=rz-k;
        k=min(k,kztm2);
        k=max(k,0);
        pshzl=(psih2[k+2-1]-psih2[k+1-1])*rdzt+psih2[k+1-1];

        simh=(pshzl-pshz+rlogt);

        ustark=ustar*vkarman;
        akms=max(ustark/simm,cxchl);
        akhs=max(ustark/simh,cxchl);

        /**--------------------------------------------------------------------
        !***  beljaars correction for ustar
        !----------------------------------------------------------------------
        */
        if(dthv<=0.)
            wstar2=wwst2*pow(abs(btgh*akhs*dthv),(2./3.));
        else
            wstar2=0.;

        ustar=max(sqrt(akms*sqrt(du2+wstar2)),epsust);
    }
}


void  Noahlsm::soil_veg_gen_parm()
{
    //CHARACTER(LEN=*), INTENT(IN) :: MMINLU, MMINSL
    int lumatch, iindex, lc, num_slope;
    int ierr;
    int open_ok = 0;

    //char mess[128] , message;
    //logical, external :: wrf_dm_on_monitor
    /**---SPECIFY VEGETATION RELATED CHARACTERISTICS :
    !             ALBBCK: SFC albedo (in percentage)
    !                 Z0: Roughness length (m)
    !             SHDFAC: Green vegetation fraction (in percentage)
    !  Note: The ALBEDO, Z0, and SHDFAC values read from the following table
    !          ALBEDO, amd Z0 are specified in LAND-USE TABLE; and SHDFAC is
    !          the monthly green vegetation data
    !             CMXTBL: MAX CNPY Capacity (m)
    !             NROTBL: Rooting depth (layer)
    !              RSMIN: Mimimum stomatal resistance (s m-1)
    !              RSMAX: Max. stomatal resistance (s m-1)
    !                RGL: Parameters used in radiation stress function
    !                 HS: Parameter used in vapor pressure deficit functio
    !               TOPT: Optimum transpiration air temperature. (K)
    !             CMCMAX: Maximum canopy water capacity
    !             CFACTR: Parameter used in the canopy inteception calculati
    !               SNUP: Threshold snow depth (in water equivalent m) that
    !                     implies 100% snow cover
    !                LAI: Leaf area index (dimensionless)
    !             MAXALB: Upper bound on maximum albedo over deep snow
    */
//!-----READ IN VEGETAION PROPERTIES FROM VEGPARM.TBL

    std::ifstream vegetbl("VEGPARM.TBL");
    if(vegetbl.good()!=true)
        throw Exception("Noahlsm",
                        "soil_veg_gen_parm",
                        "can not open vegetation parameter table file.");

    std::istringstream iss, iss0;
    int i;
    char buf[256];

    // line 1
    vegetbl.getline(buf, 256);

    // line 2
    vegetbl.getline(buf, 256);
    modi_buf(buf);
    iss.str(buf);
    iss>>lutype;
    iss.clear();

    // line 3
    vegetbl.getline(buf, 256);
    modi_buf(buf);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_lucats;
    iss.clear();

/*    __module_sf_noahlsm_MOD_shdtbl=new float[__module_sf_noahlsm_MOD_lucats];
    __module_sf_noahlsm_MOD_nrotbl=new int[__module_sf_noahlsm_MOD_lucats];
    __module_sf_noahlsm_MOD_rstbl=new float[__module_sf_noahlsm_MOD_lucats];
    __module_sf_noahlsm_MOD_rgltbl=new float[__module_sf_noahlsm_MOD_lucats];
    __module_sf_noahlsm_MOD_hstbl=new float[__module_sf_noahlsm_MOD_lucats];
    __module_sf_noahlsm_MOD_snuptbl=new float[__module_sf_noahlsm_MOD_lucats];
    __module_sf_noahlsm_MOD_maxalb=new float[__module_sf_noahlsm_MOD_lucats];
    __module_sf_noahlsm_MOD_laimintbl=new float[__module_sf_noahlsm_MOD_lucats];
    __module_sf_noahlsm_MOD_laimaxtbl=new float[__module_sf_noahlsm_MOD_lucats];
    __module_sf_noahlsm_MOD_emissmintbl=new float[__module_sf_noahlsm_MOD_lucats];
    __module_sf_noahlsm_MOD_emissmaxtbl=new float[__module_sf_noahlsm_MOD_lucats];
    __module_sf_noahlsm_MOD_albedomintbl=new float[__module_sf_noahlsm_MOD_lucats];
    __module_sf_noahlsm_MOD_albedomaxtbl=new float[__module_sf_noahlsm_MOD_lucats];
    __module_sf_noahlsm_MOD_z0mintbl=new float[__module_sf_noahlsm_MOD_lucats];
    __module_sf_noahlsm_MOD_z0maxtbl=new float[__module_sf_noahlsm_MOD_lucats];*/
    // line 4-30
    for(int n=0; n<__module_sf_noahlsm_MOD_lucats; n++)
    {
        vegetbl.getline(buf, 256);
        modi_buf(buf);
        iss.str(buf);
        //std::cout<<iss.str()<<std::endl;
        iss>>i; // index
        iss>>__module_sf_noahlsm_MOD_shdtbl[n];
        iss>>__module_sf_noahlsm_MOD_nrotbl[n];
        iss>>__module_sf_noahlsm_MOD_rstbl[n];
        iss>>__module_sf_noahlsm_MOD_rgltbl[n];
        iss>>__module_sf_noahlsm_MOD_hstbl[n];
        iss>>__module_sf_noahlsm_MOD_snuptbl[n];
        iss>>__module_sf_noahlsm_MOD_maxalb[n];
        iss>>__module_sf_noahlsm_MOD_laimintbl[n];
        iss>>__module_sf_noahlsm_MOD_laimaxtbl[n];
        iss>>__module_sf_noahlsm_MOD_emissmintbl[n];
        iss>>__module_sf_noahlsm_MOD_emissmaxtbl[n];
        iss>>__module_sf_noahlsm_MOD_albedomintbl[n];
        iss>>__module_sf_noahlsm_MOD_albedomaxtbl[n];
        iss>>__module_sf_noahlsm_MOD_z0mintbl[n];
        iss>>__module_sf_noahlsm_MOD_z0maxtbl[n];
        iss.clear();
    }

    // line 31
    vegetbl.getline(buf, 256);
    // line 32
    vegetbl.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_topt_data;
    iss.clear();

    // line 33
    vegetbl.getline(buf, 256);
    // line 34
    vegetbl.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_cmcmax_data;
    iss.clear();

    // line 35
    vegetbl.getline(buf, 256);
    // line 36
    vegetbl.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_cfactr_data;
    iss.clear();

    // line 37
    vegetbl.getline(buf, 256);
    // line 38
    vegetbl.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_rsmax_data;
    iss.clear();

    // line 39
    vegetbl.getline(buf, 256);
    // line 40
    vegetbl.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_bare;
	iss.clear();

    // line 41
    vegetbl.getline(buf, 256);
    // line 42
    vegetbl.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_natural;
    vegetbl.close();

//!-----READ IN SOIL PROPERTIES FROM SOILPARM.TBL
    std::ifstream soiltbl("SOILPARM.TBL");
    if(soiltbl.good()!=true)
        throw Exception("Noahlsm",
                        "soil_veg_gen_parm",
                        "can not open soil parameter table file.");
    iss.clear();

    // line 1
    soiltbl.getline(buf, 256);

    // line 2
    soiltbl.getline(buf, 256);
    iss.str(buf);
    iss>>sltype;
    iss.clear();

    // line 3
    soiltbl.getline(buf, 256);
    modi_buf(buf);
    iss.str(buf);
    //std::cout<<iss.str()<<std::endl;
    iss>>__module_sf_noahlsm_MOD_slcats;
    iss.clear();

    /*__module_sf_noahlsm_MOD_bb=new float[__module_sf_noahlsm_MOD_slcats];
    __module_sf_noahlsm_MOD_drysmc=new float[__module_sf_noahlsm_MOD_slcats];
    __module_sf_noahlsm_MOD_f11=new float[__module_sf_noahlsm_MOD_slcats];
    __module_sf_noahlsm_MOD_maxsmc=new float[__module_sf_noahlsm_MOD_slcats];
    __module_sf_noahlsm_MOD_refsmc=new float[__module_sf_noahlsm_MOD_slcats];
    __module_sf_noahlsm_MOD_satpsi=new float[__module_sf_noahlsm_MOD_slcats];
    __module_sf_noahlsm_MOD_satdk=new float[__module_sf_noahlsm_MOD_slcats];
    __module_sf_noahlsm_MOD_satdw=new float[__module_sf_noahlsm_MOD_slcats];
    __module_sf_noahlsm_MOD_wltsmc=new float[__module_sf_noahlsm_MOD_slcats];
    __module_sf_noahlsm_MOD_qtz=new float[__module_sf_noahlsm_MOD_slcats];*/
    // line 4 - 22
    for(int n=0; n<__module_sf_noahlsm_MOD_slcats; n++)
    {
        soiltbl.getline(buf, 256);
        modi_buf(buf);
        iss.str(buf);
        //std::cout<<iss.str()<<std::endl;
        iss>>i;
        iss>>__module_sf_noahlsm_MOD_bb[n];
        iss>>__module_sf_noahlsm_MOD_drysmc[n];
        iss>>__module_sf_noahlsm_MOD_f11[n];
        iss>>__module_sf_noahlsm_MOD_maxsmc[n];
        iss>>__module_sf_noahlsm_MOD_refsmc[n];
        iss>>__module_sf_noahlsm_MOD_satpsi[n];
        iss>>__module_sf_noahlsm_MOD_satdk[n];
        iss>>__module_sf_noahlsm_MOD_satdw[n];
        iss>>__module_sf_noahlsm_MOD_wltsmc[n];
        iss>>__module_sf_noahlsm_MOD_qtz[n];
        iss.clear();
    }
    soiltbl.close();

//!-----READ IN GENERAL PARAMETERS FROM GENPARM.TBL
    std::ifstream gen_param("GENPARM.TBL");
    if(gen_param.good()!=true)
        throw Exception("Noahlsm",
                        "soil_veg_gen_parm",
                        "can not open soil parameter table file.");
    iss.clear();
    // line 1
    gen_param.getline(buf, 256);
    // line 2
    gen_param.getline(buf, 256);
    // line 3
    gen_param.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_slpcats;
    iss.clear();

    /*__module_sf_noahlsm_MOD_slope_data=new float[__module_sf_noahlsm_MOD_slpcats];*/
    // line 4 - 12
    for(int n=0; n<__module_sf_noahlsm_MOD_slpcats; n++)
    {
        gen_param.getline(buf, 256);
        iss.str(buf);

        iss>>__module_sf_noahlsm_MOD_slope_data[n];
        iss.clear();
    }
    // line 13
    gen_param.getline(buf, 256);
    // line 14
    gen_param.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_sbeta_data;
    iss.clear();

    // line 15
    gen_param.getline(buf, 256);
    // line 16
    gen_param.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_fxexp_data;
    iss.clear();

    // line 17
    gen_param.getline(buf, 256);
    // line 18
    gen_param.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_csoil_data;
    iss.clear();

    // line 19
    gen_param.getline(buf, 256);
    // line 20
    gen_param.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_salp_data;
    iss.clear();

    // line 21
    gen_param.getline(buf, 256);
    // line 22
    gen_param.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_refdk_data;
    iss.clear();

    // line 23
    gen_param.getline(buf, 256);
    // line 24
    gen_param.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_refkdt_data;
    iss.clear();

    // line 25
    gen_param.getline(buf, 256);
    // line 26
    gen_param.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_frzk_data;
    iss.clear();

    // line 27
    gen_param.getline(buf, 256);
    // line 28
    gen_param.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_zbot_data;
    iss.clear();

    // line 29
    gen_param.getline(buf, 256);
    // line 30
    gen_param.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_czil_data;
    iss.clear();

    // line 31
    gen_param.getline(buf, 256);
    // line 32
    gen_param.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_smlow_data;
    iss.clear();

    // line 33
    gen_param.getline(buf, 256);
    // line 34
    gen_param.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_smhigh_data;
    iss.clear();

    // line 35
    gen_param.getline(buf, 256);
    // line 36
    gen_param.getline(buf, 256);
    iss.str(buf);
    iss>>__module_sf_noahlsm_MOD_lvcoef_data;
    iss.clear();

    gen_param.close();
}

void Noahlsm::modi_buf(char *buf)
{
    int nlen = std::strlen(buf);
    for(int n=0; n<nlen; n++)
    {
        if(buf[n]==',')
        {
            buf[n]=' ';
        }
        if(buf[n]=='\'' || buf[n]=='\n')
        {
            buf[n]='\0';
            break;
        }
    }
}
Time Noahlsm::time_convert(const char startdate[12]) const
{
    std::string tmp;
    tmp=startdate[0];
    tmp+=startdate[1];
    tmp+=startdate[2];
    tmp+=startdate[3];
    tmp+='-';
    tmp+=startdate[4];
    tmp+=startdate[5];
    tmp+="-";
    tmp+=startdate[6];
    tmp+=startdate[7];
    tmp+=" ";
    tmp+=startdate[8];
    tmp+=startdate[9];
    tmp+=":";
    tmp+=startdate[10];
    tmp+=startdate[11];
    return Time(tmp);
}
std::string Noahlsm::convert_time(const Time& t) const
{
    std::ostringstream os;

    os<<t.getYear();
    os.width(2);
    os.fill('0');
    os<<t.getMonth();
    os.width(2);
    os.fill('0');
    os<<t.getDay();
    os.width(2);
    os.fill('0');
    os<<t.getHour();
    os.width(2);
    os.fill('0');
    os<<t.getMinute();
    return os.str();
}
void Noahlsm::sfcdif_off(float zlm,float z0,float thz0,float thlm,float sfcspd,float czil,float& akms,float& akhs,
                         int vegtyp, int isurban, int iz0tlnd)
{
    /**---------------------------------------------------------------------
    ! calculate surface layer exchange coefficients via iterative process.
    ! see chen et al (1997, blm)
    ! ----------------------------------------------------------------------
    */
    float wwst, wwst2, vkrm, excm, beta, btg, elfc, wold, wnew;
    int  itrmx, ilech, itr;
    float pihf, epsu2, epsust, epsit, epsa, ztmin, ztmax, hpbl, sqvisc;
    float ric, rric, fhneu, rfc, rfac, zz;
    float xx, yy;
    float zilfc, zu, zt, rdz, cxch;
    float dthv, du2, btgh, wstar2, ustar, zslu, zslt, rlogu, rlogt;
    float rlmo, zetalt, zetalu, zetau, zetat, xlu4, xlt4, xu4, xt4;
    float xlu, xlt, xu, xt, psmz, simm, pshz, simh, ustark, rlmn,rlma;
    wwst = 1.2;
    wwst2 = wwst * wwst;
    vkrm = 0.40;
    excm = 0.001;
    beta = 1./270.;
    btg = beta * g;
    elfc = vkrm * btg;
    wold =.15;
    wnew = 1. - wold;
    itrmx = 5;
    pihf = 3.14159265/2.;
    epsu2 = 1.e-4;
    epsust = 0.07;
    epsit = 1.e-4;
    epsa = 1.e-8;
    ztmin = -5.;
    ztmax = 1.;
    hpbl = 1000.0;
    sqvisc = 258.2;
    ric = 0.183;
    rric = 1.0/ ric;
    fhneu = 0.8;
    rfc = 0.191;
    rfac = ric / (fhneu * rfc * rfc);


    /** ----------------------------------------------------------------------
    !     ztfc: ratio of zoh/zom  less or equal than 1
    !     c......ztfc=0.1
    !     czil: constant c in zilitinkevich, s. s.1995,:note about zt
    ! ----------------------------------------------------------------------
    */
    ilech = 0;

    if ( (iz0tlnd==0) || (vegtyp == isurban) )
        //! just use the original czil value.
        zilfc = - czil * vkrm * sqvisc;
    else
    {
        //! modify czil according to chen & zhang, 2009
        //! czil = 10 ** -0.40 h, ( where h = 10*zo )
        czil = pow(10.0 ,( -0.40 * ( z0 / 0.07 ) ));
        zilfc = - czil * vkrm * sqvisc;
    }
    zu = z0;
    rdz = 1./ zlm;
    cxch = excm * rdz;
    dthv = thlm - thz0;

//! beljars correction of ustar
    du2 = max (sfcspd * sfcspd,epsu2);
//!   if statements to avoid tangent linear problems near zero
    btgh = btg * hpbl;
    if (btgh * akhs * dthv != 0.0)
        wstar2 = wwst2* pow(abs (btgh * akhs * dthv), (2./3.));
    else
        wstar2 = 0.0;

//! zilitinkevitch approach for zt
    ustar = max (sqrt (akms * sqrt (du2+ wstar2)),epsust);

    zt = exp (zilfc * sqrt (ustar * z0))* z0;
    zslu = zlm + zu;


    zslt = zlm + zt;
    rlogu = log (zslu / zu);

    rlogt = log (zslt / zt);

    rlmo = elfc * akhs * dthv / pow(ustar,3);
//! 1./monin-obukkhov length-scale
    //do itr = 1,itrmx
    for(int itr=1; itr<=itrmx; itr++)
    {
        zetalt = max (zslt * rlmo,ztmin);
        rlmo = zetalt / zslt;
        zetalu = zslu * rlmo;
        zetau = zu * rlmo;

        zetat = zt * rlmo;
        if (ilech == 0)
            if (rlmo < 0.)
            {
                xlu4 = 1. -16.* zetalu;
                xlt4 = 1. -16.* zetalt;
                xu4 = 1. -16.* zetau;

                xt4 = 1. -16.* zetat;
                xlu = sqrt (sqrt (xlu4));
                xlt = sqrt (sqrt (xlt4));
                xu = sqrt (sqrt (xu4));

                xt = sqrt (sqrt (xt4));

                psmz = pspmu (xu);
                simm = pspmu (xlu) - psmz + rlogu;
                pshz = psphu (xt);
                simh = psphu (xlt) - pshz + rlogt;
            }
            else
            {
                zetalu = min (zetalu,ztmax);
                zetalt = min (zetalt,ztmax);

                psmz = pspms (zetau);
                simm = pspms (zetalu) - psmz + rlogu;
                pshz = psphs (zetat);
                simh = psphs (zetalt) - pshz + rlogt;
            }
//! lech's functions
        else if (rlmo < 0.)
        {
            psmz = pslmu (zetau);
            simm = pslmu (zetalu) - psmz + rlogu;
            pshz = pslhu (zetat);
            simh = pslhu (zetalt) - pshz + rlogt;
        }
        else
        {
            zetalu = min (zetalu,ztmax);

            zetalt = min (zetalt,ztmax);

            psmz = pslms (zetau);
            simm = pslms (zetalu) - psmz + rlogu;
            pshz = pslhs (zetat);
            simh = pslhs (zetalt) - pshz + rlogt;
        }

//! beljaars correction for ustar
//! zilitinkevitch fix for zt

        ustar = max (sqrt (akms * sqrt (du2+ wstar2)),epsust);

        zt = exp (zilfc * sqrt (ustar * z0))* z0;
        zslt = zlm + zt;

        rlogt = log (zslt / zt);
        ustark = ustar * vkrm;
        akms = max (ustark / simm,cxch);

//! if statements to avoid tangent linear problems near zero

        akhs = max (ustark / simh,cxch);
        if (btgh * akhs * dthv != 0.0)
            wstar2 = wwst2* pow(abs (btgh * akhs * dthv), (2./3.));
        else
            wstar2 = 0.0;

        rlmn = elfc * akhs * dthv / pow(ustar,3);
        rlma = rlmo * wold+ rlmn * wnew;
        rlmo = rlma;
    }
}
/**---------------------------------------------------------------------
! note: the two code blocks below define functions
! ----------------------------------------------------------------------
! lech's surface functions
! ----------------------------------------------------------------------
*/
float Noahlsm::pslmu(const float zz) const
{
    return -0.96* log (1.0-4.5* zz);
}
float Noahlsm::pslms(const float zz) const
{
    float ric = 0.183;
    float rric = 1.0/ ric;
    return zz * rric -2.076* (1. -1./ (zz +1.));
}
float Noahlsm::pslhu(const float zz) const
{
    return -0.96* log (1.0-4.5* zz);
}
/**---------------------------------------------------------------------
! paulson's surface functions
! ----------------------------------------------------------------------
*/
float Noahlsm::pslhs(const float zz) const
{
    float ric = 0.183;
    float fhneu = 0.8;
    float rfc = 0.191;
    float rfac = ric / (fhneu * rfc * rfc);
    return zz * rfac -2.076* (1. -1./ (zz +1.));
}
float Noahlsm::pspmu(const float xx) const
{
    float pihf = 3.14159265/2.;
    return -2.* log ( (xx +1.)*0.5) - log ( (xx * xx +1.)*0.5)+2.* atan (xx)- pihf;
}
float Noahlsm::pspms(const float yy) const
{
    return 5.* yy;
}
float Noahlsm::psphu(const float xx) const
{
    return -2.* log ( (xx * xx +1.)*0.5);
}

/**---------------------------------------------------------------------
! this routine sfcdif can handle both over open water (sea, ocean) and
! over solid surface (land, sea-ice).
! ----------------------------------------------------------------------
*/
float Noahlsm::psphs(const float yy) const
{
    return 5.* yy;
}
void Noahlsm::open_forcing_file()
{
    std::ifstream forcing(forcing_filename);
    bool sea_ice_point=false;
    int glacial_veg_category=0;

    if(forcing.good()!=true)
        throw Exception("Noahlsm",
                        "open_forcing_file",
                        "can not open config&forcing file.");
    std::istringstream iss;
    char buf[256];
    std::string tmp;
    std::string::iterator it;
    float sl_thickness[100],sl_temperature[100],sl_moisture[100],sl_liquid[100];
    for(int i=0; i<100; i++)
    {
        sl_thickness[i]=badval;
        sl_temperature[i]=badval;
        sl_moisture[i]=badval;
        sl_liquid[i]=badval;
    }
    // line 1
    forcing.getline(buf, 256);
    iss.str(buf);
    iss>>tmp;
    if (tmp=="&METADATA_NAMELIST")
    {
        //处理配置信息
        for(;;)
        {
            forcing.getline(buf,256);
            modi_buf(buf);
            iss.clear();
            iss.str(buf);
            iss>>tmp;
            if (tmp=="/")
                break;
            else if (tmp=="startdate")
            {
                iss>>tmp;//=
                iss>>tmp;//"199801010630"
                it=tmp.begin();
                tmp.erase(it);
                it=tmp.end()-1;
                tmp.erase(it);
                strcpy(startdate,tmp.c_str());
            }
            else if (tmp=="enddate")
            {
                iss>>tmp;//=
                iss>>tmp;//"199801010630"
                it=tmp.begin();
                tmp.erase(it);
                it=tmp.end()-1;
                tmp.erase(it);
                strcpy(enddate,tmp.c_str());
            }
            else if (tmp=="loop_for_a_while")
            {
                iss>>tmp;//=
                iss>>tmp;//"FALSE"
                if (tmp=="FALSE" || tmp=="false")
                    loop_for_a_while=false;
                else
                    loop_for_a_while=true;
            }
            else if (tmp=="Latitude")
            {
                iss>>tmp>>latitude;
            }
            else if (tmp=="Longitude")
                iss>>tmp>>longitude;
            else if (tmp=="Forcing_Timestep")
                iss>>tmp>>forcing_timestep;
            else if (tmp=="Noahlsm_Timestep")
                iss>>tmp>>noahlsm_timestep;
            else if (tmp=="Sea_ice_point")
            {
                iss>>tmp;//=
                iss>>tmp;//".FALSE."
                if (tmp=="FALSE" || tmp=="false"||tmp==".FALSE."||tmp==".false.")
                    sea_ice_point=false;
                else
                    sea_ice_point=true;
            }
            else if (tmp=="Soil_layer_thickness")
            {
                iss>>tmp;//=sldpth
                for(int i=0; i<100; i++)
                {
                    if (!iss.eof())
                        iss>>sl_thickness[i];
                    else
                    {
                        nsoil=i;
                        break;
                    }
                }
            }
            else if (tmp=="Soil_Temperature")
            {
                //stc
                iss>>tmp;//=
                for(int i=0; i<100; i++)
                {
                    if (!iss.eof())
                        iss>>sl_temperature[i];
                    else
                    {
                        nsoil=i;
                        break;
                    }
                }
            }
            else if (tmp=="Soil_Moisture")
            {
                //smc
                iss>>tmp;//=
                for(int i=0; i<100; i++)
                {
                    if (!iss.eof())
                        iss>>sl_moisture[i];
                    else
                    {
                        nsoil=i;
                        break;
                    }
                }
            }
            else if (tmp=="Soil_Liquid")
            {
                //sh2o
                iss>>tmp;//=
                for(int i=0; i<100; i++)
                {
                    if (!iss.eof())
                        iss>>sl_liquid[i];
                    else
                    {
                        nsoil=i;
                        break;
                    }
                }
            }
            else if (tmp=="Skin_Temperature")
            {
                iss>>tmp>>t1;
            }
            else if (tmp=="Canopy_water")
            {
                iss>>tmp>>cmc;
            }
            else if (tmp=="Snow_depth")
            {
                iss>>tmp>>snowh;
            }
            else if (tmp=="Snow_equivalent")
            {
                iss>>tmp>>sneqv;
            }
            else if (tmp=="Deep_Soil_Temperature")
            {
                iss>>tmp>>tbot;
            }
            else if (tmp=="Landuse_dataset")
            {
                iss>>tmp;//=
                iss>>tmp;//"199801010630"
                it=tmp.begin();
                tmp.erase(it);
                it=tmp.end()-1;
                tmp.erase(it);
                //lutype=it;
                strcpy(llanduse,tmp.c_str());
            }
            else if (tmp=="Soil_type_index")
            {
                iss>>tmp>>soiltyp;
            }
            else if (tmp=="Vegetation_type_index")
            {
                iss>>tmp>>vegtyp;
            }
            else if (tmp=="Urban_veg_category")
            {
                iss>>tmp>>isurban;
            }
            else if (tmp=="glacial_veg_category")
            {
                iss>>tmp>>glacial_veg_category;
            }
            else if (tmp=="Slope_type_index")
            {
                iss>>tmp>>slopetyp;
            }
            else if (tmp=="Max_snow_albedo")
            {
                iss>>tmp>>snoalb;
            }
            else if (tmp=="Air_temperature_level")
            {
                iss>>tmp>>zlvl;
            }
            else if (tmp=="Wind_level")
            {
                iss>>tmp>>zlvl_wind;
            }
            else if (tmp=="Green_Vegetation_Min")
            {
                iss>>tmp>>shdmin;
            }
            else if (tmp=="Green_Vegetation_Max")
            {
                iss>>tmp>>shdmax;
            }
            else if (tmp=="Usemonalb")
            {
                iss>>tmp;//=
                iss>>tmp;//".FALSE."
                if (tmp=="FALSE" || tmp=="false"||tmp==".FALSE."||tmp==".false.")
                    usemonalb=false;
                else
                    usemonalb=true;
            }
            else if (tmp=="Rdlai2d")
            {
                iss>>tmp;//=
                iss>>tmp;//".FALSE."
                if (tmp=="FALSE" || tmp=="false"||tmp==".FALSE."||tmp==".false.")
                    rdlai2d=false;
                else
                    rdlai2d=true;
            }
            else if (tmp=="sfcdif_option")
            {
                iss>>tmp>>sfcdif_option;
            }
            else if (tmp=="iz0tlnd")
            {
                iss>>tmp>>iz0tlnd;
            }
            else if (tmp=="Albedo_monthly")
            {
                iss>>tmp;
                for(int i=0; i<12; i++)
                    iss>>albedo_monthly[i];
            }
            else if (tmp=="Shdfac_monthly")
            {
                iss>>tmp;
                for(int i=0; i<12; i++)
                    iss>>shdfac_monthly[i];
            }
            else if (tmp=="lai_monthly")
            {
                iss>>tmp;
                for(int i=0; i<12; i++)
                    iss>>lai_monthly[i];
            }
            else if (tmp=="Z0brd_monthly")
            {
                iss>>tmp;
                for(int i=0; i<12; i++)
                    iss>>z0brd_monthly[i];
            }
        }
        if (sea_ice_point)
        {
            ice=1;
        }
        else
        {
            if ( vegtyp == glacial_veg_category )
                ice = -1;
            else
                ice = 0;
        }
        if (nsoil>0)
        {
            sldpth=new float[nsoil];
            stc=new float[nsoil];
            smc=new float[nsoil];
            sh2o=new float[nsoil];
            for(int i=0; i<nsoil; i++)
            {
                sldpth[i]=sl_thickness[i];
                stc[i]=sl_temperature[i];
                smc[i]=sl_moisture[i];
                sh2o[i]=sl_liquid[i];
            }
        }

    }
    forcing.close();
}

void Noahlsm::read_forcing_text()
{
    const float svp1=611.2;
    const float svp2=17.67;
    const float svp3=29.65;
    const float svpt0=273.15;
    const float eps=0.622;


    char buf[256];
    std::string tmp;
    std::istringstream iss;
    int year,month,day,hour,minute;
    float windspeed,winddir,temperature,humidity,pressure,swrad,lwrad,rain;
    //first time?
    Time t=time_convert(startdate);
    //t.addSeconds(noahlsm_timestep);
    if (m_time==t)
    {

        for(;;)
        {
            forcing.getline(buf,256);
            iss.clear();
            iss.str(buf);
            iss>>tmp;
            if (tmp=="<Forcing>")
                break;
        }
    }
    for(;;)
    {
        //
        forcing.getline(buf,256);
        iss.clear();
        iss.str(buf);
        iss>>year>>month>>day>>hour>>minute>>windspeed>>winddir>>temperature>>humidity>>pressure>>swrad>>lwrad>>rain;
        Time current(year,month,day,hour,minute);
		//current.addSeconds(noahlsm_timestep);
        if (m_time==current)
        {
            //return by the value?
            sfcspd=windspeed;
            sfcu=-windspeed*sin(winddir*3.14159265/180);
            sfcv=-windspeed*cos(winddir*3.14159265/180);
            sfctmp=temperature;
			sfcprs=pressure*100;

            float svp=svp1*exp((2.501e6-(4187.0-1870.0)*(temperature-svpt0))*(1.0/svpt0-1.0/temperature)/461.5);
            float qs=eps*svp/(sfcprs-(1-eps)*svp);
            float e=(sfcprs*svp*humidity*0.01)/(sfcprs-svp*(1-humidity*0.01));
            float spechumd=(eps*e)/(sfcprs-(1.0-eps)*e);
            if (spechumd<0.1e-5)
                spechumd=0.1e-5;
            if (spechumd>=qs)
                spechumd=qs*0.99;
            q2=spechumd;

            soldn=swrad;
            longwave=lwrad;
            prcp=rain;
            break;
        }
        else if (m_time>current)
        {

        }
        else
        {

        }
    }


}
