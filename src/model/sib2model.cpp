#include "sib2model.h"

using namespace ldas;
using namespace std;

SiB2Config& SiB2Config::operator=(const SiB2Config& other)
{
    if (this == &other) return *this; // handle self assignment
    //assignment operator
    forcing_path=other.forcing_path;
    vegetation_path=other.vegetation_path;
    output_path=other.output_path;
    return *this;
}

SiB2Model::SiB2Model()
{
    init();
}

SiB2Model::~SiB2Model()
{
    //dtor
}

SiB2Model::SiB2Model(const SiB2Model& other)
{
    //copy ctor
}

SiB2Model& SiB2Model::operator=(const SiB2Model& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

void SiB2Model::config(const SiB2Config& cfg)
{
	m_config=cfg;
}

SiB2Config SiB2Model::config() const
{
	return m_config;
}
/* ****************************************************************************
   This function is the constructor of the class sib.
   Functions called: veginc, cntrol, const2, driver, balan,
   inter2, rada2, begtem, endtem, updat2, balan

   ipbl = 0 : dtc, dtg only     ; as per SE-86
   ipbl = 1 : dtc, dtg, dth, dqm; as per Randall-Sellers, SA-89b

   isnow= 0 : conventional run using amazon forcings.
   isnow= 1 : snow shock test, warming to normal over 5 days.
**************************************************************************** */
void SiB2Model::run()
{

    const int ipbl	= 1;
    const int isnow	= 0;
    const int maxit	= 360;
    int iplace;

    init();
    // open general input file
    initParam(m_config.vegetation_path.c_str());

    // open the file containing forcing data
    ifstream iu(m_config.forcing_path.c_str());
    if (!iu)
        throw Exception("SiB2Model","run()","Can not open the forcing data file");

    // open/create the output file
    std::string fn=m_config.output_path+"out.txt";
    ofstream out(fn.c_str());
    fn=m_config.output_path+"out1.txt";
    ofstream out1(fn.c_str());
    fn=m_config.output_path+"out2.txt";
    ofstream out2(fn.c_str());
    fn=m_config.output_path+"out3.txt";
    ofstream out3(fn.c_str());
    fn=m_config.output_path+"out4.txt";
    ofstream out4(fn.c_str());

    for (iter=1; iter<=maxit; iter++)
    {
        const2();

        int nymd;
        driver(iu, isnow, nymd);
        balan(iplace = 1, nymd);

        inter2();
        rada2();
        begtem();
        endtem(ipbl);
        updat2();

        balan(iplace = 2, nymd );
        outer(out, out1, out2, out3, out4, nymd);
    }

    iu.close();
    out.close();
    out1.close();
    out2.close();
    out3.close();
    out4.close();
}

void SiB2Model::initParam(const char* filename)
{
    std::ifstream ichi;
    ichi.open (filename);
    if (!ichi)
        throw Exception("SiB2Model","initParam(const char* filename)","Can not open general input data file");

    /* ****************************************************************************
       read vegetation and soil parameters for SIB2.

       function called:		vegpar
    						soipar
    						dynveg
                            varcal
           subscripts:		(iv, iw, il):

    						iv  : surface layer ;
    						1 = canopy
    						2 = ground
    						iw  : radiation wavelength;
    						1 = visible
    						2 = near infrared
    						il  : vegetation state;
    						1 = live (green) and
    						2 = dead (stems and trunk)

       -----------------------------------------------------------------------
                     input parameter set
       -----------------------------------------------------------------------

       ivtype        : vegetation type

            static parameters associated with vegetation type
            -------------------------------------------------

       z2            : canopy top height
       z1            : canopy base height
       vcover        : vegetation cover fraction
       chil          : leaf angle distribution factor
       rootd         : rooting depth
       phc           : 1/2 critical leaf water potential limit
       tran(iw,il)   : leaf transmittance
       ref (iw,il)   : leaf reflectance
       effcon        : quantum efficiency
       gradm         : conductance-photosynthesis slope parameter
       binter        : conductance-photosynthesis intercept
       respcp        : respiration fraction of vmax
       atheta        : wc, we coupling parameter
       btheta        : wc & we, ws coupling parameter
       trda          : temperature coefficient in gs-a model
       trdm          : temperature coefficient in gs-a model
       trop          : temperature coefficient in gs-a model
       slti          : slope of low temperature inhibition function
       hlti          : 1/2 point of low temperature inhibition function
       shti          : slope of high temperature inhibition function
       hhti          : 1/2 point of high temperature inhibition function

       istype        : soil type

            static parameters associated with soil type
            -------------------------------------------

       sodep         : total depth of 3 soil moisture layers
       soref(iw)     : soil reflectance
       bee           : soil wetness exponent
       phsat         : soil tension at saturation
       satco         : hydraulic conductivity at saturation
       poros         : soil porosity
       slope         : cosine of mean slope

            time-space varying vegetation parameters, from
            spectral vegetation indice (svi).
            --------------------------------------------------

       zlt           : leaf area index
       green         : green leaf fraction
       fparc         : canopy absorbed fraction of photosynthetically
                     : active radiation (par)  from svi

            parameters derived from the above
            ---------------------------------

       vmax0         : rubisco velocity of sun-leaf
       gmudmu        : time-mean leaf projection ( g(mu)/ mu )
       z0d           : roughness length
       dd            : zero plane displacement
       cc1           : rb coefficient (c1)
       cc2           : rd coefficient (c2)
       zdepth        : individual depths of 3 soil moisture layers

       -----------------------------------------------------------------------

            other variables
            ---------------

          g1, g2, g3, ztz0, corb1, corb2, ha, zwind, zmet

       g1            : ratio of km(actual) to km(log-linear) at z = z2
       g2            : ratio of ra(actual) to ra(log-linear) for momentum
                       between: z = z2 and z = zx, where zx = min(zl,zwind)
       g3            : ratio of ra(actual) to ra(log-linear) for heat
                       between: z = z2 and z = zx, where zx = min(zl,zmet)
       ztz0          : parameter to determine depth of transition layer
                       above canopy, zl. zl = z2 + ztz0 * z0
       corb1         : non-neutral correction for calculation of aerodynami
                       resistance between ha and z2. when multiplied by
                       h*rbb/tm gives bulk estimate of local richardson
                       number
                       rbb = ra for heat between ha and z2.
                       corb2 = 9*g/( rhoair*cpair* (du/dz)**2 )
       corb2         : neutral value of rbb*u2 ( squared ), equivalent to
                       rdc**2 for upper canopy
       ha            : canopy source height for heat
       zwind         : reference height for wind measurement
       zmet          : reference height for temperature, humidity
                       measurement

            the above are generated from sibx + momopt output

    **************************************************************************** */
    char str[255];

    ichi.getline(str, 255);
    ichi.getline(str, 255);

    ichi >> ivtype;
    /* ****************************************************************************
       reading/setting of vegetation-type dependent static parameters
    **************************************************************************** */
    ichi.getline(str, 255);
    ichi.getline(str, 255);

    ichi >> z2 >> z1 >> vcover >> chil;
    ichi >> rootd >> phc;

    ichi >> tran[0][0] >> tran[1][0] >> tran[0][1] >> tran[1][1]
         >> ref[0][0] >> ref[1][0] >> ref[0][1] >> ref[1][1];
    ichi >> effcon >> gradm >> binter >> respcp >> atheta >> btheta;
    ichi >> trda >> trdm >> trop >> slti >> hlti >> shti >> hhti;
    ichi.getline(str, 255);
    ichi.getline(str, 255);
    ichi >> istype >> sodep >> soref[0] >> soref[1];
    /* ****************************************************************************
       reading/setting of soil-type dependent static parameters
    **************************************************************************** */
    ichi.getline(str, 255);
    ichi.getline(str, 255);
    ichi >> bee >> phsat >> satco >> poros >> slope;
    /* ****************************************************************************
       reading of phenological vegetation parameters
    **************************************************************************** */
    ichi.getline(str, 255);
    ichi.getline(str, 255);
    ichi >> fparc;
    /* ****************************************************************************
       calculation of secondary parameters from input data
    **************************************************************************** */
    ichi.getline(str, 255);
    ichi.getline(str, 255);
    ichi >> vmax0 >> gmudmu >> green >> zlt;
    ichi >> z0d >> dd >> cc1 >> cc2;
    ichi.getline(str, 255);
    ichi.getline(str, 255);
    ichi >> corb1 >> corb2 >> ha >> g1 >> g2 >> g3 >> ztz0 >> zwind >> zmet;

    rootd = amin1( rootd, sodep*0.75 );
    zdepth[0] = 0.02;
    zdepth[1] = rootd - 0.02;
    zdepth[2] = sodep - zdepth[0] - zdepth[1];

    double scatp = green * (tran[0][0] + ref[0][0]) +
                   (1.0-green) * (tran[0][1] + ref[0][1]);
    double park = sqrt(1.-scatp) * gmudmu;
    zlt = -1.0/park*log( 1.0-fparc );
    /* ****************************************************************************
       initialisation and switches
    **************************************************************************** */
    ichi.getline(str, 255);
    ichi.getline(str, 255);
    ichi >> dtt >> itrunk >> ilw;
    ichi >> zlat >> zlong >> time >> month >> day >> year >> niter;
    ichi >> tc >> tg >> td >> ta >> tm >> ht >> qa;
    ichi >> www[0] >> www[1] >> www[2];

    //if (iopt == 1) return;

    /*	cout << setfill(' ') << setw(10) << setfill('*') << setw(32) << endl;
    	cout << setfill(' ') << setw(10) << "* sib off-line test simulation *"
    		<< endl;
    	cout << setfill(' ') << setw(10) << setfill('*') << setw(32) << endl;
    	cout << "     " << setw(6) << setprecison(2) << latitude << setw(7)
    		<< setprecison(2) << longitude << endl;
    	cout << zlat << zlong;

    	if (itrunk == 1) cout << "resistances calculated from initial fluxes" << endl;
        if (itrunk == 2) cout << "resistances calculated by iteration, itrunk= "
    		<< setw(4) << itrunk << endl;
    	if (ilw == 1) cout << "downward longwave radiation read in as data" << endl;
    	if (ilw == 2) cout << "downward longwave radiation computed from brunts'
    		equation" << endl;
    	if (ilw == 3) cout << "downward longwave radiation computed as residual
    		in energy balance" << endl; */

    ichi.close();
}
/* ****************************************************************************
   init arrays
**************************************************************************** */
void SiB2Model::init()
{
    www[0] = www[1] = www[2] = 0.0;
    capac[0] = capac[1] = 0.0;
    snoww[0] = snoww[1] = 0.0;
    satcap[0] = satcap[1] = 0.0;
    rstfac[0] = rstfac[1] = rstfac[2] = rstfac[3] = 0.0;
}

/* ****************************************************************************
   initialization of physical constants

   const2 is called at every time step because in some applications
   constants may depend on environmental conditions.
**************************************************************************** */
void SiB2Model::const2()
{

    asnow    = 13.2;
    bps      = 1.0;
    clai     = 4.2 * 1000. * 0.2;
    cpair    = 1010.0;
    cw       = 4.2 * 1000.0 * 1000.0;
    epsfac   = 0.622;
    //g        = 9.81;
    kappa    = 0.286;
    //pie      = 3.14159265;
    po2m     = 20900.0;
    pco2m    = 34.0;
    psur     = 1000.0;
    rhoair   = 1.225;
    snomel   = 370518.5 * 1000.0;
    stefan   = 5.6697E-8;
    tf       = 273.16;
    vkc      = 0.41;
    rcp      = rhoair * cpair;
    timcon   = PI_NUMBER/86400.0;

    // n.b. :  snomel is expressed in j m-1
}


/* ****************************************************************************
   Reading of meteorological data. preparation of forcing variables
   data required to run sib :
   nymd     : date and time (yymmddhh)
   swdown   : shortwave downward radiation.
   zlwd     : longwave  downward radiation.
   em,tm,um : vapor pressure,temperature and wind speed
              at measurement height ( or lowest model level.)
   tprec    : precipitation.

   function called : radc2

   ****************************************************************************
   access meteorological data : manaus run measured data for comparison
   mevap          : latent heat flux
   msensh         : sensible heat flux
   mustar         : friction velocity
   iqcaws ,iqchyd : data quality check indices
**************************************************************************** */
void SiB2Model::driver(ifstream& iu, const int isnow, int& nymd)
{
    int iqcaws ,iqchyd;
    float zlwd, tprec, mevap, msensh, mustar;

    iu >> nymd >> iqcaws >> swdown >> rnetm >> zlwd >> em >> tm >> um
       >> tprec >> iqchyd >> mevap >> msensh >> mustar;

    // isnow = 0 : conventional run with received met data.
    // isnow = 1 : shock test with initial snow dump and freezing
    //			   temperatures at beginning of run warming to
    //             normal over  5-day period.

    if (isnow)
    {
        if (iter == 1)
        {
            tc = 270.0;
            tg = 270.0;
            snoww[0] = 0.1;
            snoww[1] = 0.1;
        }
        else
        {
            double cold = amax1 ( 0.0, (120.0 - (1.0*iter)) / 120.0 );
            double rhair = em/e(tm);
            tm = tm * ( 1. - cold ) + (tm - 30.) * cold;
            em = e(tm)*rhair;
            if (em < 0.0) em = 0.1;
        }
    }

    /* if(nymd < 83090300) goto 100;
    if(nymd > 85082523)
    {
    	write(icho, 900)iu, nymd, iout;
    	900   format(5x,'eof encountered for unit= ',i2,' eof date = ',i8,1x,i2);
    	return;
    } */

    um = amax1(um,0.25);
    // ustarm = mustar/100.0; //
    swdown = amax1(swdown,0.1);
    // ihour = mod(nymd,100); //
    // ptot = ptot + tprec; //
    ppl = tprec;
    ppc = tprec-ppl;

    radc2();

    cloud = (1160.*sunang - swdown) / (963. * sunang);
    cloud = amax1(cloud,0.);
    cloud = amin1(cloud,1.);
    cloud = amax1(0.58,cloud);

    double difrat = 0.0604 / ( sunang-0.0223 ) + 0.0683;
    if ( difrat < 0. ) difrat = 0.0;
    if ( difrat > 1. ) difrat = 1.0;

    difrat = difrat + ( 1. - difrat ) * cloud;
    double vnrat = ( 580. - cloud*464. ) / ( ( 580. - cloud*499. )
                   + ( 580. - cloud*464. ) );

    radn[0][0] = (1.-difrat)*vnrat*swdown;
    radn[0][1] = difrat*vnrat*swdown;
    radn[1][0] = (1.-difrat)*(1.-vnrat)*swdown;
    radn[1][1] = difrat*(1.-vnrat)*swdown;
    radn[2][0] = 0.0;
    radn[2][1] = zlwd;

    return;
}


/* ****************************************************************************
   solar zenith angle computation; downcoming radiation at bottom
**************************************************************************** */
void SiB2Model::radc2()
{
    double dayspy = 365.0;
    if (((int)year % 4) == 0) dayspy = 366.0;

    // julian day and time update; skip on 1st time step (initialized)
    if (iter > 1)
    {
        time = time + dtt / 3600.0;
        if ( time >= 23.99 ) time = 0.0;
        day = day +  dtt / 86400.0;
    }

    if ( day > dayspy ) year = year + 1.0;
    if ( day > dayspy ) day = day - dayspy;

    // solar declination calculation
    double decmax = PI_NUMBER * ( 23.5 / 180.);
    double sols   = ( 4141./24. ) + ( (int)(year+3) % 4 ) * 0.25;

    double season = ( day - sols ) / 365.2;
    double dec    = decmax * cos ( 2. * PI_NUMBER * season );

    double rfd  = PI_NUMBER / 180.0;
    double sind = sin( dec );
    double cosd = cos( dec );
    double hac  = -tan( zlat * rfd )*tan( dec );
    hac  = amin1(hac,1.0);
    hac  = amax1(hac,-1.0);

    // h is the half-day length (in radians)

    double h   = acos(hac);
    double dawn= -h;
    double dusk= +h;
    double sr  = 12.-(h/(15.*rfd));
    double ss  = 12.+(h/(15.*rfd));
    double coshr = cos( - PI_NUMBER + (time + 0.5*dtt/3600.) / 24. * 2. * PI_NUMBER );
    sunang = sin( zlat*rfd ) * sind + cos ( zlat*rfd ) * cosd * coshr;
    sunang = amax1( 0.01, sunang );

    return;
}


/* ****************************************************************************
   energy and water balance check
**************************************************************************** */
void SiB2Model::balan (int iplace, const int nymd)
{
    static double totwb;
    double endwb, errorw, pmeter, emeter;

    if (iplace == 1)
    {
        etmass = 0;
        roff = 0;
        totwb = www[0] * poros * zdepth[0]
                + www[1] * poros * zdepth[1]
                + www[2] * poros * zdepth[2]
                + capac[0] + capac[1] + snoww[0] + snoww[1];
    }

    else if (iplace == 2)
    {
        endwb = www[0] * poros * zdepth[0]
                + www[1] * poros * zdepth[1]
                + www[2] * poros * zdepth[2]
                + capac[0] + capac[1] + snoww[0] + snoww[1]
                - (ppl+ppc)/1000. + etmass/1000. + roff;

        errorw = totwb - endwb;
        errorw = totwb - endwb;
        pmeter = (ppl+ppc)/1000.0;
        emeter = etmass/1000.0;

        if(fabs(errorw) > 0.0001)
        {
            printf("\n");
            /* .&write(icho,900) nymd, totwb, endwb, errorw,.
            .&. www(1), www(2), www(3),.
            .&. capac(1), capac(2), snoww(1), snoww(2),.
            .&. pmeter, emeter, roff.
            900 format(//,10x,'** warning: water balance violation **',//,.
            .&/,1x,'date ', i8,.
            .&/,1x,'begin, end, diff ', 3(f10.7,1x),.
            .&/,1x,'www,1-3. ', 3(f10.8,1x),.
            .&/,1x,'capac,1-2. ', 2(f10.8,1x),.
            .&/,1x,'snoww,1-2. ', 2(f10.8,1x),.
            .&/,1x,'p, et, roff. ', 3(f10.8,1x) ). */
        }
        double cbal = radt[0] - chf - (ect+hc+eci)/dtt;
        double gbal = radt[1] - shf - (egs+hg+egi)/dtt - heaten/dtt;
        double zlhs = radt[0] + radt[1] - chf - shf;
        double zrhs = hflux + (ect + eci + egi + egs)/dtt + heaten/dtt;

        double errore= zlhs - zrhs;

        if(fabs(errore) > 1.0)
        {
            printf("\n");


            /* .& write(icho,910) nymd, zlhs, zrhs, radt(1), radt(2), chf, shf,.
            .& hflux, ect, eci, egi, egs, hc, hg, heaten, cbal, gbal.
            910 format(//,10x,'** warning: energy balance violation **',//,.
            .& /,1x,'date ', i8,.
            .& /,1x,'rhs, lhs.', 2g12.5,.
            .& /,1x,'rn1, rn2, chf, shf, h ', 5g12.5,.
            .& /,1x,'ect, eci, egi, egs', 4g12.5,.
            .& /,1x,'hc. hg. ',g12.5, 12x, g12.5,.
            .& /,1x,'heaten, c-bal, g-bal', 3g12.5 ).; */
        }
    }

    return;
}


/* ****************************************************************************
   calculation of  interception and drainage of rainfall and snow
   incorporating effects of patchy snow cover and temperature
   adjustments.

   (1) non-uniform precipitation convective ppn. is described
   by area-intensity relationship :-

                f(x) = a*exp(-b*x)+c

   throughfall, interception and infiltration excess are functional on this
   relationship and proportion of large-scale ppn.
   reference: sato et al.(1989b), appendix.

   (2) reorganisation of snowmelt and rain-freeze procedures.
               function adjust

   (3) additional calculation for partial snow-cover case.
               function patchs

   (4) reorganisation of overland flow.
         reference: SA-89B, appendix.

   (5) modified calaculation of soil heat capacity and conductivity.

   *************************************************************************
   subroutines in this block :	snow1
								adjust
								patchs
								snow1

   *************************************************************************
    roff           runoff (mm)
	tc             canopy temperature (K)
	tg             ground surface temperature (K)
	www[0]         ground wetness of surface layer
	capac[1]       canopy/ground liquid interception store (m)
	snoww[1]       canopy/ground snow interception store (m)

**************************************************************************** */
void SiB2Model::inter2()
{
    const double pcoefs[2][2] = { 20.0, 0.206E-8, 0.0001, 0.9999 };
    const double bp = 20.0;

    snow1();

    // prec ( pi-x )   : equation (c.3), SA-89B
    double ap = pcoefs[1][0];
    double cp = pcoefs[1][1];
    double totalp = ppc + ppl;
    if( snoww[0] > 0.0 || snoww[1] > 0.0 || tm < tf )
        ppc = 0.0;
    ppl = totalp - ppc;
    if(totalp >= 1.0E-8)
    {
        ap = ppc/totalp * pcoefs[0][0] + ppl/totalp * pcoefs[1][0];
        cp = ppc/totalp * pcoefs[0][1] + ppl/totalp * pcoefs[1][1];
    }

    roff = 0.0;
    double thru = 0.0;
    double fpi  = 0.0;

    // heat capacity of the soil, as used in force-restore heat flux
    // description. dependence of csoil on porosity and wetness is
    // based on CS-81.
    double slamda = ( 1.5*(1.-poros) + 1.3*www[0]*poros ) /
                    ( 0.75 + 0.65*poros - 0.4*www[0]*poros ) * 0.4186;
    double shcap  = ( 0.5*(1.-poros) + www[0]*poros ) * 4.186 * 1.e6;
    csoil  = sqrt( slamda * shcap * 86400./PI_NUMBER ) / 2.0;

    // input precipitation is given in mm, converted to m to give p0.
    double p0 = totalp * 0.001;

    for (int iveg = 0; iveg<2; iveg++)
    {
        double realc = (double)(1.0 - iveg);
        double realg = (double)iveg;

        double capacp = capac[iveg];
        double snowwp = snoww[iveg];

        double xsc = amax1(0.0, capac[iveg] - satcap[iveg]);
        capac[iveg] = capac[iveg] - xsc;
        double xss = amax1(0.0, snoww[iveg] - satcap[iveg]) * realc;
        snoww[iveg] = snoww[iveg] - xss;
        roff = roff + xsc + xss;

        double spechc = amin1( 0.05, (capac[iveg] + snoww[iveg])) * cw
                        + realc * zlt * clai + realg * csoil;
        double ts = tc * realc + tg * realg;

        // proportional saturated area (xs) and leaf drainage(tex)
        // tex ( d-c )     : equation (c.8), SA-89B
        // xs  ( x-s )     : equation (c.7), SA-89B
        // tex ( d-c )     : equation (c.8), SA-89B

        double chiv = chil;
        if ( fabs(chiv) <= 0.01 ) chiv = 0.01;
        double aa = 0.5 - 0.633 * chiv - 0.33 * chiv * chiv;
        double bb = 0.877 * ( 1. - 2. * aa );
        exrain = aa + bb;

        double zload = capac[iveg] + snoww[iveg];
        fpi = ( 1.0-exp( - exrain*zlt/vcover ) )*vcover*realc + realg;
        double tti = p0 * ( 1.0-fpi );
        double xs = 1.0;
        if ( p0 >= 1.e-9 )
        {
            double arg =  ( satcap[iveg]-zload )/( p0*fpi*ap ) -cp/ap;
            if ( arg >= 1.e-9 )
            {
                xs = -1./bp * log( arg );
                xs = amin1( xs, 1. );
                xs = amax1( xs, 0. );
            }
        }

        double tex = p0*fpi * ( ap/bp*( 1.- exp( -bp*xs )) + cp*xs ) -
                     ( satcap[iveg] - zload ) * xs;
        tex = amax1( tex, 0. );

        // total throughfall (thru) and store augmentation

        if ( iveg == 0 )
        {
            thru = tti + tex;
            double pinf = p0 - thru;
            if ( tm > tf ) capac[iveg] = capac[iveg] + pinf;
            if ( tm <= tf ) snoww[iveg] = snoww[iveg] + pinf;

            adjust ( tc, spechc, capacp, snowwp, iveg );

            p0 = thru;
            continue;
        }

        if ( tg > tf && snoww[1] > 0.0 )
        {
            patchs ( p0 );
            continue;
        }

        thru = tti + tex;
        if ( tg <= tf || tm <= tf ) thru = 0.0;
        double pinf = p0 - thru;
        if ( tm > tf ) capac[iveg] = capac[iveg] + pinf;
        if ( tm <= tf ) snoww[iveg] = snoww[iveg] + pinf;
        if ( tm > tf )
        {
            // instantaneous overland flow contribution ( roff )
            // roff( r-i )     : equation (c.13), SA-89B

            double equdep = satco * dtt;

            xs = 1.0;
            if ( thru >= 1.e-9 )
            {
                double arg = equdep / ( thru * ap ) -cp/ap;
                if ( arg >= 1.e-9 )
                {
                    xs = -1./bp * log( arg );
                    xs = amin1( xs, 1. );
                    xs = amax1( xs, 0. );
                }
            }
            double roffo = thru * ( ap/bp * ( 1.-exp( -bp*xs )) + cp*xs )
                           -equdep*xs;
            roffo = amax1 ( roffo, 0. );
            roff = roff + roffo;
            www[0] = www[0] + (thru - roffo) / (poros*zdepth[0]);
        }

        adjust ( tg, spechc, capacp, snowwp, iveg );
    }

    return;
}


/* ****************************************************************************
   calculation of effects of snow cover on surface morphology
   and maximum water storage values.

   z0             roughness length (m)
   d              zero plane displacement (m)
   rbc            rb coefficient (c1) (s m-1)**1/2
   rdc            rd coefficient (c2)
   satcap(2)      interception capacities (m)
   canex          fraction of exposed canopy (snow-free)
   areas          ground snow cover fraction

   *************************************************************************
   alteration of aerodynamic transfer properties in case of snow
   accumulation. calculation of maximum water storage values.

   canex       (fraction of canopy not covered by snow)
   d           (snow-modified value of dd, used in all calculations)
   z0          (snow-modified value of z0d, used in all calculations)
   rbc         (snow-modified value of cc1, used in all calculations)
   rdc         (snow-modified value of cc2, used in all calculations)
   areas       (fraction of ground covered by snow)
   satcap(1)   (s-c)   : equation (56) , SE-86, page 741 se-89
   satcap(2)   (s-g)   : 0.002, surface interception store
**************************************************************************** */
void SiB2Model::snow1()
{
    canex  = 1.-( snoww[1]*5.-z1)/(z2-z1);
    canex  = amax1( 0.1, canex );
    canex  = amin1( 1.0, canex );

    d      = z2 - ( z2-dd ) * canex;
    z0     = z0d/( z2-dd ) * ( z2-d );
    rbc    = cc1/canex;
    rdc    = cc2*canex;
    areas    = amin1(1., asnow*snoww[1]);
    satcap[0] = zlt*0.0001 * canex;
    // 3/96 changes
    // old     satcap(2) = 0.002
    satcap[1] = 0.0002;

    return;
}


/* ****************************************************************************
   temperature change due to addition of precipitation

   roff           runoff (mm)
   tc             canopy temperature (K)
   tg             ground surface temperature (K)
   www[0]         ground wetness of surface layer
   capac[1]       canopy/ground liquid interception store (m)
   snoww[1]       canopy/ground snow interception store (m)
**************************************************************************** */
void SiB2Model::adjust(double& ts, const double spechc, const double capacp,
                       const double snowwp, const int iveg)
{
    double freeze = 0.0;
    double diff = ( capac[iveg]+snoww[iveg] - capacp-snowwp )*cw;
    double ccp = spechc;
    double cct = spechc + diff;

    double tsd = ( ts * ccp + tm * diff ) / cct;
    if (!(ts > tf && tm > tf))
    {
        if (!(ts <= tf && tm <= tf))
        {
            double tta = ts;
            double ttb = tm;
            double cca = ccp;
            double ccb = diff;
            if ( tsd <= tf )
            {
                // freezing of water on canopy or ground
                double ccc = capacp * snomel;
                if ( ts < tm ) ccc = diff * snomel / cw;
                tsd = ( tta * cca + ttb * ccb + ccc ) / cct;

                freeze = ( tf * cct - ( tta * cca + ttb * ccb ) );
                freeze = (amin1 ( ccc, freeze )) / snomel;
                if(tsd > tf) tsd = tf - 0.01;
            }
            else
            {
                // melting of snow on canopy or ground, water infiltrates.
                double ccc = - snoww[iveg] * snomel;
                if ( ts > tm ) ccc = - diff * snomel / cw;

                tsd = ( tta * cca + ttb * ccb + ccc ) / cct;

                freeze = ( tf * cct - ( tta * cca + ttb * ccb ) );
                freeze = (amax1( ccc, freeze )) / snomel;
                if(tsd <= tf) tsd = tf - 0.01;
            }
        }
    }

    snoww[iveg] = snoww[iveg] + freeze;
    capac[iveg] = capac[iveg] - freeze;

    double xs = amax1( 0.0, ( capac[iveg] - satcap[iveg] ) );
    if( snoww[iveg] >= 0.0000001 ) xs = capac[iveg];
    www[0] = www[0] + xs / ( poros * zdepth[0]);
    capac[iveg] = capac[iveg] - xs;
    ts = tsd;

    return;
}


/* ****************************************************************************
   marginal situation: snow exists in patches at temperature tf
   with remaining area at temperature tg > tf.

   calculation of effect of intercepted snow and rainfall on ground.
   patchy snowcover situation involves complex treatment to keep
   energy conserved.
**************************************************************************** */
void SiB2Model::patchs(const double p0)
{

    double ex, tsd;

    double pinf = p0;
    double thru = 0.0;
    double snowhc = amin1( 0.05, snoww[1] ) * cw;
    areas = amin1( 1.0,(asnow*snoww[1]) );

    if( tm <= tf )			// snow falling onto area
    {
        double rhs = tm*pinf*cw + tf*(snowhc + csoil*areas)
                     + tg*csoil*(1.-areas);
        double dareas = amin1( asnow*pinf, ( 1.-areas ) );
        ex = rhs - tf*pinf*cw - tf*(snowhc + csoil*(areas + dareas))
             - tg*csoil*(1.-areas-dareas);
        if( (areas+dareas) >= 0.999 ) tg = tf - 0.01;

        if( ex >= 0.0 )
        {
            // excess energy is positive, some snow melts and infiltrates.
            double zmelt = ex/snomel;
            if( asnow*(snoww[1] + pinf - zmelt) <= 1.0 )
            {
                zmelt = 0.0;
                if( asnow*(snoww[1] + pinf) >= 1.0 )
                    zmelt = ( asnow*(snoww[1] + pinf) - 1. ) / asnow;
                zmelt = ( ex - zmelt*snomel )/( snomel + asnow*csoil*(tg-tf) ) + zmelt;
            }
            snoww[1] =  snoww[1] + pinf - zmelt;
            www[0] = www[0] + zmelt/(poros*zdepth[0]);
        }

        else
        {
            // excess energy is negative, bare ground cools to tf, then whole
            // area cools together to lower temperature.
            tsd = 0.0;
            if( (areas+dareas) <= 0.999 )
                tsd = ex/(csoil*( 1.-areas-dareas)) + tg;
            if( tsd <= tf )
                tsd = tf + (ex - (tf-tg)*csoil*(1.-areas-dareas))
                      / (snowhc+pinf*cw+csoil);
            tg = tsd;
            snoww[1] = snoww[1] + pinf;
        }
    }

    else			// rain falling onto area
    {
        // rain falls onto snow-free sector first.
        tsd = tf - 0.01;
        if ( areas < 0.999 ) tsd = (tm*pinf*cw + tg*csoil) / (pinf*cw + csoil);
        tg = tsd;
        www[0]= www[0]+pinf*(1.-areas)/(poros*zdepth[0]);

        // rain falls onto snow-covered sector next.
        ex = ( tm - tf )*pinf*cw*areas;
        double dcap = -ex / ( snomel + ( tg-tf )*csoil*asnow );
        if( (snoww[1] + dcap) >=  0.0 )
        {
            www[0] = www[0]+(pinf*areas-dcap)/(poros*zdepth[0]);
            snoww[1] = snoww[1] + dcap;
        }
        else
        {
            tg = ( ex - snomel*snoww[1] - ( tg-tf )*csoil*areas ) / csoil + tg;
            www[0]=www[0]+(snoww[1]+pinf*areas)/(poros*zdepth[0]);
            capac[1] = 0.0;
            snoww[1] = 0.0;
        }
    }

    return;
}


/* ****************************************************************************
   calculation of albedos via two stream approximation (direct and diffuse )
   and partition of radiant energy

   functions  called  :	snow1
						longrn

   ********************   output    ***********************
   salb[2][2]		surface albedos
   tgeff			effective surface radiative temperature (k)
   radfac[2][2][2]	radiation absorption factors
   thermk			canopy gap fraction for tir radiation

   ********************  diagnostics  **********************
   albedo[2][2][2]  component reflectances
   closs			tir emission from the canopy (w m-2)
   gloss			tir emission from the ground (w m-2)
**************************************************************************** */
void SiB2Model::rada2()
{
    double tranc1[2], tranc2[2], tranc3[2];
    double zmew;

    double f = sunang;

    // modification for effect of snow on upper story albedo
    // snow reflectance   = 0.80, 0.40 . multiply by 0.6 if melting
    // snow transmittance = 0.20, 0.54
    snow1();
    double facs  = ( tg-tf ) * 0.04;
    facs  = amax1( 0.0, facs);
    facs  = amin1( 0.4, facs);
    double fmelt = 1.0 - facs;

    for (int iwave = 0; iwave <2; iwave++)
    {

        double scov =  amin1( 0.5, snoww[0]/satcap[0]);
        double reff1 = (1.0 - scov) * ref[iwave][0] + scov * ( 1.2 -
                       (iwave+1) * 0.4 ) * fmelt;
        double reff2 = (1.0 - scov) * ref[iwave][1] + scov * ( 1.2 -
                       (iwave+1) * 0.4 ) * fmelt;
        double tran1 = tran[iwave][0] * (1.0 - scov) +
                       scov * ( 1.0- ( 1.2 - (iwave+1) * 0.4 ) * fmelt ) * tran[iwave][0];
        double tran2 = tran[iwave][1] * (1.0 - scov)  +
                       scov * ( 1.0- ( 1.2 - (iwave+1) * 0.4 ) * fmelt ) * 0.9
                       * tran[iwave][1];

        /* ***************************************************************
        calculate average scattering coefficient, leaf projection and
        other coefficients for two-stream model.

          scat  (omega)        : equation (1,2) , SE-85
          proj  (g(mu))        : equation (13)  , SE-85
          extkb (k, g(mu)/mu)  : equation (1,2) , SE-85
          zmew  (int(mu/g(mu)) : equation (1,2) , SE-85
          acss  (a-s(mu))      : equation (5)   , SE-85
          extk  (k, various)   : equation (13)  , SE-85
          upscat(omega-beta)   : equation (3)   , SE-85
          betao (beta-0)       : equation (4)   , SE-85
          psi   (h)            : appendix       , SE-85
        *************************************************************** */

        double scat = green*( tran1 + reff1 ) +( 1. - green ) * ( tran2 + reff2);
        double chiv = chil;

        if ( fabs(chiv) <= 0.01 ) chiv = 0.01;
        double aa = 0.5 - 0.633 * chiv - 0.33 * chiv * chiv;
        double bb = 0.877 * ( 1. - 2. * aa );

        double proj = aa + bb * f;
        double extkb = ( aa + bb * f ) / f;
        zmew = 1. / bb * ( 1. - aa / bb * log ( ( aa + bb ) / aa ) );
        double acss = scat / 2. * proj / ( proj + f * bb );
        acss = acss * ( 1. - f * aa / ( proj + f * bb ) *
                        log ( (proj + f * bb + f * aa) / (f * aa) ) );

        double upscat = green * tran1 + ( 1. - green ) * tran2;
        upscat = 0.5 * ( scat + ( scat - 2.0 * upscat )
                         * pow(( 1.0 - chiv ) / 2.0, 2 ));
        double betao = ( 1. + zmew * extkb ) / ( scat * zmew * extkb ) * acss;

        /* ***************************************************************
        intermediate variables identified in appendix of SE-85.

          be          (b)     : appendix      , SE-85
          ce          (c)     : appendix      , SE-85
          bot         (sigma) : appendix      , SE-85
          hh1         (h1)    : appendix      , SE-85
          hh2         (h2)    : appendix      , SE-85
          hh3         (h3)    : appendix      , SE-85
          hh4         (h4)    : appendix      , SE-85
          hh5         (h5)    : appendix      , SE-85
          hh6         (h6)    : appendix      , SE-85
          hh7         (h7)    : appendix      , SE-85
          hh8         (h8)    : appendix      , SE-85
          hh9         (h9)    : appendix      , SE-85
          hh10        (h10)   : appendix      , SE-85
          psi         (h)     : appendix      , SE-85
          zat         (l-t)   : appendix      , SE-85
          epsi        (s1)    : appendix      , SE-85
          ek          (s2)    : appendix      , SE-85
        *************************************************************** */

        double be = 1. - scat + upscat;
        double ce = upscat;
        double bot = ( zmew * extkb ) * ( zmew * extkb ) + ( ce * ce - be * be );
        if ( fabs(bot) <= 1.e-10)
        {
            scat = scat* 0.98;
            be = 1. - scat + upscat;
            bot = ( zmew * extkb ) * ( zmew * extkb ) + ( ce * ce - be * be );
        }

        double de = scat * zmew * extkb * betao;
        double fe = scat * zmew * extkb * ( 1. - betao );
        double hh1 = - de * be + zmew * de * extkb - ce * fe;
        double hh4 = - be * fe - zmew * fe * extkb - ce * de;
        double psi = sqrt(be * be - ce * ce)/zmew;
        double zat = zlt/vcover*canex;

        double power1 = amin1( psi*zat, 50.0 );
        double power2 = amin1( extkb*zat, 50.0 );
        double epsi = exp( - power1 );
        double ek = exp ( - power2 );

        albedo[1][iwave][0] = soref[iwave]*(1.0-areas)
                              + ( 1.2-(iwave+1)*0.4 )*fmelt * areas;
        albedo[1][iwave][1] = soref[iwave]*(1.0-areas)
                              + ( 1.2-(iwave+1)*0.4 )*fmelt * areas;
        double ge = albedo[1][iwave][0]/albedo[1][iwave][1];

        // calculation of diffuse albedos
        // albedo[0][ir][1] ( i-up ) : appendix , SE-85

        double f1 = be - ce / albedo[1][iwave][1];
        double zp = zmew * psi;
        double den = (be + zp) * (f1 - zp) / epsi - (be - zp) * (f1 + zp) * epsi;
        double hh7   = ce * ( f1 - zp ) / epsi / den;
        double hh8  = -ce * ( f1 + zp ) * epsi / den;

        f1 = be - ce * albedo[1][iwave][1];
        den = ( f1 + zp ) / epsi - ( f1 - zp ) * epsi;
        double hh9   = ( f1 + zp ) / epsi / den;
        double hh10  = - ( f1 - zp ) * epsi / den;
        tranc2[iwave] = hh9 * epsi + hh10 / epsi;

        albedo[0][iwave][1] =  hh7 + hh8;

        // calculation of direct albedos and canopy transmittances.
        // albedo[0][iw][0] ( i-up ) : equation(11)   , SE-85
        // tranc [iw]   ( i-down ) : equation(10)   , SE-85

        f1 = be - ce / albedo[1][iwave][1];
        double zmk = zmew * extkb;
        den = ( be + zp ) * ( f1 - zp ) / epsi - ( be - zp ) * ( f1 + zp ) * epsi;
        double hh2 = ( de - hh1/bot * ( be + zmk ) ) * ( f1 - zp ) / epsi -
                     ( be - zp ) * ( de - ce*ge - hh1/bot * ( f1 + zmk ) ) * ek;
        hh2 = hh2 / den;
        double hh3 = ( be + zp ) * (de - ce*ge - hh1/bot * ( f1 + zmk ))* ek -
                     ( de - hh1/bot * ( be + zmk ) ) * ( f1 + zp ) * epsi;
        hh3 = hh3 / den;

        f1 = be - ce * albedo[1][iwave][1];
        den = ( f1 + zp ) / epsi - ( f1 - zp ) * epsi;
        double hh5 = - hh4/bot * ( f1 + zp ) / epsi -
                     ( fe + ce*ge*albedo[1][iwave][1] + hh4/bot*( zmk-f1 ) ) * ek;
        hh5 = hh5 / den;
        double hh6 = hh4/bot * ( f1 - zp ) * epsi +
                     ( fe + ce*ge*albedo[1][iwave][1] + hh4/bot*( zmk-f1 ) ) * ek;
        hh6 = hh6 / den;
        tranc1[iwave] = ek;
        tranc3[iwave] = hh4/bot * ek + hh5 * epsi + hh6 / epsi;

        albedo[0][iwave][0] = hh1/bot + hh2 + hh3;

        // calculation of terms which multiply incoming short wave fluxes
        // to give absorption of radiation by canopy and ground
        // radfac (f(il,imu,iv)) : equation (19,20) , SE-86

        radfac[1][iwave][0] = ( 1.-vcover ) * ( 1.-albedo[1][iwave][0] )
                              + vcover * ( tranc1[iwave] * ( 1.-albedo[1][iwave][0] )
                                           + tranc3[iwave] * ( 1.-albedo[1][iwave][1] ) );

        radfac[1][iwave][1] = ( 1.-vcover ) * ( 1.-albedo[1][iwave][1] )
                              + vcover *  tranc2[iwave] * ( 1.-albedo[1][iwave][1] );

        radfac[0][iwave][0] = vcover * ( ( 1.-albedo[0][iwave][0] )
                                         - tranc1[iwave] * ( 1.-albedo[1][iwave][0] )
                                         - tranc3[iwave] * ( 1.-albedo[1][iwave][1] ) );

        radfac[0][iwave][1] = vcover * ( ( 1.-albedo[0][iwave][1] )
                                         - tranc2[iwave] * ( 1.-albedo[1][iwave][1] ) );

        // calculation of total surface albedos ( salb ) with weighting
        // for cover fractions.

        for (int irad = 0; irad <2; irad++)
        {
            salb[iwave][irad] = ( 1.-vcover ) * albedo[1][iwave][irad] +
                                vcover * albedo[0][iwave][irad];
        }
    }

    //  calculation of long-wave flux terms from canopy and ground
    //  closs ( fc - rnc )     : equation (21),  SE-86
    //  gloss ( fg - rng )     : equation (22),  SE-86

    tgs = amin1(tf,tg)*areas + tg*(1.-areas);
    double tc4 = tc  * tc  * tc  * tc;
    double tg4 = tgs * tgs * tgs * tgs;

    double zkat = 1.0/zmew * zlt / vcover;
    zkat = amin1( 50. , zkat );
    zkat = amax1( 1.e-5, zkat );
    thermk = exp(-zkat);

    double fac1 =  vcover * ( 1.-thermk );
    double fac2 =  1.0;
    double closs =  2. * fac1 * stefan * tc4;
    closs =  closs - fac2 * fac1 * stefan * tg4;
    double gloss =  fac2 * stefan * tg4;
    gloss =  gloss - fac1 * fac2 * stefan * tc4;
    zlwup =  fac1 * stefan * tc4 + (1. - fac1 ) * fac2 * stefan * tg4;

    longrn(tranc1, tranc2, tranc3);

    // calculation of absorption of radiation by surface

    radt[0] = 0.0;
    radt[1] = 0.0;

    for (int iveg = 0; iveg < 2; iveg++)
        for (int iwave = 0; iwave < 2; iwave++)
            for (int irad = 0; irad < 2; irad++)
                radt[iveg] += radfac[iveg][iwave][irad]*radn[iwave][irad];
    radt[0] += radn[2][1]*vcover*(1.0- thermk) - closs;
    radt[1] += radn[2][1]*(1.0-vcover*(1.0-thermk) ) - gloss;

    return;
}


/* ****************************************************************************
   calculation of downward longwave. this is not required in gcm if
   downward longwave is provided by gcm-radiation code as radn[2][1].
**************************************************************************** */
void SiB2Model::longrn(double tranc1[2], double tranc2[2], double tranc3[2])
{
    // downward long-wave assumed to be provided as radn[2][1]
    if (ilw == 1) return;

    // downward long-wave from brunt's equation, Monteith(1973), p37.
    if (ilw == 2)
    {
        double esky = 0.53 + 0.06*sqrt(em);
        radn[2][1]  =  esky*(1.+0.2*(cloud*cloud))*stefan*pow(tm,4);
        return;
    }

    // downward long-wave flux calculated as residual from measured
    // net radiation and outgoing longwave radiation.
    // calculation of absorbed fractions of radiation ( expendable )
    if (ilw == 3)
    {
        for (int iwave = 0; iwave < 2; iwave++)
        {
            rab[1][iwave][0] =  (1.0 - vcover)
                                * ( radn[iwave][0] * (1.0 - albedo[1][iwave][0] ) );
            rab[1][iwave][1] =  ( 1. - vcover ) *
                                radn[iwave][1] * ( 1. - albedo[1][iwave][1] );

            rab[1][iwave][0] = rab[1][iwave][0] + vcover *
                               (radn[iwave][0] * (tranc1[iwave] * (1.0 - albedo[1][iwave][0])
                                                  + tranc3[iwave] * ( 1. - albedo[1][iwave][1]) ) );
            rab[1][iwave][1] = rab[1][iwave][1] + vcover *
                               radn[iwave][1] * tranc2[iwave] * ( 1. - albedo[1][iwave][1] );

            rab[0][iwave][0] =  vcover *
                                radn[iwave][0] * ( ( 1. - albedo[0][iwave][0] ) -
                                                   tranc1[iwave] * ( 1. - albedo[1][iwave][0] ) -
                                                   tranc3[iwave] * ( 1. - albedo[1][iwave][1] ) );
            rab[0][iwave][1] =  vcover *
                                radn[iwave][1] * ( ( 1. - albedo[0][iwave][1] ) -
                                                   tranc2[iwave] * ( 1. - albedo[1][iwave][1] ) );
        }

        double swab = rab[0][0][0] + rab[0][0][1] + rab[0][1][0] + rab[0][1][1] +
                      rab[1][0][0] + rab[1][0][1] + rab[1][1][0] + rab[1][1][1];
        double swup = swdown - swab;
        radn[2][1] = rnetm - swab + zlwup;

        return;
    }
}


/* ****************************************************************************
   core routine: calculation of canopy and ground temperature increments
   over time step, fluxes derived.

   function called : snow1

   *************************************************************************
	rstfac(2)      soil moisture stress factor
	rsoil          soil surface resistance (s m-1)
	hr             soil surface relative humidity
	wc             canopy wetness fraction
	wg             ground wetness fraction
	ccx            canopy heat capacity (j m-2 k-1)
	cg             ground heat capacity (j m-2 k-1)
**************************************************************************** */
void SiB2Model::begtem()
{
    snow1();

    hlat     = ( 3150.19 - 2.378 * tm ) * 1000.0;
    psy      = cpair / hlat * psur / .622;
    snofac   = hlat/ ( hlat + snomel / 1.e3 );

    // calculation of canopy and ground heat capacities.
    // n.b. this specification does not necessarily conserve energy when
    // dealing with very large snowpacks.

    ccx = zlt * clai + (snoww[0]+capac[0])*cw;
    cg  = csoil + amin1 ( 0.05, (snoww[1]+capac[1]) ) * cw;

    // calculation of ground surface temperature and wetness fractions

    tsnow = amin1 ( tf-0.01, tg );
    rsnow = snoww[1] / (snoww[1]+capac[1]+1.e-10);

    tgs = tsnow*areas + tg*(1.-areas);

    etc   = e(tc);
    etgs  = e(tgs);
    getc  = ge(tc);
    getgs = ge(tgs);

    wc = amin1( 1.,( capac[0] + snoww[0])/satcap[0] );
    wg = amax1( 0.,  capac[1]/satcap[1] )*0.25;

    /* ***************************************************************
    calculation of soil moisture stress factor.
    average soil moisture potential in root zone (layer-2) used as
    source for transpiration.

      phroot      (psi-r) : equation (47) , SE-86
      rstfac(2)  f(psi-l) :    "     (12) , SE-89
    *************************************************************** */

    double phroot = phsat * pow(amax1( 0.02, www[1] ),( - bee ));
    phroot = amax1 ( phroot, -2.e3 );
    rstfac[1] = 1./( 1 + exp( 0.02*( phc-phroot) ) );
    rstfac[1] = amax1( 0.0001, rstfac[1] );
    rstfac[1] = amin1( 1.,     rstfac[1] );

    /* ***************************************************************
    rsoil function from fit to FIFE-87 data.  soil surface layer
    relative humidity.

      rsoil      (rsoil) : SE-92B (personal communication)
      hr         (fh)    : equation (66) , SE-86
    *************************************************************** */

    double fac = amin1( www[0], 1. );
    fac = amax1( fac, 0.02  );

    // 3/96 changes
    // old      rsoil =  amax1 (0.1, 694. - fac*1500.) + 23.6
    rsoil =  exp (8.206 - 4.255*fac );

    double psit = phsat * pow(fac, - bee);
    double argg = amax1(-10.,(psit*g/461.5/tgs));
    hr = exp(argg);

    return;
}


/* ****************************************************************************
   calculation of ea, ta, ra, rb, rd and soil moisture stress for the beginning
   of the time step

                        modifications

     (1)     : change in cog1,cog2,cogs1 to allow soil evaporation
	           from beneath ground cover canopy.

     (2)     : change in temperature tgs to account for patchy snow.
               a bulk (area-weighted) temperature is substituted for
               all energy budget calculations. energy used to warm
               exposed patches is subtracted from snow evaporation.

     (3)     : inclusion of randall-sellers backward implicit scheme
               for calculating dtc, dtg, dth, dqm. option remains to
               use original dtc, dtg scheme only using parameter ipbl.

     (4)     : inclusion of integrated canopy  photosynthesis -
               conductance model. note that soil moisture stress is
               modelled as a chronic condition, while the relative
               humidity term is solved within a feedback loop.
               reference : SE-92A

	***********************************************************************

     functions called :    rasite --> unstab,stab,rafcal
     -------------------   rbrd
                           phosib --> cycalc-sortin
                           delrn
                           delhf
                           delef
                           sibslv --> gauss
                           dtcdtg
                           newton
	***********************************************************************
	output:

       ect            canopy transpiration (j m-2)
       eci            canopy interception loss (j m-2)
       egs            ground evaporation (j m-2)
       egi            ground interception loss (j m-2)
       ec             ect + eci
       eg             egs + egi
       hc             canopy sensible heat flux (j m-2)
       hg             ground sensible heat flux (j m-2)
       chf            canopy heat storage flux (j m-2)
       shf            soil heat storage flux (j m-2)
       fc             canopy dew indicator
       fg             ground dew indicator
       heaten         heat loss to snow melt process (j m-2)
	***********************************************************************
	diagnostics:
       tsnow          snow temperature (k)
**************************************************************************** */
void SiB2Model::endtem (const int ipbl)
{
    double ecpot, ecidif, egpot, egidif;
    double taen;

    int ifirst = 1;

    fc = 1.0;
    fg = 1.0;
    ta = (tgs+tc)/2.0;
    ea = em;
    ht = 0.0;

    for (int count = 0; count <= 4; count++)
    {
        rasite();
        rbrd();
        phosib();
    }

    ifirst = 0;
    delrn();

    // dew calculation : dew condition is set at beginning of time step.
    // if surface changes state during time step, latent heat flux is set
    // to zero.
    if ( ea > etc ) fc = 0.0;
    if ( ea > etgs) fg = 0.0;

    // start of non-neutral resistance calculation loop
    // initialize newton-raphson iterative routine for ra calculation
    int nox = 0;
    int nonpos = 0;
    int iwalk = 0;
    int lx = 2;
    double finc = 50.0;

    for (int i=0; i<=itrunk; i++)
    {
        rasite();
        delhf();
        delef();

        if (ipbl == 0) dtcdtg();
        if (ipbl == 1) sibslv();

        /* *******************************************************************
        calculation of evaporative potentials (ecpot, egpot) and interception
        losses; fluxes in j m**-2.  ecidif and egidif hold the excess energy
        if all intercepted water is evaporated during the timestep.  these
        energies are later added to the sensible heat fluxes.

          eci         (e-wc)  : equation (59) , SE-86
          egi         (e-wg)  : equation (60) , SE-86
        ******************************************************************** */

        // check if interception loss term has exceeded canopy storage
        ecpot = ( (etc-ea) + (  getc-deadtc)*dtc-deadtg*dtg-deadqm*dqm );
        eci	= ecpot * wc /(2.*rb) * rcp/psy * dtt;
        ecidif = amax1(0.0,(eci-(snoww[0]+capac[0])*1.e3*hlat));
        eci	= amin1(eci,((snoww[0]+capac[0])*1.e3*hlat));

        egpot = ( (etgs-ea) + (getgs-deadtg)*dtg-deadtc*dtc-deadqm*dqm );
        egi	= egpot/rd*rcp/psy*dtt * ( wg*(1.-areas) + areas );
        egidif = amax1(0.0, egi-(snoww[1]+capac[1])*1.e3*hlat )*(1.-rsnow);
        double egit = amin1(egi, (snoww[1]+capac[1])*1.e3*hlat )*(1.-rsnow);

        // calculation of interception loss from ground-snow. if snow patch
        // shrinks, energy is taken from egi to warm exposed soil to tgs.

        double t1 = snoww[1] - 1./asnow;
        double t2 = amax1( 0.0, t1 );
        double aven = egi - t2*hlat*1.e3/snofac;
        if ( (t1-t2)*egi > 0. ) aven = egi;
        double darea = aven/( (tsnow-tg)*csoil - 1./asnow*hlat*1.e3/snofac);
        double arean = areas + darea;
        egidif = egidif - amin1( 0., arean )*asnow*hlat*1.e3/snofac*rsnow;
        darea = amax1( darea, -areas );
        darea = amin1( 1.-areas, darea );
        heaten = (tsnow-tg)*csoil*darea*rsnow;
        egit = egit + ( egi - heaten - egidif )*rsnow;
        egi = egit;

        double d1 = 1./ra + 1./rb + 1./rd;
        taen = ( (tgs+dtg)/rd + (tc+dtc)/rb + tm/ra ) / d1;
        double hend = ( taen - tm ) * rcp / ra + (ecidif + egidif)/dtt;
        double y = ht - hend;

        newton(ht, y, finc, nox, nonpos, iwalk, lx);
        if(nox != 0) break;
    }
    // exit from non-neutral calculation
    // evapotranspiration fluxes calculated first ( j m-2 )

    /* *******************************************************************
     calculation of transpiration and soil evaporation fluxes for the
     end of the timestep. see figure (2) of se-86.

      ect         (e-c)   : equation (64) , SE-86
      egs         (e-s)   : equation (66) , SE-86
    ******************************************************************* */

    double hrr = hr;
    if ( fg < 0.5 ) hrr = 1.0;

    double coct = (1.0-wc)/ ( rst*fc + 2.*rb );
    double cogs1 = (1.-areas)/(rd+rsoil*fg)*(1.-wg)*hrr;
    double cogs2 = cogs1 / hrr;

    ect = ecpot * coct * rcp/psy * dtt;
    egs = (etgs + getgs*dtg ) * cogs1
          - ( ea + deadtg*dtg + deadtc*dtc + deadqm*dqm ) * cogs2;
    egs = egs * rcp/psy * dtt;
    double egsmax = www[0] / 2. * zdepth[0] * poros * hlat * 1000.0;
    egidif = egidif + amax1( 0., egs - egsmax );
    egs = amin1 ( egs, egsmax );

    // sensible heat flux calculated with latent heat flux correction
    /* *******************************************************************
    calculation of sensible heat fluxes for the end of the timestep.
    see figure (2) of se-86.  note that interception loss excess
    energies (ecidif, egidif) are added.

      hc          (hc)    : equation (63) , SE-86
      hg          (hgs)   : equation (65) , SE-86
    ******************************************************************* */

    hc = hc + (hcdtc*dtc + hcdtg*dtg + hcdth*dth)*dtt + ecidif;
    hg = hg + (hgdtc*dtc + hgdtg*dtg + hgdth*dth)*dtt + egidif;

    /* *******************************************************************
    test of dew condition. latent heat fluxes set to zero if sign of flux
    changes over time step.excess of energy donated to sensible heat flux.
    calculation of total latent heat fluxes,  see figure (2), se-86.

      ec          (ec)    : equation (63) , SE-86
      eg          (eg)    : equation (65) , SE-86
    ******************************************************************* */

    double ecf = sign( 1.0, ecpot );
    double egf = sign( 1.0, egpot );
    double dewc = fc * 2. - 1.0;
    double dewg = fg * 2. - 1.0;

    if(dewc*ecf <= 0.0)
    {
        hc = hc + eci + ect;
        eci = 0.0;
        ect = 0.0;
    }
    if(dewg*egf <= 0.0)
    {
        hg = hg + egs + egi;
        egs = 0.0;
        egi = 0.0;
    }


    ec = eci + ect;
    eg = egi + egs;

    /* *******************************************************************
    adjustment of : temperatures and vapor pressure
    net radiation terms
    storage heat fluxes
    longwave loss and effective surface temperature
    ******************************************************************* */

    ta = taen;
    ea = ea + deadtc*dtc + deadtg*dtg;

    radt[0] = radt[0] + rncdtc*dtc + rncdtg*dtg;
    radt[1] = radt[1] + rngdtc*dtc + rngdtg*dtg;

    // calculation of storage heat fluxes
    chf = ccx / dtt * dtc;
    // 3/96 changes
    // old:  shf = cg / dtt * dtg + timcon*cg*2. * ( tgs+dtg - td )
    shf = cg / dtt * dtg + timcon*csoil*2. * ( tgs+dtg - td );

    zlwup = zlwup - rncdtc * dtc / 2.0
            - rngdtg * dtg * (1.-vcover*(1.-thermk) );

    return;
}


/* ****************************************************************************
   calculation of ustar, u2, ra and drag using Paulson's method.

     (1) site parameters derived from momopt program suite.

     (2) routine is not suitable for gcm applications; designed for
         use with forcing variables measured at a field site.

     (3) paulson psi-coefficients are constrained under unsatble
         conditions to prevent unrealistically low ra values.

     (4) wind speed (um) must be greater than or equal to 0.1 m/s

	***********************************************************************

     variables that defined in the class declaration

      tm     : air temperature at zmet
      um     : wind speed at zwind, um .ge. 0.1
      ht     : sensible heat flux from surface

     parameters that must enter by other functions

      z2     : height of canopy top
      z0     : roughness length
      d      : zero plane displacement
      vkc    : von karmans constant = 0.41
      rhoair : air density
      cpair  : air specific heat

     other parameters

      g1, g2, g3, ztz0, corb1, corb2, ha, zwind, zmet

      g1     : ratio of km(actual) to km(log-linear) at z = z2
      g2     : ratio of ra(actual) to ra(log-linear) for momentum
               between: z = z2 and z = zx, where zx = min(zl,zwind)
      g3     : ratio of ra(actual) to ra(log-linear) for heat
               between: z = z2 and z = zx, where zx = min(zl,zmet)
      ztz0   : parameter to determine depth of transition layer above
               canopy, zl. zl = z2 + ztz0 * z0
      corb1  : non-neutral correction for calculation of aerodynamic
               resistance between ha and z2. when multiplied by
               h*rbb/tm gives bulk estimate of local richardson number.
               rbb = ra for heat between ha and z2.
               corb2 = 9*g/( rhoair*cpair* (du/dz)**2 )
      corb2  : neutral value of rbb*u2 ( squared ), equivalent to
               rdc**2 for upper canopy
      ha     : canopy source height for heat
      zwind  : reference height for wind measurement
      zmet   : reference height for temperature, humidity measurement

        the above are generated from sibx + momopt output

	***********************************************************************

     variables returned from this function

      ustar  : friction velocity
      u2     : wind speed at canopy top
      ra     : aerodynamic resistance for heat flux between ha and zmet
      drag   : shear stress at canopy top

	***********************************************************************

     references
	  1. Paulson C.A. (1970) ' Mathematical representation of wind and
	  temperature profiles in the unstable atmospheric surface layer',
	  J. Appl. Met., 9, 129-861.
	  2. SE-89

**************************************************************************** */
void SiB2Model::rasite()
{
    double zx1, zx2, top;
    static double us1;
    double raf, rafmax;
    double ps1, ps2;
    int nox, nonpos, iwalk, lx;
    double finc, y;

    double hress = ht;
    double zl    = z2 + ztz0 * z0;
    double uest  = vkc*um / log((zwind-d)/z0);

    // calculation of u2 assuming neutral conditions
    if ( zwind <= zl )
    {
        top = 0.0;
        zx1 = zwind - d;
        zx2 = z2 - d;
    }
    else
    {
        zx1 = zwind - d;
        zx2 = zl - d;
        top = log( zx1 / zx2 );
        zx1 = zl - d;
        zx2 = z2 - d;
    }
    double bot = log( zx1 / zx2 );
    double ram = 1. / ( vkc * uest ) * ( top + g2 * bot );
    u2 = um - ram * uest * uest;

    // calculation of ra for heat follows : non-neutrality assumed
    zx1 = zwind - d;
    zx2 = 0.0;
    double arg1 = log ( zx1 / z0 );

    // initialize newton-raphson iterative routine
    nox = 0;
    nonpos = 1;
    iwalk = 0;
    lx = 1;
    finc = 0.2;

    // unstable case : calculation of ustar followed by ra
    if(ht > 0.0)
    {
        do
        {
            unstab ( uest, zx1, zx2, arg1, ht, ps1, ps2);

            y = um - uest/vkc * ( arg1 - ps1 );
            newton ( uest, y, finc, nox, nonpos, iwalk, lx );
        }
        while (nox == 0);

        if(nox == 2)
            cout << "convergence failure in rasite - unstable case" << endl;

        raf = rafcal ( zl, uest, ht );
    }

    // stable case : calculation of ustar
    else
    {
        /* ********************************************************************
        interpolation zone is defined: this spans negative sensible heat fluxes
        from positive side ( hence factor 0.95 ) of the following conditions:
        y = 0,    dy/du* = 0.0;
        see notes for details
        ******************************************************************** */
        double gfac = log((zwind - d)/z0);
        double hm1  = -0.95*tm*rhoair*cpair/(2.0*4.7*g*(zwind-d))*
                      pow(2.0*um/3.0,3)*(vkc/gfac)*(vkc/gfac);
        double hm2  = 5.0*hm1;
        double us2  = vkc*um/(gfac+4.7);
        if( ht >= hm2 )
        {
            ht = amax1 ( hm1 , ht );

            // ustar calculated for slightly stable conditions : ht >= hm1
            do
            {
                stab ( uest, zx1, zx2, arg1, ht, ps1, ps2);

                y = um - uest/vkc * ( arg1 - ps1 );
                newton ( uest, y, finc, nox, nonpos, iwalk, lx );
            }
            while ( nox == 0 );

            if( nox ==  2 )
                cout << "convergence failure in rasite - stable case" << endl;
            ht = hress;

            // ustar calculation in interpolation zone
            if ( ht <= hm1 )
            {
                us1 = uest;
                uest = ( ht-hm2 ) / ( hm1-hm2 ) * ( us1-us2 ) + us2;
            }
        }

        // ustar calculation for collapsed profiles
        else
            uest = us2;

        // calculation of ra for heat transfer between z2 and zmet
        raf = 1.0e5;
        rafmax = rafcal ( zl, us2, hm2 );

        if ( ht >= hm2 )
        {
            double hss = amax1 ( hm1, ht );
            double uss = amax1 ( us1, uest );

            raf = rafcal ( zl, uss, hss );

            if ( ht <= hm1 )
            {
                double raf1 = raf;
                raf  = ( ht-hm2 ) / ( hm1-hm2 ) * ( raf1 - rafmax ) + rafmax;
            }
        }
        raf = amin1 ( raf , rafmax );
    }

    // above canopy variables calculated
    double hrb = ( ht + sqrt(ht * ht) ) / 2.0 + 0.1;

    // corb1 and corb2 are calculated for between ha and z2 only.
    double rbbest = sqrt(corb2)/u2;

    // initialize newton-raphson iterative routine
    nox = 0;
    nonpos = 1;
    iwalk = 0;
    lx = 1;
    finc = 0.2;

    do
    {
        double coef3 = corb1 * hrb / tm / ( z2-ha );
        y = coef3 * pow(rbbest,3) + (u2*rbbest) * (u2*rbbest) - corb2;
        newton( rbbest, y, finc, nox, nonpos, iwalk, lx);
    }
    while ( nox != 1 );

    ra  = raf + rbbest;

    ustar = uest;
    drag = rhoair * uest*uest;

    return;
}


/* ****************************************************************************
   calculation of Paulson psi-function for unstable condition
**************************************************************************** */
void SiB2Model::unstab (const double uest, const double a, const double b,
                        const double argz, const double heat,
                        double& psione, double& psitwo)
{
    double x[2];

    double zin = a;

    for (int i=0; i<2; i++)
    {
        double zml = pow(-uest,3) * rhoair * cpair * tm;
        zml = zml / ( vkc*g*heat );
        double fac = 16.0 * zin/zml;
        x[i] = pow((1.0 - fac), 0.25);
        zin = b;
    }

    psione = 2.*log((1.+x[0])/(1.+x[1]))+log((1.+x[0]*x[0])/
             (1.+x[1]*x[1]))-2.*atan(x[0])+2.*atan(x[1]);
    psione = amin1 ( argz * 0.75, psione );

    psitwo = 2.0*log((1.+x[0]*x[0])/(1.+x[1]*x[1]));
    psitwo = amin1 ( argz * 0.75, psitwo );

    return;
}


/* ****************************************************************************
   calculation of Paulson psi-function for stable condition
**************************************************************************** */
void SiB2Model::stab (const double uest, const double a, const double b,
                      const double argz, const double heat,
                      double& psione, double& psitwo)
{
    psione = 0.0;
    psitwo = 0.0;

    if ( fabs(heat) > 1.e-4 )
    {
        double zml = pow(-uest,3) * rhoair * cpair * tm;
        zml = zml / ( vkc*g*heat );

        psione = -4.7 * ( a-b ) / zml;
        psione = amax1( -4.7, psione );

        psitwo = psione;
    }

    return;
}


/* ****************************************************************************
   calculation of ra for heat between z2 and zmet
**************************************************************************** */
double SiB2Model::rafcal (const double zl, const double uest, const double heat)
{
    double top, zx1, zx2, arg;
    double ps1, ps2;

    if ( zmet <= zl )
    {
        top = 0.0;
        zx1 = zmet - d;
        zx2 = z2 - d;
    }
    else
    {
        zx1 = zmet - d;
        zx2 = zl - d;
        arg = log( zx1 / zx2 );
        if (heat > 0.0) unstab ( uest, zx1, zx2, arg, heat, ps1, ps2);
        if (heat <= 0.0) stab ( uest, zx1, zx2, arg, heat, ps1, ps2);
        top = arg - ps2;

        zx1 = zl - d;
        zx2 = z2 - d;
    }

    arg = log ( zx1 / zx2 );
    if (heat > 0.0) unstab ( uest, zx1, zx2, arg, heat, ps1, ps2);
    if (heat <= 0.0) stab ( uest, zx1, zx2, arg, heat, ps1, ps2);
    double bot = arg - ps2;

    return ( 1.0 / ( vkc*uest ) * ( top + g3 * bot ) );
}


/* ****************************************************************************
   the newton raphson iterative routine will be used to generate new values of
   a1 if dabsolute value of y is greater than ertol; a1 is estimate, y is
   resultant error nex is exit condition  (0=no exit) or (1 when abs(y)<ertol),
   ertol is the dabsolute value of y necessary to obtain an exit, finc is
   initial increment size for second estimate of a1, nonpos=0 if quantity to be
   minimized can be less than zero; nonpos=1 if quantity can only be positive,
   l identifies which quantity is being calculated.

   control values: finc,ertol,nox,nonpos,l:must be set by user
**************************************************************************** */
void SiB2Model::newton(double& a1, double& y, const double finc, int& nox,
                       const int nonpos, int& iwolk, const int l)
{
    static int iter[3], iwalk[3], nex[3];
    static double zinc[3], a2[3], y1[3];

    const double cons	= 1.0;

    double ertol = 0.05 * finc;
    iwalk[l] = iwolk;
    nex[l]=nox;

    if (ertol <= 0.000001) ertol = 0.000001;

    if (fabs(y) <= ertol)
    {
        nex[l]	= 1;
        if (iter[l] >= 490) iter[l] = 2;
        zinc[l]	= 0.0;
        iter[l]	= 0;
        iwalk[l]=0;
        y1[l]	= 0.0;
        y		= 0.0;
        a2[l]	= 0.0;
    }

    else if ((fabs(y-y1[l])) <= 0.01*ertol && iwalk[l] == 0 )
    {
        a1 = a1 + finc*2.0;
        nex[l]=0;
    }

    else if (fabs(y1[l]) > ertol)
    {
        iter[l]=iter[l]+1;
        if(iter[l] == 20) iwalk[l]=1;
        if(iwalk[l] != 0)
        {
L2:
            if (iwalk[l] == 2)
            {
                if(sign(cons,y) == sign(cons,y1[l]))
                {
                    a2[l]=a1;
                    a1=a1+zinc[l];
                    y1[l]=y;
                    nex[l]=0;
                }
                else
                {
                    zinc[l]=-zinc[l]/4.0;
                    a2[l]=a1;
                    a1=a1+zinc[l];
                    nex[l]=0;
                    y1[l]=y;
                }
            }

            else if (iwalk[l] == 3)
            {
                if(sign(cons,y) == sign(cons,y1[l]))
                {
                    a2[l] = a1;
                    a1 = a1+finc;
                    y1[l]=y;
                    nex[l] = 0;
                }
                else
                {
                    iwalk[l] = 1;
                    goto L2;
                }
            }

            else if (sign(cons,y) == sign(cons,y1[l]))
            {
                double a=a1-y*(a1-a2[l])/(y-y1[l]);
                if(fabs(a-a1) > (10.0*finc)) a = a1+10.0*finc*sign(cons,(a-a1));
                a2[l]=a1;
                a1=a;
                y1[l]=y;
            }

            else
            {
                zinc[l]		= (a1-a2[l])/4.0;
                a1			= a2[l]+zinc[l];
                iwalk[l]	= 2;
                nex[l]		= 0;
            }
        }

        else if (fabs(y) > ertol)
        {
            double a=a1-y*(a1-a2[l])/(y-y1[l]);
            if(fabs(a-a1) > (10.0*finc)) a=a1+10.0*finc*sign(cons,(a-a1));
            a2[l]=a1;
            a1=a;
            y1[l]=y;
        }

        else
        {
            nex[l]	= 1;
            if (iter[l] >= 490) iter[l] = 2;
            zinc[l]	= 0.0;
            iter[l]	= 0;
            iwalk[l]=0;
            y1[l]	= 0.0;
            y		= 0.0;
            a2[l]	= 0.0;
        }

    }

    else
    {
        a2[l] = a1;
        double step = amin1( fabs(y), fabs(10.0 * finc) ) * sign(cons,y);
        a1 = a1-step;
        nex[l]=0;
        y1[l]=y;
        iter[l]=1;
        if (iwalk[l] != 3) iwalk[l]=0;
    }

    if(nonpos == 1 && a1 < 0.0) a1=a2[l]/2.0;
    nox = nex[l];
    iwolk = iwalk[l];

    return;
}


/* ****************************************************************************
   calculation of rb and rd as functions of u2 and temperatures
	rb (grb)       canopy to cas aerodynamic resistance (s m-1)
	rd (grd)       ground to cas aerodynamic resistance (s m-1)
	ta (gta)       cas temperature (k)
		cas : canopy air space
**************************************************************************** */
void SiB2Model::rbrd()
{
    // rb : equation (a9), SE-86
    double temdif = amax1( 0.1,  tc-tm );
    double fac = zlt/890.* sqrt( sqrt(temdif/0.05));
    rb  = 1.0/(sqrt(u2)/rbc+fac);

    // rd : equation (a15), se-86
    temdif = amax1( 0.1, tgs-tm );
    double fih = sqrt( 1. + 9. * g * temdif * z2 / tgs / ( u2*u2) );
    rd  = rdc / u2 / fih;

    // calculation of ta, ht and effective soil aero-surface conductance,
    // see equation (66) of SE-86 with snow area term added
    double d1 = 1./ra + 1./rb + 1./rd;
    ta = ( tgs/rd + tc/rb + tm/ra ) / d1;
    ht = ( ta - tm ) * rcp / ra;

    double cogr = (1.-wg)/(rsoil*fg+rd);
    double cogs =  wg/rd;
    cog1 = (cogs + cogr*hr) * (1.-areas) + areas/rd;
    cog2 = (cogs + cogr   ) * (1.-areas) + areas/rd;

    return;
}


/* ****************************************************************************
   calculation of canopy photosynthetic rate using the integrated model
   relating assimilation and stomatal conductance.
   method uses numerical solution based on extrapolation from error versus
   co2i line.
   units are converted from mks to biological units in this routine.
   base reference :  SE-92A

                        units
                       -------

    pco2m, pco2a, pco2i, po2m                : pascals
    co2a, co2s, co2i, h2oa, h2os, h2oa       : mol mol-1
    vmax0, respn, assim, gs, gb, ga, pfd     : mol m-2 s-1
    effcon                                   : mol co2 mol quanta-1
    gcan, 1/rb, 1/ra, 1/rst                  : m s-1
    evapkg                                   : kg m-2 s-1
    q                                        : kg kg-1

                     conversions
                    -------------

    1 mol h2o           = 0.018 kg
    1 mol co2           = 0.044 kg
    h2o (mol mol-1)     = ea / psur ( mb mb-1 )
    h2o (mol mol-1)     = q*mm/(q*mm + 1)
    gs  (co2)           = gs (h2o) * 1./1.6
    gs  (mol m-2 s-1 )  = gs (m s-1) * 44.6*tf/t*p/po
    par (mol m-2 s-1 )  = par(w m-2) * 4.6*1.e-6
    mm  (molair/molh2o) = 1.611


                       output
                    -------------

    assimn              = canopy net assimilation rate
    ea                  = canopy air space vapor pressure
    1/rst               = canopy conductance
    pco2i               = internal co2 concentration
    respc               = canopy respiration
    respg               = ground respiration
    rstfac[4]      canopy resistance stress factors

	**********************************************************************

	   rstfac[0] ( f(h-s) )               : equation (17,18), SE-92A
       rstfac[1] ( f(soil) )              : equation (12 mod), SE-89
       rstfac[2] ( f(temp) )              : equation (5b)   , CO-92
       rstfac[3] ( f(h-s)*f(soil)*f(temp))

**************************************************************************** */
void SiB2Model::phosib()
{
    double pco2y[6], eyy[6];
    const double respg = 0.0;

    double c3     = 0.0;
    if( effcon > 0.07 ) c3 = 1.0;
    double c4     = 1.0 - c3;

    // calculation of canopy par use parameter.
    // fparkk      (pi)     : equation (31) , SE-92A
    double scatp    =     green   * ( tran[0][0] + ref[0][0] )
                          + ( 1.-green ) * ( tran[0][1] + ref[0][1] );
    double scatg    = tran[0][0] + ref[0][0];
    double park = sqrt(1.-scatp) * gmudmu;
    fparc = 1. - exp ( -park*zlt );
    double fparkk   = fparc / park * green;

    // q-10 temperature effects :
    // qt          (qt)    : table (2)     , SE-92A
    // qt for vm changed to 2.1
    double qt = 0.1*( tc - trop );
    double respn = respcp * vmax0 * rstfac[1];
    respc = respn * pow(2.0,qt) / ( 1. + exp( trda*(tc-trdm )) );
    double vm = vmax0 * pow(2.1,qt);

    double templ = 1. + exp(slti*(hlti-tc));
    double temph = 1. + exp(shti*(tc-hhti));
    rstfac[2] = 1./( templ*temph);
    vm    = vm/temph * rstfac[1] * c3  + vm * rstfac[1]*rstfac[2] * c4;

    /* *********************************************************************
    Michaelis-Menten constants for co2 and o2, co2/o2 specificity,
    compensation point

      zkc          (kc)     : table (2)     , SE-92A
      zko          (ko)     : table (2)     , SE-92A
      spfy         (s)      : table (2)     , SE-92A
      gammas       (gamma-*): table (2)     , SE-92A
      omss         (omega-s): equation (13) , SE-92A
      bintc        (b*zlt)  : equation (35) , SE-92A
    ********************************************************************* */

    double zkc = 30. * pow(2.1,qt);
    double zko = 30000. * pow(1.2,qt);
    double spfy = 2600. * pow(0.57,qt);
    double gammas = 0.5 * po2m/spfy * c3;
    double pfd    = 4.6*1.e-6 * gmudmu * ( radn[0][0] + radn[0][1] );

    double h2oi   = etc/psur;
    double h2oa   =  ea/psur;
    double h2os	  = 0.0;
    double h2om   =  em/psur;
    double h2osl  = etgs/psur;

    double tprcor = tf*psur*100./1.013e5;

    double gbh2o  = 0.5/rb * 44.6*tprcor/tc;
    double gah2o  = 1.0/ra * 44.6*tprcor/tm;
    double gog1   = cog1   * 44.6*tprcor/tgs;
    double gog2   = cog2   * 44.6*tprcor/tgs;

    double rrkk   = zkc*( 1. + po2m/zko ) * c3 + vmax0/5.* pow(1.8,qt) * c4;
    double par    = pfd*effcon*( 1.-scatg );
    double bintc  = binter*zlt*green*amax1(0.1,rstfac[1]);

    double omss  = ( vmax0/2.0 ) * pow(1.8,qt)/templ * rstfac[1] * c3
                   + rrkk * rstfac[1] * c4;

    // first guess is midway between compensation point and maximum
    // assimilation rate.
    double range    = pco2m * ( 1.0 - 1.6/gradm ) - gammas;

    for (int ic=0; ic<6; ic++)
    {
        pco2y[ic]	= 0.0;
        eyy[ic]		= 0.0;
    }

    int ic2;
    for (int ic=0; ic<6; ic++)
    {
        ic2 = ic;

        sortin( eyy, pco2y, range, gammas, ic );

        cycalc( fparkk, vm, gradm, bintc, atheta, btheta,
                gah2o, gbh2o, gog1, gog2, wc,
                h2oi, h2om, h2osl, par, pco2m, psur,
                gammas, respc, respg, rrkk, omss, c3, c4,
                pco2y[ic], eyy[ic], gsh2o, assimn, h2os, h2oa );

        if( fabs(eyy[ic]) < 0.1 ) break;
    }

    pco2i = pco2y[ic2];

    rstfac[0] = h2os/h2oi;
    rstfac[3] = rstfac[0]*rstfac[1]* rstfac[2];
    rst   = amin1( 1.e6, 1./( gsh2o*tc/( 44.6*tprcor) ) );
    ea    = h2oa*psur;

    return;
}


/* ****************************************************************************
   arranges successive pco2/error pairs in order of increasing pco2.
   estimates next guess for pco2 using combination of linear and quadratic fits
**************************************************************************** */
void SiB2Model::sortin(double eyy[6], double pco2y[6], const double range,
                       const double gammas, const int ic)
{
    if ( ic < 3 )
    {
        pco2y[0] = gammas + 0.5*range;
        pco2y[1] = gammas + range*( 0.5 - 0.3*sign(1.0,eyy[0]) );
        pco2y[2] = pco2y[0] - (pco2y[0]-pco2y[1])/(eyy[0]-eyy[1]+1.e-10)*eyy[0];

        double pmin = amin1( pco2y[0], pco2y[1] );
        double emin = amin1(   eyy[0],   eyy[1] );
        if ( emin > 0.0 && pco2y[2] > pmin ) pco2y[2] = gammas;
    }

    else
    {
        int n = ic - 1;
        for (int j=1; j<=n; j++)
        {
            double a = eyy[j];
            double b = pco2y[j];
            int i;
            for (i = j-1; i>=0; i--)
            {
                if(eyy[i] <= a ) goto L200;
                eyy[i+1] = eyy[i];
                pco2y[i+1] = pco2y[i];
            }
            i = -1;
L200:
            eyy[i+1] = a;
            pco2y[i+1] = b;
        }

        double pco2b = 0.0;
        int is    = 0;
        for (int ix = 0; ix<=n; ix++)
        {
            if( eyy[ix] < 0.0 ) pco2b = pco2y[ix];
            if( eyy[ix] < 0.0 ) is = ix;
        }

        int i1 = is-1;
        i1 = max0(0, i1);
        i1 = min0(n-2, i1);
        int i2 = i1 + 1;
        int i3 = i1 + 2;
        int isp   = is + 1;
        isp = min0( isp, n );
        is = isp - 1;

        double pco2yl = pco2y[is] - (pco2y[is]-pco2y[isp])/(eyy[is]-eyy[isp])*eyy[is];

        // method using a quadratic fit
        double ac1 = eyy[i1]*eyy[i1] - eyy[i2]*eyy[i2];
        double ac2 = eyy[i2]*eyy[i2] - eyy[i3]*eyy[i3];
        double bc1 = eyy[i1] - eyy[i2];
        double bc2 = eyy[i2] - eyy[i3];
        double cc1 = pco2y[i1] - pco2y[i2];
        double cc2 = pco2y[i2] - pco2y[i3];
        double bterm = (cc1*ac2-cc2*ac1)/(bc1*ac2-ac1*bc2);
        double aterm = (cc1-bc1*bterm)/ac1;
        double cterm = pco2y[i2] - aterm*eyy[i2]*eyy[i2] - bterm*eyy[i2];
        double pco2yq = cterm;
        pco2yq = amax1( pco2yq, pco2b );
        pco2y[ic] = ( pco2yl+pco2yq)/2.0;
    }

    pco2y[ic] = amax1 ( pco2y[ic], 0.01 );

    return;
}


/* ****************************************************************************
   calculation equivalent to steps in figure 4 of SE-92A
   c4 calculation based on CO-92.

	******************************   output   ******************************

       pco2i          canopy internal co2 concentration (mol mol-1)
       gsh2o          canopy conductance (mol m-2 s-1)
       h2os           canopy surface h2o concentration (mol mol-1)

  	****************************  diagnostics  *****************************

       omc            rubisco limited assimilation (mol m-2 s-1)
                        (omega-c): equation (11) , SE-92A
       ome            light limited assimilation (mol m-2 s-1)
                        (omega-e): equation (12) , SE-92A
       oms            sink limited assimilation (mol m-2 s-1)
       co2s           canopy surface co2 concentration (mol mol-1)
                        equation (18c) , SE-92
       assimn         (a-n)    :  equation (14,15), SE-92A

**************************************************************************** */
void SiB2Model::cycalc(const double fparkk, const double vm, const double gradm,
                       const double bintc, const double atheta, const double btheta,
                       const double gah2o, const double gbh2o, const double gog1,
                       const double gog2, const double wc, const double h2oi,
                       const double h2om, const double h2osl, const double par,
                       const double pco2m, const double psur, const double gammas,
                       const double respc, const double respg, const double rrkk,
                       const double omss, const double c3, const double c4,
                       const double pco2i, double& eyy, double& gsh2o,
                       double& assimn, double& h2os, double& h2oa)
{
    double omc = vm  * ( pco2i-gammas )/( pco2i + rrkk ) * c3 + vm * c4;
    double ome = par * ( pco2i-gammas )/( pco2i+2.*gammas ) * c3 + par * c4;
    double sqrtin= amax1( 0.0, ((ome+omc)*(ome+omc) - 4.0*atheta*ome*omc ) );
    double omp   = ( ( ome+omc ) - sqrt( sqrtin ) ) / ( 2.0*atheta );
    double oms   = omss * c3 + omss*pco2i * c4;
    sqrtin= amax1( 0., ( (omp+oms)*(omp+oms) - 4.*btheta*omp*oms ) );
    double assim = ( ( oms+omp ) - sqrt( sqrtin ) ) / ( 2.*btheta );
    assimn= ( assim - respc) * fparkk;

    // gah2o bottom stopped to prevent negative values of pco2a
    double pco2a = pco2m - (1.4/amax1(0.446,gah2o)
                            * (assimn - respg)* psur*100.);
    double co2s  = pco2a/(psur*100.) - 1.4*assimn/gbh2o;

    double assmt = amax1( 1.e-12, assimn );
    double co2st = amin1( co2s, pco2a/(psur*100.) );
    co2st = amax1( co2st,1./(psur*100.) );

    double div2  = gah2o + gbh2o + gog2;
    double hcdma = h2oi*co2st / ( gradm*assmt );
    double alpha = hcdma*gbh2o*(1.-wc) / div2;
    double beta  = ( -hcdma*bintc*gbh2o*(1.-wc) + h2oi*gbh2o*wc
                     + h2osl*gog1 + h2om*gah2o ) / div2;
    double aquad = hcdma;
    double bquad = gbh2o*( hcdma-alpha ) - h2oi - bintc*hcdma;
    double cquad = -gbh2o*( hcdma*bintc + beta );

    sqrtin= amax1( 0., ( bquad * bquad - 4.*aquad*cquad ) );
    gsh2o = ( -bquad + sqrt ( sqrtin ) ) / (2.*aquad);
    h2os  = ( gsh2o-bintc ) * hcdma;
    h2os  = amin1( h2os, h2oi );
    h2os  = amax1( h2os, 1.e-7);
    gsh2o = h2os/hcdma + bintc;
    double div3  = gbh2o*gsh2o/(gbh2o+gsh2o)*(1.-wc) + gbh2o*wc + gog2 + gah2o;
    h2oa  = ( ( h2oi - h2oi*gsh2o/(gbh2o+gsh2o) )*gsh2o*(1.-wc)
              + h2oi*gbh2o*wc + h2osl*gog1 + h2om*gah2o ) / div3;
    h2os  = ( h2oa*gbh2o + h2oi*gsh2o )/( gbh2o + gsh2o );

    // mplied value of co2i derived from assimilation rate and stomatal
    // conductance.
    double pco2in = ( co2s - 1.6 * assimn / gsh2o )*psur*100.0;
    eyy = pco2i - pco2in;

    return;
}


/* ****************************************************************************
   partial derivatives of radiative and sensible heat fluxes
**************************************************************************** */
void SiB2Model::delrn()
{
    double tc3 = tc * tc * tc;
    double tg3 = tgs * tgs * tgs;
    double fac1 = ( 1.-thermk ) * vcover;
    double fac2 = 1.0;

    rncdtc = - 2. * 4. * fac1 * stefan * tc3;
    rncdtg = 4. * fac1 * fac2 * stefan * tg3;

    rngdtg = - 4. * fac2 * stefan * tg3;
    rngdtc = 4. * fac1 * fac2 * stefan * tc3;

    return;
}


/* ****************************************************************************
   calculation of partial derivatives of canopy and ground sensible heat
   fluxes with respect to tc, tgs, and theta-m.
   calculation of initial sensible heat fluxes.

      hc             canopy sensible heat flux (j m-2)
      hg             ground sensible heat flux (j m-2)
      hcdtc          dhc/dtc
      hcdtg          dhc/dtgs
      hcdth          dhc/dth
      hgdtc          dhg/dtc
      hgdtg          dhg/dtgs
      hgdth          dhg/dth
      aac            dh/dtc
      aag            dh/dtgs
      aam            dh/dth
**************************************************************************** */
void SiB2Model::delhf()
{

    // fluxes expressed in joules m-2
    // hc (h-c) : equation(63), SE-86; hg  (h-g) : equation(65), SE-86
    double d1 = 1./ra + 1./rb + 1./rd;
    ta = ( tgs/rd + tc/rb + tm/ra ) / d1;

    hc = rcp * ( tc - ta ) / rb * dtt;
    hg = rcp * ( tgs - ta ) / rd * dtt;

    /* *******************************************************************
      n.b.      fluxes expressed in joules m-2

      hcdtc     (dhc/dtc) : equation (14) , SA-89B
      hcdtg     (dhc/dtgs): equation (14) , SA-89B
      hcdth     (dhc/dth) : equation (14) , SA-89B
      hgdtc     (dhg/dtc) : equation (15) , SA-89B
      hgdtg     (dhg/dtgs): equation (15) , SA-89B
      hgdth     (dhg/dth) : equation (15) , SA-89B
      aac       (dh/dtc)  : equation (12) , SA-89B
      aag       (dh/dtgs) : equation (12) , SA-89B
      aam       (dh/dth)  : equation (12) , SA-89B
    ******************************************************************* */

    hcdtc = rcp / rb * ( 1./ra + 1./rd ) / d1;
    hcdtg = - rcp / ( rb * rd ) / d1;

    hgdtg = rcp / rd * ( 1./ra + 1./rb ) / d1;
    hgdtc = - rcp / ( rd * rb ) / d1;

    hcdth = - rcp/( rb*ra )/d1*bps;
    double hcdqm = 0.0;

    hgdth = -rcp/(rd*ra)/d1*bps;
    double hgdqm = 0.0;

    aag = 1.0/(rd*d1);
    aac = 1.0/(rb*d1);
    aam = 1.0/(ra*d1)*bps;

    return;
}


/* ****************************************************************************
   calculation of partial derivatives of canopy and ground latent heat fluxes
   with respect to tc, tgs, theta-m, and qm.
   calculation of initial latent heat fluxes.

      ec  ( e-c ) : equation(64), SE-86
      eg  ( e-gs) : equation(66), SE-86

       ec             ect + eci
       eg             egs + egi
       ecdtc          dec/dtc
       ecdtg          dec/dtgs
       ecdqm          dec/dqm
       egdtc          deg/dtc
       egdtg          deg/dtgs
       egdqm          deg/dqm
       bbc            de/dtc
       bbg            de/dtgs
       bbm            de/dqm
**************************************************************************** */
void SiB2Model::delef()
{
    // modification for soil dryness : hr = rel. humidity in top layer
    double hrr = hr;
    if ( fg < 0.5 ) hrr = 1.0;

    // calculation of surface resistance components, see equations (64,66)
    // of SE-86
    double rcc = rst*fc + 2. * rb;
    double coc = (1.-wc)/rcc + wc/(2.*rb);

    double cogr = (1.-wg)/(rsoil*fg+rd);
    double cogs =  wg/rd;
    cog1 = (cogs + cogr*hrr) * (1.-areas) + areas/rd;
    cog2 = (cogs + cogr    ) * (1.-areas) + areas/rd;

    double d2 = 1./ra + coc + cog2;
    double top = coc * etc + cog1 * etgs + em/ra;
    ea = top / d2;

    ec = ( etc - ea ) * coc * rcp/psy * dtt;
    eg = ( etgs*cog1 - ea*cog2 ) * rcp/psy * dtt;

    deadtc = getc * coc / d2;
    deadtg = getgs * cog1 / d2;

    /* *******************************************************************
    	ecdtc     (dec/dtc) : equation (14) , SA-89B
    	ecdtg     (dec/dtgs): equation (14) , SA-89B
    	ecdqm     (dec/dqm) : equation (14) , SA-89B
    	egdtc     (deg/dtc) : equation (15) , SA-89B
    	egdtg     (deg/dtgs): equation (15) , SA-89B
    	egdqm     (deg/dqm) : equation (15) , SA-89B
    	bbc       (de/dtc)  : equation (13) , SA-89B
    	bbg       (de/dtgs) : equation (13) , SA-89B
    	bbm       (de/dqm)  : equation (13) , SA-89B
    ******************************************************************* */

    ecdtc = ( getc - deadtc ) * coc * rcp / psy;
    ecdtg = - deadtg * coc * rcp / psy;

    egdtg = ( getgs*cog1 - deadtg*cog2 ) * rcp / psy;
    egdtc = - deadtc * cog2 * rcp / psy;

    double sh	    = epsfac * em / (psur -em );
    double deadem	= 1.0/( ra*d2 );
    demdqm	= epsfac * psur/((epsfac+sh)*(epsfac+sh));
    deadqm	= deadem*demdqm;

    ecdqm	= - deadqm*coc*rcp/psy;
    egdqm	= -deadqm*cog2*rcp/psy;

    bbg		= (cog1/d2) * getgs*epsfac * psur/((psur-etgs)*(psur-etgs));
    bbc		= (coc /d2) * getc *epsfac * psur/((psur-etc)*(psur-etc));
    bbm		= 1.0/(ra*d2);

    return;
}


/* ****************************************************************************
   solve for time changes of pbl and sib variables using a semi-implicit scheme

     dtc, dtg, dth, dqm  : equations(12-15) , SA-89B + radiation terms

	*****************************  output  ******************************
       dtc            canopy temperature increment (K)
       dtg            ground surface temperature increment (K)
       dth            mixed layer potential temperature increment (K)
       dqm            mixed layer mixing ratio increment (kg kg-1)
**************************************************************************** */
void SiB2Model::sibslv()
{
    double pblsib[4][4], cvec[4], solvec[4], chin[4][5];

    double grav2 = 0.01*g;
    double gv    = grav2*rhoair/ra;
    double psb   = 100.0;

    double etem = 0.0;
    double dthb = 0.0;
    double dwb	= 0.0;

    // cvec uses provisional values of the fluxes
    double fths = (hc+hg) / (dtt*cpair*bps);
    double fws  = (ec+eg) / (dtt*hlat     );

    // tg equation
    // 3/5/96 changes
    // old   pblsib[0][0] = cg/dtt + hgdtg+egdtg-rngdtg + timcon*cg*2.0
    pblsib[0][0] = cg/dtt + hgdtg+egdtg-rngdtg + timcon*csoil*2.0;
    pblsib[0][1] =        + hgdtc+egdtc-rngdtc;
    pblsib[0][2] =        + hgdth;
    pblsib[0][3] =				egdqm;

    // tc equation
    pblsib[1][0] =         + hcdtg+ecdtg-rncdtg;
    pblsib[1][1] =	ccx/dtt+ hcdtc+ecdtc-rncdtc;
    pblsib[1][2] =         + hcdth;
    pblsib[1][3] =                 ecdqm;

    // theta equation
    pblsib[2][0] = -gv* aag;
    pblsib[2][1] = -gv* aac;
    pblsib[2][2] = -gv*(aam-1.0) + etem + psb/dtt;
    pblsib[2][3] = 0.0;

    // sh equation
    pblsib[3][0] = -gv*bbg;
    pblsib[3][1] = -gv*bbc;
    pblsib[3][2] = 0.0;
    pblsib[3][3] = -gv*(bbm-1.0) + etem + psb/dtt;

    // 3/96 changes
    // old  cvec[0] = radt[1] - hg/dtt - eg/dtt - timcon*(tgs-td)*cg*2.
    cvec[0] = radt[1] - hg/dtt - eg/dtt - timcon*(tgs-td)*csoil*2.0;
    cvec[1] = radt[0] - hc/dtt - ec/dtt;
    cvec[2] = grav2*fths + etem*dthb;
    cvec[3] = grav2*fws  + etem*dwb;

    // solve 4 x 4 matrix equation

    for (int j=0; j<4; j++)
        for (int i=0; i<4; i++)
            chin[i][j] = pblsib[i][j];
    for (int i=0; i<4; i++)
        chin[i][4] = cvec[i];

    gauss(chin, 4, 5, solvec);

    dtg=solvec[0];
    dtc=solvec[1];
    dth=solvec[2];
    dqm=solvec[3];

    return;
}


/* ****************************************************************************
   solve a linear system by gaussian elimination.  developed
   by Dr. chin-hoh moeng.  a is the matrix of coefficients, with the vector of
   constants appended as an extra column.  x is the vector containing the
   results. the input matrix is not destroyed.
**************************************************************************** */
void SiB2Model::gauss(const double a[4][5], const int n, const int np1, double x[4])
{
    double work[4][5];
    double r;

    for (int i=0; i<n; i++)
        for (int j=0; j<np1; j++)
            work[i][j] = a[i][j];

    for (int i=1; i<n; i++)
    {
        for (int j=i; j<n; j++)
        {
            r = work[j][i-1] / work[i-1][i-1];
            for (int k=0; k<np1; k++)
                work[j][k] = work[j][k] - r * work[i-1][k];
        }
    }

    for (int i=1; i<n; i++)
    {
        int k = n-i;
        r = work[k][np1-1] / work[k][k];
        for (int j=i; j<n; j++)
        {
            int l = n-j-1;
            work[l][np1-1] = work[l][np1-1] - r * work[l][k];
        }
    }

    for (int i=0; i<n; i++)
        x[i] = work[i][np1-1] / work[i][i];

    return;
}


/* ****************************************************************************
   calculation of temperature tendencies assuming no interaction with the pbl:
   equations(69,70), SE-86
**************************************************************************** */
void SiB2Model::dtcdtg()
{
    double ccodtc = ccx / dtt - rncdtc + hcdtc + ecdtc;
    double ccodtg = - rncdtg + hcdtg + ecdtg;
    double ccorhs = radt[0] - ( hc + ec ) / dtt;

    // 3/96 changes
    //old  gcodtg = cg / dtt + timcon*cg*2. - rngdtg + hgdtg + egdtg
    double gcodtg = cg / dtt + timcon*csoil*2. - rngdtg + hgdtg + egdtg;
    double gcodtc = - rngdtc + hgdtc + egdtc;

    // 3/96 changes
    // old  gcorhs = radt[1] - timcon*cg*2. * (tgs -td) - (hg + eg) / dtt
    double gcorhs = radt[1] - timcon*csoil*2. * (tgs -td) - (hg + eg) / dtt;

    double denom = ccodtc * gcodtg - ccodtg * gcodtc;

    dtc = ( ccorhs * gcodtg - ccodtg * gcorhs ) / denom;
    dtg = ( ccodtc * gcorhs - ccorhs * gcodtc ) / denom;

    return;
}


/* ****************************************************************************
   updating of all prognostic variables.

     subroutines called   : updat2
     ------------------     snow2
                            run2

	************************** output from this block **********************

		dtc            canopy temperature increment (K)
		dtd            deep soil temperature increment (K)
		dtg            ground surface temperature increment (K)
		www[2]         ground wetness
		capac[1]       canopy/ground liquid interception store (m)
		snoww[1]       canopy/ground snow interception store (m)
		roff           runoff (mm)
		etmass (fws)   evapotranspiration (mm)
		hflux (fss)    sensible heat flux (w m-2)

	******************************* diagnostics ****************************

		ecmass         canopy evapotranspiration (mm)
		egmass         ground evapotranspiration (mm)

**************************************************************************** */
void SiB2Model::updat2()
{
    snow1();

    /* *********************************************************************
    interception losses.
    evaporation losses are expressed in j m-2 : when divided by ( hlat*1000.)
    loss is in m m-2. mass terms are in kg m-2 dt-1 interception and drainage
    treated in inter2.
    ********************************************************************* */

    rsnow = snoww[0]/(snoww[0]+capac[0]+1.e-10);
    double facks = 1. + rsnow * ( snofac-1. );
    if ( (ect+eci) <= 0.0)
    {
        eci = ect+eci;
        ect = 0.0;
        facks = 1.0 / facks;
    }
    capac[0] = capac[0] - ( 1.-rsnow )*eci*facks/hlat/1.e3;
    snoww[0] = snoww[0] -      rsnow  *eci*facks/hlat/1.e3;
    ecmass = eci*facks / hlat;

    rsnow = snoww[1]/(snoww[1]+capac[1]+1.e-10);
    facks = 1. + rsnow * ( snofac-1. );
    if ( (egs+egi) <= 0. )
    {
        egi = egs+egi;
        egs = 0.0;
        facks = 1. / facks;
    }
    capac[1] = capac[1] - ( 1.-rsnow )*egi*facks/hlat/1.e3;
    snoww[1] = snoww[1] -      rsnow  *egi*facks/hlat/1.e3;
    egmass = egi*facks / hlat;

    // dumping of small capac values onto soil surface store
    for (int iveg = 0; iveg < 2; iveg++)
    {
        if ( (snoww[iveg]+capac[iveg]) > 0.00001 ) continue;
        www[0] = www[0] + (snoww[iveg]+capac[iveg]) / ( poros*zdepth[0] );
        capac[iveg] = 0.0;
        snoww[iveg] = 0.0;
    }

    // snowmelt / refreeze calculation
    snow2();

    /* *********************************************************************
    evapotranspiration losses,
    extraction of transpiration loss from root zone, soil evaporation.

    	ect         (e-dc)  : equation (5,6), SE-86
    	egs         (e-s)   : equation (5)  , SE-86
    ********************************************************************* */

    double facl   = 1./hlat/1.e3/(poros*zdepth[1]);
    double extrak = ect*facl;
    extrak = amin1( extrak, www[1] );
    double ectdif = ect - extrak/facl;
    ect    = extrak/facl;
    hc     = hc + ectdif;
    ecmass = ecmass + ect/hlat;
    www[1] = www[1] - ect*facl;

    facl   = 1./hlat/1.e3/(poros*zdepth[0]);
    extrak = egs*facl;
    extrak = amin1( extrak, www[0] );
    double egsdif = egs - extrak/facl;
    egs    = extrak/facl;
    hg     = hg + egsdif;
    egmass = egmass + egs/hlat;
    www[0] = www[0] - egs*facl;

    // calculation of total moisture and sensible heat fluxes from surface.
    etmass = ecmass + egmass;
    hflux  = (hc+hg) / dtt;

    // calculation of interflow, infiltration excess and loss to groundwater.
    // all losses are assigned to variable 'roff'.

    run2();

    /* *********************************************************************
    update of temperatures and pbl variables. note that tc and tg are
    modified in interc as a result of precipitation inputs.
    update of interception stores.
    ********************************************************************* */

    tc  = tc + dtc;
    tg  = tg + dtg;
    td  = td + dtd;
    // double th  = th + dth;
    // double qm  = qm + dqm;

    return;
}


/* ****************************************************************************
   snowmelt / refreeze calculation

   calculation of snowmelt and modification of temperatures
     modification deals with snow patches:
          ts < tf, tsnow = ts
          ts > tf, tsnow = tf
**************************************************************************** */
void SiB2Model::snow2()
{
    for (int iveg = 0; iveg < 2; iveg++)
    {
        double realc = (1 - iveg) * 1.0;
        double realg = iveg * 1.0;

        double cctt = realc*ccx  +  realg*cg;
        double cct  = realc*ccx  +  realg*csoil;
        double ts   = realc*tc   +  realg*tg;
        double dts  = realc*dtc  +  realg*dtg;
        double flux = realc*chf  +  realg*dtg/dtt*cg;

        tsnow = amin1 ( tf-0.01, ts );
        double snowhc = amin1( 0.05, snoww[iveg] ) * cw * realg;
        double zmelt = 0.0;

        // no snow  present
        if ( snoww[iveg] <= 0.0 )
        {
            // simple thermal balance with possible freezing
            if ( ( ts+dts) <= tf )
            {
                double freeze = amin1 ( 0., (flux*dtt - ( tf-0.01 - ts )*cctt ) );
                snoww[iveg] = amin1( capac[iveg], - freeze/snomel );
                zmelt = capac[iveg] - snoww[iveg];
                capac[iveg] = 0.0;
                dts = dts + snoww[iveg]*snomel/cctt;
            }
        }

        // snow present
        else
        {
            if (!( ts < tf && (ts+dts) < tf ))
            {
                // snow present : ts < tf,  ts+dts > tf
                if ( ts <= tf )
                {
                    double avex = flux - ( tf-0.01 - ts ) * cctt/dtt;
                    double avmelt = ( avex/snomel * (areas*realg + realc))*dtt;
                    zmelt = amin1( avmelt, snoww[iveg] );
                    snoww[iveg] = snoww[iveg] - zmelt;
                    double avheat = avex*( 1.-areas )*realg + ( avmelt-zmelt )
                                    * snomel/dtt;
                    double safe = amax1( ( 1.-areas*realg ), 1.e-8 );
                    dts = tf-0.01 - ts + avheat / ( cctt*safe )*dtt;
                }

                // snow present and ts > tf : ground only
                else
                {
                    double tbulk = tsnow*areas + ts*( 1. - areas );
                    double tn = tbulk + dts;
                    double exheat = cct*( 1.001-amax1(0.1,areas)) * dts;
                    double exmelt = flux*dtt - exheat;
                    double heat = exheat;
                    double dtsg = exheat / ( cct*(1.001-areas ));
                    if ( (ts+dtsg) <= tf )
                    {
                        heat = ( tf-0.01 - ts ) * ( cct*(1.-areas) );
                        dtsg = tf-0.01 - ts;
                    }

                    exmelt = exmelt + exheat - heat;

                    if( exmelt >= 0.0 )
                    {
                        zmelt = exmelt/snomel;
                        if( asnow*(snoww[iveg]-zmelt) < 1.0 )
                            zmelt = amax1( 0., snoww[iveg] - 1./asnow );
                        snoww[iveg] = snoww[iveg] - zmelt;
                        exmelt = exmelt - zmelt*snomel;
                        double zmelt2 = exmelt/ (cct*( ts-tf)*asnow + snomel);
                        zmelt2 = amin1( zmelt2, snoww[iveg] );
                        zmelt = zmelt + zmelt2;
                        snoww[iveg] = snoww[iveg] - zmelt2;
                        exmelt = exmelt - zmelt2*( cct*( ts-tf )
                                                   * asnow + snomel );
                        dts  = dtsg + exmelt/cct;
                    }

                    else
                    {
                        double cool = amin1(0.0, tf-0.01 - (ts+dtsg))
                                      * cct * (1.0-areas);
                        double dtsg2 = amax1 ( cool, exmelt ) / ( cct*(1.001-areas));
                        exmelt = exmelt - dtsg2*cct*(1.-areas);
                        double dtsg3 =exmelt/cctt;
                        dts = dtsg + dtsg2 + dtsg3;
                    }
                }
            }
        }

        www[0] = www[0] + zmelt / ( poros * zdepth[0] );

        dtc = dtc*realg + dts*realc;
        dtg = dtg*realc + dts*realg;
    }

    double fluxef = shf - cg*dtg/dtt;

    // 3/96 changes
    // old  dtd = fluxef / ( cg * 2. * sqrt ( pie*365. ) ) * dtt;
    dtd = fluxef / ( csoil * 2. * sqrt ( PI_NUMBER*365. ) ) * dtt;

    return;
}


/* ****************************************************************************
   calculation of interflow, infiltration excess and loss to
   groundwater .  all losses are assigned to variable 'roff'
**************************************************************************** */
void SiB2Model::run2()
{
    double temw[3], temwp[3], temwpp[3];
    double aaa[2], bbb[2], ccc[2], qqq[2];

    for (int i=0; i<3; i++)
    {
        temw[i]   = amax1( 0.03, www[i] );
        temwp[i]  = pow(temw[i], -bee );
        temwpp[i] = pow((amin1( 1., temw[i])),(2.0*bee + 3.0 ));
    }

    /* *********************************************************************
    calculation of gravitationally driven drainage from w[2] : taken as an
    integral of time varying conductivity.addition of liston baseflow term
    to original q3g to insure flow in dry season. modified liston baseflow
    constant scaled by available water.

      q3g (q3) : equation (62) , SE-86
    ********************************************************************* */

    double pows = 2.*bee+2.0;
    double q3g = pow(temw[2],-pows) + satco/zdepth[2]/poros*slope*pows*dtt;
    q3g = pow(q3g, ( 1. / pows ));
    q3g = - ( 1. / q3g - www[2] ) * poros * zdepth[2] / dtt;
    q3g = amax1( 0., q3g );
    q3g = amin1( q3g, www[2]*poros*zdepth[2]/dtt );

    q3g = q3g + 0.002*poros*zdepth[2]*0.5 / 86400. * www[2];

    /* *********************************************************************
    calculation of inter-layer exchanges of water due to gravitation and
    hydraulic gradient. the values of w(x) + dw(x) are used to calculate
    the potential gradients between layers.
    modified calculation of mean conductivities follows ME-82 ),
    reduces recharge flux to top layer.

      dpdw           : estimated derivative of soil moisture potential
                       with respect to soil wetness. assumption of
                       gravitational drainage used to estimate likely
                       minimum wetness over the time step.

      qqq  (q     )  : equation (61) , SE-86;
             i,i+1

      avk  (k     )  : equation (4.14) , ME-82
             i,i+1
    ********************************************************************* */

    double wmax = amax1( www[0], amax1(www[1], amax1(www[2], 0.05 )));
    wmax = amin1( wmax, 1.0 );
    double pmax = pow(wmax,-bee);
    double wmin = pow((pmax-2./(phsat*(zdepth[0]+2.*zdepth[1]+zdepth[2]))),
                      (-1./bee));
    wmin = amin1( www[0], amin1(www[1], amin1(www[2], wmin)));
    wmin = amax1( wmin, 0.02 );
    double pmin = pow(wmin,-bee);
    double dpdw = phsat*( pmax-pmin )/( wmax-wmin );

    for (int i = 0; i < 2; i++)
    {
        double rsame = 0.0;
        double avk  = temwp[i]*temwpp[i] - temwp[i+1]*temwpp[i+1];
        double div  = temwp[i+1] - temwp[i];
        if ( fabs(div) < 1.e-6 ) rsame = 1.0;
        avk = satco*avk / ( ( 1. + 3./bee ) * div + rsame );
        double avkmin = satco * amin1( temwpp[i], temwpp[i+1] );
        double avkmax = satco * amax1( temwpp[i], temwpp[i+1] )*1.01;
        avk = amax1( avk, avkmin );
        avk = amin1( avk, avkmax );

        // 3/96 changes
        if (www[i] < www[i+1]) avk = 0.1*avk;

        // conductivities and base flow reduced when temperature drops
        // below freezing
        tsnow = amin1 ( tf-0.01, tg );
        tgs = tsnow*areas + tg*(1.-areas);
        double ts    = tgs*(1-i) + td*i;
        double props = ( ts-(tf-10.) ) / 10.0;
        props = amax1( 0.05, amin1( 1.0, props ) );
        avk  = avk * props;
        q3g  = q3g * props;

        // backward implicit calculation of flows between soil layers
        double dpdwdz = dpdw * 2./( zdepth[i] + zdepth[i+1] );
        aaa[i] = 1. + avk*dpdwdz*( 1./zdepth[i]+1./zdepth[i+1] )*dtt/poros;
        bbb[i] =-avk *   dpdwdz * 1./zdepth[1]*dtt/poros;
        ccc[i] = avk * ( dpdwdz * ( www[i]-www[i+1] ) + 1. +
                         i * dpdwdz*q3g*1./zdepth[2]*dtt/poros );
    }

    double denom  = ( aaa[0]*aaa[1] - bbb[0]*bbb[1] );
    double rdenom = 0.0;
    if ( fabs(denom) < 1.e-6 ) rdenom = 1.0;
    rdenom = ( 1.-rdenom)/( denom + rdenom );
    qqq[0]   = ( aaa[1]*ccc[0] - bbb[0]*ccc[1] ) * rdenom;
    qqq[1]   = ( aaa[0]*ccc[1] - bbb[1]*ccc[0] ) * rdenom;

    // update wetness of each soil moisture layer due to layer interflow
    // and base flow.

    www[2] = www[2] - q3g*dtt/(poros*zdepth[2]);
    roff = roff + q3g * dtt;

    for (int i = 0; i < 2; i++)
    {
        double qmax   =  www[i]   * (poros*zdepth[i]  /dtt);
        double qmin   = -www[i+1] * (poros*zdepth[i+1]/dtt);
        qqq[i] = amin1( qqq[i],qmax);
        qqq[i] = amax1( qqq[i],qmin);
        www[i]   =   www[i]   - qqq[i]/(poros*zdepth[i]  /dtt);
        www[i+1] =   www[i+1] + qqq[i]/(poros*zdepth[i+1]/dtt);
    }

    for (int i = 0; i < 3; i++)
    {
        double excess = amax1(0.0,(www[i] - 1.0));
        www[i] = www[i] - excess;
        roff   = roff   + excess * poros*zdepth[i];
    }

    // prevent negative values of www[i]
    for (int i = 0; i < 2; i++)
    {
        double deficit   = amax1 (0.,(1.e-12 - www[i]));
        www [i]   = www[i] + deficit;
        www [i+1] = www[i+1] - deficit * zdepth[i] / zdepth [i+1];
    }
    www[2]    = amax1 (www[2],1.e-12);

    return;
}


/* ****************************************************************************
   output of results to files
**************************************************************************** */
void SiB2Model::outer (ofstream& out, ofstream& out1, ofstream& out2, ofstream& out3,
                       ofstream& out4, const int nymd)
{
    double trant  = ect/hlat;
    double canil  = eci/hlat;
    double evapg  = (egs+egi)/hlat;
    double radswd = radn[0][0] + radn[0][1] + radn[1][0] + radn[1][1];
    double radswa = (1.-salb[0][0])*radn[0][0] + (1.-salb[0][1])*radn[0][1]
                    + (1.-salb[1][0])*radn[1][0] + (1.-salb[1][1])*radn[1][1];
    double radtot = radt[0] + radt[1];
    double elat   = etmass / dtt * hlat;
    double gcstor = chf + shf;

    tgs  = amin1(tf,tg)*areas + tg*(1.-areas);
    double tc4  = tc  * tc  * tc  * tc;
    double tg4  = tgs * tgs * tgs * tgs;
    double zlwf =  tc4 * stefan * vcover * ( 1. - thermk )
                   + (1. - vcover * ( 1. - thermk ) ) * stefan * tg4;
    tgeff= sqrt ( sqrt ( zlwf/stefan ) );


    out << nymd << "  " << www[0] << " " << www[1] << " " << www[2] << " "
        << ppl << " " << endl;
    out1 << nymd << "  " << roff << "  " << trant << "  " << canil << "  "
         << evapg << "  " << etmass << "  " << ppl << endl;

    out2 << nymd  << "  " << tc  << "  " << tg  << "  " << td  << "  "
         << tm  << "  " << tgeff  << "  " << capac[0]  << "  " << capac[1]
         << "  " << snoww[0]  << "  " << snoww[1]  << endl;

    out3 << nymd  << "  " << radn[2][1]  << "  " << zlwup  << "  "
         << radswd  << "  " << radswa  << "  " << radtot  << "  "
         << elat  << "  " << hflux  << "  " << gcstor << endl;
}


double SiB2Model::amin1(double a, double b)
{
    if (a<b) return a;
    else return b;
}

double SiB2Model::amax1(double a, double b)
{
    if (a>=b) return a;
    else return b;
}

int SiB2Model::min0(int a, int b)
{
    if (a<b) return a;
    else return b;
}

int SiB2Model::max0(int a, int b)
{
    if (a>=b) return a;
    else return b;
}

double SiB2Model::sign(double a, double b)
{
    if (b >= 0) return fabs(a);
    else return -fabs(a);
}

double SiB2Model::e(double x)
{
    return exp( 21.18123 - 5418.0 / x ) / 0.622;
}

double SiB2Model::ge(double x)
{
    return exp( 21.18123 - 5418.0 / x ) * 5418.0 / (x * x) / 0.622;
}
