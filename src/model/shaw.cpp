#include "shaw.h"

using namespace ldas;
using namespace shaw;
using namespace std;
SHAW::SHAW()
{
    //ctor
}

SHAW::~SHAW()
{
    //dtor
}

SHAW::SHAW(const SHAW& other)
{
    //copy ctor
}

SHAW& SHAW::operator=(const SHAW& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

void SHAW::run()
{
    INITAL=0;
    input();

    if (INPH2O!=1)
        VLCDAY[HOUR-1]=VLCDT[NS-1] + VICDT[NS-1]*RHOI/RHOL;
    else
        VLCDAY[HOUR-1]=MATDT[NS-1];
    SOITMP[HOUR-1]=TSDT[NS-1];

    dayinput();
    cloudy();

    while((JULIAN <= JEND && YEAR == YREND)||(YEAR < YREND))
    {

        DT = NHRPDT*3600.0;                // seconds in one model step
        HOUR = HOUR+NHRPDT;           // next hour (NHRPDT=1 hour)
        if(HOUR > 24)                           // next day
        {
            HOUR = HOUR-24;
            JULIAN = JULIAN + 1;
            if (JULIAN > MAXJUL)
            {
                JULIAN = JULIAN-MAXJUL;
                ++YEAR;
                MAXJUL = 365;
                if (m_time.is_leap_year(YEAR))
                    MAXJUL=366;
            }
            if((YEAR == YREND && JULIAN > JEND)||(YEAR > YREND))
                break;
            dayinput();
            cloudy();
        }
        if(JULIAN == LEVEL[1] && HOUR == LEVEL[2])
        {
            // change LEVEL of output
            int LVL = LEVEL[0];
            LEVEL[0] = LEVEL[3];
            LEVEL[1] = LEVEL[4];
            LEVEL[2] = LEVEL[5];
            LEVEL[3] = LVL;
            LEVEL[4] = JULIAN;
            LEVEL[5] = HOUR;
        }
        step();
        if(INITAL == 0)
            INITAL = 1;
    }
}

void SHAW::input()
{
    input_(&NC,&NSP,&NR,&NS,&TOLER,LEVEL,&MTSTEP,&MZCINP,&MPLTGRO,
           &INPH2O,&MWATRXT,LVLOUT,&IVLCBC,&ITMPBC,&NPLANT,PLTHGT,PLTWGT,PLTLAI,
           ROOTDP,DCHAR,TCCRIT,RSTOM0,RSTEXP,PLEAF0,RLEAF0,RROOT0,
           CANALB,&CANMA,&CANMB,&WCMAX,LANGLE,ITYPE,ZC,WCANDT,
           ZSP,DZSP,RHOSP,TSPDT,DLWDT,ICESPT,&SNOTMP,
           GMCDT,&GMCMAX,&RLOAD,&ZRTHIK,&COVER,&ALBRES,&RESCOF,
           ZS,TSDT,VLCDT,VICDT,MATDT,CONCDT,ICESDT,SALTDT,&ALBDRY,&ALBEXP,
           DGRADE,SLTDIF,ASALT,DISPER,
           &ZMSRF,&ZHSRF,&ZERSRF,&ZMSP,&ZHSP,&HEIGHT,&PONDMX,
           &JSTART,&YRSTAR,&HRSTAR,&JEND,&YREND,&NHRPDT,&WDT,
           &ALATUD,&SLOPE,&ASPECT,&HRNOON);

    JULIAN=JSTART;           // current julian day of model
    HOUR=HRSTAR;             // current hour of model
    YEAR=YRSTAR;             // current year of model
    MAXJUL=365;
    if (m_time.is_leap_year(YEAR))
        MAXJUL=366;
}
void SHAW::cloudy()
{
    cloudy_(&CLOUDS,&ALATUD,&DECLIN,&HAFDAY,SUNHOR,&JULIAN,&NHRPDT);
}
void SHAW::step()
{
    goshaw_(&JULIAN,&HOUR,&YEAR,&NHRPDT,&WDT,&DT,&INITAL,&NC,&NSP,&NR,&NS,&TOLER,LEVEL,&MZCINP,
            &INPH2O,&MWATRXT,LVLOUT,&IVLCBC,&ITMPBC,&NPLANT,PLTHGT,PLTWGT,PLTLAI,
            ROOTDP,DCHAR,TCCRIT,RSTOM0,RSTEXP,PLEAF0,RLEAF0,RROOT0,PCANDT,
            CANALB,&CANMA,&CANMB,&WCMAX,LANGLE,ITYPE,ZC,TCDT,VAPCDT,WCANDT,
            ZSP,DZSP,RHOSP,TSPDT,DLWDT,ICESPT,WLAG,&STORE,&SNOWEX,&SNOTMP,
            ZR,RHOR,TRDT,VAPRDT,GMCDT,&GMCMAX,&RLOAD,&ZRTHIK,&COVER,&ALBRES,&RESCOF,
            &DIRRES,ZS,TSDT,VLCDT,VICDT,MATDT,CONCDT,ICESDT,SALTDT,&ALBDRY,&ALBEXP,
            DGRADE,SLTDIF,ASALT,DISPER,&ZMSRF,&ZHSRF,&ZERSRF,&ZMSP,&ZHSP,&HEIGHT,&POND,&PONDMX,
            &ALATUD,&SLOPE,&ASPECT,&HRNOON,&CLOUDS,&DECLIN,&HAFDAY,
            &SUNHOR[HOUR-1],&TMPDAY[HOUR-1],&WINDAY[HOUR-1],&HUMDAY[HOUR-1],&PRECIP[HOUR-1],
            &SNODEN[HOUR-1],&SOITMP[HOUR-1],&VLCDAY[HOUR-1],&SOILXT[HOUR-1][0]);
}
void SHAW::output()
{

}
void SHAW::dayinput()
{
    dayinp_(&JULIAN,&YEAR,&MAXJUL,&HOUR,&NHRPDT,&MTSTEP,&MPLTGRO,&MWATRXT,&NS,&ALATUD,&HRNOON,SUNHOR,TMPDAY,
            WINDAY,HUMDAY,PRECIP,SNODEN,SOITMP,VLCDAY,SOILXT,&NPLANT,PLTHGT,DCHAR,PLTWGT,PLTLAI,ROOTDP);
}

void SHAWConfig::open(const char* fn)
{
    char buf[MAXBUFFER];
    ifstream in(fn);
    if(!in)
    {
        string error_msg="can't open config data file: ";
        error_msg+=fn;
        throw Exception("SHAWConfig","open(fn)",error_msg.c_str());
    }
    in.getline(buf,MAXBUFFER,' ');
    std::istringstream iss;
    iss.str(buf);
    iss>>mtstep;
    //cout<<mtstep<<endl;
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>inph2o;
    //cout<<total_steps<<endl;
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>mwatrxt;

    in.getline(buf,MAXBUFFER);
    site_path=buf;
    in.getline(buf,MAXBUFFER);
    weather_path=buf;
    in.getline(buf,MAXBUFFER);
    moisture_path=buf;
    in.getline(buf,MAXBUFFER);
    temperature_path=buf;

    if (mwatrxt!=0)
    {
        in.getline(buf,MAXBUFFER);
        soilsink_path=buf;
    }

    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>output_frequency;
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>comparison_flag;
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>temperature_frequency;
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>watercontent_frequency;
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>matricpotential_frequency;
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>energy_frequency;
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>waterbalance_frequency;
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>rootwater_frequency;
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>frostdepth_frequency;
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>salt_frequency;
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>soilsolution_frequency;
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>updatescreen_frequency;

    in.getline(buf,MAXBUFFER);
    outputgeneral_path=buf;
    in.getline(buf,MAXBUFFER);
    profile_path=buf;
    in.getline(buf,MAXBUFFER);
    simulated_temperature_path=buf;
    in.getline(buf,MAXBUFFER);
    watercontent_path=buf;
    in.getline(buf,MAXBUFFER);
    waterpotential_path=buf;
    in.getline(buf,MAXBUFFER);
    energyflux_path=buf;
    in.getline(buf,MAXBUFFER);
    waterbalance_path=buf;
    in.getline(buf,MAXBUFFER);
    waterflow_path=buf;
    in.getline(buf,MAXBUFFER);
    rootwater_path=buf;
    in.getline(buf,MAXBUFFER);
    frostdepth_path=buf;
    in.getline(buf,MAXBUFFER);
    salt_path=buf;
    in.getline(buf,MAXBUFFER);
    soilsolution_path=buf;

    in.close();
}
