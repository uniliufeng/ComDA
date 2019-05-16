#include "jmasib.h"

JMASiB::JMASiB()
{
    //ctor
}

JMASiB::~JMASiB()
{
    //dtor
}

JMASiB::JMASiB(const JMASiB& other)
{
    //copy ctor
}

JMASiB& JMASiB::operator=(const JMASiB& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

void JMASiB::init()
{
    SIB0109_mp_TMP_CNP_NOS_ALL[0][0]		= state.CanopyTemp;
    SIB0109_mp_FR_WTR_CNP_NOS_ALL[0][0]		= state.CanopyWater;
    SIB0109_mp_FR_GLA_CNP_NOS_ALL[0][0]		= state.CanopyIce;
    SIB0109_mp_TMP_CNP_SNW_ALL[0][0]		= state.CanopyTempSnow;
    SIB0109_mp_FR_WTR_CNP_SNW_ALL[0][0]		= state.CanopyWaterSnow;
    SIB0109_mp_FR_GLA_CNP_SNW_ALL[0][0]		= state.CanopyIceSnow;
    SIB0109_mp_TMP_GRSK_ALL[0][0]			= state.GroundSurfaceTemp;
    SIB0109_mp_FR_WTR_GRS_ALL[0][0]			= state.GroundSurfaceWater;
    SIB0109_mp_FR_GLA_GRS_ALL[0][0]			= state.GroundSurfaceIce;
    SIB0109_mp_TMP_SNSK_ALL[0][0]			= state.SnowSurfaceTemp;
    for (int m=0; m<3; m++)
    {
        SIB0109_mp_TMP_SOIL_ALL[m][0][0]	= state.SoilTemp[m];
        SIB0109_mp_FR_WTR_SOIL_ALL[m][0][0]	= state.SoilWater[m];
        SIB0109_mp_FR_GLA_SOIL_ALL[m][0][0]	= state.SoilIce[m];
        SIB0109_mp_TMP_SNSL_ALL[m][0][0]	= state.SoilTempSnow[m];
        SIB0109_mp_FR_WTR_SNSL_ALL[m][0][0]	= state.SoilWaterSnow[m];
        SIB0109_mp_FR_GLA_SNSL_ALL[m][0][0]	= state.SoilIceSnow[m];
    }
    for (int n=0; n<4; n++)
    {
        SIB0109_mp_TMP_SNOW_ALL[n][0][0]	= state.SnowTemp[n];
        SIB0109_mp_WTR_SNOW_ALL[n][0][0]	= state.SnowWater[n];
        SIB0109_mp_GLA_SNOW_ALL[n][0][0]	= state.SnowIce[n];
        SIB0109_mp_RHO_SNOW_INV_ALL[n][0][0]= 1.0/state.SnowDensity[n];
    }
    SIB0109_mp_AGE_SNOW_ALL[0][0]			= state.SnowAge;
    SIB0109_mp_INFO_SNOW_ALL[0][0]			= state.SnowInfo;
    SIB0109_mp_ENG_SNOW_BUCKET_ALL[0][0]	= state.SnowHeat;
    SIB0109_mp_H2O_SNOW_BUCKET_ALL[0][0]	= state.SWE;
}

void JMASiB::forcing(const int current)
{

    GU_PHY[0][0]	= force[current].WindSpeedU;
    GV_PHY[0][0]	= force[current].WindSpeedV;
    GT_PHY[0][0]	= force[current].AirTemperature;
    GQ_PHY[0][0]	= force[current].Humidity;
    GPHAF_PHY[0][0]	= force[current].AirPressureSurface;
    GPFUL_PHY[0][0]	= force[current].AirPressure;
    GPDEL_PHY[0][0]	= force[current].RefLevel;
    ZTEMP_PHY[0][0]	= force[current].SolarZenith;
    ZMEAN_PHY[0][0]	= force[current].SolarZenithPrevious;
    PPLI_PHY[0][0]	= force[current].PrecipLargescale;
    PPCI_PHY[0][0]	= force[current].PrecipGridscale;

    RVISB[0][0]		= force[current].VisDir;
    RVISD[0][0]		= force[current].VisDif;
    RNIRB[0][0]		= force[current].NearInfraredDir;
    RNIRD[0][0]		= force[current].NearInfraredDif;
    DLWB[0][0]		= force[current].LongDown;

}

void JMASiB::run()
{
    double TOTMON = 0.0;
    int ISTEP = 0;
    static double DUMMY_SCALAR = 0.0;
    static double RSEC;						// current location in 1 day
    static double RDAY;						// current location in 1 year
    static int KTM;
    static double	FSEC0;

    static int	ID_PRE[5];
    static int	ID_NOW[5];
    static int	ID_NEXT[5];

    if (AssimilationMode==1)
    {
        TimeSpan span( (long int)(DeltaModel) );
        CurrentTime	+= span;
        CurrentStep ++;
    }

    MON = CurrentTime.GetMonth();

    TS_mp_TIME_STEP(&DeltaModel, &DUMMY_SCALAR, IDSTAR, IDEND, &FSEC0, &FSECP,
                    &KTM, &KT0, &KTP, &ISTEP, IDATE, ID_NEXT, ID_NOW, ID_PRE, &RDAY, &RSEC, &TOTMON) ;

    SIB0109_mp_SIB0109_RUN_STEPINI(&FSECP);

    int IIJPHY	= IJPHY;
    int IIJPHY2	= 2*IJPHY;
    double CON = 1.0E6;
    CON_SET(&GMT_PHY[0][0],  &IIJPHY,   &CON);
    CON_SET(&GMQ_PHY[0][0],  &IIJPHY,   &CON);
    CON_SET(&GMUV_PHY[0][0], &IIJPHY,   &CON);
    CON = 1.0E-6;
    CON_SET(&GMT_PHY[1][0],  &IIJPHY,   &CON);
    CON_SET(&GMQ_PHY[1][0],  &IIJPHY,   &CON);
    CON_SET(&GMUV_PHY[1][0], &IIJPHY2,  &CON);

    // Loop initialization
    double	DAYTIME_1HR_PHY[JLPHY][IJPHY];
    DAYTIME_1HR_PHY[0][0] = DeltaModel;

    SIB0109_mp_SIB0109_RUN_LOOPINI(&JL, &MON, &KT0,
                                   GPDEL_PHY, GPHAF_PHY, GPFUL_PHY,
                                   GT_PHY, GQ_PHY, GU_PHY, GV_PHY,
                                   ZMEAN_PHY, ZTEMP_PHY, PPLI_PHY, PPCI_PHY,
                                   DAYTIME_1HR_PHY,  GMT_PHY, GMQ_PHY, GMUV_PHY);

    // radiation calculation
    // albedo
    SIB0109_mp_SIB0109_RUN_ALBEDO(&MON, &JL, AVISB, ANIRB, AVISD, ANIRD);

    // shortwave radiation
    SIB0109_mp_SIB0109_RUN_SR_CALC(&MON, &JL, RVISB, RVISD, RNIRB, RNIRD);
    SIB0109_mp_SIB0109_RUN_SR_RESTORE(&JL);

    // longwave radiation
    SIB0109_mp_SIB0109_RUN_LR_SET(DLWB);

    // land calculation
    SIB0109_mp_SIB0109_RUN_SIBMAIN(&MON, &JL, ID_NEXT, ID_NOW);

    SIB0109_mp_SIB0109_RUN_LOOPEND(&JL, TMTX2L_PHY, QMTX2L_PHY, UMTX2L_PHY,
                                   VMTX2L_PHY, RAD_LONG_SIB_2_ATM_PHY);

    // OUTPUT(IDATE, &KTP, &TOTMON, &FSECP, IMASK);


}
