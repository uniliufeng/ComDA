#ifndef JMASIB_H
#define JMASIB_H

#include <BaseModel.h>
/* ************************************************************************************
// jmasib.hpp
// This is a program to run the JMA-new-SiB in one dimension for a LDAS
// original author: Dr. Hirai at JMA
// Modified for LDAS purpose by <Li Xin @ CAREERI/CAS> in 2005
// Major modifications:
// 1. main_newsib.F-->newSib.cpp, a main program in c++
// 2. c++ APIs for data assimilation purpose are developed
************************************************************************************ */
namespace ldas
{
// Global definition for running JMA-new-SiB
// the namelists referenced by main_newsib.F
#define IMAX	1				// grid columns max
#define JMAX	1				// grid rows max
#define IDIM	IMAX			// grid columns
#define JDIM	JMAX			// grid rows
#define IJPHY	IMAX * JMAX		//
#define JLPHY	1				//
#define ISPT	1				// vegetation types in the specific grid
#define IDP		3				// soil layers
#define ISN		4				// snow layers
#define IVN		2				// canopy layers

// definition of state vector (initial variables)
extern "C" double SIB0109_mp_TMP_CNP_NOS_ALL	[JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_FR_WTR_CNP_NOS_ALL	[JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_FR_GLA_CNP_NOS_ALL	[JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_TMP_CNP_SNW_ALL	[JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_FR_WTR_CNP_SNW_ALL	[JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_FR_GLA_CNP_SNW_ALL	[JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_TMP_GRSK_ALL		[JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_FR_WTR_GRS_ALL		[JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_FR_GLA_GRS_ALL		[JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_TMP_SNSK_ALL		[JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_TMP_SOIL_ALL		[IDP][JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_FR_WTR_SOIL_ALL	[IDP][JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_FR_GLA_SOIL_ALL	[IDP][JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_TMP_SNSL_ALL		[IDP][JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_FR_WTR_SNSL_ALL	[IDP][JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_FR_GLA_SNSL_ALL	[IDP][JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_TMP_SNOW_ALL		[ISN][JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_WTR_SNOW_ALL		[ISN][JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_GLA_SNOW_ALL		[ISN][JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_RHO_SNOW_INV_ALL	[ISN][JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_AGE_SNOW_ALL		[JLPHY][ISPT*IJPHY];
extern "C" int	  SIB0109_mp_INFO_SNOW_ALL		[JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_ENG_SNOW_BUCKET_ALL[JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_H2O_SNOW_BUCKET_ALL[JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_TMP_SOIL4_ALL		[JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_ALB_SNOW_SNW_ALL	[IVN][JLPHY][ISPT*IJPHY];
extern "C" double SIB0109_mp_RAIN_1HOUR_ALL		[JLPHY][ISPT*IJPHY];

// Declaration of fortran subrountines in main_newsib.F
extern "C"
{
    void __stdcall	COM_STDINOUT_UNIT_mp_COM_STDINOUT_UNIT_INI(int* in, int* out);
    void __stdcall  COM_DEBUG_mp_COM_DEBUG_INI();
    void __stdcall  MESSAGE_mp_MESSAGE_INI();
    void __stdcall  COM_RUNCONF_mp_COM_RUNCONF_INI();
    void __stdcall  COM_RUNCONF_SIB0109_mp_COM_RUNCONF_SIB0109_INI();
    void __stdcall  COM_JOBINFO_SIB0109_mp_COM_JOBINFO_SIB0109_INI();
    void __stdcall  COM_TETEN_SIB_0109_mp_COM_TETEN_SIB_0109_INI();
    void __stdcall  CALENDAR_mp_CALENDAR_INI();
    void __stdcall  GET_KTINFO(int* IDSTAR, int* IDATE, int* KTP, int* KT0, double* FSECP);
    void __stdcall  CALENDAR_mp_CALENDAR_RUN_GETKT(int* IDATE, int* IDSTAR, int* IFLAG, int* KT);
    void __stdcall  DELT_CHECK(double* DELT_ATM, double* DELT_CNP,
                               double* DELT_SNOW, double* DELT_SOIL);
    void __stdcall  SIB0109_mp_SIB0109_INI(int IMASK[JLPHY][ISPT*IJPHY], int* IDSTAR);
    void __stdcall  MONIT(double GLON[JDIM][IDIM], double GLAT[JDIM][IDIM],
                          double* DELT_MODEL, int* MON, int* KT0, int* KTP, int* KTSTAR, int* KTEND,
                          int* IDATE, double* FSECP, int IMASK[JLPHY][ISPT*IJPHY]);
    void __stdcall  TS_mp_TIME_STEP_RESET();
    void __stdcall  TS_mp_TIME_STEP(double*, double*, int* IDSTAR, int* IDEND,
                                    double* FSEC0, double* FSECP,
                                    int* KTM, int* KT0, int* KTP, int* ISTEP, int* IDATE,
                                    int* ID_NEXT, int* ID_NOW, int* ID_PRE, double* RDAY, double* RSEC, double* TOTMON);
    void __stdcall  SIB0109_mp_SIB0109_RUN_STEPINI(double* FSECP);
    void __stdcall  CON_SET(double* GMT_PHY, int* IIJPHY, double* CONS_MTX);
    void __stdcall  SIB0109_mp_SIB0109_RUN_LOOPINI(int* JL, int* MON, int* KT0,
            double GPDEL_PHY[JLPHY][IJPHY], double GPHAF_PHY[JLPHY][IJPHY],
            double GPFUL_PHY[JLPHY][IJPHY],	double GT_PHY[JLPHY][IJPHY],
            double GQ_PHY[JLPHY][IJPHY], 	double GU_PHY[JLPHY][IJPHY],
            double GV_PHY[JLPHY][IJPHY], 	double ZMEAN_PHY[JLPHY][IJPHY],
            double ZTEMP_PHY[JLPHY][IJPHY],	double PPLI_PHY[JLPHY][IJPHY],
            double PPCI_PHY[JLPHY][IJPHY],  double DAYTIME_1HR_PHY[JLPHY][IJPHY],
            double GMT_PHY[2][IJPHY],       double GMQ_PHY[2][IJPHY],
            double GMUV_PHY[3][IJPHY]);
    void __stdcall  SIB0109_mp_SIB0109_RUN_ALBEDO(int* MON, int* JL,
            double AVISB[IJPHY], double ANIRB[IJPHY],
            double AVISD[IJPHY], double ANIRD[IJPHY]);
    void __stdcall	SIB0109_mp_SIB0109_RUN_SR_CALC(int* MON, int* JL,
            double RVISB[JLPHY][IJPHY], double RVISD[JLPHY][IJPHY],
            double RNIRB[JLPHY][IJPHY], double RNIRD[JLPHY][IJPHY]);
    void __stdcall	SIB0109_mp_SIB0109_RUN_SR_RESTORE(int* JL);
    void __stdcall	SIB0109_mp_SIB0109_RUN_LR_SET(double DLWB[JLPHY][IJPHY]);
    void __stdcall	SIB0109_mp_SIB0109_RUN_SIBMAIN(int* MON, int* JL, int* ID_NEXT, int* ID_NOW);
    void __stdcall	SIB0109_mp_SIB0109_RUN_LOOPEND(int* JL,
            double* TMTX2L_PHY, double* QMTX2L_PHY,
            double* UMTX2L_PHY, double* VMTX2L_PHY, double* RAD_LONG_SIB_2_ATM_PHY);
    void __stdcall	OUTPUT(int* IDATE, int* KTP, double* TOTMON, double* FSECP, int IMASK[JLPHY][ISPT*IJPHY]);
}
class JMASiBForcing
{
public:
    double	RefLevel;				// reference level (m)
    double	SolarZenith;			// cosine solar zenith angle temporal (-)
    double	SolarZenithPrevious;	// cosine solar zenith angle of averaged pre 1 hour (-)
    double	AirTemperature;			// air temperature (K)
    double	AirPressure;			// air pressure (mb)
    double	AirPressureSurface;		// air pressure at surface (mb)
    double	Humidity;				// specific humidity (kg kg^-1)
    double	PrecipLargescale;		// rainfall due to large scale condensation (mm)
    double	PrecipGridscale;		// rainfall due to sub-grid scale convection (mm)
    double	WindSpeedU;				// wind speed U (m s^-1)
    double	WindSpeedV;				// wind speed V (m s^-1)
    double	VisDir;					// direct solar radiation visible (W m^-2)
    double	VisDif;					// diffuse solar radiation visible (W m^-2)
    double	NearInfraredDir;		// direct solar radiation infrared (W m^-2)
    double	NearInfraredDif;		// diffuse solar radiation infrared (W m^-2)
    double	LongDown;				// longwave radiation downward (W m^-2)
};
class JMASiB : public BaseModel
{
public:
    /** Default constructor */
    JMASiB();
    /** Default destructor */
    virtual ~JMASiB();
    /** Copy constructor
     *  \param other Object to copy from
     */
    JMASiB(const JMASiB& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    JMASiB& operator=(const JMASiB& other);
    virtual run();
    virtual void init();
    virtual void forcing();
protected:
private:
    int		IMASK[JLPHY][ISPT*IJPHY];	// land type
    double	GLON[JDIM][IDIM];			// longitude
    double	GLAT[JDIM][IDIM];			// latitude
    int		JL;							// set to 1
    int		MON;						// month

    // atmospheric radiation related
    double	RNIRB[JLPHY][IJPHY];		// near-infrared beam
    double	RNIRD[JLPHY][IJPHY];		// near-infrared diffuse
    double	RVISB[JLPHY][IJPHY];		// shortwave beam
    double	RVISD[JLPHY][IJPHY];		// shortwave diffuse
    double	DLWB[JLPHY][IJPHY];			// longwave downward

    // atmosphere related
    // input
    double	GMT_PHY[2][IJPHY];
    double	GMQ_PHY[2][IJPHY];
    double	GMUV_PHY[3][IJPHY];

    // output
    double	TMTX2L_PHY[IJPHY];
    double	QMTX2L_PHY[IJPHY];
    double	UMTX2L_PHY[IJPHY];
    double	VMTX2L_PHY[IJPHY];
    double	AVISB[IJPHY];
    double	ANIRB[IJPHY];
    double	AVISD[IJPHY];
    double	ANIRD[IJPHY];
    double	RAD_LONG_SIB_2_ATM_PHY[IJPHY];

    // atmosphere
    double	ZMEAN_PHY[JLPHY][IJPHY];
    double	ZTEMP_PHY[JLPHY][IJPHY];
    double	PPLI_PHY[JLPHY][IJPHY];
    double	PPCI_PHY[JLPHY][IJPHY];
    double	GT_PHY[JLPHY][IJPHY];
    double	GQ_PHY[JLPHY][IJPHY];
    double	GU_PHY[JLPHY][IJPHY];
    double	GV_PHY[JLPHY][IJPHY];
    double	GPFUL_PHY[JLPHY][IJPHY];
    double	GPHAF_PHY[JLPHY][IJPHY];
    double	GPDEL_PHY[JLPHY][IJPHY];

    // control variables
    int		IDSTAR[5];					// initial time of integration
    int		IDEND[5];					// end time of integration
    int		IDATE[5];					// valid time
    int		KTSTAR;						// model integration start time
    int		KTEND;						// model integration end time
    double	FSECP;
    int		KTP;						// forecast time + delta t
    int		KT0;						// forecast time
};
}
#endif // JMASIB_H
