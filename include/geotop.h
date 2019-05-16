#ifndef __LDAS_GEOTOP_H
#define __LDAS_GEOTOP_H

#include "model.h"

extern "C" {
#include <sys/stat.h>
#include "struct.geotop.h"
#include "input.h"
#include "output.h"
#include "times.h"
#include "constant.h"
#include "keywords_file.h"
#include "energy.balance.h"
#include "meteo.h"
#include "water.balance.h"

#include "geomorphology.0875.h"
#include "pedo.funct.h"
#include "geo_statistic.h"
#include "networks.h"

#include "dtm_resolution.h"
#include "rw_maps.h"
#include "extensions.h"
#include "tabs.h"
#include "snow.h"
#include "micromet.h"
#include "vegetation.h"
#include "get_filenames.h"
}
extern T_INIT *UV;
extern 	STRINGBIN *files;
extern long Nl,Nr,Nc;
extern double NoV;
extern DOUBLEMATRIX *outdata_point;
extern DOUBLEVECTOR *outdata_basin;
extern char *MISSING_FILE;
namespace ldas
{
	namespace geotop
	{
const int jdz=1;						//layer thickness [mm]
const float z_evap=40.; //soil depth responsable for soil evaporation [mm]
const int jdvegprop=8;
const int nmet=12;
	}
class GeotopParameter
{
public:
    GeotopParameter();
    virtual ~GeotopParameter();
    /** Copy constructor
     *  \param other Object to copy from
     */
    GeotopParameter(const GeotopParameter& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    GeotopParameter& operator=(const GeotopParameter& other);
    void open(const char*);
    void save();
    ALLDATA* alldata();
    SOIL *S;
    WATER *W;
    LAND *L;
    PAR *P;
    TOPO *T;
    CHANNEL *C;
    ENERGY *E;
    SNOW *N;
    GLACIER *G;
    METEO *M;
    TIMES *I;
private:
    //std::string workdir;
    ALLDATA* adt;
};
class Geotop : public BaseModel
{
public:
    Geotop();
    Geotop(int Dt, double JD0, int year0, double TH,double ST, double Dt_output_discharge,
           double Dt_output_pixel, double Dt_output_basin, short nsky, double channel_thres,
           short format_out, short point_sim,short recover);
    void parameters(int Dt, double JD0, int year0, double TH,double ST, double Dt_output_discharge,
                    double Dt_output_pixel, double Dt_output_basin, short nsky, double channel_thres,
                    short format_out, short point_sim,short recover);

    void init();
    void config(const char*);

    virtual ~Geotop();
    Geotop(const Geotop& other);
    Geotop& operator=(const Geotop& other);

    void step();
    virtual void run();
    virtual void output();

    GeotopParameter *adt;

private:

// base parameters in geotop
    int 		Dt;					//THE INTEGRATION INTERVAL [s]
    double 		JD0;				//Decimal julian day at the beginning of simulation (0.0 - 365.99)
    int 		year0;				//Year at the beginning of simulation
    double 		TH;					//THE NUMBER OF DAYS OF SIMULATION
    double 		ST;					//standard time		Standard time to which all the output data are referred (difference respect UMT, in hour)
    double 		Dt_output_discharge;//Dt_output_discharge	Delta TIME (in hour) with which THE WATER DISCHARGE AT THE BASIN OUTLET IS PRINTED
    double 		Dt_output_pixel;	//Delta TIME (in hour) with which THE OUTPUT FOR SPECIFIED PIXELS IS PRINTED
    double 		Dt_output_basin;	//Delta TIME (in hour) with which THE OUTPUT FOR BASIN-ALTIMETRIC RANK VALUES IS PRINTED (if applicable)
    short 		nsky;				//Multiplying factor decreasing the dem remaps/solution for the calculation of the sky view factor
    double 		channel_thres;		//Value of the threshold for definition of pixel channel [1.0-1000.0] (used if the network map is not provided)
    short 		format_out;			//OUTPUT MAPS in fluidturtle format (=1), GRASS ASCII (=2), ESRI ASCII (=3)
    short 		point_sim;			//=0 distributed simulation, =1 point simulation (the parameter files are different in the two cases)
    short 		recover;			//=1 if you want to recover a simulation, 0 otherwise


    //T_INIT *UV;
};
}
#endif // __LDAS_GEOTOP_H
