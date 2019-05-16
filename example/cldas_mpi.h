#ifndef __LDAS_CLDAS_H
#define __LDAS_CLDAS_H

#include "commonlandmodel.h"

namespace ldas
{
enum LAND_TYPE{LAND_UNUSED,LAND_SOIL,LAND_WATER,LAND_SNOW,LAND_FROSEN_SOIL};
//define AMSR brightness temperature
const double  amsr_freq_6 = 6.925;
const double  amsr_freq_10 = 10.65;
const double  amsr_freq_18 = 18.7;
const double  amsr_freq_23 = 23.8;
const double  amsr_freq_37 = 36.5;
const double  amsr_freq_89 = 89.0;
//define SSMI brightness temperature
const double  ssmi_freq_19 = 19.35;
const double  ssmi_freq_22 = 22.235;
const double  ssmi_freq_37 = 37;
const double  ssmi_freq_85 = 85.5;

// incident angle
const double  ssmi_theta  = 53.1;
const double  amsr_theta  = 55;

// observation variacne
const double R_ssmi_19v = 9;
const double R_ssmi_19h = 16;
const double R_ssmi_37v = 16;
const double R_ssmi_37h = 25;

const double R_amsr_6v  = 4;
const double R_amsr_6h  = 8;
const double R_amsr_10v = 4;
const double R_amsr_10h = 8;
const double R_amsr_18v = 9;
const double R_amsr_18h = 12;
const double R_amsr_37v = 16;
const double R_amsr_37h = 20;
class Cldas_Observation
{
public:
	Cldas_Observation(const int r,const int c);
	virtual ~Cldas_Observation();
	virtual void open(const std::string& basename,const DATA_TYPE dt);
    short **ssmi_data_19v,**ssmi_data_19h,**ssmi_data_37v,**ssmi_data_37h,
          **amsr_data_6v,**amsr_data_6h,**amsr_data_10v,**amsr_data_10h,**amsr_data_18v,
          **amsr_data_18h,**amsr_data_37v,**amsr_data_37h;
private:
	int m_row,m_col;
};
class Cldas : public CommonLandModel
{
	public:
		/** Default constructor */
		Cldas(int,char**);
		/** Default destructor */
		virtual ~Cldas();
		/** Copy constructor
		 *  \param other Object to copy from
		 */
		Cldas(const Cldas& other);
		/** Assignment operator
		 *  \param other Object to assign from
		 *  \return A reference to this
		 */
		Cldas& operator=(const Cldas& other);

		void ensemble(const int en);
		int ensemble() const;
		virtual void run();
		virtual void run(const int ipatch);
		virtual void config(const char* fn);
	protected:
		Cldas_Observation* m_observe;
		int m_ensemble;
	private:
		double*** fvar_en;
		double*** flux_en;
		double** fvar_mean;
		double** flux_mean;
		int** snl_en;
		float** sand_grid;
		float** clay_grid;
		int sensor;
		bool obs_flag;
		int numprocs, myrank;
};
}
#endif // __LDAS_CLDAS_H
