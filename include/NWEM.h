#ifndef __LDAS_NWEM_H
#define __LDAS_NWEM_H

#include "model.h"
#include "exception.h"
#include "model.h"
#include "constant.h"
//#include "Percolation_Sim.h"
#include "ChainBinomial_Sim.h"
#include "random.hpp"
#include <iostream>
#include <string>
#include "CommonMems.h"

namespace ldas
{
///Novel coronavirus (2019-nCoV) transcription Model
/**
	NetWork Epidemiological Model
Reference:Liu F, Li X, Zhu GF. Using the contact network model and Metropolis-Hastings sampling to reconstruct the COVID-19 spread on the "Diamond Princess". Science Bulletin, 2020, 10.1016/j.scib.2020.04.043.
*/
class NWEM : public BaseModel
{
public:
    /** Default constructor */
    NWEM();
    /** Default destructor */
    virtual ~NWEM();
    /** Copy constructor
     *  \param other Object to copy from
     */
    NWEM(const NWEM& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    NWEM& operator=(const NWEM& other);
    virtual void run();
	void Init(int TotalN, double IniInfect, int psc, int infectperiod, string nfn);
	void parameterSW(double SemiDegree, double Dp, double transmissibility);
	double GetInfect();
	double GetRecovered();
	double GetR0();
	int Run2End();
	void SetIsolated(double iso);
protected:
	int PSCount;
	Network nets[4];
	//Percolation_Sim sim[4];
	ChainBinomial_Sim sim[4];
	string NetFN;
	int md;	//4 smallworld
	double mInf;
	double Rec;
	double R0;
private:
	Network net;
	int N;
	double SD;
	double p;
	double ii;
	double t;
	double o;
	int ip;
	
	int CS;
};
}
#endif // __LDAS_NWEM_H
