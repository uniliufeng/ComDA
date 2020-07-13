#ifndef __LDAS_NCSI_H
#define __LDAS_NCSI_H

#include "model.h"

#include "exception.h"
#include "model.h"
#include "constant.h"
#include "matrix.hpp"
#include <armadillo>

namespace ldas
{
///Novel coronavirus (2019-nCoV) transcription Model
/**
1、	No-Contact SI model
Reference: Liu F, Li X, Zhu GF. Using the contact network model and Metropolis-Hastings sampling to reconstruct the COVID-19 spread on the "Diamond Princess". Science Bulletin, 2020, 10.1016/j.scib.2020.04.043.
*/
class NCSI : public BaseModel
{
public:
    /** Default constructor */
    NCSI();
    /** Default destructor */
    virtual ~NCSI();
    /** Copy constructor
     *  \param other Object to copy from
     */
    NCSI(const NCSI& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    NCSI& operator=(const NCSI& other);
    virtual void run();
    void parameter(const double delta);
    void init(const double N, double S, double I, double O, int m, int step);
    void updateO(double o);
    void output() const;
    const double* result() const;
protected:
	//gross population
	double m_N;
	//Infected
	double m_I;
	//Susceptible
	double m_S;
	//isolated
	double m_O;
	//Delta
	double m_d;
	//survival days
	double m_m;
	//total step
	double m_t;
	//current step
	int cs;
	arma::vec SigmaI;
	arma::vec CumI;
	double m_state[2];
	
	void updateSigmaI();
	void updateCumI();
private:

};
}
#endif // __LDAS_NCSI_H
