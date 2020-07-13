#ifndef __LDAS_SEIR_H
#define __LDAS_SEIR_H

#include "model.h"

#include "exception.h"
#include "model.h"
#include "constant.h"

namespace ldas
{
///Novel coronavirus (2019-nCoV) transcription Model
/**
* This header file describes function type of the SEIR (Susceptible-Exposed-Infected-Recovered) model for 2019-nCoV 
* The differential equations is:
* dS/dt = -gamma*beta*I*S/N-gamma2*beta2*E*S/N
* dE/dt = gamma*beta*I*S/N+alpha*E+gamma2*beta2*E*S/N
* dI/dt = -alpha*E-gamma*I
* dR/dt = gamma*I
*/
class SEIR : public BaseModel
{
public:
    /** Default constructor */
    SEIR();
    //SEIR(const int alpha, int beta, int gamma, int beta2, int gamma2, long N, long E, long I);
    /** Default destructor */
    virtual ~SEIR();
    /** Copy constructor
     *  \param other Object to copy from
     */
    SEIR(const SEIR& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    SEIR& operator=(const SEIR& other);
    virtual void run();
    void parameter(const double alpha, double beta1, double gamma, double beta2, int r1, int r2);
    void init(const double N, double S, double E, double I, double R);
    void output() const;
    void update(const double S, double E, double I, double R);
    void updateEIR(const double E, double I, double R);
    void updateIR(const double I, double R);
    void updateI(const double I);
    const double* result() const;
protected:
	//gross population
	double m_N;
	//Exposed
	double m_E;
	//Infected
	double m_I;
	//Susceptible
	double m_S;
	//Recovered
	double m_R;
	double m_gamma, m_beta1, m_alpha, m_beta2;
	int  m_r1, m_r2;
	double m_state[4];
private:

};
}
#endif // __LDAS_SEIR_H
