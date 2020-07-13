#ifndef __LDAS_SEIQ_H
#define __LDAS_SEIQ_H

#include "model.h"

#include "exception.h"
#include "model.h"
#include "constant.h"

namespace ldas
{
///Novel coronavirus (COVID-19) transcription Model
/**
* This header file describes function type of the SEIQP (Susceptible-Exposed-Infected-Quarantined-insuscePtible) model for COVID-19
* Reference: Peng et al., 2020. Epidemic analysis of COVID-19 in China by dynamical modeling.
* The differential equations (a Generalized SEIR model) is:
* dS/dt = -beta*I*S/N-alpha*S
* dE/dt = beta*I*S/N-gamma*E
* dI/dt = gamma*E-delta*I
* dQ/dt = delta*I
* dP/dt = alpha*S
*/
class SEIQ : public BaseModel
{
public:
    /** Default constructor */
    SEIQ();
    //SEIQ(const int alpha, int beta, int gamma, int beta2, int gamma2, long N, long E, long I);
    /** Default destructor */
    virtual ~SEIQ();
    /** Copy constructor
     *  \param other Object to copy from
     */
    SEIQ(const SEIQ& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    SEIQ& operator=(const SEIQ& other);
    virtual void run();
    void parameter(const double alpha, double beta, double gamma, double delta);
    void init(const double N, double S, double E, double I, double Q, double P);
    void output() const;
    void update(const double S, double E, double I, double Q, double P);
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
	//Quarantined
	double m_Q;
	//insuscePtible
	double m_P;
	double m_gamma, m_beta, m_alpha, m_delta;
	double m_state[5];
private:

};
}
#endif // __LDAS_SEIQ_H
