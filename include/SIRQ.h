#ifndef __LDAS_SIRQ_H
#define __LDAS_SIRQ_H

#include "model.h"

#include "exception.h"
#include "model.h"
#include "constant.h"

namespace ldas
{
///Novel coronavirus (COVID-19) transcription Model
/**
* This header file describes function type of the SIRQ (Susceptible-Infected-Recovered with Quarantined) model for COVID-19
* The differential equations (a Generalized SEIR model) is:
* dS/dt = -beta*I*S/N-alpha*S
* dI/dt = beta*I*S/N
* dR/dt = gamma*I
*/
class SIRQ : public BaseModel
{
public:
    /** Default constructor */
    SIRQ();
    //SIRQ(const int alpha, int beta, int gamma, int beta2, int gamma2, long N, long E, long I);
    /** Default destructor */
    virtual ~SIRQ();
    /** Copy constructor
     *  \param other Object to copy from
     */
    SIRQ(const SIRQ& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    SIRQ& operator=(const SIRQ& other);
    virtual void run();
    void parameter(const double alpha, double beta, double gamma);
    void init(const double N, double S, double I, double R);
    void output() const;
    void update(const double S, double I, double R);
    void updateI(const double I);
    const double* result() const;
protected:
	//gross population
	double m_N;
	//Infected
	double m_I;
	//Susceptible
	double m_S;
	//Recovered
	double m_R;
	double m_gamma, m_beta, m_alpha;
	double m_state[3];
private:

};
}
#endif // __LDAS_SIRQ_H
