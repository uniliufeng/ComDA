#ifndef __LDAS_SIR_H
#define __LDAS_SIR_H

#include "model.h"

#include "exception.h"
#include "model.h"
#include "constant.h"

namespace ldas
{
///Novel coronavirus (COVID-19) transcription Model
/**
* This header file describes function type of the SIR (Susceptible-Infected-Recovered) model for COVID-19

* The differential equations (a Generalized SEIR model) is:
* dS/dt = -beta*I*S/N
* dI/dt = beta*I*S/N
* dR/dt = gamma*I
*/
class SIR : public BaseModel
{
public:
    /** Default constructor */
    SIR();
    //SIR(const int alpha, int beta, int gamma, int beta2, int gamma2, long N, long E, long I);
    /** Default destructor */
    virtual ~SIR();
    /** Copy constructor
     *  \param other Object to copy from
     */
    SIR(const SIR& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    SIR& operator=(const SIR& other);
    virtual void run();
    void parameter(const double beta, double gamma);
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
	double m_gamma, m_beta;
	double m_state[3];
private:

};
}
#endif // __LDAS_SIR_H
