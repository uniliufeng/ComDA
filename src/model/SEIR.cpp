#include "SEIR.h"

using namespace ldas;
SEIR::SEIR()
{
    //ctor
}

SEIR::~SEIR()
{
    //dtor
}

SEIR::SEIR(const SEIR& other)
{
    //copy ctor
}

SEIR& SEIR::operator=(const SEIR& other)
{
    if (this == &other) return *this; // handle self assignment

    return *this;
}

void SEIR::parameter(const double alpha, double beta1, double gamma, double beta2, int r1, int r2)
{
	m_alpha = alpha;
	m_beta1 = beta1;
	m_gamma = gamma;
	m_beta2 = beta2;
	m_r1 = r1;
	m_r2 = r2;
}

void SEIR::init(const double N, double S, double E, double I, double R)
{
	m_N = N;
	m_E = E;
	m_I = I;
	m_R = R;
	if (S==0)
	{m_S = m_N - m_I - m_E - m_R;}
	else
	{m_S = S;}
	
	
	m_state[0] = m_S;
	m_state[1] = m_E;
	m_state[2] = m_I;
	m_state[3] = m_R;
}

void SEIR::updateEIR(const double E, double I, double R)
{
	m_E = E;
	m_I = I;
	m_R = R;
	m_S = m_N - m_I - m_E - m_R;
	
	m_state[0] = m_S;
	m_state[1] = m_E;
	m_state[2] = m_I;
	m_state[3] = m_R;
}

void SEIR::update(const double S, double E, double I, double R)
{
	m_E = E;
	m_I = I;
	m_R = R;
	m_S = S;
	
	m_state[0] = m_S;
	m_state[1] = m_E;
	m_state[2] = m_I;
	m_state[3] = m_R;
}

void SEIR::updateIR(const double I, double R)
{
	m_I = I;
	m_R = R;
	m_S = m_N - m_I - m_E - m_R;
	
	m_state[0] = m_S;
	m_state[1] = m_E;
	m_state[2] = m_I;
	m_state[3] = m_R;
}

void SEIR::updateI(const double I)
{
	m_I = I;
	//m_S = m_N - m_I - m_E - m_R;
	
	//m_state[0] = m_S;
	//m_state[1] = m_E;
	m_state[2] = m_I;
	//m_state[3] = m_R;
}

void SEIR::run()
{
	//m_I = m_I+m_alpha*m_E-m_gamma*m_I;
	//double rS = m_N - m_I - m_E - m_R;
	double rS = m_S;
	double rE = m_E;
	double rI = m_I;
	double rR = m_R;
	m_S = rS-m_r1*m_beta1*rI*rS/m_N-m_r2*m_beta2*rE*rS/m_N;
	m_E = rE+m_r1*m_beta1*rI*rS/m_N-m_alpha*rE+m_r2*m_beta2*rE*rS/m_N;
	m_I = rI+m_alpha*rE-m_gamma*rI;
	m_R = rR+m_gamma*rI;
	
	m_state[0] = m_S;
	m_state[1] = m_E;
	m_state[2] = m_I;
	m_state[3] = m_R;

}

void SEIR::output() const
{
    std::cout<<m_S<<" "<<m_E<<" "<<m_I<<" "<<m_R<<std::endl;
}

const double* SEIR::result() const
{
    return m_state;
}
