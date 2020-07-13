#include "SEIQ.h"

using namespace ldas;
SEIQ::SEIQ()
{
    //ctor
}

SEIQ::~SEIQ()
{
    //dtor
}

SEIQ::SEIQ(const SEIQ& other)
{
    //copy ctor
}

SEIQ& SEIQ::operator=(const SEIQ& other)
{
    if (this == &other) return *this; // handle self assignment

    return *this;
}

void SEIQ::parameter(const double alpha, double beta, double gamma, double delta)
{
	m_alpha = alpha;
	m_beta = beta;
	m_gamma = gamma;
	m_delta = delta;
}

void SEIQ::init(const double N, double S, double E, double I, double Q, double P)
{
	m_N = N;
	m_E = E;
	m_I = I;
	m_Q = Q;
	m_P = P;
	if (S==0)
	{m_S = m_N - m_I - m_E - m_Q;}
	else
	{m_S = S;}
	
	
	m_state[0] = m_S;
	m_state[1] = m_E;
	m_state[2] = m_I;
	m_state[3] = m_Q;
	m_state[4] = m_P;
}

void SEIQ::update(const double S, double E, double I, double Q, double P)
{
	m_E = E;
	m_I = I;
	m_Q = Q;
	m_S = S;
	m_P = P;
	
	m_state[0] = m_S;
	m_state[1] = m_E;
	m_state[2] = m_I;
	m_state[3] = m_Q;
	m_state[4] = m_P;
}

void SEIQ::updateI(const double I)
{
	m_I = I;
	//m_S = m_N - m_I - m_E - m_R;
	
	//m_state[0] = m_S;
	//m_state[1] = m_E;
	m_state[2] = m_I;
	//m_state[3] = m_R;
}

void SEIQ::run()
{
	double rS = m_S;
	double rE = m_E;
	double rI = m_I;
	double rQ = m_Q;
	double rP = m_P;
	m_S = rS-m_beta*rI*rS/m_N-m_alpha*m_S;
	m_E = rE+m_beta*rI*rS/m_N-m_gamma*rE;
	m_I = rI+m_gamma*rE-m_delta*rI;
	m_Q = rQ+m_delta*rI;
	m_P = m_P + m_alpha*m_S;
	
	
	m_state[0] = m_S;
	m_state[1] = m_E;
	m_state[2] = m_I;
	m_state[3] = m_Q;
	m_state[4] = m_P;
}

void SEIQ::output() const
{
    std::cout<<m_S<<" "<<m_E<<" "<<m_I<<" "<<m_Q<<" "<<m_P<<std::endl;
}

const double* SEIQ::result() const
{
    return m_state;
}
