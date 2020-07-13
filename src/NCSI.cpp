#include "NCSI.h"

using namespace ldas;
NCSI::NCSI()
{
    //ctor
}

NCSI::~NCSI()
{
    //dtor
}

NCSI::NCSI(const NCSI& other)
{
    //copy ctor
}

NCSI& NCSI::operator=(const NCSI& other)
{
    if (this == &other) return *this; // handle self assignment

    return *this;
}

void NCSI::parameter(const double delta)
{
	m_d = delta;
}

void NCSI::init(double N, double S, double I, double O, int m, int step)
{
	m_m=m;
	m_t=step;
	m_N = N;
	m_I = I;
	m_O = O;
	if (S==0)
	{m_S = m_N - m_I - m_O;}
	else
	{m_S = S;}
	
	SigmaI.set_size(m_t);
	SigmaI.zeros();
	CumI.set_size(m_t);
	CumI.zeros();
	m_state[0] = m_S;
	m_state[1] = m_I;
	cs=0;
}

void NCSI::run()
{
	double rS = m_S;
	double rI = m_I;
	double rO = m_O;
	double sumI=0.0;
	updateSigmaI();
	//updateCumI();
	int sd=0;
	if(cs>=m_m)
	{
		sd = cs-m_m+1;
	}
	for(int i=sd;i<=cs;i++)
	{
		sumI=sumI+SigmaI(i);
	}
	
	m_S = rS-m_d*(sumI)*rS;
	m_I = rI+m_d*(sumI)*rS-rO;
	
	m_state[0] = m_S;
	m_state[1] = m_I;
	
	cs++;
}

void NCSI::updateO(double o)
{
	m_O = o;
}

void NCSI::updateSigmaI()
{
	if(SigmaI.n_rows<2)
	{
		cout<<SigmaI.n_rows<<endl;
		return;
	}
	if(cs==0)
	{
		CumI(cs)=m_I;
	}
	else
	{
		CumI(cs)=CumI(cs-1)+m_I;
	}
	
	if(cs>=m_m)
	{
		SigmaI(cs)=CumI(cs)-CumI(cs-m_m);
	}
	else
	{
		SigmaI(cs)=CumI(cs);
	}
	/*if(cs==m_t-1)
	{
		arma::mat SC = join_rows(SigmaI, CumI);
		cout<<"***";
		SC.print();
	}
	*/
}

void NCSI::updateCumI()
{
	if(CumI.n_rows<2)
	{
		cout<<CumI.n_rows<<endl;
		return;
	}
	if(cs==0)
	{
		CumI(cs)=m_I;
	}
	else
	{
		CumI(cs)=CumI(cs-1)+m_I;
	}
	if(cs==m_t-1)
	{
		cout<<"***";
		CumI.print();
	}
}

void NCSI::output() const
{
    std::cout<<m_S<<" "<<m_I<<" "<<std::endl;
}

const double* NCSI::result() const
{
    return m_state;
}
