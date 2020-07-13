#include "NWEM.h"

using namespace ldas;
NWEM::NWEM(): net("NWEMnetwork", Network::Undirected)
{
    //ctor
}

NWEM::~NWEM()
{
    //dtor
}

NWEM::NWEM(const NWEM& other)
{
    //copy ctor
}

NWEM& NWEM::operator=(const NWEM& other)
{
    if (this == &other) return *this; // handle self assignment

    return *this;
}

void NWEM::parameterSW(double SemiDegree, double Dp, double transmissibility)
{
	SD=SemiDegree;
	p=Dp;
	t=transmissibility;
	md=4;
}

void NWEM::Init(int TotalN, double IniInfect, int psc, int infectperiod, string nfn)
{
	CS = 0;
	N=TotalN;
	ii=IniInfect;
	NetFN = nfn;
	PSCount=psc;
	ip=infectperiod;
	o=0;
	
	char sep = ',';
	if(net.size()==0)
	{
		//net.read_edgelist(NetFN,sep);
		//if(net.size()==0)
		{
			net.populate(N);

			if(md==4)
			{
				net.small_world(N, SD, p);
			}
			else
			{}
			net.write_edgelist(NetFN,Network::NodeIDs,sep);
		}
/*
		Network net("Nnetwork", Network::Undirected);
		net.populate(4000);
		net.small_world(4000, 12, 0.2675);
		ChainBinomial_Sim cs(&net,16,0.007);
		cout<<" R0:"<<cs.expected_R0();
*/
	}
	
	string fn;
	for(int i=0;i<PSCount;i++)
	{
		nets[i].InitNetwork("net_"+to_string(i), Network::Undirected);
		fn = "./EpiFire/net"+to_string(N)+"_"+to_string(i)+".txt";
		//if (access(fn.c_str(),0)!=0)
		{
			string cmd = "cp "+NetFN+" "+fn;
			system(cmd.c_str());
		}
		nets[i].read_edgelist(fn,',');
		sim[i].Simulator::set_network(&nets[i]);
		sim[i].set_infectious_period(ip);
		sim[i].set_transmissibility(t);
		sim[i].define_time_dist();
		sim[i].rand_infect(ii);
	}
}

void NWEM::run()
{
	CS++;
	double cinf = 0;
	double inf = 0;
	int ee;
#pragma omp parallel for
	for(int i=0;i<PSCount;i++)
	{
		inf = sim[i].count_infected()-o;
		if (inf<1)
			continue;
		//if(o>0)
		{
			sim[i].reset();
			sim[i].rand_infect(inf);
		}
		
		sim[i].set_transmissibility(t);
		sim[i].step_simulation();
		cinf = cinf + inf;
		//crec = crec + sim[i].epidemic_size();
	}
	mInf = cinf/PSCount;
	//Rec = crec/PSCount;
	o=0;
}

double NWEM::GetInfect()
{
	return mInf;
}

double NWEM::GetRecovered()
{
	double crec = 0;
	for(int i=0;i<PSCount;i++)
	{
		crec = crec + sim[i].epidemic_size();
	}
	Rec = crec/PSCount;
	return Rec;
}

void NWEM::SetIsolated(double iso)
{
	o = iso;
}

double NWEM::GetR0()
{
	R0=0;
	for(int i=0;i<PSCount;i++)
	{
		R0 = R0 + sim[i].expected_R0();
	}
	R0 = R0/PSCount;
	return R0;
}

int NWEM::Run2End()
{
return 0;
}

