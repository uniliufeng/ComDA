#include "gausshermiteparticlefilter.h"

using namespace ldas;
GaussHerimteParticleFilter::GaussHerimteParticleFilter():ParticleFilter(),GaussHermiteFilter()
{
    //ctor
}

GaussHerimteParticleFilter::~GaussHerimteParticleFilter()
{
    //dtor
}

GaussHerimteParticleFilter::GaussHerimteParticleFilter(const GaussHerimteParticleFilter& other)
{
    //copy ctor
}

GaussHerimteParticleFilter& GaussHerimteParticleFilter::operator=(const GaussHerimteParticleFilter& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

GaussHerimteParticleFilter::GaussHerimteParticleFilter(const unsigned int state, const unsigned int obs, const unsigned int en):ParticleFilter(state,obs,en),GaussHermiteFilter(state,obs)
{
    //ctor
}

void GaussHerimteParticleFilter::update(const mat& state_en_gauss, const mat& obs_forecast_en_gauss, const mat& obs, const mat& R,mat* P)
{
	for(int k=0;k<particle_num;k++)
	{
        mat Xf=zeros<mat>(state_en_gauss.n_rows,gauss_points_num);
        mat Yf=zeros<mat>(obs_forecast_en_gauss.n_rows,gauss_points_num);
        for(int i=0; i<gauss_points_num; i++)
        {
            Xf.col(i)=state_en_gauss.col(k*gauss_points_num+i);
            Yf.col(i)=obs_forecast_en_gauss.col(k*gauss_points_num+i);
        }
		GaussHermiteFilter::update(Xf,particles.col(k),Yf,obs,R,P[k]);
        particles.col(k)=GaussHermiteFilter::Xa.col(0);
	}
	//ParticleFilter::update(particles,particles,obs);
}

void GaussHerimteParticleFilter::update(const mat& state_en_gauss, const mat& Q, mat* P)
{
	mat particles_tmp=zeros<mat>(state_en_gauss.n_rows,particle_num);
	for(int k=0;k<particle_num;k++)
	{
        mat Xf=zeros<mat>(state_en_gauss.n_rows,gauss_points_num);
        for(int i=0; i<gauss_points_num; i++)
        {
            Xf.col(i)=state_en_gauss.col(k*gauss_points_num+i);
        }
		GaussHermiteFilter::update(Xf,Q,P[k]);
        particles_tmp.col(k)=GaussHermiteFilter::Xa.col(0);
	}
	ParticleFilter::update(particles_tmp);
}
