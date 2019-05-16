#include "particlefilter.h"

using namespace ldas;

ParticleFilter::ParticleFilter():particle_num(0),m_threshold(0)
{

}

ParticleFilter::~ParticleFilter()
{

}

ParticleFilter::ParticleFilter(const unsigned int StateNum, const unsigned int ObsNum, const unsigned int EnNum):Filter(StateNum,ObsNum),particle_num(EnNum),m_threshold(0)
{
    init();
}

void ParticleFilter::particleNum(const unsigned int n)
{
    particle_num=n;
    init();
}

unsigned int ParticleFilter::particleNum() const
{
    return particle_num;
}

void ParticleFilter::threshold(const double th)
{
    m_threshold=th;
}

double ParticleFilter::threshold() const
{
    return m_threshold;
}

mat ParticleFilter::getXa() const
{
    return Xa;
}
mat ParticleFilter::getXaEn() const
{
    return particles;
}

void ParticleFilter::init()
{
    weights=1.0/particle_num*ones<mat>(1,particle_num);
    if (m_threshold==0 && particle_num>0)
        m_threshold=2*particle_num/3.0;
}

void ParticleFilter::update(const mat& state_en, const mat& Yf, const mat& obs)
{
    mat obs_en=zeros<mat>(observe_num,particle_num);
    for(int i=0; i<particle_num; i++)
        obs_en.col(i)=obs.col(0);
    mat llh=likelihood(obs_en,Yf)+1e-99;
    weights=weights%llh;
    weights=weights/sum(sum(weights));

    // calculate effective particle set size
    double Neff = 1.0 / sum(sum(pow(weights, 2)));
    if (floor(Neff) <= m_threshold)
    {
        particles=residual_resample(weights,state_en);
        weights = 1.0/particle_num*ones<mat>(1,particle_num);
    }
    else
    {
        particles=state_en;
        //std::cerr<<"Neff:"<<Neff<<", threshold:"<<m_threshold<<std::endl;
    }
    mat weights_en=zeros<mat>(state_num,particle_num);
    for(int i=0;i<state_num;i++)
		weights_en.row(i)=weights.row(0);
    Xa=sum(weights_en%particles,1);
}

void ParticleFilter::update(const mat& state_en)
{
    particles=state_en;
    //Xa=mean(particles,1);
    mat weights_en=zeros<mat>(state_num,particle_num);
    for(int i=0;i<state_num;i++)
		weights_en.row(i)=weights.row(0);
    Xa=sum(weights_en%particles,1);
}

mat ParticleFilter::likelihood(const mat &obs, const mat &state)
{
	//assumption: state=observe forcast?
    mat X=zeros<mat>(obs.n_rows,obs.n_cols);
    for(int i=0;i<obs.n_rows;i++)
		X.row(i)=obs.row(i)-state.row(i);

    mat mu = zeros<mat>(obs.n_rows,1);
    //mu=mean(X,1);
    mat coveriance=zeros<mat>(obs.n_rows,obs.n_rows);
    mat cov1=cov(trans(X))+1e-99;
    for(int i=0; i<obs.n_rows; i++)
    {
        //coveriance(i,i)=(cov1(i,i)<1e-4)?1e-4:cov1(i,i);
        /** FIXME **/
        coveriance(i,i)=1;  //why?
    }
    //gauss assumption
    return gauseval(obs.n_rows, mu, coveriance, X);

    /*mat S=trans(chol(coveriance));
    mat foo=solve(S,X);
    mat T=abs(prod(diagvec(S)));
    return exp(-0.5 * sum(foo%foo))/std::sqrt(2*PI_NUMBER)/T(0,0);*/
}

mat ParticleFilter::gauseval(const int dim, const mat& mu, const mat& cov, const mat& X)
{
    mat XX, S, foo;
    double normfact = std::pow(2 *PI_NUMBER, dim / 2.0);
    XX.set_size(X.n_rows,X.n_cols);
    for(int i=0; i<X.n_cols; i++)
        XX.col(i)=X.col(i)-mu.col(0);
    S = trans(chol(cov));
    foo = solve(S,XX);
    mat T=prod(diagvec(S));
    double t=std::abs(T(0,0));
    mat res=exp(-0.5 * sum(foo%foo)) / (normfact*t);
    return res;
}

mat ParticleFilter::residual_resample(const mat& weights,const mat& particles)
{
    ivec outIndex(particles.n_cols);	// setup output index buffer
    //FIXME:change the random seed
    std::srand(time(NULL));
    mat resample_particles=zeros<mat>(particles.n_rows,particles.n_cols);

    //=== RESIDUAL RESAMPLING  ==========================================================
    imat N_kind=zeros<imat>(1,particles.n_cols);

    // first integer part
    mat weights_res = particles.n_cols * weights;
    //N_kind = fix(weights_res);
    for(int i=0; i<particles.n_cols; i++)
    {
        if(weights_res(i)>0)
            N_kind(i)=floor(weights_res(i));
        else
            N_kind(i)=ceil(weights_res(i));
    }
    // residual number of particles to sample
    int N_res = particles.n_cols - sum(sum(N_kind));

    if(N_res)
    {
        weights_res = (weights_res - N_kind) / N_res;
        mat cumDist = cumsum(weights_res,1);

        // generate N_res ordered random variables uniformly distributed in [0,1]
        vec divec=zeros<vec>(N_res);
        for(int i=0; i<N_res; i++)
        {
            divec(i) = 1.0/(N_res - i);
        }
        vec randvec = randu<vec>(N_res);
        //u = reverse();
        vec powervec=zeros<vec>(N_res);
        vec u1=zeros<vec>(N_res);
        vec u=zeros<vec>(N_res);
        //powervec = power(randvec,divec);
        for(int i=0; i<N_res; i++)
        {
            powervec(i)=pow(randvec(i),divec(i));
        }
        //cumprod(powervec)
        u1(0)=powervec(0);
        for(int i=1; i<N_res; i++)
        {
            u1(i)=u1(i-1)*powervec(i);
        }
        //fliplr
        for(int i=0; i<N_res; i++)
        {
            u(i)=u1(N_res-1-i);
        }
        
        int j=0;
        for(int i=0; i<N_res; i++)
        {
            while(u(i)>cumDist(j))
            {
                j = j+1;
            }
            N_kind(j) = N_kind(j) + 1;
        }
    }
    
    //COPY RESAMPLED TRAJECTORIES
    int index = 0;
    for(int i=0; i<particles.n_cols; i++)
    {
        if(N_kind(i)>0)
        {
            for(int j=index; j<=index+N_kind(i)-1; j++)
            {
                outIndex(j) = i;
            }
        }
        index += N_kind(i);
    }
    //std::cout<<weights;
    //std::cout<<N_kind;
    for(int i=0; i<particles.n_cols; i++)
    {
        resample_particles.col(i)=particles.col(outIndex(i));
    }
    return resample_particles;
}
