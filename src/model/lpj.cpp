#include "lpj.h"
using namespace ldas;
LPJ::LPJ()
{
    //ctor
}

LPJ::~LPJ()
{
    //dtor
}

LPJ::LPJ(const LPJ& other)
{
    //copy ctor
}

LPJ& LPJ::operator=(const LPJ& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}
void LPJ::run()
{
	output=fopenoutput(m_config);
    year=iterate(output,grid,climate,pftpar,npft,NTYPES,m_config);
    if(output[GRID]!=NULL)
        writecoords(output[GRID],grid,m_config.ngridcell);
    fcloseoutput(output,m_config.n_out);
    printf("Simulation ended.\n");
    /* free memory */
    freeclimate(climate);
    freegrid(grid,m_config.ngridcell);
}
void LPJ::config(const std::string fn)
{
    int argc=2;
    char** argv;
    argv[0]="not_used";
    //argv[1]=fn.c_str();
    strcpy(argv[1], fn.c_str());
    if(fscanconfig(&m_config,&argc,&argv))
        throw Exception("LPJ","config()","Error reading config");
    if((npft=fscanpftpar(&pftpar,m_config.pftpar_filename,scanfcn,NTYPES))==0)
        throw Exception("LPJ","config()","Error reading pftpar file");
    if((nsoil=fscansoilpar(&soilpar,m_config.soilpar_filename))==0)
        throw Exception("LPJ","config()","Error reading soilpar file");
    if(isreadrestart(m_config))
        printf("Starting from restart file '%s'.\n",m_config.restart_filename);
    if((grid=newgrid(&m_config,pftpar,npft,soilpar,nsoil))==NULL)
        throw Exception("LPJ","config()","Error initializing grid");
    if((climate=initclimate(m_config))==NULL)
        throw Exception("LPJ","config()","Error initializing climate");
    /* Initialize random seed */
#ifdef USE_RAND48
    srand48((m_config.seed==RANDOM_SEED) ? std::time(NULL) : m_config.seed);
#else
    setseed(m_config.seed==RANDOM_SEED ? std::time(NULL) : m_config.seed);
#endif
    printf("Simulation begins...\n");

}
