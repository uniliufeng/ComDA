#ifndef __LDAS_LPJ_H
#define __LDAS_LPJ_H

#include "model.h"
#include "exception.h"
namespace ldas
{
extern "C"
{
	#include "lpj/include/lpj.h"
	#include "lpj/include/grass.h"
	#include "lpj/include/tree.h"
}
#define NTYPES 2
const Fscanpftparfcn scanfcn[NTYPES]={fscanpft_grass,fscanpft_tree};
class LPJ : public BaseModel
{
public:
    /** Default constructor */
    LPJ();
    /** Default destructor */
    virtual ~LPJ();
    /** Copy constructor
     *  \param other Object to copy from
     */
    LPJ(const LPJ& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    LPJ& operator=(const LPJ& other);
    virtual void run();
    virtual void config(const std::string);
protected:
    Config m_config;
    Climate *climate;
    Pftpar *pftpar;
    Soilpar *soilpar;
    Cell *grid;
    FILE **output;

    int npft,nsoil,year;
private:
};
}
#endif // __LDAS_LPJ_H
