#ifndef __LDAS_FILTER_H
#define __LDAS_FILTER_H

#include <armadillo>

using namespace arma;
namespace ldas
{
class Filter
{
public:
    /** Default constructor */
    Filter();
    /** Default destructor */
    virtual ~Filter();
    Filter(const unsigned int state, const unsigned int obs);
    /** Copy constructor
     *  \param other Object to copy from
     */
    Filter(const Filter& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    Filter& operator=(const Filter& other);
    /** set state variable number
    * \param n state number
    */
    void stateNum(const unsigned int n);
    unsigned int stateNum() const;

    void observeNum(const unsigned int n);
    unsigned int observeNum() const;
protected:
    unsigned int state_num;
    unsigned int observe_num;
private:
};
}
#endif // __LDAS_FILTER_H
