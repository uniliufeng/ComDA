/**
\file The SiB2 Model example
*/
#include "memls.h"
#include "radiancetransfermodel.h"
#include <iostream>
int main()
{
    ldas::SnowParameter sp;
    sp.theta=10;
    sp.frequency=10.87;
    sp.col=1;
    sp.temperature_sky=0;
    sp.temperature_ground=273;
    sp.depth.setDim(1,sp.col);
    sp.depth.set(1,sp.col,66);
    sp.density.setDim(1,sp.col);
    sp.density.set(1,sp.col,0.3);
    sp.temperature.setDim(1,sp.col);
    sp.temperature.set(1,sp.col,273);
    sp.wetness.setDim(1,sp.col);
    sp.wetness.set(1,sp.col,0.06);
    sp.pci.setDim(1,sp.col);
    sp.pci.set(1,sp.col,0.12);
    ldas::MEMLS memls(sp);
    memls.run();
    //double tbv,tbh;
    //memls.run(tbv,tbh,sp.frequency,sp.theta,sp.temperature_sky,sp.temperature_ground);
    //std::cout<<tbv<<std::endl;
    std::cout<<memls.brightness()<<std::endl;
    return 0;
}
