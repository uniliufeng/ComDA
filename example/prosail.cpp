/**
\file The ProSail Model example
*/
#include "prosail.h"
#include "radiancetransfermodel.h"
#include <iostream>
int main()
{
	//数据来源于FORTRAN源程序
    ldas::CanopyParameter cp;
    cp.observer_zenith=0;
    cp.solar_zenith=30;
    cp.azimuth=0;

    cp.chlorophyll_content=30;
    cp.carotenoid_content=10;
    cp.brownpigment_content=0;
    cp.cw=0.015;
    cp.cm=0.009;
    cp.structure_coefficient=1.2;
    cp.lai=2;
    cp.angle=50;
    cp.soil_coefficient=1;
    cp.skyl=70;
    cp.hot_spot=0.2;

    ldas::ProSail sail(cp);
    sail.run();
    sail.output();
    return 0;
}
