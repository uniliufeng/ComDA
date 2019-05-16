/**
\file The SiB2 Model example
*/
#include "sib2model.h"
using namespace std;
using namespace ldas;

int main()
{
	SiB2Config cfg;
	cfg.vegetation_path="data1";
	cfg.forcing_path="data2";
	cfg.output_path="./";
    SiB2Model sib2;
    sib2.config(cfg);
    sib2.run();
    return 0;
}
