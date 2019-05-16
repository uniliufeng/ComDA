/**
\file The LPJ Model example
*/
#include "lpj.h"
#include <iostream>
int main() {
	//char* control_file = "global_param.Arou";
	ldas::LPJ lpj;
	lpj.config("lpj.conf");
	lpj.run();
	return 0;
}
