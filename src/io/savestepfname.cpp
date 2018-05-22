#include "includes.h"


std::string kolmofname(const std::string base, const Real kf)
{

	return base+"_kf_"+std::to_string((Int)kf);
}

std::string step_fname(const std::string prefix, const std::string suffix,
		       const Int N, const Int stepNo)
{

	
	return prefix+"_N_"+std::to_string(N)+"_stepNo_"+std::to_string(stepNo)+suffix;
}

std::string stats_fname(const std::string prefix, const std::string suffix,const Int N)
{

	
	return prefix+"_N_"+std::to_string(N)+suffix;
}
