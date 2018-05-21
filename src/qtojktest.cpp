//#include "includes.h"
#include "qjkmap.h"

Int main()
{
	/*
	for (Int q=0;q<9;q++)
	{
		Intpair jkk = qtojk(q);
		std::cout << "(" << jkk.first << "," << jkk.second << ")" << "\t";
	}
  
  

	std::cout << std::endl;
	std::cout << "-----------------------------" << std::endl;
	*/
	
	Int jk[2];
	for (Int q=0;q<9;q++)
	{
		qtojk(jk,q);
		std::cout << "(" << jk[0] << "," << jk[1] << ")" << "\t";
	}

	return 0;
	
}
