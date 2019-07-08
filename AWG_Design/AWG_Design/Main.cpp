#ifndef ATTACH_H
#include "Attach.h"
#endif

int main() 
{
	testing::compute_AWG_Params(); 

	//testing::REDFINCH_AWG_Params(3); 

	std::cout << "Press enter to close\n"; 
	std::cin.get(); 

	return 0; 
}