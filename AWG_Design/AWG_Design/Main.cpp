#ifndef ATTACH_H
#include "Attach.h"
#endif

int main() 
{
	testing::compute_AWG_Params(); 

	std::cout << "Press enter to close\n"; 
	std::cin.get(); 

	return 0; 
}