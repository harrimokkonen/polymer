#include<iostream>
#include<string>

#include "polymer.hpp"

int main(int argc, char *argv[])
{
	int time = 0;
	if(argc = 2) {
		time = std::stoi(argv[1]);
	}
	Polymer pol(Harmonic, 128, 15.0, 15.0, 2.0, 1.0, 1.0,
										 0.01, 0.7, 0.0, 0.005);

	pol.ConfigureSAW(true);
	
	for(int t = 0;  t < time; ++t) {
		pol.UpdateBBK();
		std::cout << pol.E_tot() << std::endl;
	}

}