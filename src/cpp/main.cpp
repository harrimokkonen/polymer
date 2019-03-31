#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "polymer.hpp"

void WriteCoordinates(std::ofstream &fp, std::vector<float> x,
											std::vector<float> y,
											std::vector<float> z)
{
	for(auto val : x)
		fp << val << " ";
	for(auto val : y)
		fp << val << " ";
	for(auto val : z)
		fp << val << " ";
	fp << std::endl;
}

int main(int argc, char *argv[])
{
	int time = 1000;
	int N = 8;	

	std::ofstream myfile;
  myfile.open("coords.txt");

	Polymer pol(SelfAvoiding, N, 15.0, 15.0, 2.0, 1.0, 1.0,
							0.01, 0.7, 0.0, 0.005);
	pol.ConfigureSAW(true);

	for (int t = 0; t < time; ++t)
	{
		WriteCoordinates(myfile, pol.x, pol.y, pol.z);
		pol.UpdateBBK();
		std::cout << pol.E_tot() << std::endl;		
	}

  myfile.close();
	return 0;
}