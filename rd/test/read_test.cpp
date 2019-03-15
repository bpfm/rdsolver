#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>


int main(){
	std::ifstream POSITIONS_FILE, TRIANGLES_FILE;
	std::string INFO;
	int N_POINTS, N_TRIANG, N_VERT;
  
	POSITIONS_FILE.open("rand_points.txt");
	TRIANGLES_FILE.open("rand_triangles.txt");

	std::getline(POSITIONS_FILE,INFO);

	POSITIONS_FILE >> N_POINTS;
	TRIANGLES_FILE >> N_TRIANG;

	double X[N_POINTS],Y[N_POINTS];
	int VERT0[N_TRIANG],VERT1[N_TRIANG],VERT2[N_TRIANG];

	for (int i = 0; i < 10; ++i){
		POSITIONS_FILE >> X[i] >> Y[i];
		// std::cout << X[i] << "\t" << Y[i] << std::endl;
	}

	for (int j = 0; j < N_TRIANG; ++j){
		TRIANGLES_FILE >> N_VERT >> VERT0[j] >> VERT1[j] >> VERT2[j];
		std::cout << VERT0[j] << "\t" << VERT1[j] << "\t" << VERT2[j] << std::endl;;
	}

	POSITIONS_FILE.close();
	TRIANGLES_FILE.close();
  
  return 0;
}