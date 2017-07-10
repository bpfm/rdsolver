using namespace std;

void setup(int n_points, int n_cells, vector<vertex> &points, int ic){
	int i;
	vertex new_vertex;					// new_vertex = temporary vertex to be added to vector of vertices
	cell new_cell;						// new_cell = temporary cell to be added to vector of cells

	if(ic==0){
		for (i=0;i<n_points;i++){
			new_vertex.set_x(float(i)*dx);

			if(i<n_points/2){
				new_vertex.set_mass_density(1.0);						// units kg/m^3
				new_vertex.set_velocity(0.0);							// units m/s
				new_vertex.set_pressure(100000.0);						// units N/m^2
			}else{
				new_vertex.set_mass_density(0.125);						// units kg/m^3
				new_vertex.set_velocity(0.0);							// units m/s
				new_vertex.set_pressure(1000.0);						// units N/m^2
			}
			new_vertex.setup_energy_density();
			new_vertex.con_to_prim();
			new_vertex.setup_f_variables();
			new_vertex.reset_du();
			points.push_back(new_vertex);			// add new vertex to vector of all vertices
		}
	}else if(ic==1){
		for (i=0;i<n_points;i++){
			new_vertex.set_x(float(i)*dx);

			if(i<n_points/2){
				new_vertex.set_mass_density(sin(float(i)/(2.0*3.1415));						// units kg/m^3
				new_vertex.set_velocity(0.1);							// units m/s
				new_vertex.set_pressure(100000.0);						// units N/m^2
			}
			new_vertex.setup_energy_density();
			new_vertex.con_to_prim();
			new_vertex.setup_f_variables();
			new_vertex.reset_du();
			points.push_back(new_vertex);			// add new vertex to vector of all vertices
		}
	}
}