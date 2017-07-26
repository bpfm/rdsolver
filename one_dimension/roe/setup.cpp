using namespace std;

vertex setup_vertex(int n_points, int i, float dx, vertex new_vertex, int ic){

	if(ic==0){
		if(i==0){cout << "Using Sod Shock Tube 1" << endl;}
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
		return new_vertex;
	}else{
		if(i==0){cout << "Using Sine Wave" << endl;}
		new_vertex.set_x(float(i)*dx);

		new_vertex.set_mass_density(2.0+sin(float(i)/float(n_points)*(2.0*3.1415)));		// units kg/m^3
		new_vertex.set_velocity(100.0);													// units m/s
		new_vertex.set_pressure(100000.0*new_vertex.get_mass_density());												// units N/m^2

		new_vertex.setup_energy_density();
		new_vertex.con_to_prim();
		new_vertex.setup_f_variables();
		new_vertex.reset_du();
		return new_vertex;
	}
}