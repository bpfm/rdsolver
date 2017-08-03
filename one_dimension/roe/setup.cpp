using namespace std;

centre setup_centre(int n_points, int i, float dx, centre new_centre, int ic){

	if(ic==0){
		if(i==0){cout << "Using Sod Shock Tube 1" << endl;}
		new_centre.set_x((float(i)+0.5)*dx);

		if(i>n_points/3 and i<2*n_points/3){
			new_centre.set_mass_density(1.0);					// units kg/m^3
			new_centre.set_velocity(0.0);						// units m/s
			new_centre.set_pressure(100000.0);					// units N/m^2
		}else{
			new_centre.set_mass_density(0.125);					// units kg/m^3
			new_centre.set_velocity(0.0);						// units m/s
			new_centre.set_pressure(10000.0);					// units N/m^2
		}
		new_centre.setup_energy_density();
		new_centre.con_to_prim();
		new_centre.setup_f_variables();
		new_centre.reset_du();
		return new_centre;
	}else{
		if(i==0){cout << "Using Sine Wave" << endl;}
		new_centre.set_x(float(i)*dx);

		new_centre.set_mass_density(2.0+sin(float(i)/float(n_points)*(2.0*3.1415)));	// units kg/m^3
		new_centre.set_velocity(100.0);							// units m/s
		new_centre.set_pressure(100000.0*new_centre.get_mass_density());		// units N/m^2

		new_centre.setup_energy_density();
		new_centre.con_to_prim();
		new_centre.setup_f_variables();
		new_centre.reset_du();
		return new_centre;
	}
}