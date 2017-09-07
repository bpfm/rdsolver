using namespace std;

centre setup_centre(int n_points, int i, float dx, centre new_centre){

	if(IC == 0){
		if(i==0){cout << "Using Sod Shock Tube 1" << endl;}
		new_centre.set_x((float(i)+0.5)*dx);

		if(i>0.3*n_points and i<0.7*n_points){
			new_centre.set_mass_density(1.0);				// units kg/m^3
			new_centre.set_velocity(0.0);					// units m/s
			new_centre.set_pressure(500.0);					// units N/m^2
		}else{
			new_centre.set_mass_density(0.125);				// units kg/m^3
			new_centre.set_velocity(0.0);					// units m/s
			new_centre.set_pressure(100.0);					// units N/m^2
		}
		new_centre.setup_energy_density();
		new_centre.con_to_prim();
		new_centre.setup_f_variables();
		new_centre.reset_du();
		return new_centre;

	}else if(IC == 1){

                double vel_phase = 2.0;
                double wavelength = 50.0;
                double freq,period;
		
		if(i==0){
                        cout << "Using Sine Wave" << endl;
                        freq = vel_phase/wavelength;
                        period = 1.0/freq;
                        //cout << "Wave properties:\t frequency =\t" << freq << " per s\tperiod =\t" << period << " s\twavelength =\t" << wavelength << " m\tphase velocity =\t" << vel_phase << " m/s" << endl;
                }
		new_centre.set_x(float(i)*dx);

		new_centre.set_mass_density(1.0);	// units kg/m^3
		new_centre.set_velocity(sin(float(i)/float(n_points)*(2.0*3.1415)));	                                        // units m/s
		new_centre.set_pressure(100.0);		        // units N/m^2

		new_centre.setup_energy_density();
		new_centre.con_to_prim();
		new_centre.setup_f_variables();
		new_centre.reset_du();

                return new_centre;


        }else{

        if(i==0){cout << "Using 1D Sedov Blast Wave" << endl;}
                double x = (float(i) + 0.5) * dx;
                new_centre.set_x(x);
    
                new_centre.set_mass_density(1);        // units kg/m^3
                new_centre.set_velocity(0);            // units m/s
                new_centre.set_pressure(1);            // units N/m^2
    
                if(x > (50/2 - dx) and x < (50/2 + dx)){ new_centre.set_pressure(100000/dx); }

                new_centre.setup_energy_density();
                new_centre.con_to_prim();
                new_centre.setup_f_variables();
                new_centre.reset_du();
                return new_centre;
        }
}