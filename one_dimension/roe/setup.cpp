using namespace std;

centre setup_centre(int n_points, int i, float dx, centre new_centre){

        double x = (double(i) + 0.5) * dx;

	if(IC == 0){
		if(i==0){cout << "Using Sod Shock Tube" << endl;}
		new_centre.set_x(x);

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
                if(i==0){cout << "Using Sine Wave" << endl;}

                double rho,rho_0 = 50.0;
                double p,p_0 = 3.0;
                double c_s = sqrt(GAMMA*p_0/rho_0);
                double k = 2.0*3.1415/50.0;
                double epsilon = 0.0001;
                double v;

                new_centre.set_x(x);

                rho = rho_0*(1.0 + epsilon * cos(k * x));
                p = p_0 + c_s * c_s * (rho - rho_0);
                v = c_s * (rho - rho_0)/rho_0;

                //cout << v << endl;

		new_centre.set_mass_density(rho);	                // units kg/m^3
		new_centre.set_velocity(v);	                // units m/s
		new_centre.set_pressure(p);		        // units N/m^2

		new_centre.setup_energy_density();
		new_centre.con_to_prim();
		new_centre.setup_f_variables();
		new_centre.reset_du();

                return new_centre;

        }else if(IC == 2){

                if(i==0){cout << "Using Sedov Blast Wave" << endl;}
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
        }else{
                if(i==0){cout << "Using Gaussian Pulse" << endl;}
                double centre = 2.5;
                double s,w,rho,rho_0 = 10.0,rho_pulse = 50.0;
                double velocity = 5.0,pressure = 3.0;

                s = abs(centre-x);       // distance from centre of pulse
                w = 1.0;                 // characteristic width

                rho = rho_pulse*exp(-s*s/(w*w)) + rho_0*(1-exp(-s*s/(w*w)));

                new_centre.set_x(x);

                new_centre.set_mass_density(rho);
                new_centre.set_velocity(velocity);
                new_centre.set_pressure(pressure);

                new_centre.setup_energy_density();
                new_centre.con_to_prim();
                new_centre.setup_f_variables();
                new_centre.reset_du();

                return new_centre;

        }
}