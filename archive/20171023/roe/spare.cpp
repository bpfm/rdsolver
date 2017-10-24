/*e_vec[0][0] = h_tot_avg - c_sound*u_avg;		// calculate eigenvectors (for x direction)
		e_vec[0][1] = u_avg - c_sound;
		e_vec[0][2] = v_avg;
		e_vec[0][3] = w_avg;
		e_vec[0][4] = 1.0;

		e_vec[1][0] = h_tot_avg+c_sound+u_avg;
		e_vec[1][1] = u_avg+c_sound;
		e_vec[1][2] = v_avg;
		e_vec[1][3] = w_avg;
		e_vec[1][4] = 1.0;

		e_vec[2][0] = 0.5*u_avg*u_avg;
		e_vec[2][1] = u_avg;
		e_vec[2][2] = v_avg;
		e_vec[2][3] = w_avg;
		e_vec[2][4] = 1.0;

		e_vec[3][0] = v_avg*v_avg;
		e_vec[3][1] = 0.0;
		e_vec[3][2] = 1.0;
		e_vec[3][3] = 0.0;
		e_vec[3][4] = 0.0;

		e_vec[4][0] = w_avg*w_avg;
		e_vec[4][1] = 0.0;
		e_vec[4][2] = 0.0;
		e_vec[4][3] = 1.0;
		e_vec[4][4] = 0.0;*/


		/*new_vertex.set_x(float(i)*dx);

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
		new_vertex.reset_du();*/