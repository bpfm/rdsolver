/* class containing values associated with the face
	*centre_0 = pointer to left centre
	*centre_1 = pointer to right centre
	*centre_00 = pointer to left centre
	*centre_11 = pointer to right centre
*/

using namespace std;

class face{

private:

	int id;
	centre *centre_0,*centre_1,*centre_00,*centre_11;
	double q[4][5];

public:

	void set_id(int new_id){
		id = new_id;
	}

	void set_centre_0(centre* new_centre){
		centre_0 = new_centre;
	}

	void set_centre_1(centre* new_centre){
		centre_1 = new_centre;
	}

	void set_centre_00(centre* new_centre){
		centre_00 = new_centre;
	}

	void set_centre_11(centre* new_centre){
		centre_11 = new_centre;
	}

	centre* get_centre_0(){
		return centre_0;
	}

	centre* get_centre_1(){
		return centre_1;
	}

	void construct_state(double dx, double &dt, double t){
		double density[4],u[4],v[4],w[4],e_tot[4],pressure[4],h_tot[4];
		double density_avg,u_avg,v_avg,w_avg,e_tot_avg,pressure_avg,h_tot_avg,e_kin_avg;
		double e_vec[5][5],e_val[5],theta[5],phi[5],epsilon[5];
		double c_sound,gamma = 1.4,zeta,sum;
		double e_delta_q[5],delta_q[5],f0[5],f1[5],f_int[5],du0[3],du1[3];
		double r_int[5];

		//if(centre_0->get_x()>9.79 and centre_0->get_x()<9.81){
		density[0] = centre_0->get_mass_density();		// contruct state on either side of the face
		density[1] = centre_1->get_mass_density();
		density[2] = centre_00->get_mass_density();
		density[3] = centre_11->get_mass_density();

		u[0] = centre_0->get_velocity();
		u[1] = centre_1->get_velocity();
		u[2] = centre_00->get_velocity();
		u[3] = centre_11->get_velocity();

		e_tot[0] = centre_0->get_energy_density();
		e_tot[1] = centre_1->get_energy_density();
		e_tot[2] = centre_00->get_energy_density();
		e_tot[3] = centre_11->get_energy_density();

		pressure[0] = centre_0->get_pressure();
		pressure[1] = centre_1->get_pressure();
		pressure[2] = centre_00->get_pressure();
		pressure[3] = centre_11->get_pressure();

		h_tot[0] = e_tot[0]+pressure[0]/density[0];
		h_tot[1] = e_tot[1]+pressure[1]/density[1];
		h_tot[2] = e_tot[2]+pressure[2]/density[2];
		h_tot[3] = e_tot[3]+pressure[3]/density[3];

		/****** Construct Roe averages ******/

		density_avg = sqrt(density[0])*sqrt(density[1]);			// find average values at boundary
		e_tot_avg = (e_tot[0]+e_tot[1])/2.0;
		pressure_avg = (pressure[0]+pressure[1])/2.0;

		u_avg = roe_avg(density[0], u[0], density[1], u[1]);		// find Roe averages for state values
		v_avg = 0.0;
		w_avg = 0.0;
		h_tot_avg = roe_avg(density[0], h_tot[0], density[1], h_tot[1]);

		//cout << "avgs " << u_avg << " " << h_tot_avg << endl;

		/****** Constuct q vector on either side of boundary ******/

		for(int j=0;j<4;j++){
			v[j] = 0.0;
			w[j] = 0.0;
			q[j][0] = density[j]*e_tot[j];
			q[j][1] = density[j]*u[j];
			q[j][2] = density[j]*v[j];
			q[j][3] = density[j]*w[j];
			q[j][4] = density[j];
		}

		c_sound = sqrt(gamma*pressure_avg/density_avg);		// calculate averge adiabatic sound speed

		e_val[0] = u_avg - c_sound;				// calculate eigen values (for x direction)
		e_val[1] = u_avg + c_sound;
		e_val[2] = u_avg;
		e_val[3] = u_avg;
		e_val[4] = u_avg;

		for(int k=0;k<5;k++){
			if(u_avg<0){
				if(q[1][k]==q[0][k]){
					r_int[k]=0.0;
				}else{
					r_int[k] = (q[3][k]-q[0][k])/(q[0][k]-q[1][k]); // using equation 4.37 from lecture notes
				}
				theta[k] = -1.0;
				phi[k] = flux_limiter(r_int[k]);
			}else{
				if(q[1][k]==q[0][k]){
					r_int[k]=0.0;
				}else{
					r_int[k] = (q[1][k]-q[2][k])/(q[0][k]-q[1][k]); // using equation 4.37 from lecture notes
				}
				theta[k] = 1.0;
				phi[k] = flux_limiter(r_int[k]);
			}
			phi[k] = 0.0;					// simple Lax-Wendroff flux limiter (1.0) or donor cell (0.0)
			epsilon[k] = e_val[k]*dt/dx;			// (above equation 6.51 in lecture notes)
		}

		e_kin_avg = (u_avg*u_avg + v_avg*v_avg + w_avg*w_avg)/2.0;

		delta_q[0] = (density[1]*e_tot[1])-(density[0]*e_tot[0]);
		delta_q[1] = (density[1]*u[1])-(density[0]*u[0]);
		delta_q[2] = (density[1]*v[1])-(density[0]*v[0]);
		delta_q[3] = (density[1]*w[1])-(density[0]*w[0]);
		delta_q[4] = density[1]-density[0];

		zeta = u_avg*delta_q[1] + v_avg*delta_q[2] + w_avg*delta_q[3] - delta_q[0];

		e_delta_q[0] = (gamma-1.0)/(2.0*c_sound*c_sound)*(e_kin_avg*delta_q[4]-zeta)-(delta_q[1]-u_avg*delta_q[4])/(2.0*c_sound);
		e_delta_q[1] = (gamma-1.0)/(2.0*c_sound*c_sound)*(e_kin_avg*delta_q[4]-zeta)+(delta_q[1]-u_avg*delta_q[4])/(2.0*c_sound);
		e_delta_q[2] = (gamma-1.0)/(2.0*c_sound*c_sound)*((h_tot_avg-2.0*e_kin_avg)*delta_q[4]+zeta);
		e_delta_q[3] = delta_q[2]-v_avg*delta_q[4];
		e_delta_q[4] = delta_q[3]-w_avg*delta_q[4];

		f0[0] = density[0]*h_tot[0]*u[0];
		f0[1] = density[0]*u[0]*u[0]+pressure[0];
		f0[2] = density[0]*u[0]*v[0];
		f0[3] = density[0]*u[0]*w[0];
		f0[4] = density[0]*u[0];

		f1[0] = density[1]*h_tot[1]*u[1];
		f1[1] = density[1]*u[1]*u[1]+pressure[1];
		f1[2] = density[1]*u[1]*v[1];
		f1[3] = density[1]*u[1]*w[1];
		f1[4] = density[1]*u[1];

		for(int k=0;k<5;k++){
			sum = 0.0;
			for(int m=0;m<5;m++){sum  = sum + e_val[k]*e_delta_q[m];}//*(theta[m]+phi[m]*(epsilon[m]-theta[m]));}
			f_int[k] = 0.5*(f1[k]+f0[k])-0.5*(sum);
		}

		du0[0] = -f_int[4]*dt/dx;
		du0[1] = -f_int[1]*dt/dx;
		du0[2] = -f_int[0]*dt/dx;

		//cout << "du0 = " << du0[0] << " " << du0[1] << " " << du0[2] << endl;

		du1[0] = f_int[4]*dt/dx;
		du1[1] = f_int[1]*dt/dx;
		du1[2] = f_int[0]*dt/dx;

		
		/*if(centre_0->get_x()>centre_1->get_x()){	// create boundary at x=0.0 and x=20.0
			du0[2] = 0.0;
			du0[1] = 0.0;
			du0[0] = 0.0;

			du1[2] = 0.0;
			du1[1] = 0.0;
			du1[0] = 0.0;
		}*/
		
		
		/*if(centre_1->get_x()>13 and centre_1
			->get_x()<14){
			cout << "*************************************" << endl;
			cout << "positions of centres = " << centre_0->get_x() << " " << centre_1->get_x() << endl;
			cout << "starting values -> " << u[0] << " " << u[1] << endl;
			cout << "steps = " << dt << " " << dx << endl;
			cout << "delta_q = " << delta_q[0] << " " << delta_q[1] << " " << delta_q[2] << " " << delta_q[3] << " " << delta_q[4]<< endl;
			cout << "eigen values = " << e_val[0] << " " << e_val[1] << " " << e_val[2] << " " << e_val[3] << " " << e_val[4] << endl;
			cout << "sound speed components = " << pressure_avg << " " << density_avg << endl;
			cout << "components = " << gamma << " " << c_sound << " " << e_kin_avg << " " << u_avg << endl;
			cout << "e_delta_q components = " << e_kin_avg << " " << zeta << endl;
			cout << "e_delta_q = " << e_delta_q[0] << " " << e_delta_q[1] << " " << e_delta_q[2] << " " << e_delta_q[3] << " " << e_delta_q[4]<< endl;
			cout << "at " << (centre_0->get_x()+centre_1->get_x())/2.0 << " du0 = " << du0[0] << " " << du0[1] << " " << du0[2] << endl;
		}
		*/
		
		centre_0->update_du(du0);
		centre_1->update_du(du1);

		//}
	}

	// returns 1.0 for positive numbers and -1.0 for negative numbers
	double sign(double x){
		if(x>0.0){
			return 1.0;
		}else if(x<0.0){
			return -1.0;
		}else{
			return 0.0;
		}
	}

	// returns Roe average of left and right states
	double roe_avg(double l1, double l2, double r1, double r2){
		double avg;
		avg = (sqrt(l1)*l2+sqrt(r1)*r2)/(sqrt(l1)+sqrt(r1));
		return avg;
	}

	double flux_limiter(double r){
		double phi,twor;

		/*
		phi = mymin(1.0,r);					// minmod
		if(phi < 0.0){
			phi=0.0;
		}
		*/

		twor = 2.0*r;
		phi = mymax(0.0,mymin(1.0,twor),mymin(2.0,r));		// superbee

		//phi = r; 						// Beam-Warming

		//phi = (r+abs(r))/(1.0+abs(r)); 			// van Leer flux limiter

		//cout << phi << " " << twor << " " << r << endl;
		//if(phi==2){exit(0);}

		return phi;
	}

	double mymax(double val0, double val1, double val2){
		double max_val;
		double val[3];

		val[0]=val0;
		val[1]=val1;
		val[2]=val2;

		max_val = val[0];

		for (int i = 0; i < 3; ++i){
			if(val[i]>max_val){
				max_val = val[i];
			}
		}
		return max_val;
	}

	double mymin(double val0, double val1){
		double min_val;
		double val[2];

		val[0]=val0;
		val[1]=val1;

		min_val = val[0];

		for (int i = 0; i < 2; ++i){
			if(val[i]<min_val){
				min_val = val[i];
			}
		}
		return min_val;
	}

};