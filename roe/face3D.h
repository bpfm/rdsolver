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
        string direction;
        centre *centre_0,*centre_1,*centre_00,*centre_11;
        double q[4][5];

public:

        void set_id(int new_id){
                id = new_id;
        }

        void set_direction(string new_direction){
                direction = new_direction;
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

        void calculate_flux(double dx, double dt, double t, ofstream &du_file){
                double density[4],x_vel[4],y_vel[4],z_vel[4];
                double spec_energy[4],pressure[4],h_tot[4];
                double density_avg,spec_energy_avg,pressure_avg,h_tot_avg,e_kin_avg;
                double x_vel_avg,y_vel_avg,z_vel_avg;
                double gamma_bar,m_bar_sq,c_sq;
                double e_val[5];
                double c_sound,correction;
                double e_delta_q[5],delta_q[5],f0[5],f1[5],f_int[5],du0[5],du1[5];

                double r_lower[5][5],r_upper[5][5],alpha[5];
                double delta_p,delta_u,delta_d;

                double r_int[5],phi[5];

                //if(centre_0->get_x() >= (N_POINTS-1.01)*dx or centre_0->get_y() >= (N_POINTS-1.01)*dx or centre_0->get_z() >= (N_POINTS-1.01)*dx){return;}

                density[0] = centre_0->get_mass_density();              // contruct state on either side of the face
                density[1] = centre_1->get_mass_density();
                density[2] = centre_00->get_mass_density();
                density[3] = centre_11->get_mass_density();

                // import velocities from centres based on rotation of axes such that face is always aligned with x-axis

                if(direction == "x"){
                        // x-face case
                        x_vel[0] = centre_0->get_x_velocity();                        // q_i
                        x_vel[1] = centre_1->get_x_velocity();                        // q_i+1
                        x_vel[2] = centre_00->get_x_velocity();                       // q_i-1
                        x_vel[3] = centre_11->get_x_velocity();                       // q_i+2

                        y_vel[0] = centre_0->get_y_velocity();                        // q_i
                        y_vel[1] = centre_1->get_y_velocity();                        // q_i+1
                        y_vel[2] = centre_00->get_y_velocity();                       // q_i-1
                        y_vel[3] = centre_11->get_y_velocity();                       // q_i+2

                        z_vel[0] = centre_0->get_z_velocity();                        // q_i
                        z_vel[1] = centre_1->get_z_velocity();                        // q_i+1
                        z_vel[2] = centre_00->get_z_velocity();                       // q_i-1
                        z_vel[3] = centre_11->get_z_velocity();                       // q_i+2
                }else if(direction == "y"){
                        // y-face case
                        x_vel[0] = centre_0->get_z_velocity();                        // q_i
                        x_vel[1] = centre_1->get_z_velocity();                        // q_i+1
                        x_vel[2] = centre_00->get_z_velocity();                       // q_i-1
                        x_vel[3] = centre_11->get_z_velocity();                       // q_i+2

                        y_vel[0] = centre_0->get_x_velocity();                        // q_i
                        y_vel[1] = centre_1->get_x_velocity();                        // q_i+1
                        y_vel[2] = centre_00->get_x_velocity();                       // q_i-1
                        y_vel[3] = centre_11->get_x_velocity();                       // q_i+2

                        z_vel[0] = centre_0->get_y_velocity();                        // q_i
                        z_vel[1] = centre_1->get_y_velocity();                        // q_i+1
                        z_vel[2] = centre_00->get_y_velocity();                       // q_i-1
                        z_vel[3] = centre_11->get_y_velocity();                       // q_i+2
                }else if(direction == "z"){
                        // z-face case
                        x_vel[0] = centre_0->get_y_velocity();                        // q_i
                        x_vel[1] = centre_1->get_y_velocity();                        // q_i+1
                        x_vel[2] = centre_00->get_y_velocity();                       // q_i-1
                        x_vel[3] = centre_11->get_y_velocity();                       // q_i+2

                        y_vel[0] = centre_0->get_z_velocity();                        // q_i
                        y_vel[1] = centre_1->get_z_velocity();                        // q_i+1
                        y_vel[2] = centre_00->get_z_velocity();                       // q_i-1
                        y_vel[3] = centre_11->get_z_velocity();                       // q_i+2

                        z_vel[0] = centre_0->get_x_velocity();                        // q_i
                        z_vel[1] = centre_1->get_x_velocity();                        // q_i+1
                        z_vel[2] = centre_00->get_x_velocity();                       // q_i-1
                        z_vel[3] = centre_11->get_x_velocity();                       // q_i+2
                }

                spec_energy[0] = centre_0->get_specific_energy();
                spec_energy[1] = centre_1->get_specific_energy();
                spec_energy[2] = centre_00->get_specific_energy();
                spec_energy[3] = centre_11->get_specific_energy();

                pressure[0] = centre_0->get_pressure();
                pressure[1] = centre_1->get_pressure();
                pressure[2] = centre_00->get_pressure();
                pressure[3] = centre_11->get_pressure();

                h_tot[0] = spec_energy[0]+pressure[0]/density[0];
                h_tot[1] = spec_energy[1]+pressure[1]/density[1];
                h_tot[2] = spec_energy[2]+pressure[2]/density[2];
                h_tot[3] = spec_energy[3]+pressure[3]/density[3];

                /****** Construct Roe averages ******/

                density_avg = sqrt(density[0])*sqrt(density[1]);                        // find average values at boundary
                spec_energy_avg = (spec_energy[0]+spec_energy[1])/2.0;
                pressure_avg = (pressure[0]+pressure[1])/2.0;

                x_vel_avg = roe_avg(density[0], x_vel[0], density[1], x_vel[1]);           // find Roe averages for state values
                y_vel_avg = roe_avg(density[0], y_vel[0], density[1], y_vel[1]);
                z_vel_avg = roe_avg(density[0], z_vel[0], density[1], z_vel[1]);
                h_tot_avg = roe_avg(density[0], h_tot[0], density[1], h_tot[1]);

                e_kin_avg = (x_vel_avg*x_vel_avg + y_vel_avg*y_vel_avg + z_vel_avg*z_vel_avg)/2.0;

                /****** Constuct q vector on either side of boundary ******/

                for(int j=0;j<4;j++){
                        q[j][0] = density[j];
                        q[j][1] = density[j]*x_vel[j];
                        q[j][2] = density[j]*y_vel[j];
                        q[j][3] = density[j]*z_vel[j];
                        q[j][4] = density[j]*spec_energy[j];
                }

                c_sound = sqrt((GAMMA-1.0)*(h_tot_avg - e_kin_avg));

                e_val[0] = x_vel_avg - c_sound;                             // calculate eigenvalues
                e_val[1] = x_vel_avg;
                e_val[2] = x_vel_avg;
                e_val[3] = x_vel_avg;
                e_val[4] = x_vel_avg + c_sound;

                if(e_val[0] > 0.0){

                        f_int[0] = density[0]*x_vel[0];
                        f_int[1] = density[0]*x_vel[0]*x_vel[0]+pressure[0];
                        f_int[2] = density[0]*x_vel[0]*y_vel[0];
                        f_int[3] = density[0]*x_vel[0]*z_vel[0];
                        f_int[4] = density[0]*h_tot[0]*x_vel[0];

                        du0[0] = -f_int[0]*dt/dx;
                        du0[1] = -f_int[1]*dt/dx;
                        du0[2] = -f_int[2]*dt/dx;
                        du0[3] = -f_int[3]*dt/dx;
                        du0[4] = -f_int[4]*dt/dx;

                        du1[0] = f_int[0]*dt/dx;
                        du1[1] = f_int[1]*dt/dx;
                        du1[2] = f_int[2]*dt/dx;
                        du1[3] = f_int[3]*dt/dx;
                        du1[4] = f_int[4]*dt/dx;

                        centre_0->update_du(du0);
                        centre_1->update_du(du1);

                        cout << "du0 =\t" << du0[0] << "\t" << du0[1] << "\t" << du0[2] << endl;
  
                        return ;
                }

                if(e_val[2] < 0.0){
                  
                        f_int[0] = density[1]*x_vel[1];
                        f_int[1] = density[1]*x_vel[1]*x_vel[1]+pressure[1];
                        f_int[2] = density[1]*x_vel[1]*y_vel[1];
                        f_int[3] = density[1]*x_vel[1]*z_vel[1];
                        f_int[4] = density[1]*h_tot[1]*x_vel[1];

                        du0[0] = -f_int[0]*dt/dx;
                        du0[1] = -f_int[1]*dt/dx;
                        du0[2] = -f_int[2]*dt/dx;
                        du0[3] = -f_int[3]*dt/dx;
                        du0[4] = -f_int[4]*dt/dx;

                        du1[0] = f_int[0]*dt/dx;
                        du1[1] = f_int[1]*dt/dx;
                        du1[2] = f_int[2]*dt/dx;
                        du1[3] = f_int[3]*dt/dx;
                        du1[4] = f_int[4]*dt/dx;

                        centre_0->update_du(du0);
                        centre_1->update_du(du1);

                        cout << "du0 =\t" << du0[0] << "\t" << du0[1] << "\t" << du0[2] << endl;

                        return ;
                }


                // in case of transonic fluxes, there may be numerical problems
                // this is a possible fix follows:

                double dlambda, eps;
                dlambda = fabs(fabs(x_vel_avg) - c_sound) / c_sound;

                if(dlambda < 1e-6){
                        eps = 1e-6 * c_sound;
                        e_val[0] = 0.5 * (e_val[0] * e_val[0] / eps + eps);
                        e_val[4] = 0.5 * (e_val[4] * e_val[4] / eps + eps);
                }

                delta_q[0] = density[1]-density[0];
                delta_q[1] = (density[1]*x_vel[1])-(density[0]*x_vel[0]);
                delta_q[2] = (density[1]*y_vel[1])-(density[0]*y_vel[0]);
                delta_q[3] = (density[1]*z_vel[1])-(density[0]*z_vel[0]);
                delta_q[4] = (density[1]*spec_energy[1])-(density[0]*spec_energy[0]);

                /****** Construct R matrix ******/

                r_lower[0][0] = 1.0;
                r_lower[0][1] = 1.0;
                r_lower[0][2] = 0.0;
                r_lower[0][3] = 0.0;
                r_lower[0][4] = 1.0;

                r_lower[1][0] = x_vel_avg - c_sound;
                r_lower[1][1] = x_vel_avg;
                r_lower[1][2] = 0.0;
                r_lower[1][3] = 0.0;
                r_lower[1][4] = x_vel_avg + c_sound;

                r_lower[2][0] = y_vel_avg;
                r_lower[2][1] = y_vel_avg;
                r_lower[2][2] = 1.0;
                r_lower[2][3] = 0.0;
                r_lower[2][4] = y_vel_avg;

                r_lower[3][0] = z_vel_avg;
                r_lower[3][1] = z_vel_avg;
                r_lower[3][2] = 0.0;
                r_lower[3][3] = 1.0;
                r_lower[3][4] = z_vel_avg;

                r_lower[4][0] = h_tot_avg - x_vel_avg*c_sound;
                r_lower[4][1] = 0.5*(x_vel_avg*x_vel_avg + y_vel_avg*y_vel_avg + z_vel_avg*z_vel_avg);
                r_lower[4][2] = y_vel_avg;
                r_lower[4][3] = z_vel_avg;
                r_lower[4][4] = h_tot_avg + x_vel_avg*c_sound;

                gamma_bar = GAMMA-1.0;
                c_sq = c_sound*c_sound;
                m_bar_sq = (x_vel_avg*x_vel_avg + y_vel_avg*y_vel_avg + z_vel_avg*z_vel_avg)/(c_sq);

                r_upper[0][0] = x_vel_avg/(2.0*c_sound) + gamma_bar*m_bar_sq/4.0;
                r_upper[0][1] = -1.0/(2.0*c_sound) - gamma_bar*x_vel_avg/(2.0*c_sq);
                r_upper[0][2] = -1.0*gamma_bar * y_vel_avg/(2.0*c_sq);
                r_upper[0][3] = -1.0*gamma_bar * z_vel_avg/(2.0*c_sq);
                r_upper[0][4] = gamma_bar/(2.0*c_sq);

                r_upper[1][0] = 1.0 - gamma_bar*m_bar_sq/2.0;
                r_upper[1][1] = gamma_bar*x_vel_avg/(c_sq);
                r_upper[1][2] = gamma_bar*y_vel_avg/(c_sq);
                r_upper[1][3] = gamma_bar*z_vel_avg/(c_sq);
                r_upper[1][4] = -1.0*gamma_bar/(c_sq);

                r_upper[2][0] = -1.0*y_vel_avg;
                r_upper[2][1] = 0.0;
                r_upper[2][2] = 1.0;
                r_upper[2][3] = 0.0;
                r_upper[2][4] = 0.0;

                r_upper[3][0] = -1.0*z_vel_avg;
                r_upper[3][1] = 0.0;
                r_upper[3][2] = 0.0;
                r_upper[3][3] = 1.0;
                r_upper[3][4] = 0.0;

                r_upper[4][0] = -1.0*x_vel_avg/(2.0*c_sound) + gamma_bar*m_bar_sq/4.0;
                r_upper[4][1] = 1.0/(2.0*c_sound) - gamma_bar*x_vel_avg/(2.0*c_sq);
                r_upper[4][2] = -1.0*gamma_bar * y_vel_avg/(2.0*c_sq);
                r_upper[4][3] = -1.0*gamma_bar * z_vel_avg/(2.0*c_sq);
                r_upper[4][4] = gamma_bar/(2.0*c_sq);

                /****** Check matrices give identity matrix when multiplied together ******/

                /*
                cout << r_lower[0][0]*r_upper[0][0] + r_lower[0][1]*r_upper[1][0] + r_lower[0][2]*r_upper[2][0] + r_lower[0][3]*r_upper[3][0] + r_lower[0][4]*r_upper[4][0] << "\t";
                cout << r_lower[0][0]*r_upper[0][1] + r_lower[0][1]*r_upper[1][1] + r_lower[0][2]*r_upper[2][1] + r_lower[0][3]*r_upper[3][1] + r_lower[0][4]*r_upper[4][1] << "\t";
                cout << r_lower[0][0]*r_upper[0][2] + r_lower[0][1]*r_upper[1][2] + r_lower[0][2]*r_upper[2][2] + r_lower[0][3]*r_upper[3][2] + r_lower[0][4]*r_upper[4][2] << "\t";
                cout << r_lower[0][0]*r_upper[0][3] + r_lower[0][1]*r_upper[1][3] + r_lower[0][2]*r_upper[2][3] + r_lower[0][3]*r_upper[3][3] + r_lower[0][4]*r_upper[4][3] << "\t";
                cout << r_lower[0][0]*r_upper[0][4] + r_lower[0][1]*r_upper[1][4] + r_lower[0][2]*r_upper[2][4] + r_lower[0][3]*r_upper[3][4] + r_lower[0][4]*r_upper[4][4] << endl;

                cout << r_lower[1][0]*r_upper[0][0] + r_lower[1][1]*r_upper[1][0] + r_lower[1][2]*r_upper[2][0] + r_lower[1][3]*r_upper[3][0] + r_lower[1][4]*r_upper[4][0] << "\t";
                cout << r_lower[1][0]*r_upper[0][1] + r_lower[1][1]*r_upper[1][1] + r_lower[1][2]*r_upper[2][1] + r_lower[1][3]*r_upper[3][1] + r_lower[1][4]*r_upper[4][1] << "\t";
                cout << r_lower[1][0]*r_upper[0][2] + r_lower[1][1]*r_upper[1][2] + r_lower[1][2]*r_upper[2][2] + r_lower[1][3]*r_upper[3][2] + r_lower[1][4]*r_upper[4][2] << "\t";
                cout << r_lower[1][0]*r_upper[0][3] + r_lower[1][1]*r_upper[1][3] + r_lower[1][2]*r_upper[2][3] + r_lower[1][3]*r_upper[3][3] + r_lower[1][4]*r_upper[4][3] << "\t";
                cout << r_lower[1][0]*r_upper[0][4] + r_lower[1][1]*r_upper[1][4] + r_lower[1][2]*r_upper[2][4] + r_lower[1][3]*r_upper[3][4] + r_lower[1][4]*r_upper[4][4] << endl;

                cout << r_lower[2][0]*r_upper[0][0] + r_lower[2][1]*r_upper[1][0] + r_lower[2][2]*r_upper[2][0] + r_lower[2][3]*r_upper[3][0] + r_lower[2][4]*r_upper[4][0] << "\t";
                cout << r_lower[2][0]*r_upper[0][1] + r_lower[2][1]*r_upper[1][1] + r_lower[2][2]*r_upper[2][1] + r_lower[2][3]*r_upper[3][1] + r_lower[2][4]*r_upper[4][1] << "\t";
                cout << r_lower[2][0]*r_upper[0][2] + r_lower[2][1]*r_upper[1][2] + r_lower[2][2]*r_upper[2][2] + r_lower[2][3]*r_upper[3][2] + r_lower[2][4]*r_upper[4][2] << "\t";
                cout << r_lower[2][0]*r_upper[0][3] + r_lower[2][1]*r_upper[1][3] + r_lower[2][2]*r_upper[2][3] + r_lower[2][3]*r_upper[3][3] + r_lower[2][4]*r_upper[4][3] << "\t";
                cout << r_lower[2][0]*r_upper[0][4] + r_lower[2][1]*r_upper[1][4] + r_lower[2][2]*r_upper[2][4] + r_lower[2][3]*r_upper[3][4] + r_lower[2][4]*r_upper[4][4] << endl;

                cout << r_lower[3][0]*r_upper[0][0] + r_lower[3][1]*r_upper[1][0] + r_lower[3][2]*r_upper[2][0] + r_lower[3][3]*r_upper[3][0] + r_lower[3][4]*r_upper[4][0] << "\t";
                cout << r_lower[3][0]*r_upper[0][1] + r_lower[3][1]*r_upper[1][1] + r_lower[3][2]*r_upper[2][1] + r_lower[3][3]*r_upper[3][1] + r_lower[3][4]*r_upper[4][1] << "\t";
                cout << r_lower[3][0]*r_upper[0][2] + r_lower[3][1]*r_upper[1][2] + r_lower[3][2]*r_upper[2][2] + r_lower[3][3]*r_upper[3][2] + r_lower[3][4]*r_upper[4][2] << "\t";
                cout << r_lower[3][0]*r_upper[0][3] + r_lower[3][1]*r_upper[1][3] + r_lower[3][2]*r_upper[2][3] + r_lower[3][3]*r_upper[3][3] + r_lower[3][4]*r_upper[4][3] << "\t";
                cout << r_lower[3][0]*r_upper[0][4] + r_lower[3][1]*r_upper[1][4] + r_lower[3][2]*r_upper[2][4] + r_lower[3][3]*r_upper[3][4] + r_lower[3][4]*r_upper[4][4] << endl;

                cout << r_lower[4][0]*r_upper[0][0] + r_lower[4][1]*r_upper[1][0] + r_lower[4][2]*r_upper[2][0] + r_lower[4][3]*r_upper[3][0] + r_lower[4][4]*r_upper[4][0] << "\t";
                cout << r_lower[4][0]*r_upper[0][1] + r_lower[4][1]*r_upper[1][1] + r_lower[4][2]*r_upper[2][1] + r_lower[4][3]*r_upper[3][1] + r_lower[4][4]*r_upper[4][1] << "\t";
                cout << r_lower[4][0]*r_upper[0][2] + r_lower[4][1]*r_upper[1][2] + r_lower[4][2]*r_upper[2][2] + r_lower[4][3]*r_upper[3][2] + r_lower[4][4]*r_upper[4][2] << "\t";
                cout << r_lower[4][0]*r_upper[0][3] + r_lower[4][1]*r_upper[1][3] + r_lower[4][2]*r_upper[2][3] + r_lower[4][3]*r_upper[3][3] + r_lower[4][4]*r_upper[4][3] << "\t";
                cout << r_lower[4][0]*r_upper[0][4] + r_lower[4][1]*r_upper[1][4] + r_lower[4][2]*r_upper[2][4] + r_lower[4][3]*r_upper[3][4] + r_lower[4][4]*r_upper[4][4] << endl;
                
                exit(0);

                */

                /****** Calculate alpha  ******/

                // brute force method - need to wokr out algebraic solution to multiplication

                alpha[0] = alpha[1] = alpha[2] = alpha[3] = alpha[4] = 0.0;

                for(int k=0;k<5;k++){alpha[0] += r_upper[k][0]*delta_q[k];}
                for(int k=0;k<5;k++){alpha[1] += r_upper[k][1]*delta_q[k];}
                for(int k=0;k<5;k++){alpha[2] += r_upper[k][2]*delta_q[k];}
                for(int k=0;k<5;k++){alpha[3] += r_upper[k][3]*delta_q[k];}
                for(int k=0;k<5;k++){alpha[4] += r_upper[k][4]*delta_q[k];}

                /****** Constuct fluxes either side of face ******/

                f0[0] = density[0]*x_vel[0];
                f0[1] = density[0]*x_vel[0]*x_vel[0]+pressure[0];
                f0[2] = density[0]*x_vel[0]*y_vel[0];
                f0[3] = density[0]*x_vel[0]*z_vel[0];
                f0[4] = density[0]*h_tot[0]*x_vel[0];

                f1[0] = density[1]*x_vel[1];
                f1[1] = density[1]*x_vel[1]*x_vel[1]+pressure[1];
                f1[2] = density[1]*x_vel[1]*y_vel[1];
                f1[3] = density[1]*x_vel[1]*z_vel[1];
                f1[4] = density[1]*h_tot[1]*x_vel[1];

                /****** Calculate flux limiter ******/

                for(int i=0;i<5;i++){
                        if(x_vel_avg > 0.0){
                                r_int[i] = (q[1][i]-q[0][i])/(q[0][i]-q[2][i]);
                        }else if(x_vel_avg < 0.0){
                                r_int[i] = (q[0][i]-q[1][i])/(q[1][i]-q[3][i]);
                        }
                }

                if(FLUX_LIMITER_OI == 0){
                        phi[0] = 1.0;
                        phi[1] = 1.0;
                        phi[2] = 1.0;
                }else{
                        for(int i=0;i<3;i++){
                                phi[i] = construct_flux_limiter(r_int[i]);
                        }
                }

                /****** Calculate interface flux ******/

                for(int k=0;k<5;k++){
                        correction = 0.0;
                        for(int m=0;m<5;m++){correction += alpha[m]*abs(e_val[m])*r_lower[k][m];}
                        f_int[k] = 0.5*phi[k]*(f1[k]+f0[k]);//-correction);
                }

                /****** Convert flux at face to change in variable in cells ******/

                du0[0] = -f_int[0]*dt/dx;
                du0[1] = -f_int[1]*dt/dx;
                du0[2] = -f_int[2]*dt/dx;
                du0[3] = -f_int[3]*dt/dx;
                du0[4] = -f_int[4]*dt/dx;

                du1[0] = f_int[0]*dt/dx;
                du1[1] = f_int[1]*dt/dx;
                du1[2] = f_int[2]*dt/dx;
                du1[3] = f_int[3]*dt/dx;
                du1[4] = f_int[4]*dt/dx;

                /****** Distribute changes to appropriate centres ******/

                centre_0->update_du(du0);
                centre_1->update_du(du1);

                /*
                if(centre_0->get_y()==22.5 and centre_0->get_z()==22.5){
                        cout << "x =\t" << centre_0->get_x() << "\ty =\t" << centre_0->get_y() << "\tz =\t" << centre_0->get_z() << "\tface =\t" << direction;
                        cout << "\tdensity =\t" << density[0] << "\t" << density[1] << "\tdu1 =\t" << du1[0] << "\t" << du1[1] << "\t" << du1[2] << "\t" << du1[3] << "\t" << du1[4] << endl;
                        //cout << "\tLHS fluxes =\t" << f0[0] << "\t" << f0[1] << "\t" << f0[2] << "\t" << f0[3] << "\t" <<f0[4];
                        //cout << "\tRHS fluxes =\t" << f1[0] << "\t" << f1[1] << "\t" << f1[2] << "\t" << f1[3] << "\t" <<f1[4]<< endl;
                }
                */

                if(isnan(du1[0]) or isnan(du0[0])){
                        cout << "time =\t" << t << endl;
                        cout << "direction =\t" << direction << endl;
                        cout << "x =\t" << centre_0->get_x() << "\ty =\t" << centre_0->get_y() << "\tz =\t" << centre_0->get_z() << endl;
                        cout << "vel_avg =\t" << x_vel_avg << "\t" << y_vel_avg << "\t" << z_vel_avg << endl;
                        cout << "spec_energy =\t" << spec_energy[0] << "\t" << spec_energy[1] << endl;
                        cout << "pressure =\t" << pressure[0] << "\t" << spec_energy[1] << endl;
                        cout << "h_tot =\t" << h_tot[0] << "\t" <<  h_tot[1] << endl;
                        cout << "gamma_bar =\t" << gamma_bar << endl;
                        cout << "m_bar_sq =\t" << m_bar_sq << endl;
                        cout << "delta_q =\t" << delta_q[0] << "\t" << delta_q[1] << "\t" << delta_q[2] << "\t" << delta_q[3] << "\t" << delta_q[4] << "\t" << endl;
                        /*
                        cout << r_upper[0][0] << "\t" << r_upper[0][1] << "\t" << r_upper[0][2] << "\t" << r_upper[0][3] << "\t" << r_upper[0][4] << endl;
                        cout << r_upper[1][0] << "\t" << r_upper[1][1] << "\t" << r_upper[1][2] << "\t" << r_upper[1][3] << "\t" << r_upper[1][4] << endl;
                        cout << r_upper[2][0] << "\t" << r_upper[2][1] << "\t" << r_upper[2][2] << "\t" << r_upper[2][3] << "\t" << r_upper[2][4] << endl;
                        cout << r_upper[3][0] << "\t" << r_upper[3][1] << "\t" << r_upper[3][2] << "\t" << r_upper[3][3] << "\t" << r_upper[3][4] << endl;
                        cout << r_upper[4][0] << "\t" << r_upper[4][1] << "\t" << r_upper[4][2] << "\t" << r_upper[4][3] << "\t" << r_upper[4][4] << endl;
                        */
                        cout << "alpha =\t" << alpha[0] << "\t" << alpha[1] << "\t" << alpha[2] << "\t" << alpha[3] << "\t" << alpha[4] << "\t" << endl;
                        cout << "correction =\t" << correction << endl;
                        cout << "c_sound =\t" << c_sound << endl;
                        cout << "h_tot_avg =\t" << h_tot_avg << endl;
                        cout << "e_kin_avg =\t" << e_kin_avg << endl;
                        cout << "GAMMA =\t" << GAMMA << endl;
                        cout << "density =\t" << density[0] << "\t" << density[1] << endl;
                        cout << "EXITING" << endl;
                        exit(0);
                }

                return ;
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

        double construct_flux_limiter(double r){
                double phi,twor;

                if(FLUX_LIMITER == "MINMOD"){
                        phi = mymin(1.0,r);                                     // minmod
                }else if(FLUX_LIMITER == "SUERBEE"){
                        twor = 2.0*r;
                        phi = mymax(0.0,mymin(1.0,twor),mymin(2.0,r));          // superbee
                }else if(FLUX_LIMITER == "BEAM-WARMING"){
                        phi = r;                                              // Beam-Warming
                }else if(FLUX_LIMITER == "VAN-LEER"){
                        phi = (r+abs(r))/(1.0+abs(r));                        // van Leer
                }

                if(phi < 0.0){
                        phi=0.0;
                }

                cout << "phi = " << phi << "\tr = " << r << endl;

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

        // calculate analytic value for du for gaussian pulse case
        double analytic_solution(double dx, double x, double dt, double t){
                double rho,rho_0,rho_1,x_0,sig,v,du,x_0t;
                rho_0 = 10.0;
                rho_1 = 50.0;
                x_0 = 10.0;
                sig = 2.0;
                v = 1.0;
                x_0t = x_0 + v*t;
                rho = rho_1 * exp(-(x-x_0t)*(x-x_0t)/(sig*sig)) + rho_0*(1.0-exp(-(x-x_0t)*(x-x_0t)/(sig*sig)));
                du = rho * v * dt/dx;
                return du;
        }

};