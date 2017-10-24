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
        double q[4][3];

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

        void calculate_flux(double dx, double dt, double t, ofstream &du_file){
                double x_face,x0,x1,density[4],u[4],e_tot[4],pressure[4],h_tot[4];
                double density_avg,u_avg,e_tot_avg,pressure_avg,h_tot_avg,e_kin_avg;
                double e_vec[3][3],e_val[3];
                double c_sound,correction;
                double e_delta_q[3],delta_q[3],f0[3],f1[3],f_int[3],du0[3],du1[3];

                double r_lower[3][3],r_upper[3][3],alpha[3];
                double delta_p,delta_u,delta_d;

                double r_int[3],phi[3];

                x_face = (centre_0->get_x()+centre_1->get_x())/2.0;
                x0 = centre_0->get_x();
                x1 = centre_1->get_x();

                density[0] = centre_0->get_mass_density();              // contruct state on either side of the face
                density[1] = centre_1->get_mass_density();
                density[2] = centre_00->get_mass_density();
                density[3] = centre_11->get_mass_density();

                u[0] = centre_0->get_velocity();                        // q_i
                u[1] = centre_1->get_velocity();                        // q_i+1
                u[2] = centre_00->get_velocity();                       // q_i-1
                u[3] = centre_11->get_velocity();                       // q_i+2

                e_tot[0] = centre_0->get_specific_energy();
                e_tot[1] = centre_1->get_specific_energy();
                e_tot[2] = centre_00->get_specific_energy();
                e_tot[3] = centre_11->get_specific_energy();

                pressure[0] = centre_0->get_pressure();
                pressure[1] = centre_1->get_pressure();
                pressure[2] = centre_00->get_pressure();
                pressure[3] = centre_11->get_pressure();

                h_tot[0] = e_tot[0]+pressure[0]/density[0];
                h_tot[1] = e_tot[1]+pressure[1]/density[1];
                h_tot[2] = e_tot[2]+pressure[2]/density[2];
                h_tot[3] = e_tot[3]+pressure[3]/density[3];

                /****** Construct Roe averages ******/

                density_avg = sqrt(density[0])*sqrt(density[1]);                        // find average values at boundary
                e_tot_avg = (e_tot[0]+e_tot[1])/2.0;
                pressure_avg = (pressure[0]+pressure[1])/2.0;

                u_avg = roe_avg(density[0], u[0], density[1], u[1]);                    // find Roe averages for state values
                h_tot_avg = roe_avg(density[0], h_tot[0], density[1], h_tot[1]);

                e_kin_avg = (u_avg*u_avg)/2.0;

                /****** Constuct q vector on either side of boundary ******/

                for(int j=0;j<4;j++){
                        q[j][0] = density[j];
                        q[j][1] = density[j]*u[j];
                        q[j][2] = density[j]*e_tot[j];
                }

                c_sound = sqrt((GAMMA-1.0)*(h_tot_avg - e_kin_avg));

                e_val[0] = u_avg - c_sound;                             // calculate eigenvalues
                e_val[1] = u_avg;
                e_val[2] = u_avg + c_sound;

                if(e_val[0] > 0.0){

                        cout << "supersonic " << c_sound << endl;

                        f_int[0] = density[0]*u[0];
                        f_int[1] = density[0]*u[0]*u[0]+pressure[0];
                        f_int[2] = density[0]*h_tot[0]*u[0];

                        du0[0] = -f_int[0]*dt/dx;
                        du0[1] = -f_int[1]*dt/dx;
                        du0[2] = -f_int[2]*dt/dx;

                        du1[0] = f_int[0]*dt/dx;
                        du1[1] = f_int[1]*dt/dx;
                        du1[2] = f_int[2]*dt/dx;

                        centre_0->update_du(du0);
                        centre_1->update_du(du1);

                        //cout << "du0 =\t" << du0[0] << "\t" << du0[1] << "\t" << du0[2] << endl;
  
                        return ;
                }

                if(e_val[2] < 0.0){
                  
                        f_int[0] = density[1]*u[1];
                        f_int[1] = density[1]*u[1]*u[1]+pressure[1];
                        f_int[2] = density[1]*h_tot[1]*u[1];

                        du0[0] = -f_int[0]*dt/dx;
                        du0[1] = -f_int[1]*dt/dx;
                        du0[2] = -f_int[2]*dt/dx;

                        du1[0] = f_int[0]*dt/dx;
                        du1[1] = f_int[1]*dt/dx;
                        du1[2] = f_int[2]*dt/dx;

                        centre_0->update_du(du0);
                        centre_1->update_du(du1);

                        //cout << "du0 =\t" << du0[0] << "\t" << du0[1] << "\t" << du0[2] << endl;

                        return ;
                }


                // in case of transonic fluxes, there may be numerical problems
                // a possible fix follows

                double dlambda, eps;
                dlambda = fabs(fabs(u_avg) - c_sound) / c_sound;

                if(dlambda < 1e-6){
                        eps = 1e-6 * c_sound;
                        e_val[0] = 0.5 * (e_val[0] * e_val[0] / eps + eps);
                        e_val[2] = 0.5 * (e_val[2] * e_val[2] / eps + eps);
                }

                delta_q[0] = density[1]-density[0];
                delta_q[1] = (density[1]*u[1])-(density[0]*u[0]);
                delta_q[2] = (density[1]*e_tot[1])-(density[0]*e_tot[0]);

                /****** Construct R matrix ******/

                r_lower[0][0] = 1.0;
                r_lower[0][1] = 1.0;
                r_lower[0][2] = 1.0;
                r_lower[1][0] = u_avg - c_sound;
                r_lower[1][1] = u_avg;
                r_lower[1][2] = u_avg + c_sound;
                r_lower[2][0] = h_tot_avg - u_avg*c_sound;
                r_lower[2][1] = 0.5*u_avg*u_avg;
                r_lower[2][2] = h_tot_avg + u_avg*c_sound;

                /****** Calculate alpha  ******/

                delta_p = pressure[1]-pressure[0];
                delta_u = u[1]-u[0];
                delta_d = density[1]-density[0];

                alpha[0] = (delta_p - c_sound*density_avg*delta_u)/(2.0*c_sound*c_sound);
                alpha[1] = (delta_d - delta_p/(c_sound*c_sound));
                alpha[2] = (delta_p + c_sound*density_avg*delta_u)/(2.0*c_sound*c_sound);

                /****** Constuct fluxes either side of face ******/

                f0[0] = density[0]*u[0];
                f0[1] = density[0]*u[0]*u[0]+pressure[0];
                f0[2] = density[0]*h_tot[0]*u[0];

                f1[0] = density[1]*u[1];
                f1[1] = density[1]*u[1]*u[1]+pressure[1];
                f1[2] = density[1]*h_tot[1]*u[1];

                /****** Calculate flux limiter ******/

                for(int i=0;i<3;i++){
                        if(u_avg > 0.0){
                                r_int[i] = (q[1][i]-q[0][i])/(q[0][i]-q[2][i]);
                        }else if(u_avg < 0.0){
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

                for(int k=0;k<3;k++){
                        correction = 0.0;
                        for(int m=0;m<3;m++){correction  = correction + alpha[m]*abs(e_val[m])*r_lower[k][m];}
                        f_int[k] = 0.5*phi[k]*(f1[k]+f0[k]-correction);
                }

                /****** Convert flux at face to change in variable ******/

                du0[0] = -f_int[0]*dt/dx;
                du0[1] = -f_int[1]*dt/dx;
                du0[2] = -f_int[2]*dt/dx;

                du1[0] = f_int[0]*dt/dx;
                du1[1] = f_int[1]*dt/dx;
                du1[2] = f_int[2]*dt/dx;

                if(FORCE_ANALYTIC_PULSE == 1){
                        //cout << "before du1[0] = " << du1[0] << endl;

                        du1[0] = analytic_solution(dx,x_face,dt,t);
                        du0[0] = -1.0*du1[0];
                        du_file << x_face << "\t" << du1[0] << endl;

                        //cout << "after du1[0] = " << du1[0] << endl;
                }

                /****** Distribute changes to appropriate cell centres ******/

                centre_0->update_du(du0);
                centre_1->update_du(du1);

                //cout << "du1 =\t" << du1[0] << "\t" << du1[1] << "\t" << du1[2] << endl;

                if(isnan(du1[0]) or isnan(du0[0])){
                        cout << "f0 =\t" << f0[0] << "\t" << f0[1] << "\t" << f0[2] << endl;
                        cout << "f1 =\t" << f1[0] << "\t" << f1[1] << "\t" << f1[2] << endl;
                        cout << "correction =\t" << correction << endl;
                        cout << "alpha =\t" << alpha[0] << "\t" << alpha[0] << "\t" << alpha[0] << endl;
                        cout << "delta_p =\t" << delta_p << endl;
                        cout << "delta_u =\t" << delta_u << endl;
                        cout << "delta_d =\t" << delta_d << endl;
                        cout << "c_sound =\t" << c_sound << endl;
                        cout << "h_tot_avg =\t" << h_tot_avg << endl;
                        cout << "e_kin_avg =\t" << e_kin_avg << endl;
                        cout << "GAMMA =\t" << GAMMA << endl; 
                        cout << "c_sound =\t" << c_sound << "\t" << GAMMA << "\th_tot_avg =\t" << h_tot_avg << endl;
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

        double analytic_solution(double dx, double x, double dt, double t){
                double rho,rho_0,rho_1,x_0,sig,v,du,x_0t;
                rho_0 = 10.0;
                rho_1 = 50.0;
                x_0 = 10.0;
                sig = 2.0;
                v = 1.0;
                x_0t = x_0 + v*t;
                //du = v*(rho_1 - rho_0)*(((-2.0)*(x-x_0t)/(sig*sig))*exp(-((x-x_0t)*(x-x_0t))/(sig*sig)))*dt/dx;
                rho = rho_1 * exp(-(x-x_0t)*(x-x_0t)/(sig*sig)) + rho_0*(1.0-exp(-(x-x_0t)*(x-x_0t)/(sig*sig)));
                du = rho * v * dt/dx;
                //cout << "du = " << du << "\t" << x_0t << endl;
                return du;
        }

};