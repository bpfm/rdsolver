using namespace std;

centre setup_centre(int N_POINTS, int i, float dx){
        centre new_centre;

        double x = (double(i) + 0.5) * dx;

        if(IC == 0){
                if(i==0){cout << "Using 1D Sod Shock Tube" << endl;}
                new_centre.set_x(x);

                if(i<0.5*N_POINTS){
                        new_centre.set_mass_density(1.0);                               // units kg/m^3
                        new_centre.set_velocity(0.00000001);                                   // units m/s
                        new_centre.set_pressure(500.0);                                 // units N/m^2
                }else{
                        new_centre.set_mass_density(0.125);                             // units kg/m^3
                        new_centre.set_velocity(0.0);                                   // units m/s
                        new_centre.set_pressure(100.0);                                 // units N/m^2
                }
                new_centre.setup_specific_energy();
                new_centre.prim_to_con();
                new_centre.setup_f_variables();
                new_centre.reset_du();

                return new_centre;

        }else if(IC == 1){
                if(i==0){cout << "Using 1D Sine Wave" << endl;}

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

                new_centre.set_mass_density(rho);                       // units kg/m^3
                new_centre.set_velocity(v);                             // units m/s
                new_centre.set_pressure(p);                             // units N/m^2

                new_centre.setup_specific_energy();
                new_centre.prim_to_con();
                new_centre.setup_f_variables();
                new_centre.reset_du();

                return new_centre;

        }else if(IC == 2){

                if(i==0){cout << "Using 1D Sedov Blast Wave" << endl;}
                new_centre.set_x(x);
    
                new_centre.set_mass_density(1);                        // units kg/m^3
                new_centre.set_velocity(0);                            // units m/s
                new_centre.set_pressure(1);                            // units N/m^2
    
                if(x > (SIDE_LENGTH/2 - dx) and x < (SIDE_LENGTH/2 + dx)){new_centre.set_pressure(100000/dx);}

                new_centre.setup_specific_energy();
                new_centre.prim_to_con();
                new_centre.setup_f_variables();
                new_centre.reset_du();

                return new_centre;

        }else if(IC==3){
                if(i==0){cout << "Using 1D Gaussian Pulse" << endl;}
                double centre = 25.0;
                double s,w,rho,rho_0 = 10.0,rho_pulse = 50.0;
                double velocity = 1.0,pressure = 1000.0;

                s = abs(centre-x);                                       // distance from centre of pulse
                w = 2.0;                                                 // characteristic width

                rho = rho_pulse*exp(-s*s/(w*w)) + rho_0*(1-exp(-s*s/(w*w)));

                new_centre.set_x(x);

                new_centre.set_mass_density(rho);
                new_centre.set_velocity(velocity);
                new_centre.set_pressure(pressure);

                new_centre.setup_specific_energy();
                new_centre.prim_to_con();
                new_centre.setup_f_variables();
                new_centre.reset_du();

                return new_centre;

        }else if(IC==4){
                if(i==0){cout << "Using 1D Square Wave" << endl;}

                new_centre.set_x(x);

                if(i>0.25*N_POINTS and i<0.75*N_POINTS){
                        new_centre.set_mass_density(0.032);                             // units kg/m^3
                        new_centre.set_velocity(1.0);                                   // units m/s
                        new_centre.set_pressure(10.0);                                 // units N/m^2
                }else{
                        new_centre.set_mass_density(0.0075);                            // units kg/m^3
                        new_centre.set_velocity(1.0);                                   // units m/s
                        new_centre.set_pressure(10.0);                                 // units N/m^2
                }
                new_centre.setup_specific_energy();
                new_centre.prim_to_con();
                new_centre.setup_f_variables();
                new_centre.reset_du();

                return new_centre;

        }else if(IC==5){
                if(i==0){cout << "Using 1D Steady Flow" << endl;}

                new_centre.set_x(x);

                new_centre.set_mass_density(1.0);                               // units kg/m^3
                new_centre.set_velocity(1.0);                                   // units m/s
                new_centre.set_pressure(1.0);                                   // units N/m^2

                new_centre.setup_specific_energy();
                new_centre.prim_to_con();
                new_centre.setup_f_variables();
                new_centre.reset_du();

                return new_centre;

        }else if(IC==6){
                if(i==0){cout << "Using 1D Linear Wave Test (ATHENA)" << endl;}

                double A,R[3],PI;

                PI = 3.141529;

                A = 1.0e-3;

                R[0] = 1.0;
                R[1] = 1.0;
                R[2] = 1.5;

                new_centre.set_x(x);

                new_centre.set_mass_density(1.0 + (A*R[0]*sin(2.0*PI*x)));                  // units kg/m^3
                new_centre.set_velocity(1.0 + (A*R[1]*sin(2.0*PI*x)));                      // units m/s
                //new_centre.set_pressure(1.0/GAMMA);                                         // units N/m^2

                new_centre.setup_specific_energy();
                new_centre.prim_to_con();
                new_centre.setup_f_variables();
                new_centre.reset_du();

                return new_centre;

        }else{
                if(i==0){cout << "Using 1D Interacting Blastwaves Test (ATHENA)" << endl;}

                new_centre.set_x(x);

                new_centre.set_mass_density(1.0);                             // units kg/m^3
                new_centre.set_velocity(0.0);                                 // units m/s

                if(x > 0.1*SIDE_LENGTH and x < 0.2*SIDE_LENGTH){
                        new_centre.set_pressure(1000);                                 // units N/m^2
                }else if(x > 0.8*SIDE_LENGTH and x < 0.9*SIDE_LENGTH){
                        new_centre.set_pressure(100);                                 // units N/m^2
                }else{
                        new_centre.set_pressure(0.1);                                 // units N/m^2
                }

                new_centre.setup_specific_energy();
                new_centre.prim_to_con();
                new_centre.setup_f_variables();
                new_centre.reset_du();

                return new_centre;
                
        }
}

face setup_face(int i, vector<centre> &points){
        int centre_id_0,centre_id_1;                                    // centre_id_0 and centre_id_1 = index number of cells on either side of face
        centre *centre_0,*centre_1,*centre_00,*centre_11;               // *centre_0 and *centre_1 = pointers to vertices
        face new_face;

        centre_id_0 = i % N_POINTS;                                     // setup preiodic boundary
        centre_id_1 = (i+1) % N_POINTS;

        centre_0 = &points[centre_id_0];                                // setup pointers to lower and upper centre
        centre_1 = &points[centre_id_1];
        centre_00 = &points[centre_id_0-1];                             // add neighbouring vertices for flux limiter
        centre_11 = &points[centre_id_1+1];

        if(centre_id_0 == 0){centre_00 = &points[N_POINTS-1];}
        if(centre_id_1 == N_POINTS-1){centre_11 = &points[0];}

        new_face.set_centre_0(centre_0);                                // pass these pointers to the face
        new_face.set_centre_1(centre_1);
        new_face.set_centre_00(centre_00);
        new_face.set_centre_11(centre_11);

        return new_face;
}