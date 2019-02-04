using namespace std;

centre setup_centre(int N_POINTS, int i, int j, int k, float dx){
        centre new_centre;

        double x = (double(i) + 0.5) * dx;
        double y = (double(j) + 0.5) * dx;
        double z = (double(k) + 0.5) * dx;

        if(IC == 0){
                if(i==0 and j==0 and k==0){cout << "Using 3D Sod Shock Tube" << endl;}

                new_centre.set_x(x);
                new_centre.set_y(y);
                new_centre.set_z(z);

                if(x >= 0.3*SIDE_LENGTH and x <= 0.7*SIDE_LENGTH){
                        new_centre.set_mass_density(1.0);                               // units kg/m^3
                        new_centre.set_pressure(500.0);                                 // units N/m^2
                }else{
                        new_centre.set_mass_density(0.125);                             // units kg/m^3
                        new_centre.set_pressure(100.0);                                 // units N/m^2
                }

                new_centre.set_x_velocity(0.0);                                   // units m/s
                new_centre.set_y_velocity(0.0);
                new_centre.set_z_velocity(0.0);

                new_centre.setup_specific_energy();
                new_centre.prim_to_con();
                new_centre.setup_f_variables();
                new_centre.reset_du();

                //cout << x << "\t" << y << "\t" << z << "\t" << new_centre.get_mass_density() << endl;

                return new_centre;
        }else if(IC == 1){
                if(i==0 and j==0 and k==0){cout << "Using 3D Sine Wave" << endl;}

                double rho,rho_0 = 50.0;
                double p,p_0 = 3.0;
                double c_s = sqrt(GAMMA*p_0/rho_0);
                double k = 2.0*3.1415/SIDE_LENGTH;
                double epsilon = 0.0001;
                double v;

                new_centre.set_x(x);
                new_centre.set_y(y);
                new_centre.set_z(z);

                rho = rho_0*(1.0 + epsilon * cos(k * x));
                p = p_0 + c_s * c_s * (rho - rho_0);
                v = c_s * (rho - rho_0)/rho_0;

                new_centre.set_mass_density(rho);                       // units kg/m^3
                new_centre.set_x_velocity(v);                           // units m/s
                new_centre.set_y_velocity(0.0);                         // units m/s
                new_centre.set_z_velocity(0.0);                         // units m/s
                new_centre.set_pressure(p);                             // units N/m^2

                new_centre.setup_specific_energy();
                new_centre.prim_to_con();
                new_centre.setup_f_variables();
                new_centre.reset_du();

                return new_centre;
        }else if(IC == 2){
                if(i==0 and j==0 and k==0){cout << "Using 3D Sedov Blast Wave" << endl;}

                new_centre.set_x(x);
                new_centre.set_y(y);
                new_centre.set_z(z);
    
                new_centre.set_mass_density(1.0);                        // units kg/m^3
                new_centre.set_x_velocity(0.0);                          // units m/s
                new_centre.set_y_velocity(0.0);                          // units m/s
                new_centre.set_z_velocity(0.0);                          // units m/s
                new_centre.set_pressure(1.0);                            // units N/m^2
    
                if(x > (SIDE_LENGTH/2.0 - dx) and x < (SIDE_LENGTH/2.0 + dx) and y > (SIDE_LENGTH/2.0 - dx) and y < (SIDE_LENGTH/2.0 + dx) and z > (SIDE_LENGTH/2.0 - dx) and z < (SIDE_LENGTH/2.0 + dx)){new_centre.set_pressure(100000/dx);new_centre.set_pressure(1000.0);}

                new_centre.setup_specific_energy();
                new_centre.prim_to_con();
                new_centre.setup_f_variables();
                new_centre.reset_du();

                return new_centre;
        }else if(IC == 3){
                if(i==0 and j==0 and k==0){cout << "Using 3D Gaussian Pulse" << endl;}
                double centre = 20.0;
                double s,w,rho,rho_0 = 10.0,rho_pulse = 50.0;
                double velocity = 5.0,pressure = 1000.0;

                s = abs(centre-x);                                       // distance from centre of pulse
                w = 2.0;                                                 // characteristic width

                rho = rho_pulse*exp(-s*s/(w*w)) + rho_0*(1-exp(-s*s/(w*w)));

                new_centre.set_x(x);
                new_centre.set_y(y);
                new_centre.set_z(z);

                new_centre.set_mass_density(rho);
                new_centre.set_x_velocity(velocity);
                new_centre.set_y_velocity(0.0);
                new_centre.set_z_velocity(0.0);
                new_centre.set_pressure(pressure);

                new_centre.setup_specific_energy();
                new_centre.prim_to_con();
                new_centre.setup_f_variables();
                new_centre.reset_du();

                return new_centre;
        }else{
                if(i==0 and j==0 and k==0){cout << "Using 3D Steady Flow" << endl;}

                new_centre.set_x(x);
                new_centre.set_y(y);
                new_centre.set_z(z);

                new_centre.set_mass_density(2.0);                       // units kg/m^3
                new_centre.set_x_velocity(20.0);                           // units m/s
                new_centre.set_y_velocity(0.0);                         // units m/s
                new_centre.set_z_velocity(0.0);                         // units m/s
                new_centre.set_pressure(2.0);                             // units N/m^2

                new_centre.setup_specific_energy();
                new_centre.prim_to_con();
                new_centre.setup_f_variables();
                new_centre.reset_du();

                return new_centre;
        }
}

face setup_face(int i, double dx, vector<centre> &points, string direction){
        int j,centre_id_0,centre_id_1;                                    // centre_id_0 and centre_id_1 = index number of cells on either side of face
        float x,y,z;
        centre *centre_0,*centre_1,*centre_00,*centre_11;               // *centre_0 and *centre_1 = pointers to vertices
        centre current_centre;
        face new_face;
        vector<centre>::iterator it_cent;

        current_centre = points[i];

        x = current_centre.get_x();
        y = current_centre.get_y();
        z = current_centre.get_z();

        if(direction == "z"){
                centre_id_0 = i;
                centre_id_1 = (i + 1);
                new_face.set_direction(direction);
                if(current_centre.get_z()+dx > SIDE_LENGTH){
                        for(it_cent=points.begin(),j=0;it_cent<points.end();it_cent++,j++){
                                if(points[j].get_z()==0.5*dx and points[j].get_y()==y and points[j].get_x()==x){centre_id_1=j;}
                        }
                }
                //cout << "\tz_face position =\t" << points[centre_id_0].get_z() << "\t" << points[centre_id_1].get_z() << endl;;
        }else if(direction == "y"){
                centre_id_0 = i;
                centre_id_1 = (i + N_POINTS);
                new_face.set_direction(direction);
                if(current_centre.get_y()+dx > SIDE_LENGTH){
                        for(it_cent=points.begin(),j=0;it_cent<points.end();it_cent++,j++){
                                if(points[j].get_z()==z and points[j].get_y()==0.5*dx and points[j].get_x()==x){centre_id_1=j;}
                        }
                }
                //cout << "\ty_face position =\t" << points[centre_id_0].get_y() << "\t" << points[centre_id_1].get_y();
        }else if(direction == "x"){
                centre_id_0 = i;
                centre_id_1 = (i + N_POINTS * N_POINTS);
                new_face.set_direction(direction);
                if(current_centre.get_x()+dx > SIDE_LENGTH){
                        for(it_cent=points.begin(),j=0;it_cent<points.end();it_cent++,j++){
                                if(points[j].get_z()==z and points[j].get_y()==y and points[j].get_x()==0.5*dx){centre_id_1=j;}
                        }
                }
                //cout << "x_face position =\t" << points[centre_id_0].get_x() << "\t" << points[centre_id_1].get_x();
        }

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