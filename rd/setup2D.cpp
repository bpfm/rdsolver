VERTEX setup_vertex(int i, int j, double &DX, double &DY){
        VERTEX NEW_VERTEX;
        double X,Y;

        DX = SIDE_LENGTH_X/double(N_POINTS_X);
        DY = SIDE_LENGTH_Y/double(N_POINTS_Y);//sqrt(3.0)*DX/2.0;

        if((j % 2)== 0){
                X=double(i)*DX;
                Y=double(j)*DY;
        }else{
                X=(double(i)+0.5)*DX;
                Y=double(j)*DY;
        }

        NEW_VERTEX.set_x(X);
        NEW_VERTEX.set_y(Y);

        NEW_VERTEX.set_dx(DX);
        NEW_VERTEX.set_dy(DY);

        NEW_VERTEX.calculate_dual();

        if(IC == 0){
                if(i==0 and j==0){std::cout << "Using 2D Sod Shock Tube" << std::endl;}

                if(i<0.5*N_POINTS_X){
                        NEW_VERTEX.set_mass_density(1.0);                               // units kg/m^3
                        NEW_VERTEX.set_x_velocity(0.00000000001);                       // units m/s
                        NEW_VERTEX.set_y_velocity(0.00000000001);                       // units m/s
                        NEW_VERTEX.set_pressure(500.0);                                 // units N/m^2
                }else{
                        NEW_VERTEX.set_mass_density(0.125);                             // units kg/m^3
                        NEW_VERTEX.set_x_velocity(0.00000000001);                       // units m/s
                        NEW_VERTEX.set_y_velocity(0.00000000001);                       // units m/s
                        NEW_VERTEX.set_pressure(100.0);                                 // units N/m^2
                }
        }else if(IC == 1){
                if(i==0 and j==0){std::cout << "Using 2D Sine Wave" << std::endl;}

                //std::cout << "x =\t" << X << "\ty =\t" << Y << std::endl;

                double RHO,RHO_0 = 50.0;
                double P,P_0 = 3.0;
                double C_S = sqrt(GAMMA*P_0/RHO_0);
                double KB = 2.0*3.1415/SIDE_LENGTH_X;
                double EPSILON = 0.1;
                double V;

                RHO = RHO_0*(1.0 + EPSILON * sin(KB * X));
                P = P_0 + C_S * C_S * (RHO - RHO_0);
                V = C_S * (RHO - RHO_0)/RHO_0 + 0.000000001;

                NEW_VERTEX.set_mass_density(RHO);                       // units kg/m^3
                NEW_VERTEX.set_x_velocity(V);                             // units m/s
                NEW_VERTEX.set_y_velocity(0.00000001);
                NEW_VERTEX.set_pressure(P);                             // units N/m^2
        }else if(IC == 2){
                if(i==0 and j==0){std::cout << "Using 2D Sedov Blast" << std::endl;}

                double RHO = 1000.0;
                double V = 0.00000001;
                double P = 100.0;

                if((X > (SIDE_LENGTH_X/2 - DX) and X < (SIDE_LENGTH_X/2 + DX)) and (Y > (SIDE_LENGTH_Y/2 - DY) and Y < (SIDE_LENGTH_Y/2 + DY))){
                        P = 100000/DX;
                        std::cout << "Setting blast pressure point" << std::endl;
                }

                NEW_VERTEX.set_mass_density(RHO);                       // units kg/m^3
                NEW_VERTEX.set_x_velocity(V);                             // units m/s
                NEW_VERTEX.set_y_velocity(V);
                NEW_VERTEX.set_pressure(P);                             // units N/m^2
        }else if(IC == 3){
                if(i==0 and j==0){std::cout << "Using 2D Gaussian pulse" << std::endl;}

                double CENTRE = 10.0;
                double S,W,RHO,RHO_0 = 10.0,RHO_PULSE = 50.0;
                double X_VELOCITY = 1.0,PRESSURE = 1000.0;

                S = abs(CENTRE - X);                                       // distance from centre of pulse
                W = 2.0;                                                 // characteristic width

                RHO = RHO_PULSE*exp(-S*S/(W*W)) + RHO_0*(1.0-exp(-S*S/(W*W)));

                NEW_VERTEX.set_mass_density(RHO);
                NEW_VERTEX.set_x_velocity(X_VELOCITY);
                NEW_VERTEX.set_y_velocity(0.00000001);
                NEW_VERTEX.set_pressure(PRESSURE);
        }else if(IC == 4){
                if(i==0 and j==0){std::cout << "Using 2D Sod Shock Tube (Varied in Y)" << std::endl;}

                if(j<0.5*N_POINTS_Y){
                        NEW_VERTEX.set_mass_density(1.0);                               // units kg/m^3
                        NEW_VERTEX.set_x_velocity(0.00000001);                                 // units m/s
                        NEW_VERTEX.set_y_velocity(0.00000001);                                 // units m/s
                        NEW_VERTEX.set_pressure(500.0);                                 // units N/m^2
                }else{
                        NEW_VERTEX.set_mass_density(0.125);                             // units kg/m^3
                        NEW_VERTEX.set_x_velocity(0.00000001);                                 // units m/s
                        NEW_VERTEX.set_y_velocity(0.00000001);                                 // units m/s
                        NEW_VERTEX.set_pressure(100.0);                                 // units N/m^2
                }
        }

        NEW_VERTEX.setup_specific_energy();
        NEW_VERTEX.prim_to_con();
        NEW_VERTEX.reset_du_half();
        NEW_VERTEX.reset_du();

        return NEW_VERTEX;
}

TRIANGLE setup_triangle(int i0, int j0, std::vector<std::vector<VERTEX> > &POINTS){
        int i,j;
        int VERTEX_I_ID_0,VERTEX_J_ID_0,VERTEX_I_ID_1,VERTEX_J_ID_1,VERTEX_I_ID_2,VERTEX_J_ID_2;            // VERTEX_id_0 and VERTEX_id_1 = index number of cells on either side of TRIANGLE
        VERTEX *VERTEX_0,*VERTEX_1,*VERTEX_2;                                                               // *VERTEX_0, *VERTEX_1 and *VERTEX_2 = pointers to vertices (labelled anticlockwise)
        TRIANGLE NEW_TRIANGLE;

        if((j0 % 2) == 0){
                i = i0;
                j = j0/2;
                VERTEX_I_ID_0 = i % N_POINTS_X;
                VERTEX_J_ID_0 = j % N_POINTS_Y;
                VERTEX_I_ID_1 = (i+1) % N_POINTS_X;
                VERTEX_J_ID_1 = j % N_POINTS_Y;
                VERTEX_I_ID_2 = i % N_POINTS_X;
                VERTEX_J_ID_2 = (j+1) % N_POINTS_Y;
        }else{
                i = i0;
                j = (j0-1)/2;
                VERTEX_I_ID_0 = (i+1) % N_POINTS_X;
                VERTEX_J_ID_0 = j % N_POINTS_Y;
                VERTEX_I_ID_1 = (i+1) % N_POINTS_X;
                VERTEX_J_ID_1 = (j+1) % N_POINTS_Y;
                VERTEX_I_ID_2 = i % N_POINTS_X;
                VERTEX_J_ID_2 = (j+1) % N_POINTS_Y;
        }

        VERTEX_0 = &POINTS[VERTEX_J_ID_0][VERTEX_I_ID_0];
        VERTEX_1 = &POINTS[VERTEX_J_ID_1][VERTEX_I_ID_1];
        VERTEX_2 = &POINTS[VERTEX_J_ID_2][VERTEX_I_ID_2];

#ifdef DEBUG
        std::cout << "Setting up\t" << i0 << "\t" << j0 << "\t(x,y) =\t" << VERTEX_0->get_x() << "\t" << VERTEX_0->get_y() << "\t" << VERTEX_1->get_x() << "\t" << VERTEX_1->get_y() << "\t" << VERTEX_2->get_x() << "\t" << VERTEX_2->get_y() << std::endl;
#endif

        NEW_TRIANGLE.set_vertex_0(VERTEX_0);                                // pass these pointers to the TRIANGLE
        NEW_TRIANGLE.set_vertex_1(VERTEX_1);
        NEW_TRIANGLE.set_vertex_2(VERTEX_2);

        return NEW_TRIANGLE;
}