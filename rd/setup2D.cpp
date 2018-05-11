using namespace std;

VERTEX setup_vertex(int N_POINTS, int i, int j, double &DX, double &DY){
        VERTEX NEW_VERTEX;
        double X,Y;

        DX = SIDE_LENGTH/double(N_POINTS);
        DY = SIDE_LENGTH/double(N_POINTS);//sqrt(3.0)*DX/2.0;

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
                if(i==0 and j==0){cout << "Using 2D Sod Shock Tube" << endl;}

                //cout << "x =\t" << X << "\ty =\t" << Y << endl;

                if(i>0.3*N_POINTS and i<0.7*N_POINTS){
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
        }else if(IC == 1){
                if(i==0 and j==0){cout << "Using 2D Sine Wave" << endl;}

                //cout << "x =\t" << X << "\ty =\t" << Y << endl;

                double RHO,RHO_0 = 50.0;
                double P,P_0 = 3.0;
                double C_S = sqrt(GAMMA*P_0/RHO_0);
                double KB = 2.0*3.1415/SIDE_LENGTH;
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
                if(i==0 and j==0){cout << "Using 2D Sedov Blast" << endl;}

                double RHO = 100.0;
                double V = 0.00000001;
                double P = 1.0;

                if((X > (50/2 - 2.0*DX) and X < (50/2 + 2.0*DX)) and (Y > (50/2 - 2.0*DY) and Y < (50/2 + 2.0*DY))){
                        P = 100000/DX;
                        cout << "Setting blast pressure point" << endl;
                }

                NEW_VERTEX.set_mass_density(RHO);                       // units kg/m^3
                NEW_VERTEX.set_x_velocity(V);                             // units m/s
                NEW_VERTEX.set_y_velocity(V);
                NEW_VERTEX.set_pressure(P);                             // units N/m^2
        }

        NEW_VERTEX.setup_specific_energy();
        NEW_VERTEX.prim_to_con();
        NEW_VERTEX.reset_du_half();
        NEW_VERTEX.reset_du();

        return NEW_VERTEX;
}

TRIANGLE setup_triangle(int i, int j0, vector<vector<VERTEX> > &POINTS){
        int j;
        int VERTEX_I_ID_0,VERTEX_J_ID_0,VERTEX_I_ID_1,VERTEX_J_ID_1,VERTEX_I_ID_2,VERTEX_J_ID_2;            // VERTEX_id_0 and VERTEX_id_1 = index number of cells on either side of TRIANGLE
        VERTEX *VERTEX_0,*VERTEX_1,*VERTEX_2;                                                               // *VERTEX_0, *VERTEX_1 and *VERTEX_2 = pointers to vertices (labelled anticlockwise)
        TRIANGLE NEW_TRIANGLE;

        if((j0 % 2) == 0){
                j = j0/2;
                VERTEX_I_ID_2 = j % N_POINTS;
                VERTEX_J_ID_2 = i % N_POINTS;
                VERTEX_I_ID_1 = (j+1) % N_POINTS;
                VERTEX_J_ID_1 = i % N_POINTS;
                VERTEX_I_ID_0 = j % N_POINTS;
                VERTEX_J_ID_0 = (i+1) % N_POINTS;
        }else{
                j = (j0-1)/2;
                VERTEX_I_ID_2 = (j+1) % N_POINTS;
                VERTEX_J_ID_2 = i % N_POINTS;
                VERTEX_I_ID_1 = (j+1) % N_POINTS;
                VERTEX_J_ID_1 = (i+1) % N_POINTS;
                VERTEX_I_ID_0 = j % N_POINTS;
                VERTEX_J_ID_0 = (i+1) % N_POINTS;
        }

        //cout << "IDs =\t" << VERTEX_I_ID_1 << "\t" << VERTEX_J_ID_2 << endl;

        VERTEX_0 = &POINTS[VERTEX_I_ID_0][VERTEX_J_ID_0];
        VERTEX_1 = &POINTS[VERTEX_I_ID_1][VERTEX_J_ID_1];
        VERTEX_2 = &POINTS[VERTEX_I_ID_2][VERTEX_J_ID_2];

        NEW_TRIANGLE.set_vertex_0(VERTEX_0);                                // pass these pointers to the TRIANGLE
        NEW_TRIANGLE.set_vertex_1(VERTEX_1);
        NEW_TRIANGLE.set_vertex_2(VERTEX_2);

        return NEW_TRIANGLE;
}