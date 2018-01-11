using namespace std;

VERTEX setup_vertex(int N_POINTS, int i, int j){
        VERTEX NEW_VERTEX;
        double X,Y,DX,DY;

        DX = SIDE_LENGTH/double(N_POINTS);
        DY = sqrt(3.0)*DX/2.0;

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
                        NEW_VERTEX.set_x_velocity(0.0);                                 // units m/s
                        NEW_VERTEX.set_y_velocity(0.0);                                 // units m/s
                        NEW_VERTEX.set_pressure(500.0);                                 // units N/m^2
                }else{
                        NEW_VERTEX.set_mass_density(0.125);                             // units kg/m^3
                        NEW_VERTEX.set_x_velocity(0.0);                                 // units m/s
                        NEW_VERTEX.set_y_velocity(0.0);                                 // units m/s
                        NEW_VERTEX.set_pressure(100.0);                                 // units N/m^2
                }
                NEW_VERTEX.setup_specific_energy();
                NEW_VERTEX.prim_to_con();
                NEW_VERTEX.reset_du();

                return NEW_VERTEX;

        }else if(IC == 1){
                if(i==0 and j==0){cout << "Using 2D Sine Wave" << endl;}

                //cout << "x =\t" << X << "\ty =\t" << Y << endl;

                double RHO,RHO_0 = 50.0;
                double P,P_0 = 3.0;
                double C_S = sqrt(GAMMA*P_0/RHO_0);
                double KB = 2.0*3.1415/SIDE_LENGTH;
                double EPSILON = 0.1;
                double V;

                NEW_VERTEX.set_x(X);
                NEW_VERTEX.set_y(Y);

                RHO = RHO_0*(1.0 + EPSILON * sin(KB * X));
                P = P_0 + C_S * C_S * (RHO - RHO_0);
                V = C_S * (RHO - RHO_0)/RHO_0;

                NEW_VERTEX.set_mass_density(RHO);                       // units kg/m^3
                NEW_VERTEX.set_x_velocity(V);                             // units m/s
                NEW_VERTEX.set_y_velocity(0.0);
                NEW_VERTEX.set_pressure(P);                             // units N/m^2

                NEW_VERTEX.setup_specific_energy();
                NEW_VERTEX.prim_to_con();
                NEW_VERTEX.reset_du();

                return NEW_VERTEX;

        }
}

TRIANGLE setup_triangle(int i, int j0, vector<vector<VERTEX> > &POINTS, ofstream &positions){
        int j;
        int VERTEX_I_ID_0,VERTEX_J_ID_0,VERTEX_I_ID_1,VERTEX_J_ID_1,VERTEX_I_ID_2,VERTEX_J_ID_2;            // VERTEX_id_0 and VERTEX_id_1 = index number of cells on either side of TRIANGLE
        VERTEX *VERTEX_0,*VERTEX_1,*VERTEX_2;                                                               // *VERTEX_0, *VERTEX_1 and *VERTEX_2 = pointers to vertices (labelled anticlockwise)
        TRIANGLE NEW_TRIANGLE;

        if((j0 % 2) == 0){
                j = j0/2;
                VERTEX_I_ID_0 = j;
                VERTEX_J_ID_0 = i;
                VERTEX_I_ID_1 = j+1;
                VERTEX_J_ID_1 = i;
                VERTEX_I_ID_2 = j;
                VERTEX_J_ID_2 = i+1;
        }else{
                j = (j0-1)/2;
                VERTEX_I_ID_0 = j+1;
                VERTEX_J_ID_0 = i;
                VERTEX_I_ID_1 = j+1;
                VERTEX_J_ID_1 = i+1;
                VERTEX_I_ID_2 = j;
                VERTEX_J_ID_2 = i+1;
        }

        VERTEX_0 = &POINTS[VERTEX_I_ID_0][VERTEX_J_ID_0];
        VERTEX_1 = &POINTS[VERTEX_I_ID_1][VERTEX_J_ID_1];
        VERTEX_2 = &POINTS[VERTEX_I_ID_2][VERTEX_J_ID_2];

        NEW_TRIANGLE.set_vertex_0(VERTEX_0);                                // pass these pointers to the TRIANGLE
        NEW_TRIANGLE.set_vertex_1(VERTEX_1);
        NEW_TRIANGLE.set_vertex_2(VERTEX_2);

        positions << j << "\t" << i << "\t" << NEW_TRIANGLE.get_vertex_0()->get_x() << "\t" << NEW_TRIANGLE.get_vertex_0()->get_y() << endl;
        positions << j << "\t" << i << "\t" << NEW_TRIANGLE.get_vertex_1()->get_x() << "\t" << NEW_TRIANGLE.get_vertex_1()->get_y() << endl;
        positions << j << "\t" << i << "\t" << NEW_TRIANGLE.get_vertex_2()->get_x() << "\t" << NEW_TRIANGLE.get_vertex_2()->get_y() << endl;
        positions << j << "\t" << i << "\t" << "\t" << endl;

        return NEW_TRIANGLE;
}