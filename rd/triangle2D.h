/* class containing values associated with the face
        *VERTEX_0 = pointer to left VERTEX
        *VERTEX_1 = pointer to right VERTEX
*/

using namespace std;

class TRIANGLE{

private:

        int ID;
        VERTEX *VERTEX_0,*VERTEX_1,*VERTEX_2;

public:

        void set_id(int NEW_ID){
                ID = NEW_ID;
        }

        void set_vertex_0(VERTEX* NEW_VERTEX){
                VERTEX_0 = NEW_VERTEX;
        }

        void set_vertex_1(VERTEX* NEW_VERTEX){
                VERTEX_1 = NEW_VERTEX;
        }

        void set_vertex_2(VERTEX* NEW_VERTEX){
                VERTEX_2 = NEW_VERTEX;
        }

        VERTEX* get_vertex_0(){
                return VERTEX_0;
        }

        VERTEX* get_vertex_1(){
                return VERTEX_1;
        }

        VERTEX* get_vertex_2(){
                return VERTEX_2;
        }

        void calculate_change(double T){
                int i,j,k;
                double X[3],Y[3],X_VEL[3],Y_VEL[3],U_N[4][3],U_HALF[4][3];
                double DU0[4],DU1[4],DU2[4];
                double NORMAL[3][2];
                double LAMBDA[4][3][2];
                double INFLOW[4][3],INFLOW_PLUS[4][3],INFLOW_MINUS[4][3];
                double U_IN[4],U_OUT[4];
                double IN_TOP,IN_BOTTOM,OUT_TOP,OUT_BOTTOM,PRESSURE,C_SOUND;
                double FLUC[4][3],BETA[4][3];
                int INOUT[4][3];

                // Import conditions and positions of vertices

                X[0] = VERTEX_0->get_x();
                X[1] = VERTEX_1->get_x();
                X[2] = VERTEX_2->get_x();

                Y[0] = VERTEX_0->get_y();
                Y[1] = VERTEX_1->get_y();
                Y[2] = VERTEX_2->get_y();

                X_VEL[0] = VERTEX_0->get_x_velocity();
                X_VEL[1] = VERTEX_1->get_x_velocity();
                X_VEL[2] = VERTEX_2->get_x_velocity();

                Y_VEL[0] = VERTEX_0->get_y_velocity();
                Y_VEL[1] = VERTEX_1->get_y_velocity();
                Y_VEL[2] = VERTEX_2->get_y_velocity();

                U_N[0][0] = VERTEX_0->get_u0();
                U_N[0][1] = VERTEX_1->get_u0();
                U_N[0][2] = VERTEX_2->get_u0();

                U_N[1][0] = VERTEX_0->get_u1();
                U_N[1][1] = VERTEX_1->get_u1();
                U_N[1][2] = VERTEX_2->get_u1();

                U_N[2][0] = VERTEX_0->get_u2();
                U_N[2][1] = VERTEX_1->get_u2();
                U_N[2][2] = VERTEX_2->get_u2();

                U_N[3][0] = VERTEX_0->get_u3();
                U_N[3][1] = VERTEX_1->get_u3();
                U_N[3][2] = VERTEX_2->get_u3();

                if(DEBUG==1){
                        cout << "---------------------------------------------------------" << endl;
                        cout << "i =\t" << X[0] << "\t" << Y[0] << "\tj =\t" << X[1] << "\t" << Y[1] << "\tk =\t" << X[2] << "\t" << Y[2] << endl;
                        cout << "State =" << "\t\t rho" << "\t\t x_mom" << "\t\t y_mom" << "\t\t energy" << endl;
                        for(i=0;i<3;i++){
                                cout << i << " =\t" << U_N[0][i] << "\t" << U_N[1][i] << "\t" << U_N[2][i] << "\t" << U_N[3][i] << endl;
                        }
                }

                // Placeholder pressure calculation

                PRESSURE = VERTEX_0->get_pressure();
                C_SOUND = sqrt(GAMMA*PRESSURE/U_N[0][0]);

                if(T==0.0){caclulate_normals(X,Y,NORMAL[0][0],NORMAL[0][1],NORMAL[1][0],NORMAL[1][1],NORMAL[2][0],NORMAL[2][1]);}

                for(i=0;i<3;i++){
                        LAMBDA[0][i][0] = X_VEL[i] - C_SOUND;
                        LAMBDA[1][i][0] = X_VEL[i];
                        LAMBDA[2][i][0] = X_VEL[i];
                        LAMBDA[3][i][0] = X_VEL[i] + C_SOUND;

                        LAMBDA[0][i][1] = Y_VEL[i] - C_SOUND;
                        LAMBDA[1][i][1] = Y_VEL[i];
                        LAMBDA[2][i][1] = Y_VEL[i];
                        LAMBDA[3][i][1] = Y_VEL[i] + C_SOUND;
                }

                // Take dot product of advection velocity and the normal of the opposite face to get the inflow parameters

                for(i=0;i<4;i++){
                        for(j=0;j<3;j++){
                                INOUT[i][j] = 0;
                                INFLOW[i][j] = 0.0;
                                for(k=0;k<2;k++){INFLOW[i][j] = INFLOW[i][j] + 0.5*LAMBDA[i][j][k]*NORMAL[j][k];}
                                INFLOW_PLUS[i][j] = max_val(0,INFLOW[i][j]);
                                INFLOW_MINUS[i][j] = min_val(0,INFLOW[i][j]);
                                if(INFLOW[i][j]>0){INOUT[i][j]=1;}                                        // Work out which are inflow vertices and which are outflow vertices
                        }
                }

                // Calculate inflow and outflow state of the element

                for(i=0;i<4;i++){
                        IN_TOP = 0.0;
                        IN_BOTTOM = 0.0;
                        OUT_TOP = 0.0;
                        OUT_BOTTOM = 0.0;
                        for(j=0;j<3;j++){
                                IN_TOP = IN_TOP + INFLOW_MINUS[i][j]*U_N[i][j];
                                IN_BOTTOM = IN_BOTTOM + INFLOW_MINUS[i][j];
                                OUT_TOP = OUT_TOP + INFLOW_PLUS[i][j]*U_N[i][j];
                                OUT_BOTTOM = OUT_BOTTOM + INFLOW_PLUS[i][j];
                        }
                        if(DEBUG==1){cout << "Sums (" << i << ") =\t" << IN_TOP << "\t" << IN_BOTTOM << "\t" << OUT_TOP << "\t" << OUT_BOTTOM << endl;}
                        
                        if(OUT_BOTTOM != 0.0){
                                U_IN[i] = IN_TOP/IN_BOTTOM;
                                U_OUT[i] = OUT_TOP/OUT_BOTTOM;
                                for(j=0;j<3;j++){BETA[i][j] = INFLOW_PLUS[i][j]/OUT_BOTTOM;}
                        }else{
                                U_IN[i] = 0.0;
                                U_OUT[i] = 0.0;
                                for(j=0;j<3;j++){BETA[i][j] = 0.0;}
                        }
                        if(DEBUG==1){cout << "k+ =\t" << INFLOW_PLUS[i][0] << "\t" << INFLOW_PLUS[i][1] << "\t" << INFLOW_PLUS[i][2] << endl;}
                }

                for(i=0;i<4;i++){
                        for(j=0;j<3;j++){
                                FLUC[i][j] = INFLOW_PLUS[i][j]*(U_OUT[i]-U_IN[i]);
                        }
                }

                for(i=0;i<4;i++){
                        DU0[i] = -1.0*FLUC[i][0];//-1.0*BETA[i][0]*FLUC[i][0];
                        DU1[i] = -1.0*FLUC[i][1];//-1.0*BETA[i][1]*FLUC[i][1];
                        DU2[i] = -1.0*FLUC[i][2];//-1.0*BETA[i][2]*FLUC[i][2];
                }

                VERTEX_0->update_du(DU0);
                VERTEX_1->update_du(DU1);
                VERTEX_2->update_du(DU2);

                if(DEBUG==1){
                        if(U_IN[0] != U_OUT[0]){
                                for(i=0;i<4;i++){cout << "u_in =\t" << U_IN[i] << "\tu_out =\t" << U_OUT[i] << endl;}
                                for(i=0;i<4;i++){cout << "Element fluctuation =\t" << FLUC[i][0] << "\t" << FLUC[i][1] << "\t" << FLUC[i][2] << endl;}
                                for(i=0;i<4;i++){cout << "Beta (" << i << ") =\t" << BETA[i][0] << "\t" << BETA[i][1] << "\t" << BETA[i][2] << "\tTotal =\t" << BETA[i][0]+BETA[i][1]+BETA[i][2] << endl;}
                                cout << "Change (rho) =\t" << DU0[0] << "\t" << DU1[0] << "\t" << DU2[0] << endl;
                                cout << "Change (x mom) =\t" << DU0[1] << "\t" << DU1[1] << "\t" << DU2[1] << endl;
                                cout << "Change (y mom) =\t" << DU0[2] << "\t" << DU1[2] << "\t" << DU2[2] << endl;
                                cout << "Change (energy) =\t" << DU0[3] << "\t" << DU1[3] << "\t" << DU2[3] << endl;
                                cout << "---------------------------------------------------------" << endl;
                                exit(0);
                        }
                }
                return ;
        }

        // returns Roe average of left and right states
        double roe_avg(double L1, double L2, double R1, double R2){
                double AVG;
                AVG = (sqrt(L1)*L2+sqrt(R1)*R2)/(sqrt(L1)+sqrt(R1));
                return AVG;
        }

        void caclulate_normals(double X[3],double Y[3], double &NORMAL00, double &NORMAL01, double &NORMAL10, double &NORMAL11, double &NORMAL20, double &NORMAL21){
                int i;
                double PERP[3][2],MAG,NORMAL[3][2];

                PERP[0][0] = (Y[1] - Y[2]);
                PERP[0][1] = (X[2] - X[1]);

                PERP[1][0] = (Y[2] - Y[0]);
                PERP[1][1] = (X[0] - X[2]);

                PERP[2][0] = (Y[0] - Y[1]);
                PERP[2][1] = (X[1] - X[0]);

                for(i=0;i<3;i++){
                        MAG = sqrt(PERP[i][0]*PERP[i][0]+PERP[i][1]*PERP[i][1]);
                        NORMAL[i][0] = -1.0*PERP[i][0]/MAG;
                        NORMAL[i][1] = -1.0*PERP[i][1]/MAG;
                }

                NORMAL00 = NORMAL[0][0];
                NORMAL01 = NORMAL[0][1];
                NORMAL10 = NORMAL[1][0];
                NORMAL11 = NORMAL[1][1];
                NORMAL20 = NORMAL[2][0];
                NORMAL21 = NORMAL[2][1];

                if(DEBUG==1){
                        cout << "X =\t" << X[0] << "\t" << X[1] << "\t" << X[2] << endl;
                        cout << "Y =\t" << Y[0] << "\t" << Y[1] << "\t" << Y[2] << endl;
                        cout << "Normal i =\t" << NORMAL00 << "\t" << NORMAL01 << endl;
                        cout << "Normal j =\t" << NORMAL10 << "\t" << NORMAL11 << endl;
                        cout << "Normal k =\t" << NORMAL20 << "\t" << NORMAL21 << endl;
                }

                return ;
        }

        double max_val(double A, double B){
                if(A>B){
                        return A;
                }else{
                        return B;
                }
        }

        double min_val(double A, double B){
                if(A<B){
                        return A;
                }else{
                        return B;
                }
        }

};