/* class containing values associated with the face
        *VERTEX_0 = pointer to left VERTEX
        *VERTEX_1 = pointer to right VERTEX
*/

using namespace std;

class TRIANGLE{

private:
        int ID;
        VERTEX *VERTEX_0,*VERTEX_1,*VERTEX_2;
        double AREA;

        double X[3],Y[3],DUAL[3];
        double NORMAL[3][2];

        double U_N[4][3];
        double U_HALF[4][3];

        double U_IN[4],U_OUT[4],BETA[4][3];
        double U_IN_HALF[4],U_OUT_HALF[4],BETA_HALF[4][3];

        double LAMBDA[4][3];
        double LAMBDA_HALF[4][3][2];

        double VEL[3][2];
        double VEL_HALF[3][2];

        double PRESSURE[3],C_SOUND[3];
        double PRESSURE_HALF[3],C_SOUND_HALF[3];

        double FLUC[4][3];
        double FLUC_HALF[4][3];

        double PHI[4];

public:

        void set_id(int NEW_ID){ID = NEW_ID;}

        void set_vertex_0(VERTEX* NEW_VERTEX){VERTEX_0 = NEW_VERTEX;}
        void set_vertex_1(VERTEX* NEW_VERTEX){VERTEX_1 = NEW_VERTEX;}
        void set_vertex_2(VERTEX* NEW_VERTEX){VERTEX_2 = NEW_VERTEX;}

        int get_id(){return ID;}

        VERTEX* get_vertex_0(){return VERTEX_0;}
        VERTEX* get_vertex_1(){return VERTEX_1;}
        VERTEX* get_vertex_2(){return VERTEX_2;}

        void setup_positions(){
                X[0] = VERTEX_0->get_x();
                X[1] = VERTEX_1->get_x();
                X[2] = VERTEX_2->get_x();

                Y[0] = VERTEX_0->get_y();
                Y[1] = VERTEX_1->get_y();
                Y[2] = VERTEX_2->get_y();
        }

        void setup_initial_state(){
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

                VEL[0][0] = VERTEX_0->get_x_velocity();
                VEL[1][0] = VERTEX_1->get_x_velocity();
                VEL[2][0] = VERTEX_2->get_x_velocity();

                VEL[0][1] = VERTEX_0->get_y_velocity();
                VEL[1][1] = VERTEX_1->get_y_velocity();
                VEL[2][1] = VERTEX_2->get_y_velocity();

                PRESSURE[0] = VERTEX_0->get_pressure();
                PRESSURE[1] = VERTEX_1->get_pressure();
                PRESSURE[2] = VERTEX_2->get_pressure();

                DUAL[0] = VERTEX_0->get_dual();
                DUAL[1] = VERTEX_1->get_dual();
                DUAL[2] = VERTEX_2->get_dual();
        }

        void setup_half_state(){
                U_HALF[0][0] = VERTEX_0->get_u0_half();
                U_HALF[0][1] = VERTEX_1->get_u0_half();
                U_HALF[0][2] = VERTEX_2->get_u0_half();

                U_HALF[1][0] = VERTEX_0->get_u1_half();
                U_HALF[1][1] = VERTEX_1->get_u1_half();
                U_HALF[1][2] = VERTEX_2->get_u1_half();

                U_HALF[2][0] = VERTEX_0->get_u2_half();
                U_HALF[2][1] = VERTEX_1->get_u2_half();
                U_HALF[2][2] = VERTEX_2->get_u2_half();

                U_HALF[3][0] = VERTEX_0->get_u3_half();
                U_HALF[3][1] = VERTEX_1->get_u3_half();
                U_HALF[3][2] = VERTEX_2->get_u3_half();

                VEL_HALF[0][0] = VERTEX_0->get_x_velocity_half();
                VEL_HALF[1][0] = VERTEX_1->get_x_velocity_half();
                VEL_HALF[2][0] = VERTEX_2->get_x_velocity_half();

                VEL_HALF[0][1] = VERTEX_0->get_y_velocity_half();
                VEL_HALF[1][1] = VERTEX_1->get_y_velocity_half();
                VEL_HALF[2][1] = VERTEX_2->get_y_velocity_half();

                PRESSURE_HALF[0] = VERTEX_0->get_pressure_half();
                PRESSURE_HALF[1] = VERTEX_1->get_pressure_half();
                PRESSURE_HALF[2] = VERTEX_2->get_pressure_half();
        }

        //**********************************************************************************************************************

        // Calculate first half timestep change, passing change to vertice
        void calculate_first_half(double T, double DT_TOT, double DX, double DY){
                int i,j,k,m,p;

                double DU0[4],DU1[4],DU2[4];

                double INFLOW[4][4][3][3];
                double INFLOW_PLUS_INVERSE[4][4],INFLOW_MINUS_INVERSE[4][4];
                double INFLOW_PLUS_SUM[4][4], INFLOW_MINUS_SUM[4][4];

                double IN_TOP[4],OUT_TOP[4];
                double DT = DT_TOT;

                // Import conditions and positions of vertices

                setup_positions();
                setup_initial_state();

                // Calculate normals (just in first timestep for static grid)

                // double X_MOD[3],Y_MOD[3];

                // if(T==0.0){
                //         for(m=0;m<3;++m){
                //                 X_MOD[m] = X[m];
                //                 Y_MOD[m] = Y[m];
                //                 if(X[m] + DX > 50.0){
                //                         X_MOD[m] = X[m] - SIDE_LENGTH;
                //                         return ;
                //                 }
                //                 if(Y[m] + DY > 50.0){
                //                         Y_MOD[m] = Y[m] - SIDE_LENGTH;
                //                         return ;
                //                 }
                //         }
                //         caclulate_normals(X_MOD,Y_MOD,NORMAL[0][0],NORMAL[0][1],NORMAL[1][0],NORMAL[1][1],NORMAL[2][0],NORMAL[2][1]);
                // }

                if(abs(X[0] - X[1]) > 10.0 or abs(X[0] - X[2]) > 10.0 or abs(X[1] - X[2]) > 10.0){
#ifdef DEBUG
                        cout << "Skipping x boundary" << endl;
#endif
                        return ;
                }else if(abs(Y[0] - Y[1]) > 10.0 or abs(Y[0] - Y[2]) > 10.0 or abs(Y[1] - Y[2]) > 10.0){
#ifdef DEBUG
                        cout << "Skipping y boundary" << endl;
#endif
                        return ;
                }

                if(T==0){caclulate_normals(X,Y,NORMAL[0][0],NORMAL[0][1],NORMAL[1][0],NORMAL[1][1],NORMAL[2][0],NORMAL[2][1]);}

                for(i=0;i<3;i++){C_SOUND[i] = sqrt(GAMMA*PRESSURE[i]/U_N[0][i]);}

#ifdef DEBUG
                cout << "-- FIRST  -------------------------------------------------------" << endl;
                cout << "Time     =\t" << T << endl;
                cout << "i =\t" << X[0] << "\t" << Y[0] << endl;
                cout << "j =\t" << X[1] << "\t" << Y[1] << endl;
                cout << "k =\t" << X[2] << "\t" << Y[2] << endl;
                cout << "State    =" << "\trho" << "\tx_mom" << "\ty_mom" << "\tenergy" << endl;
                for(i=0;i<3;i++){cout << i << " =\t" << U_N[0][i] << "\t" << U_N[1][i] << "\t" << U_N[2][i] << "\t" << U_N[3][i] << endl;}
                cout << "Pressure =\t" << PRESSURE[0] << "\t" << PRESSURE[1] << "\t" << PRESSURE[2] << endl;
                cout << "Sound Speed =\t" << C_SOUND[0] << "\t" << C_SOUND[1] << "\t" << C_SOUND[2] << endl;
#endif

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                // Calculate inflow parameters

                double U_AVG[4];
                double C,U,U_C,V,V_C,H,H_C,GAMMA_1,GAMMA_2,ALPHA,ALPHA_C,W;
                double Z[4][3];
                double L_1,L_2,L_3,L_4,L_12,L_123;
                double PRESSURE_AVG,C_SOUND_AVG;
                double LAMBDA_PLUS[4][3],LAMBDA_MINUS[4][3];
                double N_X[3],N_Y[3];
                double BETA[4][4][3];

                for(i=0;i<4;++i){U_AVG[i] = (U_N[i][0] + U_N[i][1] + U_N[i][2])/3.0;}

                PRESSURE_AVG = (PRESSURE[0] + PRESSURE[1] + PRESSURE[2])/3.0;
                C_SOUND_AVG  = (C_SOUND[0]  + C_SOUND[1]  + C_SOUND[2])/3.0;

#ifdef DEBUG
                cout << "U_AVG =\t" << U_AVG[0] << "\t" << U_AVG[1] << "\t" << U_AVG[2] << "\t" << U_AVG[3] << endl;
                cout << "PRESSURE_AVG =\t" << PRESSURE_AVG << endl;
                cout << "C_SOUND_AVG =\t" << C_SOUND_AVG << endl;
#endif

                for(m=0;m<3;++m){
                        Z[0][m] = sqrt(U_N[0][m]);
                        Z[1][m] = U_N[1][m]/Z[0][m];
                        Z[2][m] = U_N[2][m]/Z[0][m];
                        Z[3][m] = (U_N[3][m] + PRESSURE[m])/Z[0][m];

                        N_X[m]  = NORMAL[m][0];
                        N_Y[m]  = NORMAL[m][1];
                }

                C   = C_SOUND_AVG;

                U   = U_AVG[1]/U_AVG[0];
                U_C = U/C;

                V   = U_AVG[2]/U_AVG[0];
                V_C = V/C;

                H   = (U_AVG[3] + PRESSURE_AVG)/U_AVG[0];
                H_C = H/C;

                GAMMA_1 = GAMMA - 1.0;
                GAMMA_2 = GAMMA - 2.0;

                ALPHA   = GAMMA_1*(U*U + V*V)/2.0;
                ALPHA_C = ALPHA/C;

#ifdef DEBUG
                cout << "U =\t" << U << "\t" << U_C << endl;
                cout << "V =\t" << V << "\t" << V_C << endl;
                cout << "ALPHA =\t" << ALPHA << "\t" << ALPHA_C << endl;
                cout << "W =\t";
#endif

                for(m=0;m<3;++m){

                        W = U*N_X[m] + V*N_Y[m];

                        LAMBDA[0][m] = W + C;
                        LAMBDA[1][m] = W - C;
                        LAMBDA[2][m] = W;
                        LAMBDA[3][m] = W;

#ifdef DEBUG
                        cout << W << "\t";
#endif

                        for(i=0;i<4;++i){
                                LAMBDA_PLUS[i][m]  = max_val(0.0,LAMBDA[i][m]);
                                LAMBDA_MINUS[i][m] = min_val(0.0,LAMBDA[i][m]);
                        }

                        for(p=0;p<3;++p){
                                if(p==0){
                                        L_1 = LAMBDA_PLUS[0][m];
                                        L_2 = LAMBDA_PLUS[1][m];
                                        L_3 = LAMBDA_PLUS[2][m];
                                        L_4 = LAMBDA_PLUS[3][m];
                                }else if(p==1){
                                        L_1 = LAMBDA_MINUS[0][m];
                                        L_2 = LAMBDA_MINUS[1][m];
                                        L_3 = LAMBDA_MINUS[2][m];
                                        L_4 = LAMBDA_MINUS[3][m];
                                }else{
                                        L_1 = LAMBDA[0][m];
                                        L_2 = LAMBDA[1][m];
                                        L_3 = LAMBDA[2][m];
                                        L_4 = LAMBDA[3][m];

                                }

                                L_12  = (L_1 - L_2)/2.0;
                                L_123 = (L_1 + L_2 - 2.0*L_3)/2.0;

                                INFLOW[0][0][m][p] = ALPHA_C*L_123/C - W*L_12/C + L_3;
                                INFLOW[0][1][m][p] = -1.0*GAMMA_1*U_C*L_123/C +  N_X[m]*L_12/C;
                                INFLOW[0][2][m][p] = -1.0*GAMMA_1*V_C*L_123/C +  N_Y[m]*L_12/C;
                                INFLOW[0][3][m][p] = GAMMA_1*L_123/(C*C);

                                INFLOW[1][0][m][p] = (ALPHA_C*U_C - W*N_X[m])*L_123 + (ALPHA_C*N_X[m] - U_C*W)*L_12;
                                INFLOW[1][1][m][p] = (N_X[m]*N_X[m] - GAMMA_1*U_C*U_C)*L_123 - (GAMMA_2*U_C*N_X[m]*L_12) + L_3;
                                INFLOW[1][2][m][p] = (N_X[m]*N_Y[m] - GAMMA_1*U_C*V_C)*L_123 + (U_C*N_Y[m] - GAMMA_1*V_C*N_X[m])*L_12;
                                INFLOW[1][3][m][p] = GAMMA_1*U_C*L_123/C + GAMMA_1*N_X[m]*L_12/C;

                                INFLOW[2][0][m][p] = (ALPHA_C*V_C - W*N_Y[m])*L_123 + (ALPHA_C*N_Y[m] - V_C*W)*L_12;
                                INFLOW[2][1][m][p] = (N_X[m]*N_Y[m] - GAMMA_1*U_C*V_C)*L_123 + (V_C*N_X[m] - GAMMA_1*U_C*N_Y[m])*L_12;
                                INFLOW[2][2][m][p] = (N_Y[m]*N_Y[m] - GAMMA_1*V_C*V_C)*L_123 - (GAMMA_2*V_C*N_Y[m]*L_12) + L_3;
                                INFLOW[2][3][m][p] = GAMMA_1*V_C*L_123/C + GAMMA_1*N_Y[m]*L_12/C;

                                INFLOW[3][0][m][p] = (ALPHA_C*H_C - W*W)*L_123 + W*(ALPHA_C - H_C)*L_12;
                                INFLOW[3][1][m][p] = (W*N_X[m] - U - ALPHA_C*U_C)*L_123 + (H_C*N_X[m] - GAMMA_1*U_C*W)*L_12;
                                INFLOW[3][2][m][p] = (W*N_Y[m] - V - ALPHA_C*V_C)*L_123 + (H_C*N_Y[m] - GAMMA_1*V_C*W)*L_12;
                                INFLOW[3][3][m][p] = GAMMA_1*H_C*L_123/C + GAMMA_1*W*L_12/C + L_3;
                        }
                }

#ifdef DEBUG
                cout << endl;;
#endif

#ifdef DEBUG
                cout << "Lambda + =\t" << LAMBDA_PLUS[0][0] << "\t" << LAMBDA_PLUS[0][1] << "\t" << LAMBDA_PLUS[0][2] << endl;
                cout << "Lambda + =\t" << LAMBDA_PLUS[1][0] << "\t" << LAMBDA_PLUS[1][1] << "\t" << LAMBDA_PLUS[1][2] << endl;
                cout << "Lambda + =\t" << LAMBDA_PLUS[2][0] << "\t" << LAMBDA_PLUS[2][1] << "\t" << LAMBDA_PLUS[2][2] << endl;
                cout << "Lambda + =\t" << LAMBDA_PLUS[3][0] << "\t" << LAMBDA_PLUS[3][1] << "\t" << LAMBDA_PLUS[3][2] << endl;

                cout << "Lambda - =\t" << LAMBDA_MINUS[0][0] << "\t" << LAMBDA_MINUS[0][1] << "\t" << LAMBDA_MINUS[0][2] << endl;
                cout << "Lambda - =\t" << LAMBDA_MINUS[1][0] << "\t" << LAMBDA_MINUS[1][1] << "\t" << LAMBDA_MINUS[1][2] << endl;
                cout << "Lambda - =\t" << LAMBDA_MINUS[2][0] << "\t" << LAMBDA_MINUS[2][1] << "\t" << LAMBDA_MINUS[2][2] << endl;
                cout << "Lambda - =\t" << LAMBDA_MINUS[3][0] << "\t" << LAMBDA_MINUS[3][1] << "\t" << LAMBDA_MINUS[3][2] << endl;

                cout << "Lambda   =\t" << LAMBDA[0][0] << "\t" << LAMBDA[0][1] << "\t" << LAMBDA[0][2] << endl;
                cout << "Lambda   =\t" << LAMBDA[1][0] << "\t" << LAMBDA[1][1] << "\t" << LAMBDA[1][2] << endl;
                cout << "Lambda   =\t" << LAMBDA[2][0] << "\t" << LAMBDA[2][1] << "\t" << LAMBDA[2][2] << endl;
                cout << "Lambda   =\t" << LAMBDA[3][0] << "\t" << LAMBDA[3][1] << "\t" << LAMBDA[3][2] << endl;
#endif

                for(i=0;i<4;++i){
                        PHI[i] = 0.0;
                        for(m=0;m<3;++m){
                                PHI[i] += INFLOW[i][0][m][2]*Z[0][m] + INFLOW[i][1][m][2]*Z[1][m] + INFLOW[i][2][m][2]*Z[2][m] + INFLOW[i][3][m][2]*Z[3][m];
                        }
                }

#ifdef DEBUG
                cout << "PHI =\t" << PHI[0] << "\t" << PHI[1] << "\t" << PHI[2] << "\t" << PHI[3] << endl;
#endif

                for(i=0;i<4;++i){
                        for(j=0;j<4;++j){
                                INFLOW_MINUS_SUM[i][j] = 0.0;
                                for(m=0;m<3;++m){
                                        INFLOW_MINUS_SUM[i][j] += INFLOW[i][j][m][1];
                                }
                        }
                }

                matInv(&INFLOW_MINUS_SUM[0][0],4);

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                BETA[i][0][m] = -1.0*INFLOW[i][0][m][0] * INFLOW_MINUS_SUM[0][0] + -1.0*INFLOW[i][1][m][0] * INFLOW_MINUS_SUM[1][0] + -1.0*INFLOW[i][2][m][0] * INFLOW_MINUS_SUM[2][0] + -1.0*INFLOW[i][3][m][0] * INFLOW_MINUS_SUM[3][0];
                                BETA[i][1][m] = -1.0*INFLOW[i][0][m][0] * INFLOW_MINUS_SUM[0][1] + -1.0*INFLOW[i][1][m][0] * INFLOW_MINUS_SUM[1][1] + -1.0*INFLOW[i][2][m][0] * INFLOW_MINUS_SUM[2][1] + -1.0*INFLOW[i][3][m][0] * INFLOW_MINUS_SUM[3][1];
                                BETA[i][2][m] = -1.0*INFLOW[i][0][m][0] * INFLOW_MINUS_SUM[0][2] + -1.0*INFLOW[i][1][m][0] * INFLOW_MINUS_SUM[1][2] + -1.0*INFLOW[i][2][m][0] * INFLOW_MINUS_SUM[2][2] + -1.0*INFLOW[i][3][m][0] * INFLOW_MINUS_SUM[3][2];
                                BETA[i][3][m] = -1.0*INFLOW[i][0][m][0] * INFLOW_MINUS_SUM[0][3] + -1.0*INFLOW[i][1][m][0] * INFLOW_MINUS_SUM[1][3] + -1.0*INFLOW[i][2][m][0] * INFLOW_MINUS_SUM[2][3] + -1.0*INFLOW[i][3][m][0] * INFLOW_MINUS_SUM[3][3];
                        }
                }

                // Calculate spatial splitting for first half timestep

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                FLUC[i][m] = BETA[i][0][m] * (PHI[0]) + BETA[i][1][m] * (PHI[1]) + BETA[i][2][m] * (PHI[2]) + BETA[i][3][m] * (PHI[3]);
                                //cout << "FLUCTUATION =\t" << i << "\t" << m << "\t" << FLUC[i][m] << endl;
                        }
                }

                // Calculate change to be distributed

                for(i=0;i<4;i++){
                        DU0[i] = DT*FLUC[i][0]/DUAL[0];
                        DU1[i] = DT*FLUC[i][1]/DUAL[1];
                        DU2[i] = DT*FLUC[i][2]/DUAL[2];
                }

                VERTEX_0->update_du_half(DU0);
                VERTEX_1->update_du_half(DU1);
                VERTEX_2->update_du_half(DU2);

#ifdef DEBUG
                        for(i=0;i<4;i++){cout << "Element fluctuation =\t" << FLUC[i][0] << "\t" << FLUC[i][1] << "\t" << FLUC[i][2] << endl;}
                        cout << "Dual =\t" << VERTEX_0->get_dual() << "\t" << VERTEX_1->get_dual() << "\t" << VERTEX_2->get_dual() << endl;
                        cout << "Change (rho) =\t"    << DU0[0] << "\t" << DU1[0] << "\t" << DU2[0] << endl;
                        cout << "Change (x mom) =\t"  << DU0[1] << "\t" << DU1[1] << "\t" << DU2[1] << endl;
                        cout << "Change (y mom) =\t"  << DU0[2] << "\t" << DU1[2] << "\t" << DU2[2] << endl;
                        cout << "Change (energy) =\t" << DU0[3] << "\t" << DU1[3] << "\t" << DU2[3] << endl;
                        cout << "-----------------------------------------------------------------" << endl;
                        //if(U_N[0][0] != U_N[0][1] or U_N[0][0] != U_N[0][2] or U_N[0][1] != U_N[0][2]){exit(0);}
#endif

                return ;
        }

        //**********************************************************************************************************************

        void calculate_second_half(double T, double DT_TOT, double DX, double DY){
                int i,j,k;
                double DU0[4],DU1[4],DU2[4];
//                 double INFLOW_HALF[4][3],INFLOW_PLUS_HALF[4][3],INFLOW_MINUS_HALF[4][3];
//                 double IN_TOP,IN_BOTTOM,OUT_TOP,OUT_BOTTOM;
//                 double MASS[4][3][3],MASS_GAL[4][3][3];
//                 double SUM_MASS[4][3];
//                 double DT = DT_TOT;

//                 setup_half_state();

//                 if(abs(X[0] - X[1]) > 10.0 or abs(X[0] - X[2]) > 10.0 or abs(X[1] - X[2]) > 10.0){
// #ifdef DEBUG
//                         cout << "Skipping x boundary" << endl;
// #endif
//                         return ;
//                 }else if(abs(Y[0] - Y[1]) > 10.0 or abs(Y[0] - Y[2]) > 10.0 or abs(Y[1] - Y[2]) > 10.0){
// #ifdef DEBUG
//                         cout << "Skipping y boundary" << endl;
// #endif
//                         return ;
//                 }

// #ifdef DEBUG
//                 cout << "-- SECOND -------------------------------------------------------" << endl;
//                 cout << "Time     =\t" << T << endl;
//                 cout << "i =\t" << X[0] << "\t" << Y[0] << endl;
//                 cout << "j =\t" << X[1] << "\t" << Y[1] << endl;
//                 cout << "k =\t" << X[2] << "\t" << Y[2] << endl;
//                 cout << "Half State    =" << "\trho" << "\tx_mom" << "\ty_mom" << "\tenergy" << endl;
//                 for(i=0;i<3;i++){cout << i << " =\t" << U_HALF[0][i] << "\t" << U_HALF[1][i] << "\t" << U_HALF[2][i] << "\t" << U_HALF[3][i] << endl;}
//                 cout << "Pressure =\t" << PRESSURE_HALF[0] << "\t" << PRESSURE_HALF[1] << "\t" << PRESSURE_HALF[2] << endl;
// #endif

//                 // calcualte sound speed for half time state

//                 for(i=0;i<3;i++){C_SOUND_HALF[i] = sqrt(GAMMA*PRESSURE_HALF[i]/U_HALF[0][i]);}

//                 for(i=0;i<3;i++){
//                         for (j=0;j<2;j++){
//                                 LAMBDA_HALF[0][i][j] = VEL_HALF[i][j] - C_SOUND_HALF[i];
//                                 LAMBDA_HALF[1][i][j] = VEL_HALF[i][j];
//                                 LAMBDA_HALF[2][i][j] = VEL_HALF[i][j];
//                                 LAMBDA_HALF[3][i][j] = VEL_HALF[i][j] + C_SOUND_HALF[i];
//                         }

// #ifdef DEBUG
//                         for(k=0;k<4;++k){cout << "lambda " << i << " =\t" << LAMBDA_HALF[k][i][0] << "\t" << LAMBDA_HALF[k][i][1] << endl;}
// #endif
//                 }

//                 // Calculate inflow parameters for half timestep state

//                 for(i=0;i<4;i++){
//                         for(j=0;j<3;j++){
//                                 INFLOW_HALF[i][j] = 0.0 ;
//                                 for(k=0;k<2;k++){INFLOW_HALF[i][j] = INFLOW_HALF[i][j] + 0.5*LAMBDA_HALF[i][j][k]*NORMAL[j][k];}
//                                 INFLOW_PLUS_HALF[i][j]  = max_val(0,INFLOW_HALF[i][j]);
//                                 INFLOW_MINUS_HALF[i][j] = min_val(0,INFLOW_HALF[i][j]);
//                         }
//                 }

//                 // Calculate inflow and outflow state of element at half timestep

//                 for(i=0;i<4;i++){
//                         IN_TOP     = 0.0;
//                         IN_BOTTOM  = 0.0;
//                         OUT_TOP    = 0.0;
//                         OUT_BOTTOM = 0.0;
//                         for(j=0;j<3;j++){
//                                 IN_TOP     = IN_TOP     + INFLOW_MINUS_HALF[i][j] * U_HALF[i][j];
//                                 IN_BOTTOM  = IN_BOTTOM  + INFLOW_MINUS_HALF[i][j];
//                                 OUT_TOP    = OUT_TOP    + INFLOW_PLUS_HALF[i][j]  * U_HALF[i][j];
//                                 OUT_BOTTOM = OUT_BOTTOM + INFLOW_PLUS_HALF[i][j];
//                         }
// #ifdef DEBUG
//                         cout << "Sums (" << i << ") =\t" << IN_TOP << "\t" << IN_BOTTOM << "\t" << OUT_TOP << "\t" << OUT_BOTTOM << endl;
// #endif
//                         if(OUT_BOTTOM != 0.0 and IN_BOTTOM != 0.0){
//                                 U_IN_HALF[i]  = IN_TOP  / IN_BOTTOM;
//                                 U_OUT_HALF[i] = OUT_TOP / OUT_BOTTOM;
//                                 for(j=0;j<3;j++){BETA_HALF[i][j] = INFLOW_PLUS_HALF[i][j] / OUT_BOTTOM;}
//                         }else{
//                                 U_IN_HALF[i]  = 0.0;
//                                 U_OUT_HALF[i] = 0.0;
//                                 for(j=0;j<3;j++){BETA_HALF[i][j] = 0.0;}
//                         }
//                 }

//                 for(i=0;i<4;i++){
//                         for(j=0;j<3;j++){
//                                 FLUC_HALF[i][j] = INFLOW_PLUS_HALF[i][j]*(U_OUT_HALF[i]-U_IN_HALF[i]);
//                         }
//                 }

//                 for(i=0;i<4;i++){
//                         for(j=0;j<3;j++){
//                                 SUM_MASS[i][j] = 0.0;
//                                 for(k=0;k<3;k++){
//                                         MASS[i][j][k] = AREA*BETA[i][j]/3.0;
//                                         if(j == k){
//                                                 MASS_GAL[i][j][k] = AREA/6.0;
//                                         }else{
//                                                 MASS_GAL[i][j][k] = AREA/12.0;
//                                         }
//                                         SUM_MASS[i][j] += (MASS[i][j][k]-MASS_GAL[i][j][k])*(U_HALF[i][j]-U_N[i][j])/DT;
//                                 }
// #ifdef DEBUG
//                                 cout << "Mass Matrix Sum: " << i <<  " =\t" << SUM_MASS[i][0] << "\t" << SUM_MASS[i][1] << "\t" << SUM_MASS[i][2] << endl;
// #endif
//                         }
//                 }



                for(i=0;i<4;i++){
                        DU0[i] = 0.0;//(DT/DUAL[0])*(SUM_MASS[i][0]+0.5*(FLUC[i][0]+FLUC_HALF[i][0]));
                        DU1[i] = 0.0;//(DT/DUAL[1])*(SUM_MASS[i][1]+0.5*(FLUC[i][1]+FLUC_HALF[i][1]));
                        DU2[i] = 0.0;//(DT/DUAL[2])*(SUM_MASS[i][2]+0.5*(FLUC[i][2]+FLUC_HALF[i][2]));
                }

                VERTEX_0->update_du(DU0);
                VERTEX_1->update_du(DU1);
                VERTEX_2->update_du(DU2);

// #ifdef DEBUG
//                 for(i=0;i<4;i++){cout << "u_in =\t" << U_IN_HALF[i] << "\tu_out =\t" << U_OUT_HALF[i] << endl;}
//                 for(i=0;i<4;i++){cout << "Element fluctuation =\t" << FLUC_HALF[i][0] << "\t" << FLUC_HALF[i][1] << "\t" << FLUC_HALF[i][2] << endl;}
//                 for(i=0;i<4;i++){cout << "Beta (" << i << ") =\t" << BETA_HALF[i][0] << "\t" << BETA_HALF[i][1] << "\t" << BETA_HALF[i][2] << "\tTotal =\t" << BETA_HALF[i][0]+BETA_HALF[i][1]+BETA_HALF[i][2] << endl;}
//                 cout << "Dual =\t" << VERTEX_0->get_dual() << "\t" << VERTEX_1->get_dual() << "\t" << VERTEX_2->get_dual() << endl;
//                 cout << "Change (rho) =\t"    << DU0[0] << "\t" << DU1[0] << "\t" << DU2[0] << endl;
//                 cout << "Change (x mom) =\t"  << DU0[1] << "\t" << DU1[1] << "\t" << DU2[1] << endl;
//                 cout << "Change (y mom) =\t"  << DU0[2] << "\t" << DU1[2] << "\t" << DU2[2] << endl;
//                 cout << "Change (energy) =\t" << DU0[3] << "\t" << DU1[3] << "\t" << DU2[3] << endl;
//                 cout << "-----------------------------------------------------------------" << endl;
// #endif

//                 return ;
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
                        NORMAL[i][0] = PERP[i][0];
                        NORMAL[i][1] = PERP[i][1];
                }

                AREA = 0.5*(X[0]*Y[1]-Y[0]*X[1]);

                NORMAL00 = NORMAL[0][0];
                NORMAL01 = NORMAL[0][1];
                NORMAL10 = NORMAL[1][0];
                NORMAL11 = NORMAL[1][1];
                NORMAL20 = NORMAL[2][0];
                NORMAL21 = NORMAL[2][1];

#ifdef DEBUG
                        cout << "X =\t" << X[0] << "\t" << X[1] << "\t" << X[2] << endl;
                        cout << "Y =\t" << Y[0] << "\t" << Y[1] << "\t" << Y[2] << endl;
                        cout << "Normal i =\t" << NORMAL00 << "\t" << NORMAL01 << endl;
                        cout << "Normal j =\t" << NORMAL10 << "\t" << NORMAL11 << endl;
                        cout << "Normal k =\t" << NORMAL20 << "\t" << NORMAL21 << endl;
#endif

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