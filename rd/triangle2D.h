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

        double LAMBDA[4][3]
        ;
        double LAMBDA_HALF[4][3][2];

        double VEL[3][2];
        double VEL_HALF[3][2];

        double PRESSURE[3],C_SOUND[3];
        double PRESSURE_HALF[3],C_SOUND_HALF[3];

        double FLUC[4][3];
        double FLUC_HALF[4][3];

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
        void calculate_first_half(double T, double DT_TOT){
                int i,j,k,m;

                double DU0[4],DU1[4],DU2[4];

                double INFLOW[4][4][3],INFLOW_PLUS[4][4][3],INFLOW_MINUS[4][4][3];
                double INFLOW_PLUS_INVERSE[4][4],INFLOW_MINUS_INVERSE[4][4];
                double INFLOW_PLUS_SUM[4][4], INFLOW_MINUS_SUM[4][4];

                double IN_TOP[4],OUT_TOP[4];
                double DT = 0.5*DT_TOT;

                // Import conditions and positions of vertices

                setup_positions();
                setup_initial_state();

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

                for(i=0;i<3;i++){C_SOUND[i] = sqrt(GAMMA*PRESSURE[i]/U_N[0][i]);}

                // Calculate normals (just in first timestep for static grid)

                if(T==0.0){caclulate_normals(X,Y,NORMAL[0][0],NORMAL[0][1],NORMAL[1][0],NORMAL[1][1],NORMAL[2][0],NORMAL[2][1]);}

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

                double A[4][4],B[4][4],U[4];
                double INFLOW_INVERSE_0[4][4],INFLOW_INVERSE_1[4][4],INFLOW_INVERSE_2[4][4];

                for(i=0;i<4;++i){U[i] = (U_N[i][0] + U_N[i][1] + U_N[i][2])/3.0;}

                A[0][0] = 0.0;
                A[0][1] = 1.0;
                A[0][2] = 0.0;
                A[0][3] = 0.0;

                A[1][0] = (GAMMA - 3.0) * (U[1] * U[1]) / (2.0 * U[0] * U[0]) + (GAMMA - 1.0) * (U[2] * U[2]) / (2.0 * U[0] * U[0]);
                A[1][1] = (3.0 - GAMMA) * U[1] / U[0];
                A[1][2] = (1.0 - GAMMA) * U[2] / U[0];
                A[1][3] = GAMMA - 1.0;

                A[2][0] = -1.0 * (U[1] * U[2])/U[0];
                A[2][1] = U[2] / U[0];
                A[2][2] = U[1] / U[0];
                A[2][3] = 0.0;

                A[3][0] = -1.0 * (U[1] * U[3] * GAMMA) / (U[0] * U[0]) + (GAMMA - 1.0) * (U[1] * U[1] * U[1] + U[1] * U[2] * U[2]) / (U[0] * U[0] * U[0]);
                A[3][1] = GAMMA * U[3] / U[0] + (1 - GAMMA) * (0.5) * ((3.0 * U[1] * U[1]) / (U[0] * U[0]) + (U[2] * U[2])/(U[0] * U[0]));
                A[3][2] = U[1] * (1 - GAMMA) * U[2] * U[0];
                A[3][3] = GAMMA * U[1] / U[0];


                B[0][0] = 0.0;
                B[0][1] = 0.0;
                B[0][2] = 1.0;
                B[0][3] = 0.0;

                B[1][0] = -1.0 * (U[1]*U[2]) / (U[0] * U[0]);
                B[1][1] = U[2] / U[0];
                B[1][2] = U[1] / U[0];
                B[1][3] = 0.0;

                B[2][0] = -1.0 * (U[2] * U[2]) / (U[0] * U[0]) + (GAMMA - 1.0)*((U[1] * U[1]) / (2.0 * U[0] * U[0]) + (U[2] * U[2]) / (2.0 * U[0] * U[0]));
                B[2][1] = (1.0 - GAMMA) * (U[1] / U[0]);
                B[2][2] = (3.0 - GAMMA) * (U[2] / U[0]);
                B[2][3] = GAMMA - 1.0;

                B[3][0] = -1.0 * (U[2] * U[3]) / (U[0] * U[0]) + U[2] * (GAMMA - 1.0) * ((-1.0 * (U[3] / (U[0] * U[0]))) + ((U[1] * U[1]) + (U[2] * U[2])) / (U[0] * U[0] * U[0]));
                B[3][1] = U[2] * (1.0 - GAMMA) * (U[1] / (U[0] * U[0]));
                B[3][2] = GAMMA * U[3] / U[0] - 0.5 * (GAMMA - 1.0) * (U[1] * U[1] + 3.0 * U[2] * U[2]) / (U[0] * U[0]);
                B[3][3] = GAMMA * U[2] / U[0];

                for(i=0;i<4;++i){
                        for(j=0;j<4;++j){
                                INFLOW_PLUS_SUM[i][j]  = 0.0;
                                INFLOW_MINUS_SUM[i][j] = 0.0;
                                for(m=0;m<3;++m){
                                        INFLOW[i][j][m] = 0.5*(A[i][j] * NORMAL[m][0] + B[i][j] * NORMAL[m][1]);
                                        if(m == 0){
                                                INFLOW_INVERSE_0[i][j] = INFLOW[i][j][m];
                                        }else if(m == 1){
                                                INFLOW_INVERSE_1[i][j] = INFLOW[i][j][m];
                                        }else{
                                                INFLOW_INVERSE_2[i][j] = INFLOW[i][j][m];
                                        }
                                }
                        }
                }

                matFac(&INFLOW_INVERSE_0[0][0],4);
                matFac(&INFLOW_INVERSE_1[0][0],4);
                matFac(&INFLOW_INVERSE_2[0][0],4);

                for(i=0;i<4;++i){
                        LAMBDA[i][0] = INFLOW_INVERSE_0[i][i];
                        LAMBDA[i][1] = INFLOW_INVERSE_1[i][i];
                        LAMBDA[i][2] = INFLOW_INVERSE_2[i][i];
#ifdef DEBUG
                        cout << "Lambda =\t" << LAMBDA[i][0] << "\t" << LAMBDA[i][1] << "\t" << LAMBDA[i][2] << endl;
#endif
                }

                for(i=0;i<4;++i){
                        for(j=0;j<4;++j){
                                for(m=0;m<4;++m){
                                        if(LAMBDA[i][m] >= 0.0){
                                                INFLOW_PLUS[i][j][m]  = INFLOW[i][j][m];
                                                INFLOW_MINUS[i][j][m] = 0.0;
                                        }else{
                                                INFLOW_PLUS[i][j][m]  = 0.0;
                                                INFLOW_MINUS[i][j][m] = INFLOW[i][j][m];
                                        }
                                        INFLOW_PLUS_INVERSE[i][j]  = INFLOW_PLUS_SUM[i][j] += INFLOW_PLUS[i][j][m];
                                        INFLOW_MINUS_INVERSE[i][j] = INFLOW_MINUS_SUM[i][j] += INFLOW_MINUS[i][j][m];
                                }
                        }
                }

                for(i=0;i<4;++i){
                        IN_TOP[i] = 0.0;
                        OUT_TOP[i] = 0.0;
                        for(m=0;m<3;++m){
                                IN_TOP[i]  += INFLOW_PLUS[i][0][m]  * U_N[0][m] + INFLOW_PLUS[i][1][m]  * U_N[1][m] + INFLOW_PLUS[i][2][m]  * U_N[2][m] + INFLOW_PLUS[i][3][m]  * U_N[3][m];
                                OUT_TOP[i] += INFLOW_MINUS[i][0][m] * U_N[0][m] + INFLOW_MINUS[i][1][m] * U_N[1][m] + INFLOW_MINUS[i][2][m] * U_N[2][m] + INFLOW_MINUS[i][3][m] * U_N[3][m];
                        }
                }

                matInv(&INFLOW_PLUS_SUM[0][0],4);
                matInv(&INFLOW_MINUS_SUM[0][0],4);

                for(i=0;i<4;++i){
                        U_OUT[i] = 0.0;
                        U_OUT[i] += INFLOW_PLUS_INVERSE[i][0] * OUT_TOP[0] + INFLOW_PLUS_INVERSE[i][1] * OUT_TOP[1] + INFLOW_PLUS_INVERSE[i][2] * OUT_TOP[2] + INFLOW_PLUS_INVERSE[i][3] * OUT_TOP[3];
                }

                for(i=0;i<4;++i){
                        U_IN[i] = 0.0;
                        U_IN[i] += INFLOW_MINUS_INVERSE[i][0] * IN_TOP[0] + INFLOW_MINUS_INVERSE[i][1] * IN_TOP[1] + INFLOW_MINUS_INVERSE[i][2] * IN_TOP[2] + INFLOW_MINUS_INVERSE[i][3] * IN_TOP[3];
                }

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                FLUC[i][m] = INFLOW_PLUS[i][0][m] * (U_OUT[0] - U_IN[0]) + INFLOW_PLUS[i][1][m] * (U_OUT[1] - U_IN[1]) + INFLOW_PLUS[i][2][m] * (U_OUT[2] - U_IN[2]) + INFLOW_PLUS[i][3][m] * (U_OUT[3] - U_IN[3]);
                                //cout << "FLUCTUATION =\t" << i << "\t" << m << "\t" << FLUC[i][m] << endl;
                        }
                }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                // Calculate spatial splitting for first half timestep

                // for(i=0;i<4;i++){
                //         for(j=0;j<3;j++){
                //                 FLUC[i][j] = INFLOW_PLUS[i][j]*(U_OUT[i]-U_IN[i]);
                //         }
                // }

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
                        for(i=0;i<4;i++){cout << "u_in =\t" << U_IN[i] << "\tu_out =\t" << U_OUT[i] << endl;}
                        for(i=0;i<4;i++){cout << "Element fluctuation =\t" << FLUC[i][0] << "\t" << FLUC[i][1] << "\t" << FLUC[i][2] << endl;}
                        for(i=0;i<4;i++){cout << "Beta (" << i << ") =\t" << BETA[i][0] << "\t" << BETA[i][1] << "\t" << BETA[i][2] << "\tTotal =\t" << BETA[i][0]+BETA[i][1]+BETA[i][2] << endl;}
                        cout << "Dual =\t" << VERTEX_0->get_dual() << "\t" << VERTEX_1->get_dual() << "\t" << VERTEX_2->get_dual() << endl;
                        cout << "Change (rho) =\t"    << DU0[0] << "\t" << DU1[0] << "\t" << DU2[0] << endl;
                        cout << "Change (x mom) =\t"  << DU0[1] << "\t" << DU1[1] << "\t" << DU2[1] << endl;
                        cout << "Change (y mom) =\t"  << DU0[2] << "\t" << DU1[2] << "\t" << DU2[2] << endl;
                        cout << "Change (energy) =\t" << DU0[3] << "\t" << DU1[3] << "\t" << DU2[3] << endl;
                        cout << "-----------------------------------------------------------------" << endl;
#endif

                //exit(0);
                return ;
        }

        //**********************************************************************************************************************

        void calculate_second_half(double T, double DT_TOT){
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
                        NORMAL[i][0] = PERP[i][0]/MAG;
                        NORMAL[i][1] = PERP[i][1]/MAG;
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