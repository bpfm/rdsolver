/* class containing values associated with the face
        *VERTEX_0 = pointer to left VERTEX
        *VERTEX_1 = pointer to right VERTEX
*/

class TRIANGLE{

private:
        int ID;
        VERTEX *VERTEX_0,*VERTEX_1,*VERTEX_2;
        double AREA;

        double X[3],Y[3],DUAL[3];
        double NORMAL[3][2];

        double U_N[4][3];
        double U_HALF[4][3];

        double FLUC[4][3];
        double FLUC_HALF[4][3];

        double PRESSURE[3];

        double PHI[4];
        double BETA[4][4][3];

        double MAG[3];

public:

        void set_id(int NEW_ID){ID = NEW_ID;}

        void set_vertex_0(VERTEX* NEW_VERTEX){VERTEX_0 = NEW_VERTEX;}
        void set_vertex_1(VERTEX* NEW_VERTEX){VERTEX_1 = NEW_VERTEX;}
        void set_vertex_2(VERTEX* NEW_VERTEX){VERTEX_2 = NEW_VERTEX;}

        int get_id(){return ID;}

        VERTEX* get_vertex_0(){return VERTEX_0;}
        VERTEX* get_vertex_1(){return VERTEX_1;}
        VERTEX* get_vertex_2(){return VERTEX_2;}

        void print_triangle_state(){
                std::cout << "U[0] =\t" << U_N[0][0] << "\t" << U_N[0][1] << "\t" << U_N[0][2] << std::endl;
                std::cout << "U[1] =\t" << U_N[1][0] << "\t" << U_N[1][1] << "\t" << U_N[1][2] << std::endl;
                std::cout << "U[2] =\t" << U_N[2][0] << "\t" << U_N[2][1] << "\t" << U_N[2][2] << std::endl;
                std::cout << "U[4] =\t" << U_N[3][0] << "\t" << U_N[3][1] << "\t" << U_N[3][2] << std::endl;

                std::cout << "U_HALF[0] =\t" << U_HALF[0][0] << "\t" << U_HALF[0][1] << "\t" << U_HALF[0][2] << std::endl;
                std::cout << "U_HALF[1] =\t" << U_HALF[1][0] << "\t" << U_HALF[1][1] << "\t" << U_HALF[1][2] << std::endl;
                std::cout << "U_HALF[2] =\t" << U_HALF[2][0] << "\t" << U_HALF[2][1] << "\t" << U_HALF[2][2] << std::endl;
                std::cout << "U_HALF[4] =\t" << U_HALF[3][0] << "\t" << U_HALF[3][1] << "\t" << U_HALF[3][2] << std::endl;
        }

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

                PRESSURE[0] = VERTEX_0->get_pressure();
                PRESSURE[1] = VERTEX_1->get_pressure();
                PRESSURE[2] = VERTEX_2->get_pressure();
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

                PRESSURE[0] = VERTEX_0->get_pressure_half();
                PRESSURE[1] = VERTEX_1->get_pressure_half();
                PRESSURE[2] = VERTEX_2->get_pressure_half();
        }

        //**********************************************************************************************************************

        // Calculate first half timestep change, passing change to vertice
        void calculate_first_half(double T, double DT_TOT, double DX, double DY, std::ofstream &TEMP_FILE){
                int i,j,m,p;

                double DU0[4],DU1[4],DU2[4];

                double INFLOW[4][4][3][3];
                double INFLOW_PLUS_SUM[4][4], INFLOW_MINUS_SUM[4][4];

                double DT = DT_TOT;

                double C_SOUND[3];

                // Import conditions and positions of vertices

                setup_positions();
                setup_initial_state();

                if(T == 0.0){setup_normals(DX,DY);}

#ifdef CLOSED
                if(std::abs(X[0] - X[1]) > 2.0*DX or std::abs(X[0] - X[2]) > 2.0*DX or std::abs(X[1] - X[2]) > 2.0*DX){
#ifdef DEBUG
                        std::cout << "Skipping x boundary\t" << (X[0]+X[1]+X[2])/3.0 << "\t" << (Y[0]+Y[1]+Y[2])/3.0 << "\t" << std::endl;
#endif
                        return ;
                }else if(std::abs(Y[0] - Y[1]) > 2.0*DY or std::abs(Y[0] - Y[2]) > 2.0*DY or std::abs(Y[1] - Y[2]) > 2.0*DY){
#ifdef DEBUG
                        std::cout << "Skipping y boundary\t" << (X[0]+X[1]+X[2])/3.0 << "\t" << (Y[0]+Y[1]+Y[2])/3.0 << "\t" << std::endl;
#endif
                        return ;
                }
#endif

                for(i=0; i<3; ++i){C_SOUND[i] = sqrt(GAMMA*PRESSURE[i]/U_N[0][i]);}

#ifdef DEBUG
                std::cout << "-- FIRST  -------------------------------------------------------" << std::endl;
                std::cout << "Time     =\t" << T << std::endl;
                std::cout << "0 =\t" << X[0] << "\t" << Y[0] << std::endl;
                std::cout << "1 =\t" << X[1] << "\t" << Y[1] << std::endl;
                std::cout << "2 =\t" << X[2] << "\t" << Y[2] << std::endl;
                std::cout << "State    =" << "\trho" << "\tx_mom" << "\ty_mom" << "\tenergy" << std::endl;
                for(i=0;i<3;i++){std::cout << i << " =\t" << U_N[0][i] << "\t" << U_N[1][i] << "\t" << U_N[2][i] << "\t" << U_N[3][i] << std::endl;}
                std::cout << "Pressure =\t" << PRESSURE[0] << "\t" << PRESSURE[1] << "\t" << PRESSURE[2] << std::endl;
                std::cout << "Sound Speed =\t" << C_SOUND[0] << "\t" << C_SOUND[1] << "\t" << C_SOUND[2] << std::endl;
#endif
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // Calculate inflow parameters

                double H[3];
                double RHO,C,U,U_C,V,V_C,H_AVG,H_C,ALPHA,ALPHA_C,W;
                double Z[4][3];
                double VALUE1,VALUE2,VALUE3,VALUE4,VALUE12,VALUE123;
                double PRESSURE_AVG,C_SOUND_AVG;
                double LAMBDA[4][3],LAMBDA_PLUS[4][3],LAMBDA_MINUS[4][3];
                double N_X[3],N_Y[3];

                double Z_BAR[4],W_HAT[4][3];

                // Construct Roe vector Z

                for(m=0;m<3;++m){
                        Z[0][m] = sqrt(U_N[0][m]);
                        Z[1][m] = U_N[1][m]/Z[0][m];
                        Z[2][m] = U_N[2][m]/Z[0][m];
                        Z[3][m] = (U_N[3][m] + PRESSURE[m])/Z[0][m];

                        N_X[m]  = NORMAL[m][0];
                        N_Y[m]  = NORMAL[m][1];

                        H[m] = (U_N[3][m] + PRESSURE[m])/U_N[0][m];
#ifdef DEBUG
                        std::cout << "Z =\t" << Z[0][m] << "\t" << Z[1][m] << "\t" << Z[2][m] << "\t" << Z[3][m] << std::endl;
#endif
                }

                for(i=0; i<4; ++i){Z_BAR[i] = (Z[i][0] + Z[i][1] + Z[i][2])/3.0;}

#ifdef DEBUG
                for(i=0; i<4; ++i){std::cout << "Z_BAR " << i << " =\t" << Z_BAR[i] << std::endl;}
                std::cout << std::endl;
#endif


                for(m=0; m<3; ++m){
                        W_HAT[0][m] =  2.0*Z_BAR[0]*Z[0][m];
                        W_HAT[1][m] =  Z_BAR[1]*Z[0][m] + Z_BAR[0]*Z[1][m];
                        W_HAT[2][m] =  Z_BAR[2]*Z[0][m] + Z_BAR[0]*Z[2][m];
                        W_HAT[3][m] = (Z_BAR[3]*Z[0][m] + GAMMA_1*Z_BAR[1]*Z[1][m] + GAMMA_1*Z_BAR[2]*Z[2][m] + Z_BAR[0]*Z[3][m])/GAMMA;

#ifdef DEBUG
                        for(i=0;i<4;++i){std::cout << "W_HAT " << i << "\t" << m << " =\t" << W_HAT[i][m] << std::endl;}
                        std::cout << std::endl;
#endif
                }

                // Construct average state for element

                RHO   = pow((sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]))/3.0, 2);
                U     = (sqrt(U_N[0][0])*U_N[1][0]/U_N[0][0] + sqrt(U_N[0][1])*U_N[1][1]/U_N[0][1] + sqrt(U_N[0][2])*U_N[1][2]/U_N[0][2])/(sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]));        // U now represents x velocity
                V     = (sqrt(U_N[0][0])*U_N[2][0]/U_N[0][0] + sqrt(U_N[0][1])*U_N[2][1]/U_N[0][1] + sqrt(U_N[0][2])*U_N[2][2]/U_N[0][2])/(sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]));        // V represents y velocity
                H_AVG = (sqrt(U_N[0][0])*H[0] + sqrt(U_N[0][1])*H[1] + sqrt(U_N[0][2])*H[2])/(sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]));
               
                // E     = (sqrt(U_N[0][0])*H[0]/U_N[0][0] + sqrt(U_N[0][1])*H[1]/U_N[0][1] + sqrt(U_N[0][2])*H[2]/U_N[0][2])/(sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]));

                PRESSURE_AVG = (PRESSURE[0] + PRESSURE[1] + PRESSURE[2])/3.0;
                C_SOUND_AVG = sqrt((GAMMA-1.0) * H_AVG - (GAMMA-1.0) * (U*U + V*V)/2.0);
                // C_SOUND_AVG = sqrt(GAMMA*PRESSURE_AVG/RHO);

#ifdef DEBUG
                std::cout << "PRESSURE_AVG =\t" << PRESSURE_AVG << std::endl;
                std::cout << "C_SOUND_AVG  =\t" << C_SOUND_AVG  << std::endl;
#endif

                // Reassign variables to local equivalents

                C   = C_SOUND_AVG;

                U_C = U/C;

                V_C = V/C;

                H_C = H_AVG/C;

                ALPHA   = GAMMA_1*(U*U + V*V)/2.0;
                ALPHA_C = ALPHA/C;

#ifdef DEBUG
                std::cout << "U =\t" << U << "\tU_C =\t" << U_C << std::endl;
                std::cout << "V =\t" << V << "\tV_C =\t" << V_C << std::endl;
                std::cout << "H =\t" << H_AVG << "\tH_C =\t" << H_C << std::endl;
                std::cout << "ALPHA =\t" << ALPHA << "\tALPHA_C =\t" << ALPHA_C << std::endl;
                std::cout << std::endl;
#endif

                // Calculate K+,K- and K matrices for each vertex i,j,k

                for(m=0;m<3;++m){

                        W = U*N_X[m] + V*N_Y[m];

#ifdef DEBUG
                        std::cout << "W =\t" << W << std::endl;
#endif

                        LAMBDA[0][m] = W + C;
                        LAMBDA[1][m] = W - C;
                        LAMBDA[2][m] = W;
                        LAMBDA[3][m] = W;

                        for(i=0;i<4;++i){
                                LAMBDA_PLUS[i][m]  = max_val(0.0,LAMBDA[i][m]);
                                LAMBDA_MINUS[i][m] = min_val(0.0,LAMBDA[i][m]);
                        } 

                        for(p=0;p<3;++p){
                                if(p==0){       // Identify and select positive eigenvalues
                                        VALUE1 = LAMBDA_PLUS[0][m];
                                        VALUE2 = LAMBDA_PLUS[1][m];
                                        VALUE3 = LAMBDA_PLUS[2][m];
                                        VALUE4 = LAMBDA_PLUS[3][m];
#ifdef DEBUG
                                        std::cout << "K+" << std::endl;
#endif
                                }else if(p==1){ // Identify and select negative eigenvalues
                                        VALUE1 = LAMBDA_MINUS[0][m];
                                        VALUE2 = LAMBDA_MINUS[1][m];
                                        VALUE3 = LAMBDA_MINUS[2][m];
                                        VALUE4 = LAMBDA_MINUS[3][m];
#ifdef DEBUG
                                        std::cout << "K-" << std::endl;
#endif
                                }else{          // Select all eigenvalues
                                        VALUE1 = LAMBDA[0][m];
                                        VALUE2 = LAMBDA[1][m];
                                        VALUE3 = LAMBDA[2][m];
                                        VALUE4 = LAMBDA[3][m];
#ifdef DEBUG
                                        std::cout << "K" << std::endl;
#endif
                                }

#ifdef ASTRIX_COMP
                                if((X[m] >= 4.921 and X[m] <= 4.922) and (Y[m] >= 4.531 and Y[m] <= 4.532) and p==1){
                                        std::cout << "\tL   =\t" << VALUE1 << "\t" << VALUE2 << "\t" << VALUE3 << "\t" << VALUE4 << std::endl;
                                }
#endif

                                VALUE12  = (VALUE1 - VALUE2)/2.0;
                                VALUE123 = (VALUE1 + VALUE2 - 2.0*VALUE3)/2.0;

                                INFLOW[0][0][m][p] = 0.5*MAG[m]*(ALPHA_C*VALUE123/C - W*VALUE12/C + VALUE3);
                                INFLOW[0][1][m][p] = 0.5*MAG[m]*(-1.0*GAMMA_1*U_C*VALUE123/C + N_X[m]*VALUE12/C);
                                INFLOW[0][2][m][p] = 0.5*MAG[m]*(-1.0*GAMMA_1*V_C*VALUE123/C + N_Y[m]*VALUE12/C);
                                INFLOW[0][3][m][p] = 0.5*MAG[m]*(GAMMA_1*VALUE123/(C*C));

                                INFLOW[1][0][m][p] = 0.5*MAG[m]*((ALPHA_C*U_C - W*N_X[m])*VALUE123 + (ALPHA_C*N_X[m] - U_C*W)*VALUE12);
                                INFLOW[1][1][m][p] = 0.5*MAG[m]*((N_X[m]*N_X[m] - GAMMA_1*U_C*U_C)*VALUE123 - (GAMMA_2*U_C*N_X[m]*VALUE12) + VALUE3);
                                INFLOW[1][2][m][p] = 0.5*MAG[m]*((N_X[m]*N_Y[m] - GAMMA_1*U_C*V_C)*VALUE123 + (U_C*N_Y[m] - GAMMA_1*V_C*N_X[m])*VALUE12);
                                INFLOW[1][3][m][p] = 0.5*MAG[m]*(GAMMA_1*U_C*VALUE123/C + GAMMA_1*N_X[m]*VALUE12/C);

                                INFLOW[2][0][m][p] = 0.5*MAG[m]*((ALPHA_C*V_C - W*N_Y[m])*VALUE123 + (ALPHA_C*N_Y[m] - V_C*W)*VALUE12);
                                INFLOW[2][1][m][p] = 0.5*MAG[m]*((N_X[m]*N_Y[m] - GAMMA_1*U_C*V_C)*VALUE123 + (V_C*N_X[m] - GAMMA_1*U_C*N_Y[m])*VALUE12);
                                INFLOW[2][2][m][p] = 0.5*MAG[m]*((N_Y[m]*N_Y[m] - GAMMA_1*V_C*V_C)*VALUE123 - (GAMMA_2*V_C*N_Y[m]*VALUE12) + VALUE3);
                                INFLOW[2][3][m][p] = 0.5*MAG[m]*(GAMMA_1*V_C*VALUE123/C + GAMMA_1*N_Y[m]*VALUE12/C);

                                INFLOW[3][0][m][p] = 0.5*MAG[m]*((ALPHA_C*H_C - W*W)*VALUE123 + W*(ALPHA_C - H_C)*VALUE12);
                                INFLOW[3][1][m][p] = 0.5*MAG[m]*((W*N_X[m] - U - ALPHA_C*U_C)*VALUE123 + (H_C*N_X[m] - GAMMA_1*U_C*W)*VALUE12);
                                INFLOW[3][2][m][p] = 0.5*MAG[m]*((W*N_Y[m] - V - ALPHA_C*V_C)*VALUE123 + (H_C*N_Y[m] - GAMMA_1*V_C*W)*VALUE12);
                                INFLOW[3][3][m][p] = 0.5*MAG[m]*(GAMMA_1*H_C*VALUE123/C + GAMMA_1*W*VALUE12/C + VALUE3);

#ifdef DEBUG
                                for(i=0; i<4; ++i){
                                        for (j=0; j<4; ++j){
                                                std::cout << INFLOW[i][j][m][p] << "\t";
                                        }
                                        std::cout << std::endl;
                                }
                                std::cout << std::endl;
#endif
                        }
                }


#ifdef DEBUG
                std::cout << "Lambda + =\t" << LAMBDA_PLUS[0][0] << "\t" << LAMBDA_PLUS[0][1] << "\t" << LAMBDA_PLUS[0][2] << std::endl;
                std::cout << "Lambda + =\t" << LAMBDA_PLUS[1][0] << "\t" << LAMBDA_PLUS[1][1] << "\t" << LAMBDA_PLUS[1][2] << std::endl;
                std::cout << "Lambda + =\t" << LAMBDA_PLUS[2][0] << "\t" << LAMBDA_PLUS[2][1] << "\t" << LAMBDA_PLUS[2][2] << std::endl;
                std::cout << "Lambda + =\t" << LAMBDA_PLUS[3][0] << "\t" << LAMBDA_PLUS[3][1] << "\t" << LAMBDA_PLUS[3][2] << std::endl;

                std::cout << "Lambda - =\t" << LAMBDA_MINUS[0][0] << "\t" << LAMBDA_MINUS[0][1] << "\t" << LAMBDA_MINUS[0][2] << std::endl;
                std::cout << "Lambda - =\t" << LAMBDA_MINUS[1][0] << "\t" << LAMBDA_MINUS[1][1] << "\t" << LAMBDA_MINUS[1][2] << std::endl;
                std::cout << "Lambda - =\t" << LAMBDA_MINUS[2][0] << "\t" << LAMBDA_MINUS[2][1] << "\t" << LAMBDA_MINUS[2][2] << std::endl;
                std::cout << "Lambda - =\t" << LAMBDA_MINUS[3][0] << "\t" << LAMBDA_MINUS[3][1] << "\t" << LAMBDA_MINUS[3][2] << std::endl;

                std::cout << "Lambda   =\t" << LAMBDA[0][0] << "\t" << LAMBDA[0][1] << "\t" << LAMBDA[0][2] << std::endl;
                std::cout << "Lambda   =\t" << LAMBDA[1][0] << "\t" << LAMBDA[1][1] << "\t" << LAMBDA[1][2] << std::endl;
                std::cout << "Lambda   =\t" << LAMBDA[2][0] << "\t" << LAMBDA[2][1] << "\t" << LAMBDA[2][2] << std::endl;
                std::cout << "Lambda   =\t" << LAMBDA[3][0] << "\t" << LAMBDA[3][1] << "\t" << LAMBDA[3][2] << std::endl;
#endif

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                for(i=0;i<4;++i){
                        PHI[i] = 0.0;
#ifdef DEBUG
                        std::cout << "Calculating Phi " << i << std::endl;
#endif
                        for(m=0;m<3;++m){
                                PHI[i] += INFLOW[i][0][m][2]*W_HAT[0][m] + INFLOW[i][1][m][2]*W_HAT[1][m] + INFLOW[i][2][m][2]*W_HAT[2][m] + INFLOW[i][3][m][2]*W_HAT[3][m];
                        }
                }

#ifdef DEBUG
                std::cout << "PHI =\t" << PHI[0] << "\t" << PHI[1] << "\t" << PHI[2] << "\t" << PHI[3] << std::endl;
#endif

                for(i=0;i<4;++i){
                        for(j=0;j<4;++j){
                                INFLOW_MINUS_SUM[i][j] = 0.0;
                                for(m=0;m<3;++m){
                                        INFLOW_PLUS_SUM[i][j]  += INFLOW[i][j][m][0];
                                        INFLOW_MINUS_SUM[i][j] += INFLOW[i][j][m][1];
                                }
                        }
                }

                matInv(&INFLOW_MINUS_SUM[0][0],4,X[0],Y[0]);

                // /std::cout << "Post-inversion =" << std::endl;

#ifdef DEBUG
                for(i=0;i<4;++i){
                        for(j=0;j<4;++j){
                                std::cout << INFLOW_MINUS_SUM[i][j] << "\t";
                        }
                        std::cout << std::endl;
                }
#endif

                // Calculate spatial splitting for first half timestep

#ifdef LDA_SCHEME

#ifdef DEBUG
                std::cout << "BETA_0 =" << std::endl;
#endif

                for(i=0;i<4;++i){
                        for(j=0;j<4;++j){
                                for(m=0;m<3;++m){
                                        BETA[i][j][m] = -1.0*(INFLOW[i][0][m][0] * INFLOW_MINUS_SUM[0][j] + INFLOW[i][1][m][0] * INFLOW_MINUS_SUM[1][j] + INFLOW[i][2][m][0] * INFLOW_MINUS_SUM[2][j] + INFLOW[i][3][m][0] * INFLOW_MINUS_SUM[3][j]);
                                }
                        }
#ifdef DEBUG
                        std::cout << BETA[i][0][0] << "\t" << BETA[i][1][0] << "\t" << BETA[i][2][0] << "\t" << BETA[i][3][0] << std::endl;
#endif
                }

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                FLUC[i][m] = BETA[i][0][m] * PHI[0] + BETA[i][1][m] * PHI[1] + BETA[i][2][m] * PHI[2] + BETA[i][3][m] * PHI[3];
                        }
                }
#endif


#ifdef N_SCHEME
                double BRACKET[4][3];
                double KZ_SUM[4];

                for(i=0;i<4;++i){
                        KZ_SUM[i] = 0.0;
                        for(m=0;m<3;++m){
                                KZ_SUM[i] += INFLOW[i][0][m][1] * W_HAT[0][m] + INFLOW[i][1][m][1] * W_HAT[1][m] + INFLOW[i][2][m][1] * W_HAT[2][m] + INFLOW[i][3][m][1] * W_HAT[3][m];
                        }
                }

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                BRACKET[i][m] = W_HAT[i][m] - (INFLOW_MINUS_SUM[i][0]*KZ_SUM[0] + INFLOW_MINUS_SUM[i][1]*KZ_SUM[1] + INFLOW_MINUS_SUM[i][2]*KZ_SUM[2] + INFLOW_MINUS_SUM[i][3]*KZ_SUM[3]);
                        }
                }

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                FLUC[i][m] = INFLOW[i][0][m][0]*BRACKET[0][m] + INFLOW[i][1][m][0]*BRACKET[1][m] + INFLOW[i][2][m][0]*BRACKET[2][m] + INFLOW[i][3][m][0]*BRACKET[3][m];
                        }
                }
#endif

                DUAL[0] = VERTEX_0->get_dual();
                DUAL[1] = VERTEX_1->get_dual();
                DUAL[2] = VERTEX_2->get_dual();

                // Calculate change to be distributed

                for(i=0;i<4;i++){
                        DU0[i] = DT*FLUC[i][0]/DUAL[0];
                        DU1[i] = DT*FLUC[i][1]/DUAL[1];
                        DU2[i] = DT*FLUC[i][2]/DUAL[2];
                }

                VERTEX_0->update_du_half(DU0);
                VERTEX_1->update_du_half(DU1);
                VERTEX_2->update_du_half(DU2);

#ifdef ASTRIX_COMP
                for(m=0;m<3;m++){
                        // std::cout << m << "\t" << X[m] << "\t" << Y[m] << std::endl;
                        if((X[m] >= 4.921 and X[m] <= 4.922) and (Y[m] >= 4.531 and Y[m] <= 4.532)){
                                int K_INDEX = 2;
                                if(m == 0){
                                        // std::cout << "+\t" << m << "\tX =\t" << X[m] << "\tY =\t" << Y[m] << "\tDU0 =\t" << DU0[0] << std::endl;
                                        std::cout << "+\t" << m << "\tX =\t" << X[m] << "\tY =\t" << Y[m] << "\tFLUC0 =\t" << FLUC[0][m] << std::endl;
                                        std::cout << "\tDX  =\t" << DX << "\tDY =\t" << DY << std::endl;
                                        std::cout << "\tN_x =\t" << N_X[0] << "\t" << N_X[1] << "\t" << N_X[2] << std::endl;
                                        std::cout << "\tN_y =\t" << N_Y[0] << "\t" << N_Y[1] << "\t" << N_Y[2] << std::endl;
                                        std::cout << "\tW_h =\t" << W_HAT[0][m] << "\t" << W_HAT[1][m] << "\t" << W_HAT[2][m] << "\t" << W_HAT[3][m] << "\t" << std::endl;
                                        std::cout << "\tPHI =\t" << PHI[0] << "\t" << PHI[1] << "\t" << PHI[2] << "\t" << PHI[3] << std::endl;
                                        std::cout << "\tDU0 =\t" << DU0[0] << "\t" << DU0[1] << "\t" << DU0[2] << "\t" << DU0[3] << std::endl;
                                        // std::cout << "\tBETA =\t" << std::endl;
                                        // for(i=0;i<4;++i){std::cout << BETA[i][0][m] << "\t" << BETA[i][1][m] << "\t" << BETA[i][2][m] << "\t" << BETA[i][3][m] << std::endl;}
                                        std::cout << "\tK =\t" << std::endl;
                                        for(i=0;i<4;++i){std::cout << INFLOW[i][0][m][K_INDEX] << "\t" << INFLOW[i][1][m][K_INDEX] << "\t" << INFLOW[i][2][m][K_INDEX] << "\t" << INFLOW[i][3][m][K_INDEX] << std::endl;}
                                        std::cout << std::endl; 
                                        // std::cout << "DT =\t" << DT << "\tDUAL =\t" << DUAL[m] << std::endl;
                                }else if(m == 1){
                                        // std::cout << "+\t" << m << "\tX =\t" << X[m] << "\tY =\t" << Y[m] << "\tDU1 =\t" << DU1[0] << std::endl;
                                        std::cout << "+\t" << m << "\tX =\t" << X[m] << "\tY =\t" << Y[m] << "\tFLUC1 =\t" << FLUC[0][m] << std::endl;
                                        std::cout << "\tN_x =\t" << N_X[0] << "\t" << N_X[1] << "\t" << N_X[2] << std::endl;
                                        std::cout << "\tN_y =\t" << N_Y[0] << "\t" << N_Y[1] << "\t" << N_Y[2] << std::endl;
                                        std::cout << "\tW_h =\t" << W_HAT[0][m] << "\t" << W_HAT[1][m] << "\t" << W_HAT[2][m] << "\t" << W_HAT[3][m] << "\t" << std::endl;
                                        std::cout << "\tPHI =\t" << PHI[0] << "\t" << PHI[1] << "\t" << PHI[2] << "\t" << PHI[3] << std::endl;
                                        std::cout << "\tDU1 =\t" << DU1[0] << "\t" << DU1[1] << "\t" << DU1[2] << "\t" << DU1[3] << std::endl;
                                        // std::cout << "\tBETA =\t" << std::endl;
                                        // for(i=0;i<4;++i){std::cout << BETA[i][0][m] << "\t" << BETA[i][1][m] << "\t" << BETA[i][2][m] << "\t" << BETA[i][3][m] << std::endl;}
                                        std::cout << "\tK =\t" << std::endl;
                                        for(i=0;i<4;++i){std::cout << INFLOW[i][0][m][K_INDEX] << "\t" << INFLOW[i][1][m][K_INDEX] << "\t" << INFLOW[i][2][m][K_INDEX] << "\t" << INFLOW[i][3][m][K_INDEX] << std::endl;}
                                        std::cout << std::endl;
                                        // std::cout << "DT =\t" << DT << "\tDUAL =\t" << DUAL[m] << std::endl;
                                }else{
                                        // std::cout << "+\t" << m << "\tX =\t" << X[m] << "\tY =\t" << Y[m] << "\tDU2 =\t" << DU2[0] << std::endl;
                                        std::cout << "+\t" << m << "\tX =\t" << X[m] << "\tY =\t" << Y[m] << "\tFLUC2 =\t" << FLUC[0][m] << std::endl;
                                        std::cout << "\tN_x =\t" << N_X[0] << "\t" << N_X[1] << "\t" << N_X[2] << std::endl;
                                        std::cout << "\tN_y =\t" << N_Y[0] << "\t" << N_Y[1] << "\t" << N_Y[2] << std::endl;
                                        std::cout << "\tW_h  =\t" << W_HAT[0][m] << "\t" << W_HAT[1][m] << "\t" << W_HAT[2][m] << "\t" << W_HAT[3][m] << "\t" << std::endl;
                                        std::cout << "\tPHI =\t" << PHI[0] << "\t" << PHI[1] << "\t" << PHI[2] << "\t" << PHI[3] << std::endl;
                                        std::cout << "\tDU2 =\t" << DU2[0] << "\t" << DU2[1] << "\t" << DU2[2] << "\t" << DU2[3] << std::endl;
                                        // std::cout << "\tBETA =\t" << std::endl;
                                        // for(i=0;i<4;++i){std::cout << BETA[i][0][m] << "\t" << BETA[i][1][m] << "\t" << BETA[i][2][m] << "\t" << BETA[i][3][m] << std::endl;}
                                        std::cout << "\tK =\t" << std::endl;
                                        for(i=0;i<4;++i){std::cout << INFLOW[i][0][m][K_INDEX] << "\t" << INFLOW[i][1][m][K_INDEX] << "\t" << INFLOW[i][2][m][K_INDEX] << "\t" << INFLOW[i][3][m][K_INDEX] << std::endl;}
                                        std::cout << std::endl;
                                        // std::cout << "DT =\t" << DT << "\tDUAL =\t" << DUAL[m] << std::endl;
                                }
                                // std::cout << "\tKZ_SUM =\t" << KZ_SUM[0] << "\t" << KZ_SUM[1] << "\t" << KZ_SUM[2] << "\t" << KZ_SUM[3] << std::endl;
                                // std::cout << "\tBRACKET =\t" << BRACKET[0][m] << "\t" << BRACKET[1][m] << "\t" << BRACKET[2][m] << "\t" << BRACKET[3][m] << std::endl;
                                // std::cout << "0\t" << BRACKET[0][0] << "\t" << BRACKET[1][0] << "\t" << BRACKET[2][0] << "\t" << BRACKET[3][0] << std::endl;
                                // std::cout << "1\t" << BRACKET[0][1] << "\t" << BRACKET[1][1] << "\t" << BRACKET[2][1] << "\t" << BRACKET[3][1] << std::endl;
                                // std::cout << "2\t" << BRACKET[0][2] << "\t" << BRACKET[1][2] << "\t" << BRACKET[2][2] << "\t" << BRACKET[3][2] << std::endl;
                                std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl; 
                                // std::cout << "\tN_m =\t" << std::endl;
                                // for(i=0;i<4;++i){std::cout << INFLOW_MINUS_SUM[i][0] << "\t" << INFLOW_MINUS_SUM[i][1] << "\t" << INFLOW_MINUS_SUM[i][2] << "\t" << INFLOW_MINUS_SUM[i][3] << std::endl;}        
                        }
                }
#endif


#ifdef DEBUG
                        for(i=0;i<4;i++){std::cout << "Element fluctuation =\t" << FLUC[i][0] << "\t" << FLUC[i][1] << "\t" << FLUC[i][2] << std::endl;}
                        std::cout << "Dual =\t" << VERTEX_0->get_dual() << "\t" << VERTEX_1->get_dual() << "\t" << VERTEX_2->get_dual() << std::endl;
                        std::cout << "Change (rho) =\t"    << DU0[0] << "\t" << DU1[0] << "\t" << DU2[0] << std::endl;
                        std::cout << "Change (x mom) =\t"  << DU0[1] << "\t" << DU1[1] << "\t" << DU2[1] << std::endl;
                        std::cout << "Change (y mom) =\t"  << DU0[2] << "\t" << DU1[2] << "\t" << DU2[2] << std::endl;
                        std::cout << "Change (energy) =\t" << DU0[3] << "\t" << DU1[3] << "\t" << DU2[3] << std::endl;
                        std::cout << "-----------------------------------------------------------------" << std::endl;
                        // if(U_N[0][0] != U_N[0][1] or U_N[0][0] != U_N[0][2] or U_N[0][1] != U_N[0][2]){exit(0);}
#endif

                return ;
        }

        //**********************************************************************************************************************

        void calculate_second_half(double T, double DT_TOT, double DX, double DY){
                int i,j,m,p;

                double DU0[4],DU1[4],DU2[4];
                double INFLOW[4][4][3][3];

                double DT = DT_TOT;

                double C_SOUND[3];

                setup_half_state();

#ifdef CLOSED
                if(std::abs(X[0] - X[1]) > 2.0*DX or std::abs(X[0] - X[2]) > 2.0*DX or std::abs(X[1] - X[2]) > 2.0*DX){
#ifdef DEBUG
                        std::cout << "Skipping x boundary\t" << (X[0]+X[1]+X[2])/3.0 << "\t" << (Y[0]+Y[1]+Y[2])/3.0 << "\t" << std::endl;
#endif
                        return ;
                }else if(std::abs(Y[0] - Y[1]) > 2.0*DY or std::abs(Y[0] - Y[2]) > 2.0*DY or std::abs(Y[1] - Y[2]) > 2.0*DY){
#ifdef DEBUG
                        std::cout << "Skipping y boundary\t" << (X[0]+X[1]+X[2])/3.0 << "\t" << (Y[0]+Y[1]+Y[2])/3.0 << "\t" << std::endl;
#endif
                        return ;
                }
#endif

                // Calcualte sound speed for half time state

                for(i=0;i<3;i++){C_SOUND[i] = sqrt(GAMMA*PRESSURE[i]/U_HALF[0][i]);}

#ifdef DEBUG
                std::cout << "-- SECOND -------------------------------------------------------" << std::endl;
                std::cout << "Time     =\t" << T << std::endl;
                std::cout << "0 =\t" << X[0] << "\t" << Y[0] << std::endl;
                std::cout << "1 =\t" << X[1] << "\t" << Y[1] << std::endl;
                std::cout << "2 =\t" << X[2] << "\t" << Y[2] << std::endl;
                std::cout << "Half State    =" << "\trho" << "\tx_mom" << "\ty_mom" << "\tenergy" << std::endl;
                for(i=0;i<3;i++){std::cout << i << " =\t" << U_HALF[0][i] << "\t" << U_HALF[1][i] << "\t" << U_HALF[2][i] << "\t" << U_HALF[3][i] << std::endl;}
                std::cout << "Pressure =\t" << PRESSURE[0] << "\t" << PRESSURE[1] << "\t" << PRESSURE[2] << std::endl;
#endif


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // Calculate inflow parameters

                double H[3];
                double RHO,C,U,U_C,V,V_C,H_AVG,H_C,ALPHA,ALPHA_C,W;
                double Z[4][3];
                double VALUE1,VALUE2,VALUE3,VALUE4,VALUE12,VALUE123;
                double PRESSURE_AVG,C_SOUND_AVG;
                double LAMBDA[4][3],LAMBDA_PLUS[4][3],LAMBDA_MINUS[4][3];
                double N_X[3],N_Y[3];

                double Z_BAR[4],W_HAT[4][3];

                // Construct Roe vector Z

                for(m=0;m<3;++m){
                        Z[0][m] = sqrt(U_N[0][m]);
                        Z[1][m] = U_N[1][m]/Z[0][m];
                        Z[2][m] = U_N[2][m]/Z[0][m];
                        Z[3][m] = (U_N[3][m] + PRESSURE[m])/Z[0][m];

                        N_X[m]  = NORMAL[m][0];
                        N_Y[m]  = NORMAL[m][1];

                        H[m] = (U_N[3][m] + PRESSURE[m])/U_N[0][m];
#ifdef DEBUG
                        std::cout << "Z =\t" << Z[0][m] << "\t" << Z[1][m] << "\t" << Z[2][m] << "\t" << Z[3][m] << std::endl;
#endif
                }

                for(i=0; i<4; ++i){Z_BAR[i] = (Z[i][0] + Z[i][1] + Z[i][2])/3.0;}

#ifdef DEBUG
                for(i=0; i<4; ++i){std::cout << "Z_BAR " << i << " =\t" << Z_BAR[i] << std::endl;}
                std::cout << std::endl;
#endif


                for(m=0; m<3; ++m){
                        W_HAT[0][m] =  2.0*Z_BAR[0]*Z[0][m];
                        W_HAT[1][m] =  Z_BAR[1]*Z[0][m] + Z_BAR[0]*Z[1][m];
                        W_HAT[2][m] =  Z_BAR[2]*Z[0][m] + Z_BAR[0]*Z[2][m];
                        W_HAT[3][m] = (Z_BAR[3]*Z[0][m] + GAMMA_1*Z_BAR[1]*Z[1][m] + GAMMA_1*Z_BAR[2]*Z[2][m] + Z_BAR[0]*Z[3][m])/GAMMA;

#ifdef DEBUG
                        for(i=0;i<4;++i){std::cout << "W_HAT " << i << "\t" << m << " =\t" << W_HAT[i][m] << std::endl;}
                        std::cout << std::endl;
#endif
                }

                // Construct average state for element

                RHO   = pow((sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]))/3.0, 2);
                U     = (sqrt(U_N[0][0])*U_N[1][0]/U_N[0][0] + sqrt(U_N[0][1])*U_N[1][1]/U_N[0][1] + sqrt(U_N[0][2])*U_N[1][2]/U_N[0][2])/(sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]));        // U now represents x velocity
                V     = (sqrt(U_N[0][0])*U_N[2][0]/U_N[0][0] + sqrt(U_N[0][1])*U_N[2][1]/U_N[0][1] + sqrt(U_N[0][2])*U_N[2][2]/U_N[0][2])/(sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]));        // V represents y velocity
                H_AVG = (sqrt(U_N[0][0])*H[0] + sqrt(U_N[0][1])*H[1] + sqrt(U_N[0][2])*H[2])/(sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]));
               
                // E     = (sqrt(U_N[0][0])*H[0]/U_N[0][0] + sqrt(U_N[0][1])*H[1]/U_N[0][1] + sqrt(U_N[0][2])*H[2]/U_N[0][2])/(sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]));

                PRESSURE_AVG = (PRESSURE[0] + PRESSURE[1] + PRESSURE[2])/3.0;
                C_SOUND_AVG = sqrt((GAMMA-1.0) * H_AVG - (GAMMA-1.0) * (U*U + V*V)/2.0);
                // C_SOUND_AVG = sqrt(GAMMA*PRESSURE_AVG/RHO);

#ifdef DEBUG
                std::cout << "PRESSURE_AVG =\t" << PRESSURE_AVG << std::endl;
                std::cout << "C_SOUND_AVG  =\t" << C_SOUND_AVG  << std::endl;
#endif

                // Reassign variables to local equivalents

                C   = C_SOUND_AVG;

                U_C = U/C;

                V_C = V/C;

                H_C = H_AVG/C;

                ALPHA   = GAMMA_1*(U*U + V*V)/2.0;
                ALPHA_C = ALPHA/C;

#ifdef DEBUG
                std::cout << "U =\t" << U << "\tU_C =\t" << U_C << std::endl;
                std::cout << "V =\t" << V << "\tV_C =\t" << V_C << std::endl;
                std::cout << "H =\t" << H_AVG << "\tH_C =\t" << H_C << std::endl;
                std::cout << "ALPHA =\t" << ALPHA << "\tALPHA_C =\t" << ALPHA_C << std::endl;
                std::cout << std::endl;
#endif

                // Calculate K+,K- and K matrices for each vertex i,j,k

                for(m=0;m<3;++m){

                        W = U*N_X[m] + V*N_Y[m];

#ifdef DEBUG
                        std::cout << "W =\t" << W << std::endl;
#endif

                        LAMBDA[0][m] = W + C;
                        LAMBDA[1][m] = W - C;
                        LAMBDA[2][m] = W;
                        LAMBDA[3][m] = W;

                        for(i=0;i<4;++i){
                                LAMBDA_PLUS[i][m]  = max_val(0.0,LAMBDA[i][m]);
                                LAMBDA_MINUS[i][m] = min_val(0.0,LAMBDA[i][m]);
                        } 

                        for(p=0;p<3;++p){
                                if(p==0){       // Identify and select positive eigenvalues
                                        VALUE1 = LAMBDA_PLUS[0][m];
                                        VALUE2 = LAMBDA_PLUS[1][m];
                                        VALUE3 = LAMBDA_PLUS[2][m];
                                        VALUE4 = LAMBDA_PLUS[3][m];
#ifdef DEBUG
                                        std::cout << "K+" << std::endl;
#endif
                                }else if(p==1){ // Identify and select negative eigenvalues
                                        VALUE1 = LAMBDA_MINUS[0][m];
                                        VALUE2 = LAMBDA_MINUS[1][m];
                                        VALUE3 = LAMBDA_MINUS[2][m];
                                        VALUE4 = LAMBDA_MINUS[3][m];
#ifdef DEBUG
                                        std::cout << "K-" << std::endl;
#endif
                                }else{          // Select all eigenvalues
                                        VALUE1 = LAMBDA[0][m];
                                        VALUE2 = LAMBDA[1][m];
                                        VALUE3 = LAMBDA[2][m];
                                        VALUE4 = LAMBDA[3][m];
#ifdef DEBUG
                                        std::cout << "K" << std::endl;
#endif
                                }

                                VALUE12  = (VALUE1 - VALUE2)/2.0;
                                VALUE123 = (VALUE1 + VALUE2 - 2.0*VALUE3)/2.0;

                                INFLOW[0][0][m][p] = 0.5*MAG[m]*(ALPHA_C*VALUE123/C - W*VALUE12/C + VALUE3);
                                INFLOW[0][1][m][p] = 0.5*MAG[m]*(-1.0*GAMMA_1*U_C*VALUE123/C + N_X[m]*VALUE12/C);
                                INFLOW[0][2][m][p] = 0.5*MAG[m]*(-1.0*GAMMA_1*V_C*VALUE123/C + N_Y[m]*VALUE12/C);
                                INFLOW[0][3][m][p] = 0.5*MAG[m]*(GAMMA_1*VALUE123/(C*C));

                                INFLOW[1][0][m][p] = 0.5*MAG[m]*((ALPHA_C*U_C - W*N_X[m])*VALUE123 + (ALPHA_C*N_X[m] - U_C*W)*VALUE12);
                                INFLOW[1][1][m][p] = 0.5*MAG[m]*((N_X[m]*N_X[m] - GAMMA_1*U_C*U_C)*VALUE123 - (GAMMA_2*U_C*N_X[m]*VALUE12) + VALUE3);
                                INFLOW[1][2][m][p] = 0.5*MAG[m]*((N_X[m]*N_Y[m] - GAMMA_1*U_C*V_C)*VALUE123 + (U_C*N_Y[m] - GAMMA_1*V_C*N_X[m])*VALUE12);
                                INFLOW[1][3][m][p] = 0.5*MAG[m]*(GAMMA_1*U_C*VALUE123/C + GAMMA_1*N_X[m]*VALUE12/C);

                                INFLOW[2][0][m][p] = 0.5*MAG[m]*((ALPHA_C*V_C - W*N_Y[m])*VALUE123 + (ALPHA_C*N_Y[m] - V_C*W)*VALUE12);
                                INFLOW[2][1][m][p] = 0.5*MAG[m]*((N_X[m]*N_Y[m] - GAMMA_1*U_C*V_C)*VALUE123 + (V_C*N_X[m] - GAMMA_1*U_C*N_Y[m])*VALUE12);
                                INFLOW[2][2][m][p] = 0.5*MAG[m]*((N_Y[m]*N_Y[m] - GAMMA_1*V_C*V_C)*VALUE123 - (GAMMA_2*V_C*N_Y[m]*VALUE12) + VALUE3);
                                INFLOW[2][3][m][p] = 0.5*MAG[m]*(GAMMA_1*V_C*VALUE123/C + GAMMA_1*N_Y[m]*VALUE12/C);

                                INFLOW[3][0][m][p] = 0.5*MAG[m]*((ALPHA_C*H_C - W*W)*VALUE123 + W*(ALPHA_C - H_C)*VALUE12);
                                INFLOW[3][1][m][p] = 0.5*MAG[m]*((W*N_X[m] - U - ALPHA_C*U_C)*VALUE123 + (H_C*N_X[m] - GAMMA_1*U_C*W)*VALUE12);
                                INFLOW[3][2][m][p] = 0.5*MAG[m]*((W*N_Y[m] - V - ALPHA_C*V_C)*VALUE123 + (H_C*N_Y[m] - GAMMA_1*V_C*W)*VALUE12);
                                INFLOW[3][3][m][p] = 0.5*MAG[m]*(GAMMA_1*H_C*VALUE123/C + GAMMA_1*W*VALUE12/C + VALUE3);

#ifdef DEBUG
                                for(i=0; i<4; ++i){
                                        for (j=0; j<4; ++j){
                                                std::cout << INFLOW[i][j][m][p] << "\t";
                                        }
                                        std::cout << std::endl;
                                }
                                std::cout << std::endl;
#endif
                        }
                }


#ifdef DEBUG
                std::cout << "Lambda + =\t" << LAMBDA_PLUS[0][0] << "\t" << LAMBDA_PLUS[0][1] << "\t" << LAMBDA_PLUS[0][2] << std::endl;
                std::cout << "Lambda + =\t" << LAMBDA_PLUS[1][0] << "\t" << LAMBDA_PLUS[1][1] << "\t" << LAMBDA_PLUS[1][2] << std::endl;
                std::cout << "Lambda + =\t" << LAMBDA_PLUS[2][0] << "\t" << LAMBDA_PLUS[2][1] << "\t" << LAMBDA_PLUS[2][2] << std::endl;
                std::cout << "Lambda + =\t" << LAMBDA_PLUS[3][0] << "\t" << LAMBDA_PLUS[3][1] << "\t" << LAMBDA_PLUS[3][2] << std::endl;

                std::cout << "Lambda - =\t" << LAMBDA_MINUS[0][0] << "\t" << LAMBDA_MINUS[0][1] << "\t" << LAMBDA_MINUS[0][2] << std::endl;
                std::cout << "Lambda - =\t" << LAMBDA_MINUS[1][0] << "\t" << LAMBDA_MINUS[1][1] << "\t" << LAMBDA_MINUS[1][2] << std::endl;
                std::cout << "Lambda - =\t" << LAMBDA_MINUS[2][0] << "\t" << LAMBDA_MINUS[2][1] << "\t" << LAMBDA_MINUS[2][2] << std::endl;
                std::cout << "Lambda - =\t" << LAMBDA_MINUS[3][0] << "\t" << LAMBDA_MINUS[3][1] << "\t" << LAMBDA_MINUS[3][2] << std::endl;

                std::cout << "Lambda   =\t" << LAMBDA[0][0] << "\t" << LAMBDA[0][1] << "\t" << LAMBDA[0][2] << std::endl;
                std::cout << "Lambda   =\t" << LAMBDA[1][0] << "\t" << LAMBDA[1][1] << "\t" << LAMBDA[1][2] << std::endl;
                std::cout << "Lambda   =\t" << LAMBDA[2][0] << "\t" << LAMBDA[2][1] << "\t" << LAMBDA[2][2] << std::endl;
                std::cout << "Lambda   =\t" << LAMBDA[3][0] << "\t" << LAMBDA[3][1] << "\t" << LAMBDA[3][2] << std::endl;
#endif

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                double PHI_HALF[4];

#ifdef LDA_SCHEME
                for(i=0;i<4;++i){
                        PHI_HALF[i] = 0.0;
                        for(m=0;m<3;++m){
                                PHI_HALF[i] += INFLOW[i][0][m][2]*W_HAT[0][m] + INFLOW[i][1][m][2]*W_HAT[1][m] + INFLOW[i][2][m][2]*W_HAT[2][m] + INFLOW[i][3][m][2]*W_HAT[3][m];
                        }
                }
§
                // Calculate spatial splitting for first half timestep

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                FLUC_HALF[i][m] = BETA[i][0][m] * (PHI_HALF[0]) + BETA[i][1][m] * (PHI_HALF[1]) + BETA[i][2][m] * (PHI_HALF[2]) + BETA[i][3][m] * (PHI_HALF[3]);
                        }
                }

#ifdef DEBUG
                for(m=0;m<3;++m){
                        for(i=0;i<4;++i){
                                std::cout << "BETA =\t" << m << "\t" <<  BETA[i][0][m] << "\t" << BETA[i][1][m] << "\t" << BETA[i][2][m] << "\t" << BETA[i][3][m] << std::endl;
                        }
                }
#endif

                double MASS[4][4][3];
                double DIFF[4][3];
                double MASS_DIFF[4][3];
                double SUM_MASS[4];

                for(i=0;i<4;++i){
                        for(j=0;j<4;++j){
                                for(m=0;m<3;++m){MASS[i][j][m]= AREA * BETA[i][j][m]/3.0;}
                        }
                }

#ifdef DEBUG
                for(m=0;m<3;++m){
                        for(i=0;i<4;++i){
                                std::cout << "MASS =\t" << i << "\t" <<  MASS[i][0][m] << "\t" << MASS[i][1][m] << "\t" << MASS[i][2][m] << "\t" << MASS[i][3][m] << std::endl;
                        }
                }
#endif

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                DIFF[i][m] = U_HALF[i][m] - U_N[i][m];
                        }
                }

#ifdef DEBUG
                for(i=0;i<4;++i){
                        std::cout << "DIFF =\t" << i << "\t" <<  DIFF[i][0] << "\t" << DIFF[i][1] << "\t" << DIFF[i][2] << std::endl;
                }
#endif

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                MASS_DIFF[i][m] = MASS[i][0][m] * DIFF[0][m] + MASS[i][1][m] * DIFF[1][m] + MASS[i][2][m] * DIFF[2][m] + MASS[i][3][m] * DIFF[3][m];
                        }
                }

#ifdef DEBUG
                for(i=0;i<4;++i){
                        std::cout << "MASS_DIFF =\t" << i << "\t" <<  MASS_DIFF[i][0] << "\t" << MASS_DIFF[i][1] << "\t" << MASS_DIFF[i][2] << std::endl;
                }
#endif

                for(i=0;i<4;i++){
                        SUM_MASS[i] = 0.0;
                        for(m=0;m<3;++m){
                                SUM_MASS[i] += MASS_DIFF[i][m]/DT;
                        }
                }

                for(i=0;i<4;i++){
                        DU0[i] = (DT/DUAL[0])*(SUM_MASS[i] + 0.5*(FLUC[i][0] + FLUC_HALF[i][0]));
                        DU1[i] = (DT/DUAL[1])*(SUM_MASS[i] + 0.5*(FLUC[i][1] + FLUC_HALF[i][1]));
                        DU2[i] = (DT/DUAL[2])*(SUM_MASS[i] + 0.5*(FLUC[i][2] + FLUC_HALF[i][2]));
                }
#endif

#ifdef N_SCHEME

                double INFLOW_PLUS_SUM[4][4], INFLOW_MINUS_SUM[4][4];

                for(i=0;i<4;++i){
                        for(j=0;j<4;++j){
                                INFLOW_MINUS_SUM[i][j] = 0.0;
                                for(m=0;m<3;++m){
                                        INFLOW_PLUS_SUM[i][j] += INFLOW[i][j][m][0];
                                        INFLOW_MINUS_SUM[i][j] += INFLOW[i][j][m][1];
                                }
                        }
                }

                matInv(&INFLOW_MINUS_SUM[0][0],4,X[0],Y[0]);

                double BRACKET[4][3];
                double KZ_SUM[4];
                double DIFF[4][3];

                for(i=0;i<4;++i){
                        KZ_SUM[i] = 0.0;
                        for(m=0;m<3;++m){
                                KZ_SUM[i] += INFLOW[i][0][m][1] * W_HAT[0][m] + INFLOW[i][1][m][1] * W_HAT[1][m] + INFLOW[i][2][m][1] * W_HAT[2][m] + INFLOW[i][3][m][1] * W_HAT[3][m];
                        }
                }

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                BRACKET[i][m] = W_HAT[i][m] - (INFLOW_MINUS_SUM[i][0]*KZ_SUM[0] + INFLOW_MINUS_SUM[i][1]*KZ_SUM[1] + INFLOW_MINUS_SUM[i][2]*KZ_SUM[2] + INFLOW_MINUS_SUM[i][3]*KZ_SUM[3]);
                        }
                }

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                FLUC_HALF[i][m] = INFLOW[i][0][m][0]*BRACKET[0][m] + INFLOW[i][1][m][0]*BRACKET[1][m] + INFLOW[i][2][m][0]*BRACKET[2][m] + INFLOW[i][3][m][0]*BRACKET[3][m];
                        }
                }

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                DIFF[i][m] = AREA*(U_HALF[i][m] - U_N[i][m])/3.0;
                        }
                }

                for(i=0;i<4;i++){
                        DU0[i] = (DT/DUAL[0])*(DIFF[i][0]/DT + 0.5*(FLUC[i][0] + FLUC_HALF[i][0]));
                        DU1[i] = (DT/DUAL[1])*(DIFF[i][1]/DT + 0.5*(FLUC[i][1] + FLUC_HALF[i][1]));
                        DU2[i] = (DT/DUAL[2])*(DIFF[i][2]/DT + 0.5*(FLUC[i][2] + FLUC_HALF[i][2]));

                        DU0[i] = 0.0;
                        DU1[i] = 0.0;
                        DU2[i] = 0.0;

                }

#endif

                VERTEX_0->update_du(DU0);
                VERTEX_1->update_du(DU1);
                VERTEX_2->update_du(DU2);

#ifdef DEBUG
                        for(i=0;i<4;i++){std::cout << "Element fluctuation =\t" << FLUC_HALF[i][0] << "\t" << FLUC_HALF[i][1] << "\t" << FLUC_HALF[i][2] << std::endl;}
                        std::cout << "Change (rho)    =\t" << DU0[0] << "\t" << DU1[0] << "\t" << DU2[0] << std::endl;
                        std::cout << "Change (x mom)  =\t" << DU0[1] << "\t" << DU1[1] << "\t" << DU2[1] << std::endl;
                        std::cout << "Change (y mom)  =\t" << DU0[2] << "\t" << DU1[2] << "\t" << DU2[2] << std::endl;
                        std::cout << "Change (energy) =\t" << DU0[3] << "\t" << DU1[3] << "\t" << DU2[3] << std::endl;
                        std::cout << "-----------------------------------------------------------------" << std::endl;
                        // if(U_N[0][0] != U_N[0][1] or U_N[0][0] != U_N[0][2] or U_N[0][1] != U_N[0][2]){exit(0);}
#endif

                 return ;
        }

        // Returns Roe average of left and right states
        double roe_avg(double L1, double L2, double R1, double R2){
                double AVG;
                AVG = (sqrt(L1)*L2+sqrt(R1)*R2)/(sqrt(L1)+sqrt(R1));
                return AVG;
        }

        void calculate_normals(double X[3],double Y[3]){
                int i;
                double PERP[3][2];

                PERP[0][0] = (Y[1] - Y[2]);
                PERP[0][1] = (X[2] - X[1]);

                PERP[1][0] = (Y[2] - Y[0]);
                PERP[1][1] = (X[0] - X[2]);

                PERP[2][0] = (Y[0] - Y[1]);
                PERP[2][1] = (X[1] - X[0]);


                // calculate area of triangle and pass 1/3 to each vertex for dual
                double THETA0 = atan(PERP[0][1]/PERP[0][0]);
                double THETA1 = atan(PERP[1][1]/PERP[1][0]);

                double THETA = std::abs(THETA0 - THETA1);

                if(THETA > 3.14159/2.0){THETA = 3.14159 - THETA;}

                AREA = 0.5*(sqrt(PERP[0][0]*PERP[0][0] + PERP[0][1]*PERP[0][1])*sqrt(PERP[1][0]*PERP[1][0] + PERP[1][1]*PERP[1][1]))*sin(THETA);

                // VERTEX_0->calculate_dual(AREA/3.0);
                // VERTEX_1->calculate_dual(AREA/3.0);
                // VERTEX_2->calculate_dual(AREA/3.0);

                for(i=0;i<3;i++){
                        MAG[i] = sqrt(PERP[i][0]*PERP[i][0]+PERP[i][1]*PERP[i][1]);
                        NORMAL[i][0] = PERP[i][0]/MAG[i];
                        NORMAL[i][1] = PERP[i][1]/MAG[i];
                }

#ifdef DEBUG
                        std::cout << "X =\t" << X[0] << "\t" << X[1] << "\t" << X[2] << std::endl;
                        std::cout << "Y =\t" << Y[0] << "\t" << Y[1] << "\t" << Y[2] << std::endl;
                        std::cout << "Normal i =\t" << NORMAL[0][0] << "\t" << NORMAL[0][1] << std::endl;
                        std::cout << "Normal j =\t" << NORMAL[1][0] << "\t" << NORMAL[1][1] << std::endl;
                        std::cout << "Normal k =\t" << NORMAL[2][0] << "\t" << NORMAL[2][1] << std::endl;
#endif

                return ;
        }

        void setup_normals(double DX, double DY){
                // Calculate normals (just in first timestep for static grid)

                double X_MOD[3],Y_MOD[3];

                setup_positions();

                for(int m=0; m<3; ++m){X_MOD[m] = X[m];Y_MOD[m] = Y[m];}
                for(int i=0; i<3; ++i){
                        for(int j=0; j<3; ++j){
                                if(X[j] - X[i] > 2.0*DX){
                                        X_MOD[i] = X[i] + SIDE_LENGTH_X;
                                }
                                if(Y[j] - Y[i] > 2.0*DY){
                                        Y_MOD[i] = Y[i] + SIDE_LENGTH_Y;
                                }
                        }
                }
                calculate_normals(X_MOD,Y_MOD);

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