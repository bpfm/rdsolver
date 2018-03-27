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

        double LAMBDA[4][3][2];
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
                int i,j,k;
                double DU0[4],DU1[4],DU2[4];
                double INFLOW[4][3],INFLOW_PLUS[4][3],INFLOW_MINUS[4][3];
                double IN_TOP,IN_BOTTOM,OUT_TOP,OUT_BOTTOM;
                double DT = DT_TOT;

                // Import conditions and positions of vertices

                setup_positions();
                setup_initial_state();

                if(abs(X[0] - X[1]) > 10.0 or abs(X[0] - X[2]) > 10.0 or abs(X[1] -X[2]) > 10.0){
#ifdef DEBUG
                        cout << "Skipping x boundary" << endl;
#endif
                        return ;
                }else if(abs(Y[0] - Y[1]) > 10.0 or abs(Y[0] - Y[2]) > 10.0 or abs(Y[1] -Y[2]) > 10.0){
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

                for(i=0;i<3;i++){
                        for (j=0;j<2;j++){
                                LAMBDA[0][i][j] = VEL[i][j] - C_SOUND[i];
                                LAMBDA[1][i][j] = VEL[i][j];
                                LAMBDA[2][i][j] = VEL[i][j];
                                LAMBDA[3][i][j] = VEL[i][j] + C_SOUND[i];
                        }

#ifdef DEBUG
                        cout << "lambda " << i << " =\t" << LAMBDA[0][i][0] << "\t" << LAMBDA[0][i][1] << endl;
                        cout << "lambda " << i << " =\t" << LAMBDA[1][i][0] << "\t" << LAMBDA[1][i][1] << endl;
                        cout << "lambda " << i << " =\t" << LAMBDA[2][i][0] << "\t" << LAMBDA[2][i][1] << endl;
                        cout << "lambda " << i << " =\t" << LAMBDA[3][i][0] << "\t" << LAMBDA[3][i][1] << endl;
#endif
                }

                // Calculate inflow parameters

                for(i=0;i<4;i++){
                        for(j=0;j<3;j++){
                                INFLOW[i][j] = 0.0 ;
                                for(k=0;k<2;k++){INFLOW[i][j] = INFLOW[i][j] + 0.5*LAMBDA[i][j][k]*NORMAL[j][k];}
                                INFLOW_PLUS[i][j]  = max_val(0,INFLOW[i][j]);
                                INFLOW_MINUS[i][j] = min_val(0,INFLOW[i][j]);
                        }
                }
                
                // Calculate inflow and outflow state of the element

                for(i=0;i<4;i++){
                        IN_TOP     = 0.0;
                        IN_BOTTOM  = 0.0;
                        OUT_TOP    = 0.0;
                        OUT_BOTTOM = 0.0;
                        for(j=0;j<3;j++){
                                IN_TOP     = IN_TOP     + INFLOW_MINUS[i][j] * U_N[i][j];
                                IN_BOTTOM  = IN_BOTTOM  + INFLOW_MINUS[i][j];
                                OUT_TOP    = OUT_TOP    + INFLOW_PLUS[i][j]  * U_N[i][j];
                                OUT_BOTTOM = OUT_BOTTOM + INFLOW_PLUS[i][j];
                        }
#ifdef DEBUG
                        cout << "Sums (" << i << ") =\t" << IN_TOP << "\t" << IN_BOTTOM << "\t" << OUT_TOP << "\t" << OUT_BOTTOM << endl;
#endif
                        if(OUT_BOTTOM != 0.0 and IN_BOTTOM != 0.0){
                                U_IN[i]  = IN_TOP  / IN_BOTTOM;
                                U_OUT[i] = OUT_TOP / OUT_BOTTOM;
                                for(j=0;j<3;j++){BETA[i][j] = INFLOW_PLUS[i][j] / OUT_BOTTOM;}
                        }else{
                                U_IN[i]  = 0.0;
                                U_OUT[i] = 0.0;
                                for(j=0;j<3;j++){BETA[i][j] = 0.0;}
                        }
                }

                // Calculate spatial splitting for first half timestep

                for(i=0;i<4;i++){
                        for(j=0;j<3;j++){
                                FLUC[i][j] = INFLOW_PLUS[i][j]*(U_OUT[i]-U_IN[i]);
                        }
                }

                // Calculate change to be distributed

                for(i=0;i<4;i++){
                        DU0[i] = DT*FLUC[i][0]/DUAL[0];//-1.0*BETA[i][0]*FLUC[i][0];
                        DU1[i] = DT*FLUC[i][1]/DUAL[1];//-1.0*BETA[i][1]*FLUC[i][1];
                        DU2[i] = DT*FLUC[i][2]/DUAL[2];//-1.0*BETA[i][2]*FLUC[i][2];
                }

                VERTEX_0->update_du_half(DU0);
                VERTEX_1->update_du_half(DU1);
                VERTEX_2->update_du_half(DU2);

                //VERTEX_0->update_du(DU0);
                //VERTEX_1->update_du(DU1);
                //VERTEX_2->update_du(DU2);

#ifdef DEBUG
                //if(U_IN[0] != U_OUT[0]){
                        for(i=0;i<4;i++){cout << "u_in =\t" << U_IN[i] << "\tu_out =\t" << U_OUT[i] << endl;}
                        for(i=0;i<4;i++){cout << "Element fluctuation =\t" << FLUC[i][0] << "\t" << FLUC[i][1] << "\t" << FLUC[i][2] << endl;}
                        for(i=0;i<4;i++){cout << "Beta (" << i << ") =\t" << BETA[i][0] << "\t" << BETA[i][1] << "\t" << BETA[i][2] << "\tTotal =\t" << BETA[i][0]+BETA[i][1]+BETA[i][2] << endl;}
                        cout << "Dual =\t" << VERTEX_0->get_dual() << "\t" << VERTEX_1->get_dual() << "\t" << VERTEX_2->get_dual() << endl;
                        cout << "Change (rho) =\t"    << DU0[0] << "\t" << DU1[0] << "\t" << DU2[0] << endl;
                        cout << "Change (x mom) =\t"  << DU0[1] << "\t" << DU1[1] << "\t" << DU2[1] << endl;
                        cout << "Change (y mom) =\t"  << DU0[2] << "\t" << DU1[2] << "\t" << DU2[2] << endl;
                        cout << "Change (energy) =\t" << DU0[3] << "\t" << DU1[3] << "\t" << DU2[3] << endl;
                        cout << "-----------------------------------------------------------------" << endl;
                        if(isnan(DU0[3])){
                                cout << "Exiting on NaN (first half)" << endl;
                                exit(0);
                        }
                        //if(T>0.0){exit(0);}
                        //exit(0);
                //}
#endif

                return ;
        }

        //**********************************************************************************************************************

        void calculate_second_half(double T, double DT_TOT){
                int i,j,k;
                double DU0[4],DU1[4],DU2[4];
                double INFLOW_HALF[4][3],INFLOW_PLUS_HALF[4][3],INFLOW_MINUS_HALF[4][3];
                double IN_TOP,IN_BOTTOM,OUT_TOP,OUT_BOTTOM;
                double MASS[4][3][3],MASS_GAL[4][3][3];
                double SUM_MASS[4][3];
                double DT = DT_TOT;

                setup_half_state();


                if(abs(X[0] - X[1]) > 10.0 or abs(X[0] - X[2]) > 10.0 or abs(X[1] -X[2]) > 10.0){
#ifdef DEBUG
                        cout << "Skipping x boundary" << endl;
#endif
                        return ;
                }else if(abs(Y[0] - Y[1]) > 10.0 or abs(Y[0] - Y[2]) > 10.0 or abs(Y[1] -Y[2]) > 10.0){
#ifdef DEBUG
                        cout << "Skipping y boundary" << endl;
#endif
                        return ;
                }

#ifdef DEBUG
                cout << "-- SECOND -------------------------------------------------------" << endl;
                cout << "Time     =\t" << T << endl;
                cout << "i =\t" << X[0] << "\t" << Y[0] << endl;
                cout << "j =\t" << X[1] << "\t" << Y[1] << endl;
                cout << "k =\t" << X[2] << "\t" << Y[2] << endl;
                cout << "Half State    =" << "\trho" << "\tx_mom" << "\ty_mom" << "\tenergy" << endl;
                for(i=0;i<3;i++){cout << i << " =\t" << U_HALF[0][i] << "\t" << U_HALF[1][i] << "\t" << U_HALF[2][i] << "\t" << U_HALF[3][i] << endl;}
                cout << "Pressure =\t" << PRESSURE_HALF[0] << "\t" << PRESSURE_HALF[1] << "\t" << PRESSURE_HALF[2] << endl;
#endif

                // calcualte sound speed for half time state

                for(i=0;i<3;i++){C_SOUND_HALF[i] = sqrt(GAMMA*PRESSURE_HALF[i]/U_HALF[0][i]);}

                for(i=0;i<3;i++){
                        for (j=0;j<2;j++){
                                LAMBDA_HALF[0][i][j] = VEL_HALF[i][j] - C_SOUND_HALF[i];
                                LAMBDA_HALF[1][i][j] = VEL_HALF[i][j];
                                LAMBDA_HALF[2][i][j] = VEL_HALF[i][j];
                                LAMBDA_HALF[3][i][j] = VEL_HALF[i][j] + C_SOUND_HALF[i];
                        }

#ifdef DEBUG
                        cout << "lambda " << i << " =\t" << LAMBDA_HALF[0][i][0] << "\t" << LAMBDA_HALF[0][i][1] << endl;
                        cout << "lambda " << i << " =\t" << LAMBDA_HALF[1][i][0] << "\t" << LAMBDA_HALF[1][i][1] << endl;
                        cout << "lambda " << i << " =\t" << LAMBDA_HALF[2][i][0] << "\t" << LAMBDA_HALF[2][i][1] << endl;
                        cout << "lambda " << i << " =\t" << LAMBDA_HALF[3][i][0] << "\t" << LAMBDA_HALF[3][i][1] << endl;
#endif
                }

                // Calculate inflow parameters for half timestep state

                for(i=0;i<4;i++){
                        for(j=0;j<3;j++){
                                INFLOW_HALF[i][j] = 0.0 ;
                                for(k=0;k<2;k++){INFLOW_HALF[i][j] = INFLOW_HALF[i][j] + 0.5*LAMBDA_HALF[i][j][k]*NORMAL[j][k];}
                                INFLOW_PLUS_HALF[i][j]  = max_val(0,INFLOW_HALF[i][j]);
                                INFLOW_MINUS_HALF[i][j] = min_val(0,INFLOW_HALF[i][j]);
                        }
                }

                // Calculate inflow and outflow state of element at half timestep

                for(i=0;i<4;i++){
                        IN_TOP     = 0.0;
                        IN_BOTTOM  = 0.0;
                        OUT_TOP    = 0.0;
                        OUT_BOTTOM = 0.0;
                        for(j=0;j<3;j++){
                                IN_TOP     = IN_TOP     + INFLOW_MINUS_HALF[i][j] * U_HALF[i][j];
                                IN_BOTTOM  = IN_BOTTOM  + INFLOW_MINUS_HALF[i][j];
                                OUT_TOP    = OUT_TOP    + INFLOW_PLUS_HALF[i][j]  * U_HALF[i][j];
                                OUT_BOTTOM = OUT_BOTTOM + INFLOW_PLUS_HALF[i][j];
                        }
#ifdef DEBUG
                        cout << "Sums (" << i << ") =\t" << IN_TOP << "\t" << IN_BOTTOM << "\t" << OUT_TOP << "\t" << OUT_BOTTOM << endl;
#endif
                        if(OUT_BOTTOM != 0.0 and IN_BOTTOM != 0.0){
                                U_IN_HALF[i]  = IN_TOP  / IN_BOTTOM;
                                U_OUT_HALF[i] = OUT_TOP / OUT_BOTTOM;
                                for(j=0;j<3;j++){BETA_HALF[i][j] = INFLOW_PLUS_HALF[i][j] / OUT_BOTTOM;}
                        }else{
                                U_IN_HALF[i]  = 0.0;
                                U_OUT_HALF[i] = 0.0;
                                for(j=0;j<3;j++){BETA_HALF[i][j] = 0.0;}
                        }
                }

                for(i=0;i<4;i++){
                        for(j=0;j<3;j++){
                                FLUC_HALF[i][j] = INFLOW_PLUS_HALF[i][j]*(U_OUT_HALF[i]-U_IN_HALF[i]);
                        }
                }

                for(i=0;i<4;i++){
                        for(j=0;j<3;j++){
                                SUM_MASS[i][j] = 0.0;
                                for(k=0;k<3;k++){
                                        MASS[i][j][k] = AREA*BETA[i][j]/3.0;
                                        if(j == k){
                                                MASS_GAL[i][j][k] = AREA/6.0;
                                        }else{
                                                MASS_GAL[i][j][k] = AREA/12.0;
                                        }
                                        SUM_MASS[i][j] += (MASS[i][j][k]-MASS_GAL[i][j][k])*(U_HALF[i][j]-U_N[i][j])/DT;
                                }
#ifdef DEBUG
                                cout << "Mass Matrix Sum: " << i <<  " =\t" << SUM_MASS[i][0] << "\t" << SUM_MASS[i][1] << "\t" << SUM_MASS[i][2] << endl;
#endif
                        }
                }

                for(i=0;i<4;i++){
                        DU0[i] = (DT/DUAL[0])*(SUM_MASS[i][0]+0.5*(FLUC[i][0]+FLUC_HALF[i][0]));
                        DU1[i] = (DT/DUAL[1])*(SUM_MASS[i][1]+0.5*(FLUC[i][1]+FLUC_HALF[i][1]));
                        DU2[i] = (DT/DUAL[2])*(SUM_MASS[i][2]+0.5*(FLUC[i][2]+FLUC_HALF[i][2]));
                }

                VERTEX_0->update_du(DU0);
                VERTEX_1->update_du(DU1);
                VERTEX_2->update_du(DU2);

#ifdef DEBUG
                //if(U_IN[0] != U_OUT[0]){
                        for(i=0;i<4;i++){cout << "u_in =\t" << U_IN_HALF[i] << "\tu_out =\t" << U_OUT_HALF[i] << endl;}
                        for(i=0;i<4;i++){cout << "Element fluctuation =\t" << FLUC_HALF[i][0] << "\t" << FLUC_HALF[i][1] << "\t" << FLUC_HALF[i][2] << endl;}
                        for(i=0;i<4;i++){cout << "Beta (" << i << ") =\t" << BETA_HALF[i][0] << "\t" << BETA_HALF[i][1] << "\t" << BETA_HALF[i][2] << "\tTotal =\t" << BETA_HALF[i][0]+BETA_HALF[i][1]+BETA_HALF[i][2] << endl;}
                        cout << "Dual =\t" << VERTEX_0->get_dual() << "\t" << VERTEX_1->get_dual() << "\t" << VERTEX_2->get_dual() << endl;
                        cout << "Change (rho) =\t"    << DU0[0] << "\t" << DU1[0] << "\t" << DU2[0] << endl;
                        cout << "Change (x mom) =\t"  << DU0[1] << "\t" << DU1[1] << "\t" << DU2[1] << endl;
                        cout << "Change (y mom) =\t"  << DU0[2] << "\t" << DU1[2] << "\t" << DU2[2] << endl;
                        cout << "Change (energy) =\t" << DU0[3] << "\t" << DU1[3] << "\t" << DU2[3] << endl;
                        cout << "-----------------------------------------------------------------" << endl;
                        if(isnan(DU0[3])){
                                cout << "Exiting on NaN (second half)" << endl;
                                exit(0);
                        }
                        //if(T>0.0){exit(0);}
                        //exit(0);
                //}
#endif

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