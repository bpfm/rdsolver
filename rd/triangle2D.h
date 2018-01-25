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

        void calculate_change(double DX, double DT, double T){
                int i,j;
                double X[3],Y[3],U_N[4][3],U_HALF[4][3];
                double DU0[4],DU1[4],DU2[4];
                double NORMAL[3][2];
                double LAMBDA[4][3][2];
                double INFLOW[4][3],INLFOW_PLUS[4][3],INFLOW_MINUS[4][3];

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
                                INFLOW[i][j] = 0.0;
                                for(k=0;k<2){
                                        INFLOW[i][j] = INFLOW[i][j] + 0.5*LAMBDA[i][j][k]*NORMAL[j][k];
                                }
                                INLFOW_PLUS[i][j] = max_val(0,INFLOW[i][j]);
                                INFLOW_MINUS[i][j] = min_Val(0,INFLOW[i][j]);

                        }
                }

                // Work out which are inflow vertices and which are outflow vertices

                // for (int i = 0; i < 3; ++i){
                //         DU0[i] = 0.02*X[i];
                //         DU1[i] = 0.01*Y[i]*sin(Y[i]);
                //         DU2[i] = 0.0;
                // }

                VERTEX_0->update_du(DU0);
                VERTEX_1->update_du(DU1);
                VERTEX_2->update_du(DU2);

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

                NORMAL00 = NORMAL[0][0];
                NORMAL01 = NORMAL[0][1];
                NORMAL10 = NORMAL[1][0];
                NORMAL11 = NORMAL[1][1];
                NORMAL20 = NORMAL[2][0];
                NORMAL21 = NORMAL[2][1];

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