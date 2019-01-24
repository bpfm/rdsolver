#include <stdio.h>
#include <iostream>
#include <cmath>

#include "cblas.h"
#include "lapacke.h"

using namespace std;

lapack_int matInv(double *A, unsigned n){
        int ipiv[n+1];
        lapack_int ret;

        ret =  LAPACKE_dgetrf(LAPACK_COL_MAJOR,n,n,A,n,ipiv);

        cout << "ret =\t" << ret << "\t(0 = done, <0 = illegal arguement, >0 = singular)" << endl;

        if(ret !=0){return ret;}

        // for(int i=0;i<16;++i){
        //         cout << A[i] << "\t";
        //         if((i+1)%4 == 0){cout << endl;}
        // }

        ret = LAPACKE_dgetri(LAPACK_COL_MAJOR,n,A,n,ipiv);

        return ret;
}

lapack_int matFac(double *A, unsigned n){
        int ipiv[n+1];
        lapack_int ret;

        ret =  LAPACKE_dgetrf(LAPACK_COL_MAJOR,n,n,A,n,ipiv);

        cout << "ret =\t" << ret << "\t(0 = done, <0 = illegal arguement, >0 = singular)" << endl;

        cout << endl;

        for(int i=0;i<16;++i){
                cout << A[i] << "\t";
                if((i+1)%4 == 0){cout << endl;}
        }

        cout << endl;

        return ret;
}

int main(){
        int i,j,m;
        double U[4],U_N[4][3];
        double LAMBDA[4][3];
        double INFLOW[4][4][3],INFLOW_PLUS[4][4][3],INFLOW_MINUS[4][4][3];
        double INFLOW_PLUS_INVERSE[4][4],INFLOW_MINUS_INVERSE[4][4];
        double A[4][4],B[4][4];
        double NORMAL[3][2];
        double PRESSURE;
        double VEL,C_SOUND;
        double GAMMA = 1.4;

        U_N[0][0] = 0.125;
        U_N[1][0] = 0.1;
        U_N[2][0] = 0.0;
        U_N[3][0] = 250.0;

        U_N[0][1] = 0.125;
        U_N[1][1] = 0.0;
        U_N[2][1] = 0.0;
        U_N[3][1] = 250.0;

        U_N[0][2] = 1.0;
        U_N[1][2] = 0.0;
        U_N[2][2] = 0.0;
        U_N[3][2] = 1250.0;

        cout << "U" << endl;
        cout << endl;

        for(i=0;i<4;++i){
                U[i] = (U_N[i][0] + U_N[i][1] + U_N[i][2])/3.0;
                cout << U[i] << endl;
        }

        NORMAL[0][0] = 0.0;
        NORMAL[0][1] = 1.0;

        NORMAL[1][0] = -0.866;
        NORMAL[1][1] = -0.5;

        NORMAL[2][0] = 0.866;
        NORMAL[2][1] = -0.5;

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

        cout << endl;
        cout << "A" << endl;
        cout << endl;

        for(i=0;i<4;++i){cout << A[i][0] << "\t" << A[i][1] << "\t" << A[i][2] << "\t" << A[i][3] << endl;}

        cout << endl;

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

        cout << "B" << endl;
        cout << endl;

        for(i=0;i<4;++i){cout << B[i][0] << "\t" << B[i][1] << "\t" << B[i][2] << "\t" << B[i][3] << endl;}

        cout << endl;
        cout << "INFLOW" << endl;
        cout << endl;

        for(i=0;i<4;++i){
                for(j=0;j<4;++j){
                        for(m=0;m<3;++m){
                                INFLOW[i][j][m] = 0.5*(A[i][j] * NORMAL[m][0] + B[i][j] * NORMAL[m][1]);
                        }
                }
        }

        for(m=0;m<3;++m){
                for(i=0;i<4;++i){
                        cout << INFLOW[i][0][m] << "\t" << INFLOW[i][1][m] << "\t" << INFLOW[i][2][m] << "\t" << INFLOW[i][3][m] << endl;
                }
                cout << endl;
        }

        double INFLOW_FLAT_0[16],INFLOW_FLAT_1[16],INFLOW_FLAT_2[16];

        INFLOW_FLAT_0[0]  = INFLOW[0][0][0];
        INFLOW_FLAT_0[1]  = INFLOW[0][1][0];
        INFLOW_FLAT_0[2]  = INFLOW[0][2][0];
        INFLOW_FLAT_0[3]  = INFLOW[0][3][0];

        INFLOW_FLAT_0[4]  = INFLOW[1][0][0];
        INFLOW_FLAT_0[5]  = INFLOW[1][1][0];
        INFLOW_FLAT_0[6]  = INFLOW[1][2][0];
        INFLOW_FLAT_0[7]  = INFLOW[1][3][0];

        INFLOW_FLAT_0[8]  = INFLOW[2][0][0];
        INFLOW_FLAT_0[9]  = INFLOW[2][1][0];
        INFLOW_FLAT_0[10] = INFLOW[2][2][0];
        INFLOW_FLAT_0[11] = INFLOW[2][3][0];

        INFLOW_FLAT_0[12] = INFLOW[3][0][0];
        INFLOW_FLAT_0[13] = INFLOW[3][1][0];
        INFLOW_FLAT_0[14] = INFLOW[3][2][0];
        INFLOW_FLAT_0[15] = INFLOW[3][3][0];


        INFLOW_FLAT_1[0]  = INFLOW[0][0][0];
        INFLOW_FLAT_1[1]  = INFLOW[0][1][0];
        INFLOW_FLAT_1[2]  = INFLOW[0][2][0];
        INFLOW_FLAT_1[3]  = INFLOW[0][3][0];

        INFLOW_FLAT_1[4]  = INFLOW[1][0][0];
        INFLOW_FLAT_1[5]  = INFLOW[1][1][0];
        INFLOW_FLAT_1[6]  = INFLOW[1][2][0];
        INFLOW_FLAT_1[7]  = INFLOW[1][3][0];

        INFLOW_FLAT_1[8]  = INFLOW[2][0][0];
        INFLOW_FLAT_1[9]  = INFLOW[2][1][0];
        INFLOW_FLAT_1[10] = INFLOW[2][2][0];
        INFLOW_FLAT_1[11] = INFLOW[2][3][0];

        INFLOW_FLAT_1[12] = INFLOW[3][0][0];
        INFLOW_FLAT_1[13] = INFLOW[3][1][0];
        INFLOW_FLAT_1[14] = INFLOW[3][2][0];
        INFLOW_FLAT_1[15] = INFLOW[3][3][0];


        INFLOW_FLAT_2[0]  = INFLOW[0][0][0];
        INFLOW_FLAT_2[1]  = INFLOW[0][1][0];
        INFLOW_FLAT_2[2]  = INFLOW[0][2][0];
        INFLOW_FLAT_2[3]  = INFLOW[0][3][0];

        INFLOW_FLAT_2[4]  = INFLOW[1][0][0];
        INFLOW_FLAT_2[5]  = INFLOW[1][1][0];
        INFLOW_FLAT_2[6]  = INFLOW[1][2][0];
        INFLOW_FLAT_2[7]  = INFLOW[1][3][0];

        INFLOW_FLAT_2[8]  = INFLOW[2][0][0];
        INFLOW_FLAT_2[9]  = INFLOW[2][1][0];
        INFLOW_FLAT_2[10] = INFLOW[2][2][0];
        INFLOW_FLAT_2[11] = INFLOW[2][3][0];

        INFLOW_FLAT_2[12] = INFLOW[3][0][0];
        INFLOW_FLAT_2[13] = INFLOW[3][1][0];
        INFLOW_FLAT_2[14] = INFLOW[3][2][0];
        INFLOW_FLAT_2[15] = INFLOW[3][3][0];

        matFac(INFLOW_FLAT_0,4);
        matFac(INFLOW_FLAT_1,4);
        matFac(INFLOW_FLAT_2,4);

        for(i=0;i<4;++i){
                LAMBDA[i][0] = INFLOW_FLAT_0[5*i];
                LAMBDA[i][1] = INFLOW_FLAT_0[5*i];
                LAMBDA[i][2] = INFLOW_FLAT_0[5*i];
        }

        for(i=0;i<4;++i){
                for(j=0;j<4;++j){
                        for(m=0;m<3;++m){
                                if(LAMBDA[i][m] >= 0.0){
                                        INFLOW_PLUS[i][j][m] = INFLOW[i][j][m];
                                        INFLOW_MINUS[i][j][m] = 0.0;
                                        INFLOW_PLUS_INVERSE[i][j] = INFLOW[i][j][m];
                                        INFLOW_MINUS_INVERSE[i][j] = 0.0;
                                }else{
                                        INFLOW_PLUS[i][j][m] = 0.0;
                                        INFLOW_MINUS[i][j][m] = INFLOW[i][j][m];
                                        INFLOW_PLUS_INVERSE[i][j] = 0.0;
                                        INFLOW_MINUS_INVERSE[i][j] = INFLOW[i][j][m];
                                }
                        }
                }
        }

        cout << endl;
        cout << "INFLOW plus" << endl;
        cout << endl;

        for(m=0;m<3;++m){
                for(i=0;i<4;++i){
                        cout << INFLOW_PLUS[i][0][m] << "\t" << INFLOW_PLUS[i][1][m] << "\t" << INFLOW_PLUS[i][2][m] << "\t" << INFLOW_PLUS[i][3][m] << endl;
                }
                cout << endl;
        }

        cout << "INFLOW minus" << endl;
        cout << endl;

        for(m=0;m<3;++m){
                for(i=0;i<4;++i){
                        cout << INFLOW_MINUS[i][0][m] << "\t" << INFLOW_MINUS[i][1][m] << "\t" << INFLOW_MINUS[i][2][m] << "\t" << INFLOW_MINUS[i][3][m] << endl;
                }
                cout << endl;
        }

        double IN_TOP[4];
        double INFLOW_PLUS_SUM[4][4];

        INFLOW_PLUS_SUM[0][0] = INFLOW_PLUS_SUM[0][1] = INFLOW_PLUS_SUM[0][2] = INFLOW_PLUS_SUM[0][3] = 0.0;
        INFLOW_PLUS_SUM[1][0] = INFLOW_PLUS_SUM[1][1] = INFLOW_PLUS_SUM[1][2] = INFLOW_PLUS_SUM[1][3] = 0.0;
        INFLOW_PLUS_SUM[2][0] = INFLOW_PLUS_SUM[2][1] = INFLOW_PLUS_SUM[2][2] = INFLOW_PLUS_SUM[2][3] = 0.0;
        INFLOW_PLUS_SUM[3][0] = INFLOW_PLUS_SUM[3][1] = INFLOW_PLUS_SUM[3][2] = INFLOW_PLUS_SUM[3][3] = 0.0;

        for(i=0;i<4;++i){
                IN_TOP[i] = 0.0;
                for(m=0;m<3;++m){
                        IN_TOP[i] += INFLOW_PLUS[i][0][m] * U_N[0][m] + INFLOW_PLUS[i][1][m] * U_N[1][m] + INFLOW_PLUS[i][2][m] * U_N[2][m] + INFLOW_PLUS[i][3][m] * U_N[3][m];
                        for(j=0;j<4;++j){
                                INFLOW_PLUS_SUM[i][j] += INFLOW_PLUS[i][j][m];
                        }
                }
        }

        cout << "INFLOW+ Sum" << endl;
        cout << endl;

        for(i=0;i<4;++i){
                cout << INFLOW_PLUS_SUM[i][0] << "\t" << INFLOW_PLUS_SUM[i][1] << "\t" << INFLOW_PLUS_SUM[i][2] << "\t" << INFLOW_PLUS_SUM[i][3] << endl;
        }

        double OUT_TOP[4];
        double INFLOW_MINUS_SUM[4][4];

        INFLOW_MINUS_SUM[0][0] = INFLOW_MINUS_SUM[0][1] = INFLOW_MINUS_SUM[0][2] = INFLOW_MINUS_SUM[0][3] = 0.0;
        INFLOW_MINUS_SUM[1][0] = INFLOW_MINUS_SUM[1][1] = INFLOW_MINUS_SUM[1][2] = INFLOW_MINUS_SUM[1][3] = 0.0;
        INFLOW_MINUS_SUM[2][0] = INFLOW_MINUS_SUM[2][1] = INFLOW_MINUS_SUM[2][2] = INFLOW_MINUS_SUM[2][3] = 0.0;
        INFLOW_MINUS_SUM[3][0] = INFLOW_MINUS_SUM[3][1] = INFLOW_MINUS_SUM[3][2] = INFLOW_MINUS_SUM[3][3] = 0.0;

        for(i=0;i<4;++i){
                OUT_TOP[i] = 0.0;
                for(m=0;m<3;++m){
                        OUT_TOP[i] += INFLOW_MINUS[i][0][m] * U_N[0][m] + INFLOW_MINUS[i][1][m] * U_N[1][m] + INFLOW_MINUS[i][2][m] * U_N[2][m] + INFLOW_MINUS[i][3][m] * U_N[3][m];
                        for(j=0;j<4;++j){
                                INFLOW_MINUS_SUM[i][j] += INFLOW_MINUS[i][j][m];
                        }
                }
        }

        cout << endl;
        cout << "INFLOW- Sum" << endl;
        cout << endl;

        for(i=0;i<4;++i){
                cout << INFLOW_MINUS_SUM[i][0] << "\t" << INFLOW_MINUS_SUM[i][1] << "\t" << INFLOW_MINUS_SUM[i][2] << "\t" << INFLOW_MINUS_SUM[i][3] << endl;
        }

        double INFLOW_PLUS_SUM_FLAT_[16],INFLOW_MINUS_SUM_FLAT_[16];

        INFLOW_PLUS_SUM_FLAT_[0]  = INFLOW_PLUS_SUM[0][0];
        INFLOW_PLUS_SUM_FLAT_[1]  = INFLOW_PLUS_SUM[0][1];
        INFLOW_PLUS_SUM_FLAT_[2]  = INFLOW_PLUS_SUM[0][2];
        INFLOW_PLUS_SUM_FLAT_[3]  = INFLOW_PLUS_SUM[0][3];

        INFLOW_PLUS_SUM_FLAT_[4]  = INFLOW_PLUS_SUM[1][0];
        INFLOW_PLUS_SUM_FLAT_[5]  = INFLOW_PLUS_SUM[1][1];
        INFLOW_PLUS_SUM_FLAT_[6]  = INFLOW_PLUS_SUM[1][2];
        INFLOW_PLUS_SUM_FLAT_[7]  = INFLOW_PLUS_SUM[1][3];

        INFLOW_PLUS_SUM_FLAT_[8]  = INFLOW_PLUS_SUM[2][0];
        INFLOW_PLUS_SUM_FLAT_[9]  = INFLOW_PLUS_SUM[2][1];
        INFLOW_PLUS_SUM_FLAT_[10] = INFLOW_PLUS_SUM[2][2];
        INFLOW_PLUS_SUM_FLAT_[11] = INFLOW_PLUS_SUM[2][3];

        INFLOW_PLUS_SUM_FLAT_[12] = INFLOW_PLUS_SUM[3][0];
        INFLOW_PLUS_SUM_FLAT_[13] = INFLOW_PLUS_SUM[3][1];
        INFLOW_PLUS_SUM_FLAT_[14] = INFLOW_PLUS_SUM[3][2];
        INFLOW_PLUS_SUM_FLAT_[15] = INFLOW_PLUS_SUM[3][3];

        INFLOW_MINUS_SUM_FLAT_[0]  = INFLOW_MINUS_SUM[0][0];
        INFLOW_MINUS_SUM_FLAT_[1]  = INFLOW_MINUS_SUM[0][1];
        INFLOW_MINUS_SUM_FLAT_[2]  = INFLOW_MINUS_SUM[0][2];
        INFLOW_MINUS_SUM_FLAT_[3]  = INFLOW_MINUS_SUM[0][3];

        INFLOW_MINUS_SUM_FLAT_[4]  = INFLOW_MINUS_SUM[1][0];
        INFLOW_MINUS_SUM_FLAT_[5]  = INFLOW_MINUS_SUM[1][1];
        INFLOW_MINUS_SUM_FLAT_[6]  = INFLOW_MINUS_SUM[1][2];
        INFLOW_MINUS_SUM_FLAT_[7]  = INFLOW_MINUS_SUM[1][3];

        INFLOW_MINUS_SUM_FLAT_[8]  = INFLOW_MINUS_SUM[2][0];
        INFLOW_MINUS_SUM_FLAT_[9]  = INFLOW_MINUS_SUM[2][1];
        INFLOW_MINUS_SUM_FLAT_[10] = INFLOW_MINUS_SUM[2][2];
        INFLOW_MINUS_SUM_FLAT_[11] = INFLOW_MINUS_SUM[2][3];

        INFLOW_MINUS_SUM_FLAT_[12] = INFLOW_MINUS_SUM[3][0];
        INFLOW_MINUS_SUM_FLAT_[13] = INFLOW_MINUS_SUM[3][1];
        INFLOW_MINUS_SUM_FLAT_[14] = INFLOW_MINUS_SUM[3][2];
        INFLOW_MINUS_SUM_FLAT_[15] = INFLOW_MINUS_SUM[3][3];

        matInv(INFLOW_PLUS_SUM_FLAT_,4);
        matInv(INFLOW_MINUS_SUM_FLAT_,4);


        INFLOW_PLUS_INVERSE[0][0] = 0.0;
        INFLOW_PLUS_INVERSE[0][1] = 0.0;
        INFLOW_PLUS_INVERSE[0][2] = 0.0;
        INFLOW_PLUS_INVERSE[0][3] = 0.0;

        INFLOW_PLUS_INVERSE[1][0] = -1.73208;
        INFLOW_PLUS_INVERSE[1][1] = 0.0;
        INFLOW_PLUS_INVERSE[1][2] = 0.0;
        INFLOW_PLUS_INVERSE[1][3] = 0.0;

        INFLOW_PLUS_INVERSE[2][0] = -1.00004;
        INFLOW_PLUS_INVERSE[2][1] = 0.0;
        INFLOW_PLUS_INVERSE[2][2] = 0.0;
        INFLOW_PLUS_INVERSE[2][3] = 0.0;

        INFLOW_PLUS_INVERSE[3][0] = 0.0;
        INFLOW_PLUS_INVERSE[3][1] = 4.33019;
        INFLOW_PLUS_INVERSE[3][2] = 2.50011;
        INFLOW_PLUS_INVERSE[3][3] = 0.0;

        INFLOW_MINUS_INVERSE[0][0] = 0.0;
        INFLOW_MINUS_INVERSE[0][1] = 0.0;
        INFLOW_MINUS_INVERSE[0][2] = 0.0;
        INFLOW_MINUS_INVERSE[0][3] = 0.0;

        INFLOW_MINUS_INVERSE[1][0] = 1.73208;
        INFLOW_MINUS_INVERSE[1][1] = 0.0;
        INFLOW_MINUS_INVERSE[1][2] = 0.0;
        INFLOW_MINUS_INVERSE[1][3] = 0.0;

        INFLOW_MINUS_INVERSE[2][0] = 1.00004;
        INFLOW_MINUS_INVERSE[2][1] = 0.0;
        INFLOW_MINUS_INVERSE[2][2] = 0.0;
        INFLOW_MINUS_INVERSE[2][3] = 0.0;

        INFLOW_MINUS_INVERSE[3][0] = 0.0;
        INFLOW_MINUS_INVERSE[3][1] = -4.33019;
        INFLOW_MINUS_INVERSE[3][2] = -2.50011;
        INFLOW_MINUS_INVERSE[3][3] = 0.0;

        double U_OUT[4],U_IN[4],PHI[4][3];

        for(i=0;i<4;++i){U_OUT[i] += INFLOW_PLUS_INVERSE[i][0] * OUT_TOP[0] + INFLOW_PLUS_INVERSE[i][1] * OUT_TOP[1] + INFLOW_PLUS_INVERSE[i][2] * OUT_TOP[2] + INFLOW_PLUS_INVERSE[i][3] * OUT_TOP[3];}
        for(i=0;i<4;++i){U_IN[i] += INFLOW_MINUS_INVERSE[i][0] * IN_TOP[0] + INFLOW_MINUS_INVERSE[i][1] * IN_TOP[1] + INFLOW_MINUS_INVERSE[i][2] * IN_TOP[2] + INFLOW_MINUS_INVERSE[i][3] * IN_TOP[3];}

        for(i=0;i<4;++i){
                cout << endl;
                for(m=0;m<3;++m){
                        PHI[i][m] = INFLOW_PLUS[i][0][m] * (U_OUT[0] - U_IN[0]) + INFLOW_PLUS[i][1][m] * (U_OUT[1] - U_IN[1]) + INFLOW_PLUS[i][2][m] * (U_OUT[2] - U_IN[2]) + INFLOW_PLUS[i][3][m] * (U_OUT[3] - U_IN[3]);
                        cout << "PHI " << i << "\t" << PHI[i][m] << endl;;
                }
        }

        return 0;

}