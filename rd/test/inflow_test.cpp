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
        double K[4][4][3],K_PLUS[4][4][3],K_MINUS[4][4][3];
        double K_PLUS_INVERSE[4][4],K_MINUS_INVERSE[4][4];
        double A[4][4],B[4][4];
        double NORMAL[3][2];
        double PRESSURE;
        double VEL,C_SOUND;
        double GAMMA = 1.4;

        U_N[0][0] = 0.125;
        U_N[1][0] = 0.0;
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
        cout << "k" << endl;
        cout << endl;

        for(i=0;i<4;++i){
                for(j=0;j<4;++j){
                        for(m=0;m<3;++m){
                                K[i][j][m] = 0.5*(A[i][j] * NORMAL[m][0] + B[i][j] * NORMAL[m][1]);
                        }
                }
        }

        for(m=0;m<3;++m){
                for(i=0;i<4;++i){
                        cout << K[i][0][m] << "\t" << K[i][1][m] << "\t" << K[i][2][m] << "\t" << K[i][3][m] << endl;
                }
                cout << endl;
        }

        double K_FLAT0[16],K_FLAT1[16],K_FLAT2[16];

        K_FLAT0[0]  = K[0][0][0];
        K_FLAT0[1]  = K[0][1][0];
        K_FLAT0[2]  = K[0][2][0];
        K_FLAT0[3]  = K[0][3][0];

        K_FLAT0[4]  = K[1][0][0];
        K_FLAT0[5]  = K[1][1][0];
        K_FLAT0[6]  = K[1][2][0];
        K_FLAT0[7]  = K[1][3][0];

        K_FLAT0[8]  = K[2][0][0];
        K_FLAT0[9]  = K[2][1][0];
        K_FLAT0[10] = K[2][2][0];
        K_FLAT0[11] = K[2][3][0];

        K_FLAT0[12] = K[3][0][0];
        K_FLAT0[13] = K[3][1][0];
        K_FLAT0[14] = K[3][2][0];
        K_FLAT0[15] = K[3][3][0];


        K_FLAT1[0]  = K[0][0][0];
        K_FLAT1[1]  = K[0][1][0];
        K_FLAT1[2]  = K[0][2][0];
        K_FLAT1[3]  = K[0][3][0];

        K_FLAT1[4]  = K[1][0][0];
        K_FLAT1[5]  = K[1][1][0];
        K_FLAT1[6]  = K[1][2][0];
        K_FLAT1[7]  = K[1][3][0];

        K_FLAT1[8]  = K[2][0][0];
        K_FLAT1[9]  = K[2][1][0];
        K_FLAT1[10] = K[2][2][0];
        K_FLAT1[11] = K[2][3][0];

        K_FLAT1[12] = K[3][0][0];
        K_FLAT1[13] = K[3][1][0];
        K_FLAT1[14] = K[3][2][0];
        K_FLAT1[15] = K[3][3][0];


        K_FLAT2[0]  = K[0][0][0];
        K_FLAT2[1]  = K[0][1][0];
        K_FLAT2[2]  = K[0][2][0];
        K_FLAT2[3]  = K[0][3][0];

        K_FLAT2[4]  = K[1][0][0];
        K_FLAT2[5]  = K[1][1][0];
        K_FLAT2[6]  = K[1][2][0];
        K_FLAT2[7]  = K[1][3][0];

        K_FLAT2[8]  = K[2][0][0];
        K_FLAT2[9]  = K[2][1][0];
        K_FLAT2[10] = K[2][2][0];
        K_FLAT2[11] = K[2][3][0];

        K_FLAT2[12] = K[3][0][0];
        K_FLAT2[13] = K[3][1][0];
        K_FLAT2[14] = K[3][2][0];
        K_FLAT2[15] = K[3][3][0];

        matFac(K_FLAT0,4);
        matFac(K_FLAT1,4);
        matFac(K_FLAT2,4);

        for(i=0;i<4;++i){
                LAMBDA[i][0] = K_FLAT0[5*i];
                LAMBDA[i][1] = K_FLAT0[5*i];
                LAMBDA[i][2] = K_FLAT0[5*i];
        }

        for(i=0;i<4;++i){
                for(j=0;j<4;++j){
                        for(m=0;m<3;++m){
                                if(LAMBDA[i][m] >= 0.0){
                                        K_PLUS[i][j][m] = K[i][j][m];
                                        K_MINUS[i][j][m] = 0.0;
                                        K_PLUS_INVERSE[i][j] = K[i][j][m];
                                        K_MINUS_INVERSE[i][j] = 0.0;
                                }else{
                                        K_PLUS[i][j][m] = 0.0;
                                        K_MINUS[i][j][m] = K[i][j][m];
                                        K_PLUS_INVERSE[i][j] = 0.0;
                                        K_MINUS_INVERSE[i][j] = K[i][j][m];
                                }
                        }
                }
        }

        cout << endl;
        cout << "k plus" << endl;
        cout << endl;

        for(m=0;m<3;++m){
                for(i=0;i<4;++i){
                        cout << K_PLUS[i][0][m] << "\t" << K_PLUS[i][1][m] << "\t" << K_PLUS[i][2][m] << "\t" << K_PLUS[i][3][m] << endl;
                }
                cout << endl;
        }

        cout << "k minus" << endl;
        cout << endl;

        for(m=0;m<3;++m){
                for(i=0;i<4;++i){
                        cout << K_MINUS[i][0][m] << "\t" << K_MINUS[i][1][m] << "\t" << K_MINUS[i][2][m] << "\t" << K_MINUS[i][3][m] << endl;
                }
                cout << endl;
        }

        double IN_TOP[4];
        double K_PLUS_SUM[4][4];

        K_PLUS_SUM[0][0] = K_PLUS_SUM[0][1] = K_PLUS_SUM[0][2] = K_PLUS_SUM[0][3] = 0.0;
        K_PLUS_SUM[1][0] = K_PLUS_SUM[1][1] = K_PLUS_SUM[1][2] = K_PLUS_SUM[1][3] = 0.0;
        K_PLUS_SUM[2][0] = K_PLUS_SUM[2][1] = K_PLUS_SUM[2][2] = K_PLUS_SUM[2][3] = 0.0;
        K_PLUS_SUM[3][0] = K_PLUS_SUM[3][1] = K_PLUS_SUM[3][2] = K_PLUS_SUM[3][3] = 0.0;

        for(i=0;i<4;++i){
                IN_TOP[i] = 0.0;
                for(m=0;m<3;++m){
                        IN_TOP[i] += K_PLUS[i][0][m] * U_N[0][m] + K_PLUS[i][1][m] * U_N[1][m] + K_PLUS[i][2][m] * U_N[2][m] + K_PLUS[i][3][m] * U_N[3][m];
                        for(j=0;j<4;++j){
                                K_PLUS_SUM[i][j] += K_PLUS[i][j][m];
                        }
                }
        }

        cout << "k+ Sum" << endl;
        cout << endl;

        for(i=0;i<4;++i){
                cout << K_PLUS_SUM[i][0] << "\t" << K_PLUS_SUM[i][1] << "\t" << K_PLUS_SUM[i][2] << "\t" << K_PLUS_SUM[i][3] << endl;
        }

        double OUT_TOP[4];
        double K_MINUS_SUM[4][4];

        K_MINUS_SUM[0][0] = K_MINUS_SUM[0][1] = K_MINUS_SUM[0][2] = K_MINUS_SUM[0][3] = 0.0;
        K_MINUS_SUM[1][0] = K_MINUS_SUM[1][1] = K_MINUS_SUM[1][2] = K_MINUS_SUM[1][3] = 0.0;
        K_MINUS_SUM[2][0] = K_MINUS_SUM[2][1] = K_MINUS_SUM[2][2] = K_MINUS_SUM[2][3] = 0.0;
        K_MINUS_SUM[3][0] = K_MINUS_SUM[3][1] = K_MINUS_SUM[3][2] = K_MINUS_SUM[3][3] = 0.0;

        for(i=0;i<4;++i){
                OUT_TOP[i] = 0.0;
                for(m=0;m<3;++m){
                        OUT_TOP[i] += K_MINUS[i][0][m] * U_N[0][m] + K_MINUS[i][1][m] * U_N[1][m] + K_MINUS[i][2][m] * U_N[2][m] + K_MINUS[i][3][m] * U_N[3][m];
                        for(j=0;j<4;++j){
                                K_MINUS_SUM[i][j] += K_MINUS[i][j][m];
                        }
                }
        }

        cout << endl;
        cout << "k- Sum" << endl;
        cout << endl;

        for(i=0;i<4;++i){
                cout << K_MINUS_SUM[i][0] << "\t" << K_MINUS_SUM[i][1] << "\t" << K_MINUS_SUM[i][2] << "\t" << K_MINUS_SUM[i][3] << endl;
        }

        double K_PLUS_SUM_FLAT[16],K_MINUS_SUM_FLAT[16];

        K_PLUS_SUM_FLAT[0]  = K_PLUS_SUM[0][0];
        K_PLUS_SUM_FLAT[1]  = K_PLUS_SUM[0][1];
        K_PLUS_SUM_FLAT[2]  = K_PLUS_SUM[0][2];
        K_PLUS_SUM_FLAT[3]  = K_PLUS_SUM[0][3];

        K_PLUS_SUM_FLAT[4]  = K_PLUS_SUM[1][0];
        K_PLUS_SUM_FLAT[5]  = K_PLUS_SUM[1][1];
        K_PLUS_SUM_FLAT[6]  = K_PLUS_SUM[1][2];
        K_PLUS_SUM_FLAT[7]  = K_PLUS_SUM[1][3];

        K_PLUS_SUM_FLAT[8]  = K_PLUS_SUM[2][0];
        K_PLUS_SUM_FLAT[9]  = K_PLUS_SUM[2][1];
        K_PLUS_SUM_FLAT[10] = K_PLUS_SUM[2][2];
        K_PLUS_SUM_FLAT[11] = K_PLUS_SUM[2][3];

        K_PLUS_SUM_FLAT[12] = K_PLUS_SUM[3][0];
        K_PLUS_SUM_FLAT[13] = K_PLUS_SUM[3][1];
        K_PLUS_SUM_FLAT[14] = K_PLUS_SUM[3][2];
        K_PLUS_SUM_FLAT[15] = K_PLUS_SUM[3][3];

        K_MINUS_SUM_FLAT[0]  = K_MINUS_SUM[0][0];
        K_MINUS_SUM_FLAT[1]  = K_MINUS_SUM[0][1];
        K_MINUS_SUM_FLAT[2]  = K_MINUS_SUM[0][2];
        K_MINUS_SUM_FLAT[3]  = K_MINUS_SUM[0][3];

        K_MINUS_SUM_FLAT[4]  = K_MINUS_SUM[1][0];
        K_MINUS_SUM_FLAT[5]  = K_MINUS_SUM[1][1];
        K_MINUS_SUM_FLAT[6]  = K_MINUS_SUM[1][2];
        K_MINUS_SUM_FLAT[7]  = K_MINUS_SUM[1][3];

        K_MINUS_SUM_FLAT[8]  = K_MINUS_SUM[2][0];
        K_MINUS_SUM_FLAT[9]  = K_MINUS_SUM[2][1];
        K_MINUS_SUM_FLAT[10] = K_MINUS_SUM[2][2];
        K_MINUS_SUM_FLAT[11] = K_MINUS_SUM[2][3];

        K_MINUS_SUM_FLAT[12] = K_MINUS_SUM[3][0];
        K_MINUS_SUM_FLAT[13] = K_MINUS_SUM[3][1];
        K_MINUS_SUM_FLAT[14] = K_MINUS_SUM[3][2];
        K_MINUS_SUM_FLAT[15] = K_MINUS_SUM[3][3];

        matInv(K_PLUS_SUM_FLAT,4);
        matInv(K_MINUS_SUM_FLAT,4);

        exit(0);

        K_PLUS_INVERSE[0][0] = 0.0;
        K_PLUS_INVERSE[0][1] = 0.0;
        K_PLUS_INVERSE[0][2] = 0.0;
        K_PLUS_INVERSE[0][3] = 0.0;

        K_PLUS_INVERSE[1][0] = -1.73208;
        K_PLUS_INVERSE[1][1] = 0.0;
        K_PLUS_INVERSE[1][2] = 0.0;
        K_PLUS_INVERSE[1][3] = 0.0;

        K_PLUS_INVERSE[2][0] = -1.00004;
        K_PLUS_INVERSE[2][1] = 0.0;
        K_PLUS_INVERSE[2][2] = 0.0;
        K_PLUS_INVERSE[2][3] = 0.0;

        K_PLUS_INVERSE[3][0] = 0.0;
        K_PLUS_INVERSE[3][1] = 4.33019;
        K_PLUS_INVERSE[3][2] = 2.50011;
        K_PLUS_INVERSE[3][3] = 0.0;

        K_MINUS_INVERSE[0][0] = 0.0;
        K_MINUS_INVERSE[0][1] = 0.0;
        K_MINUS_INVERSE[0][2] = 0.0;
        K_MINUS_INVERSE[0][3] = 0.0;

        K_MINUS_INVERSE[1][0] = 1.73208;
        K_MINUS_INVERSE[1][1] = 0.0;
        K_MINUS_INVERSE[1][2] = 0.0;
        K_MINUS_INVERSE[1][3] = 0.0;

        K_MINUS_INVERSE[2][0] = 1.00004;
        K_MINUS_INVERSE[2][1] = 0.0;
        K_MINUS_INVERSE[2][2] = 0.0;
        K_MINUS_INVERSE[2][3] = 0.0;

        K_MINUS_INVERSE[3][0] = 0.0;
        K_MINUS_INVERSE[3][1] = -4.33019;
        K_MINUS_INVERSE[3][2] = -2.50011;
        K_MINUS_INVERSE[3][3] = 0.0;

        double U_OUT[4],U_IN[4],PHI[4][3];

        for(i=0;i<4;++i){U_OUT[i] += K_PLUS_INVERSE[i][0] * OUT_TOP[0] + K_PLUS_INVERSE[i][1] * OUT_TOP[1] + K_PLUS_INVERSE[i][2] * OUT_TOP[2] + K_PLUS_INVERSE[i][3] * OUT_TOP[3];}
        for(i=0;i<4;++i){U_IN[i] += K_MINUS_INVERSE[i][0] * IN_TOP[0] + K_MINUS_INVERSE[i][1] * IN_TOP[1] + K_MINUS_INVERSE[i][2] * IN_TOP[2] + K_MINUS_INVERSE[i][3] * IN_TOP[3];}

        for(i=0;i<4;++i){
                cout << endl;
                for(m=0;m<3;++m){
                        PHI[i][m] = K_PLUS[i][0][m] * (U_OUT[0] - U_IN[0]) + K_PLUS[i][1][m] * (U_OUT[1] - U_IN[1]) + K_PLUS[i][2][m] * (U_OUT[2] - U_IN[2]) + K_PLUS[i][3][m] * (U_OUT[3] - U_IN[3]);
                        cout << "PHI " << i << "\t" << PHI[i][m] << endl;;
                }
        }

        return 0;

}