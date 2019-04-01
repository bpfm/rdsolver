#include <stdio.h>
#include "cblas.h"
#include "lapacke.h"

lapack_int matSchur(double *A, double *WR, double *WI, double *VS, unsigned n){
        lapack_int ret;

        char jobvs = 'V';

        LAPACK_D_SELECT2 select;

        ret =  LAPACKE_dgees(LAPACK_COL_MAJOR,
                                'V',
                                'N',
                                select,
                                n,
                                A,
                                n,
                                0,
                                WR,
                                WI,
                                VS,
                                n);

        return ret;
}

int main(){

        double A[] = {
            0.378589,   0.971711,   0.016087,   0.037668,   0.312398,
            0.756377,   0.345708,   0.922947,   0.846671,   0.856103,
            0.732510,   0.108942,   0.476969,   0.398254,   0.507045,
            0.162608,   0.227770,   0.533074,   0.807075,   0.180335,
            0.517006,   0.315992,   0.914848,   0.460825,   0.731980
        };

        for (int i=0; i<25; i++) {
                if ((i%5) == 0){putchar('\n');}
                printf("%+12.8f ",A[i]);
        }
        putchar('\n');

        //matSchur(A,5);

        for (int i=0; i<25; i++) {
                if ((i%5) == 0){putchar('\n');}
                printf("%+12.8f ",A[i]);
        }
        putchar('\n');

        return 0;
}