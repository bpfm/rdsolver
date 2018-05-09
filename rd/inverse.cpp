using namespace std;

lapack_int matInv(double *A, unsigned n){
        int ipiv[n+1];
        lapack_int ret;

        ret =  LAPACKE_dgetrf(LAPACK_COL_MAJOR,n,n,A,n,ipiv);

#ifdef DEBUG
        cout << "ret =\t" << ret << "\t(0 = done, <0 = illegal arguement, >0 = singular)" << endl;
#endif

        if(ret !=0){return ret;}

        ret = LAPACKE_dgetri(LAPACK_COL_MAJOR,n,A,n,ipiv);

        return ret;
}

lapack_int matFac(double *A, unsigned n){
        int ipiv[n+1];
        lapack_int ret;

        ret =  LAPACKE_dgetrf(LAPACK_COL_MAJOR,n,n,A,n,ipiv);

#ifdef DEBUG
        cout << "ret =\t" << ret << "\t(0 = done, <0 = illegal arguement, >0 = singular)" << endl;
#endif

        // for(int i=0;i<16;++i){
        //         cout << A[i] << "\t";
        //         if((i+1)%4 == 0){cout << endl;}
        // }

        return ret;
}