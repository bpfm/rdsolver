#ifdef ANALYTIC_GRAVITY
void plummer_gravity(std::vector<VERTEX> &MY_POINTS, double DT, int N_POINTS){
        // Fixed Plummer potential at (XC,YC)
        int i;
        double X, Y, GM, AX, AY, DU[4];
        double XC = 5.0, YC = 5.0, MPERT = 3.28E+05, EPS = 0.145;
        double MASS_DENSITY, DELTAX, DELTAY, RAD2;

        VERTEX MY_VERTEX;

        DU[0] = DU[3] = 0.0;

        for(i=0;i<N_POINTS;++i){
                MY_VERTEX = MY_POINTS[i];
                X = MY_VERTEX.get_x();
                Y = MY_VERTEX.get_y();
                MASS_DENSITY = MY_VERTEX.get_mass_density();
                DELTAX = X - XC, DELTAY = Y - YC;
                RAD2 = DELTAX*DELTAX + DELTAY*DELTAY;
                GM = GRAV * MPERT / (sqrt(RAD2 + EPS*EPS)*sqrt(RAD2 + EPS*EPS)*sqrt(RAD2 + EPS*EPS));
                AX = DELTAX * GM;
                AY = DELTAY * GM;
                DU[1] = -1.0*AX * DT * MASS_DENSITY;
                DU[2] = -1.0*AY * DT * MASS_DENSITY;
                MY_POINTS[i].update_du(DU);
        }
        return ;
}
#endif

void sources(std::vector<VERTEX> &MY_POINTS, double DT, int N_POINTS){
#ifdef ANALYTIC_GRAVITY
        plummer_gravity(MY_POINTS, DT, N_POINTS);
#endif
}
