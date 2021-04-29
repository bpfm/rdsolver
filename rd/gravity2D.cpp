#ifdef SELF_GRAVITY
void direct_gravity(std::vector<VERTEX> &MY_POINTS, int N_POINTS, double DT){
        double X0, X1, Y0, Y1, X_DIFF, Y_DIFF, AX, AY, MASS, MASS_DENSITY, DU[4];
        DU[0] = DU[3] = 0.0;
        for(int i=0; i<N_POINTS; ++i){
                X0 = MY_POINTS[i].get_x();
                Y0 = MY_POINTS[i].get_y();
                MASS = MY_POINTS[i].get_mass();
                MASS_DENSITY = MY_POINTS[i].get_mass_density();
                for(int j=0; j<N_POINTS; ++j){
                        if(j!=i){
                                X1 = MY_POINTS[j].get_x();
                                Y1 = MY_POINTS[j].get_y();
                                X_DIFF = X1 - X0;
                                Y_DIFF = Y1 - Y0;
                                AX = AX - GRAV * MASS * X_DIFF/pow((sqrt(X_DIFF*X_DIFF + Y_DIFF*Y_DIFF)),3);
                                AY = AY - GRAV * MASS * Y_DIFF/pow((sqrt(X_DIFF*X_DIFF + Y_DIFF*Y_DIFF)),3);
                                DU[1] = -1.0*AX*DT*MASS_DENSITY;
                                DU[2] = -1.0*AY*DT*MASS_DENSITY;
                                MY_POINTS[i].update_du(DU);
                        }
                }
        }
}
#endif