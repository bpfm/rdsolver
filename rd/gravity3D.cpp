#ifdef SELF_GRAVITY
void direct_gravity(std::vector<VERTEX> MY_POINTS, int N_POINTS, double DT){
        double X0, X1, Y0, Y1, Z0, Z1, X_DIFF, Y_DIFF, Z_DIFF, AX, AY, AZ, MASS, MASS_DENSITY, DU[5];
        DU[0] = DU[4] = 0.0;
        for(int i=0; i<N_POINTS; ++i){
                X0 = MY_POINTS[i].get_x();
                Y0 = MY_POINTS[i].get_y();
                Z0 = MY_POINTS[i].get_z();
                MASS = MY_POINTS[i].get_mass();
                for(int j=0; j<N_POINTS; ++j){
                        if(j!=i){
                                X1 = MY_POINTS[j].get_x();
                                Y1 = MY_POINTS[j].get_y();
                                Z1 = MY_POINTS[j].get_z();
                                X_DIFF = X1 - X0;
                                Y_DIFF = Y1 - Y0;
                                Z_DIFF = Z1 - Z0;
                                AX = AX - GRAV * MASS * X_DIFF/pow((sqrt(X_DIFF*X_DIFF + Y_DIFF*Y_DIFF + Z_DIFF*Z_DIFF)),3);
                                AY = AY - GRAV * MASS * Y_DIFF/pow((sqrt(X_DIFF*X_DIFF + Y_DIFF*Y_DIFF + Z_DIFF*Z_DIFF)),3);
                                AZ = AZ - GRAV * MASS * Z_DIFF/pow((sqrt(X_DIFF*X_DIFF + Y_DIFF*Y_DIFF + Z_DIFF*Z_DIFF)),3);
                                DU[1] = -1.0*AX*DT*MASS_DENSITY;
                                DU[2] = -1.0*AY*DT*MASS_DENSITY;
                                DU[3] = -1.0*AZ*DT*MASS_DENSITY;
                                MY_POINTS[i].update_du(DU);
                        }
                }
        }
}
#endif