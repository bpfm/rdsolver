#ifdef SELF_GRAVITY
void direct_gravity(std::vector<VERTEX> POINTS, int N_POINTS, double DT){
        double X0, X1, Y0, Y1, X_DIFF, Y_DIFF, AX, AY, MASS;
        for(int i=0; i<N_POINTS; ++i){
                X0 = POINTS[i].get_x();
                Y0 = POINTS[i].get_y();
                MASS = POINTS[i].get_mass();
                for(int j=0; j<N_POINTS; ++j){
                        if(j!=i){
                                X1 = POINTS[j].get_x();
                                Y1 = POINTS[j].get_y();
                                X_DIFF = X1 - X0;
                                Y_DIFF = Y1 - Y0;
                                AX = AX - GRAV * MASS * X_DIFF/pow((sqrt(X_DIFF*X_DIFF + Y_DIFF*Y_DIFF)),3);
                                AY = AY - GRAV * MASS * Y_DIFF/pow((sqrt(X_DIFF*X_DIFF + Y_DIFF*Y_DIFF)),3);
                                POINTS[i].accelerate(0,AX,AY,DT);
                        }
                }
        }
}
#endif