#include "constants.h"
#include "vertex2D.h"
#include <iostream>
#include <cmath>

extern int POINT_CHECK;

double F(double X){
        double FUNC = exp(-1.0/X);
        if(X<0.0){FUNC = 0.0;}
        return FUNC;
}

double G(double X){
        double FUNC = F(X)/(F(X) + F(1.0-X));
        return FUNC;
}


VERTEX setup_vertex(double X, double Y){
        VERTEX NEW_VERTEX;

        NEW_VERTEX.set_x(X);
        NEW_VERTEX.set_y(Y);

        double DX = SIDE_LENGTH_X/1000.0;               // hot garbage

        NEW_VERTEX.set_dx(SIDE_LENGTH_X/1000.0);
        NEW_VERTEX.set_dy(SIDE_LENGTH_Y/1000.0);

        NEW_VERTEX.set_dual(0.0);

#ifdef SODX
        // std::cout << "Using 1D Sod Shock Tube (Varied in X)" << std::endl;}

        if(X>0.25*SIDE_LENGTH_X and X<0.75*SIDE_LENGTH_X){
                NEW_VERTEX.set_mass_density(1.0);                               // units kg/m^3
                NEW_VERTEX.set_x_velocity(0.000000001);                       // units m/s
                NEW_VERTEX.set_y_velocity(0.000000001);                       // units m/s
                NEW_VERTEX.set_pressure(1.0);                                 // units N/m^2
        }else{
                NEW_VERTEX.set_mass_density(0.125);                             // units kg/m^3
                NEW_VERTEX.set_x_velocity(0.000000001);                       // units m/s
                NEW_VERTEX.set_y_velocity(0.000000001);                       // units m/s
                NEW_VERTEX.set_pressure(0.1);                                 // units N/m^2
        }

#endif
#ifdef SODY
        // if(i==0 and j==0){std::cout << "Using 1D Sod Shock Tube (Varied in Y)" << std::endl;}

        if(Y<0.5*SIDE_LENGTH_Y){
                NEW_VERTEX.set_mass_density(1.0);                               // units kg/m^3
                NEW_VERTEX.set_x_velocity(0.00000001);                                 // units m/s
                NEW_VERTEX.set_y_velocity(0.00000001);                                 // units m/s
                NEW_VERTEX.set_pressure(1.0);                                 // units N/m^2
        }else{
                NEW_VERTEX.set_mass_density(0.125);                             // units kg/m^3
                NEW_VERTEX.set_x_velocity(0.00000001);                                 // units m/s
                NEW_VERTEX.set_y_velocity(0.00000001);                                 // units m/s
                NEW_VERTEX.set_pressure(0.1);                                 // units N/m^2
        }

#endif
#ifdef SINEX
        // if(i==0 and j==0){std::cout << "Using 1D Sine Wave" << std::endl;}

        double RHO,RHO_0 = 50.0;
        double P,P_0 = 3.0;
        double C_S = sqrt(GAMMA*P_0/RHO_0);
        double KB = 2.0*3.1415/SIDE_LENGTH_X;
        double EPSILON = 0.1;
        double V;

        RHO = RHO_0*(1.0 + EPSILON * sin(KB * X));
        P = P_0 + C_S * C_S * (RHO - RHO_0);
        V = C_S * (RHO - RHO_0)/RHO_0 + 0.000000001;

        NEW_VERTEX.set_mass_density(RHO);                       // units kg/m^3
        NEW_VERTEX.set_x_velocity(V);                           // units m/s
        NEW_VERTEX.set_y_velocity(0.00000001);
        NEW_VERTEX.set_pressure(P);                             // units N/m^2

#endif
#ifdef SEDOV
    // if(i==0 and j==0){std::cout << "Using 2D Sedov Blast" << std::endl;}
        double RHO = 1.0;
        double V = 0.00000001;
        double P = 100.0;

        double R = sqrt((X - 5.0)*(X - 5.0) + (Y - 5.0)*(Y - 5.0));

        if(R < R_BLAST){
                P = 1000000.0;
                // P = 1000.0;
                std::cout << POINT_CHECK << "\tSetting blast pressure point at\t" << X << "\t" << Y << "\t" << P << std::endl;
                POINT_CHECK ++;
        }

        NEW_VERTEX.set_mass_density(RHO);                       // units kg/m^3
        NEW_VERTEX.set_x_velocity(V);                             // units m/s
        NEW_VERTEX.set_y_velocity(V);
        NEW_VERTEX.set_pressure(P);                             // units N/m^2

#endif
#ifdef GAUSSX
        // if(i==0 and j==0){std::cout << "Using 1D Gaussian pulse" << std::endl;}

        double CENTRE = 0.5;
        double S,W,RHO,RHO_0 = 10.0,RHO_PULSE = 50.0;
        double X_VELOCITY = 2.0,PRESSURE = 10.0;

        S = std::abs(CENTRE - X);                                       // distance from centre of pulse
        W = 0.1;                                                 // characteristic width

        RHO = RHO_PULSE*exp(-S*S/(W*W)) + RHO_0*(1.0-exp(-S*S/(W*W)));

        NEW_VERTEX.set_mass_density(RHO);
        NEW_VERTEX.set_x_velocity(X_VELOCITY);
        NEW_VERTEX.set_y_velocity(0.00000001);
        //NEW_VERTEX.set_y_velocity(0.0);
        NEW_VERTEX.set_pressure(PRESSURE);

#endif
#ifdef GAUSSY
        // if(i==0 and j==0){std::cout << "Using 1D Gaussian pulse (y-direction)" << std::endl;}

        double CENTRE = 0.3;
        double S,W,RHO,RHO_0 = 10.0,RHO_PULSE = 50.0;
        double Y_VELOCITY = 2.0,PRESSURE = 100.0;

        S = std::abs(CENTRE - Y);                                       // distance from centre of pulse
        W = 0.1;                                                 // characteristic width

        RHO = RHO_PULSE*exp(-S*S/(W*W)) + RHO_0*(1.0-exp(-S*S/(W*W)));

        NEW_VERTEX.set_mass_density(RHO);
        NEW_VERTEX.set_x_velocity(0.00000001);
        NEW_VERTEX.set_y_velocity(Y_VELOCITY);
        NEW_VERTEX.set_pressure(PRESSURE);

#endif
#ifdef UNIFORM
        // if(i==0 and j==0){std::cout << "Using Flat Start" << std::endl;}

        NEW_VERTEX.set_mass_density(1.0);                               // units kg/m^3
        NEW_VERTEX.set_x_velocity(10.0);                                 // units m/s
        NEW_VERTEX.set_y_velocity(0.00000001);                                 // units m/s
        NEW_VERTEX.set_pressure(5.0);                                 // units N/m^2

#endif
#ifdef NOH
        // if(i==0 and j==0){std::cout << "Using 2D Noh Problem" << std::endl;}

        double X_VEL,Y_VEL,X_VEL_NORM,Y_VEL_NORM;
        double X_C = SIDE_LENGTH_X/2.0;
        double Y_C = SIDE_LENGTH_Y/2.0;
        double RHO0 = 1.0;
        double PRESSURE = 0.000001;

        X_VEL = (X_C - X);
        Y_VEL = (Y_C - Y);

        X_VEL_NORM = X_VEL/sqrt(X_VEL*X_VEL + Y_VEL*Y_VEL);
        Y_VEL_NORM = Y_VEL/sqrt(X_VEL*X_VEL + Y_VEL*Y_VEL);

        if(X == X_C and Y == Y_C){X_VEL = Y_VEL = 0.00000001;}

        // std::cout << X << "\t" << Y << "\t" << X_VEL_NORM << "\t" << Y_VEL_NORM << "\t" << sqrt(X_VEL_NORM*X_VEL_NORM + Y_VEL_NORM*Y_VEL_NORM) << std::endl;

        NEW_VERTEX.set_mass_density(RHO0);
        NEW_VERTEX.set_x_velocity(X_VEL_NORM);
        NEW_VERTEX.set_y_velocity(Y_VEL_NORM);
        NEW_VERTEX.set_pressure(PRESSURE);

#endif
#ifdef KHX
        // if(i==0 and j==0){std::cout << "Using KH instability test (x flow)" << std::endl;}

        if(Y < 0.25*SIDE_LENGTH_Y or Y > 0.75*SIDE_LENGTH_Y){
                NEW_VERTEX.set_x_velocity(0.5);
                NEW_VERTEX.set_mass_density(1.0);
        }else{
                NEW_VERTEX.set_x_velocity(-0.5);
                NEW_VERTEX.set_mass_density(2.0);
        }

        NEW_VERTEX.set_y_velocity(0.05*sin((2.0*3.1415/SIDE_LENGTH_X)*X));
        NEW_VERTEX.set_pressure(2.5);

#endif
#ifdef KHY
        // if(i==0 and j==0){std::cout << "Using KH instability test (y flow)" << std::endl;}

        if(X < 0.25*SIDE_LENGTH_X or X > 0.75*SIDE_LENGTH_X){
                NEW_VERTEX.set_y_velocity(0.5);
                NEW_VERTEX.set_mass_density(1.0);
        }else{
                NEW_VERTEX.set_y_velocity(-0.5);
                NEW_VERTEX.set_mass_density(2.0);
        }

        NEW_VERTEX.set_x_velocity(0.05*sin((2.0*3.1415/SIDE_LENGTH_Y)*Y));
        NEW_VERTEX.set_pressure(2.5);

#endif
#ifdef KHXSMOOTH
        // if(i==0 and j==0){std::cout << "Using KH instability test (x flow)" << std::endl;}

        double VEL0  = -0.5, RHOL = 1.0, RHOH = 2.0;
        double WIDTH = 0.1;
        double DENS  = (RHOH - RHOL) * G(0.5+(Y-0.25)/WIDTH) * G(0.5-(Y-0.75)/WIDTH) + RHOL;
        double VEL   = 2.0 * VEL0    * G(0.5+(Y-0.25)/WIDTH) * G(0.5-(Y-0.75)/WIDTH) - VEL0;

        NEW_VERTEX.set_x_velocity(VEL);
        NEW_VERTEX.set_mass_density(DENS);

        NEW_VERTEX.set_y_velocity(0.05*sin((2.0*3.1415/SIDE_LENGTH_X)*X));
        NEW_VERTEX.set_pressure(2.5);

#endif
#ifdef KHYSMOOTH
        // if(i==0 and j==0){std::cout << "Using KH instability test (y flow)" << std::endl;}

        double VEL0  = -0.5, RHOL = 1.0, RHOH = 2.0;
        double WIDTH = 0.1;
        double DENS  = (RHOH - RHOL) * G(0.5+(X-0.25)/WIDTH) * G(0.5-(X-0.75)/WIDTH) + RHOL;
        double VEL   = 2.0 * VEL0    * G(0.5+(X-0.25)/WIDTH) * G(0.5-(X-0.75)/WIDTH) - VEL0;

        NEW_VERTEX.set_y_velocity(VEL);
        NEW_VERTEX.set_mass_density(DENS);

        NEW_VERTEX.set_x_velocity(0.05*sin((2.0*3.1415/SIDE_LENGTH_Y)*Y));
        NEW_VERTEX.set_pressure(2.5);

#endif
#ifdef BLOB
        // if(i==0 and j==0){std::cout << "Using Blob test" << std::endl;}
        double CENTRE_X = 0.25*SIDE_LENGTH_X;
        double CENTRE_Y = 0.50*SIDE_LENGTH_Y;
        double MACH = 1.5;

        double RADIUS = sqrt((X - CENTRE_X)*(X - CENTRE_X) + (Y - CENTRE_Y)*(Y - CENTRE_Y));

        if(RADIUS < 1.0){
                NEW_VERTEX.set_mass_density(100.0);
                NEW_VERTEX.set_x_velocity(0.00000001);
        }else{
                NEW_VERTEX.set_mass_density(10.0);
                NEW_VERTEX.set_x_velocity(MACH*4.1);
        }

        NEW_VERTEX.set_y_velocity(0.00000001);
        NEW_VERTEX.set_pressure(100.0);

#endif
#ifdef DF
        // if(i==0 and j==0){std::cout << "Grav Test" << std::endl;}
        double RHO0 = 1000.0;
        double X_VEL = 0.00000001;//-1.0*MACH*0.0387;
        double Y_VEL = 0.00000001;
        double PRESSURE = 1.0;

        NEW_VERTEX.set_mass_density(RHO0);
        NEW_VERTEX.set_x_velocity(X_VEL);
        NEW_VERTEX.set_y_velocity(Y_VEL);
        NEW_VERTEX.set_pressure(PRESSURE);
#endif
#ifdef GRESHO
     // if(i==0 and j==0){std::cout << "Using 1D Gresho pulse" << std::endl;}

        double X_CENTER = 1.0;
        double Y_CENTER = 1.0;
//        double X_CENTER = 0.5;
//       double Y_CENTER = 0.5;
        double R = sqrt((X - X_CENTER)*(X - X_CENTER) + (Y - Y_CENTER)*(Y - Y_CENTER));
        double RHO_0 = 1.0;

        double COS_PHI = (X-X_CENTER)/R;
        double SIN_PHI = (Y-Y_CENTER)/R;

        double VEL_PHI, PRESSURE;
        if (0.0 <= R and R < 0.2) {
            VEL_PHI = 5.0*R; PRESSURE = p0 + 25.0/2.0*R*R;
        }else if (0.2 <= R and R < 0.4){
            VEL_PHI = 2.0 - 5.0*R;
            PRESSURE = p0 + 4.0 - 4.0*log(0.2) + 25.0/2.0*R*R - 20.0*R + 4.0*log(R);
        }else{
            VEL_PHI = 0.0; PRESSURE = p0 - 2.0 + 4.0*log(2.0);
        }

        //smooth
//        double xsmooth = 0.02;
//        if (0.2-xsmooth <= R and R < 0.2+xsmooth){
//            VEL_PHI = -2.5/xsmooth*R*R + R/xsmooth + 1.0 - 0.1/xsmooth - 2.5*xsmooth;
//        }

        //VEL_PHI = -VEL_PHI; //clockwise
        //double X_VELOCITY = -VEL_PHI*SIN_PHI - 0.00000001;

        double X_VELOCITY = -VEL_PHI*SIN_PHI;
        double Y_VELOCITY = VEL_PHI*COS_PHI;

        if (R == 0.0){
            X_VELOCITY = 0.0; Y_VELOCITY = 0.0;
        }

        X_VELOCITY += 0.00000001;
        //X_VELOCITY += 1.0;

        NEW_VERTEX.set_mass_density(RHO_0);
        NEW_VERTEX.set_x_velocity(X_VELOCITY);
        NEW_VERTEX.set_y_velocity(Y_VELOCITY);
        NEW_VERTEX.set_pressure(PRESSURE);

#endif


        NEW_VERTEX.setup_specific_energy();
        NEW_VERTEX.prim_to_con();
        NEW_VERTEX.reset_du_half();
        NEW_VERTEX.reset_du();

        //setup initial mesh velocity
        NEW_VERTEX.set_sigma_x(1.5);
        NEW_VERTEX.set_sigma_y(0.0);

        return NEW_VERTEX;

}
