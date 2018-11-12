/*      class containing conserved and primative variables of fluid at the position of the TRIANGLE
                X = x position
                Y = y position
                DX = change in x between vertices
                DY = change in y between vertices
                DUAL = dual cell area
                MASS_DENSITY = mass_density of material in cell
                VELOCITY = velocity of material in cell
                PRESSURE = pressure in cell
                SPECIFIC_ENERGY = SPECIFIC_ENERGY of cell (?) 
                U_VARIABLES = array conatining values of std::vector U (see README)
*/

class VERTEX{

private:

        double X, Y, DX, DY;
        double DUAL;
        double U_VARIABLES[4], DU[4];
        double MASS_DENSITY, X_VELOCITY, Y_VELOCITY;
        double PRESSURE, SPECIFIC_ENERGY;
        double U_HALF[4], DU_HALF[4];
        double MASS_DENSITY_HALF, X_VELOCITY_HALF, Y_VELOCITY_HALF;
        double PRESSURE_HALF, SPECIFIC_ENERGY_HALF;
#ifdef UPDATE_TEST
        double UPDATE;
        double UPDATE_HALF;
#endif
        std::vector<int> ASSOC_TRIANG;

public:

        // setter functions preventing varaibles being changed accidentally
        // (no setter functions for U and F(U) as these are set by the other variables)
        void set_x( double NEW_X){X   = NEW_X;}
        void set_y( double NEW_Y){Y   = NEW_Y;}
        void set_dx(double NEW_DX){DX = NEW_DX;}
        void set_dy(double NEW_DY){DY = NEW_DY;}
        void set_mass_density(double NEW_MASS_DENSITY){MASS_DENSITY = NEW_MASS_DENSITY;}
        void set_x_velocity(  double NEW_X_VELOCITY){  X_VELOCITY   = NEW_X_VELOCITY;}
        void set_y_velocity(  double NEW_Y_VELOCITY){  Y_VELOCITY   = NEW_Y_VELOCITY;}
        void set_pressure(    double NEW_PRESSURE){    PRESSURE     = NEW_PRESSURE;}

        // add triangle to list of those associated with this vertex

        void add_triang(int NEW_TRIANGLE){ASSOC_TRIANG.push_back(NEW_TRIANGLE);}


        // getter functions for eXtracting values of variables
        double get_x(){   return X;}
        double get_y(){   return Y;}
        double get_dx(){  return DX;}
        double get_dy(){  return DY;}
        double get_dual(){return DUAL;}

        double get_specific_energy(){return SPECIFIC_ENERGY;}
        double get_mass_density(){   return MASS_DENSITY;}
        double get_x_velocity(){     return X_VELOCITY;}
        double get_y_velocity(){     return Y_VELOCITY;}
        double get_pressure(){       return PRESSURE;}
        double get_u0(){return U_VARIABLES[0];}
        double get_u1(){return U_VARIABLES[1];}
        double get_u2(){return U_VARIABLES[2];}
        double get_u3(){return U_VARIABLES[3];}

        double get_specific_energy_half(){return SPECIFIC_ENERGY_HALF;}
        double get_mass_density_half(){   return MASS_DENSITY_HALF;}
        double get_x_velocity_half(){     return X_VELOCITY_HALF;}
        double get_y_velocity_half(){     return Y_VELOCITY_HALF;}
        double get_pressure_half(){       return PRESSURE_HALF;}
        double get_u0_half(){return U_HALF[0];}
        double get_u1_half(){return U_HALF[1];}
        double get_u2_half(){return U_HALF[2];}
        double get_u3_half(){return U_HALF[3];}


        // functions to set up the specific energy varaible, as well as u and f arrays
        void setup_specific_energy(){
                double VEL_SQ_SUM = X_VELOCITY*X_VELOCITY + Y_VELOCITY*Y_VELOCITY;
                SPECIFIC_ENERGY = PRESSURE/((GAMMA-1.0)*MASS_DENSITY) + VEL_SQ_SUM/2.0; // calculate specific energy
        }

        void calculate_dual(){
                DUAL = 6.0/(2.0*3.0)*DX*DY;
        }

        void prim_to_con(){
                U_VARIABLES[0] = MASS_DENSITY;                          // U0 = mass density
                U_VARIABLES[1] = MASS_DENSITY * X_VELOCITY;             // U1 = x momentum
                U_VARIABLES[2] = MASS_DENSITY * Y_VELOCITY;             // U2 = y momentum
                U_VARIABLES[3] = MASS_DENSITY * SPECIFIC_ENERGY;        // U3 = energy density
        }

        void reset_u_half(){
                U_HALF[0] = U_VARIABLES[0];
                U_HALF[1] = U_VARIABLES[1];
                U_HALF[2] = U_VARIABLES[2];
                U_HALF[3] = U_VARIABLES[3];
        }

        // convert conserved variables to primitive variables
        void con_to_prim(){
                MASS_DENSITY    = U_VARIABLES[0];
                X_VELOCITY      = U_VARIABLES[1]/MASS_DENSITY;
                Y_VELOCITY      = U_VARIABLES[2]/MASS_DENSITY;
                SPECIFIC_ENERGY = U_VARIABLES[3]/MASS_DENSITY;
                check_values();
        }

        void con_to_prim_half(){
                MASS_DENSITY_HALF    = U_HALF[0];
                X_VELOCITY_HALF      = U_HALF[1]/MASS_DENSITY;
                Y_VELOCITY_HALF      = U_HALF[2]/MASS_DENSITY;
                SPECIFIC_ENERGY_HALF = U_HALF[3]/MASS_DENSITY;
                check_values();
        }

        // recacluate PRESSURE based on updated primitive varaibles
        void recalculate_pressure(){
                double VEL_SQ_SUM = X_VELOCITY*X_VELOCITY + Y_VELOCITY*Y_VELOCITY;
                PRESSURE = (GAMMA-1.0) * MASS_DENSITY * (SPECIFIC_ENERGY - VEL_SQ_SUM/2.0);
        }

        void recalculate_pressure_half(){
                double VEL_SQ_SUM = X_VELOCITY_HALF*X_VELOCITY_HALF + Y_VELOCITY_HALF*Y_VELOCITY_HALF;
                PRESSURE_HALF = (GAMMA-1.0) * MASS_DENSITY_HALF * (SPECIFIC_ENERGY_HALF - VEL_SQ_SUM/2.0);
        }

        // reset the changes in primative variables
        void reset_du(){     DU[0]      = DU[1]      = DU[2]      = DU[3]      = 0.0;}// UPDATE = 0;}
        void reset_du_half(){DU_HALF[0] = DU_HALF[1] = DU_HALF[2] = DU_HALF[3] = 0.0;}// UPDATE_HALF = 0;}

        // Update DU with value from face
        void update_du(double NEW_DU[4]){
#ifdef UPDATE_TEST
                UPDATE_HALF += 1;
#endif
                DU[0] = DU[0] + NEW_DU[0];
                DU[1] = DU[1] + NEW_DU[1];
                DU[2] = DU[2] + NEW_DU[2];
                DU[3] = DU[3] + NEW_DU[3];
        }

        void update_du_half(double NEW_DU[4]){
#ifdef UPDATE_TEST
                UPDATE_HALF += 1;
#endif
                DU_HALF[0] = DU_HALF[0] + NEW_DU[0];
                DU_HALF[1] = DU_HALF[1] + NEW_DU[1];
                DU_HALF[2] = DU_HALF[2] + NEW_DU[2];
                DU_HALF[3] = DU_HALF[3] + NEW_DU[3];
        }

#ifdef UPDATE_TEST
        void output_update_count(){std::cout << X << "\t" << Y << "\t" << UPDATE_HALF << std::endl;}
#endif

        // Change conserved variables based on accumulated DU from all faces of cell
        void update_u_variables(){
                //std::cout << DU[0] << "\t" << DU[1] << "\t" << DU[2] << "\t" << DU[3] << std::endl;
                U_VARIABLES[0] = U_HALF[0] - DU[0];
                U_VARIABLES[1] = U_HALF[1] - DU[1];
                U_VARIABLES[2] = U_HALF[2] - DU[2];
                U_VARIABLES[3] = U_HALF[3] - DU[3];
        }

        void update_u_half(){
                //std::cout << DU_HALF[0] << "\t" << DU_HALF[1] << "\t" << DU_HALF[2] << "\t" << DU_HALF[3] << std::endl;
                U_HALF[0] = U_VARIABLES[0] - DU_HALF[0];
                U_HALF[1] = U_VARIABLES[1] - DU_HALF[1];
                U_HALF[2] = U_VARIABLES[2] - DU_HALF[2];
                U_HALF[3] = U_VARIABLES[3] - DU_HALF[3];
        }

        void check_values(){
#ifdef DEBUG
                std::cout << "Checking vertex state at " << X << "\t" << Y << std::endl;
#endif
                if (MASS_DENSITY < 0.0){
                        std::cout << "B WARNING: Exiting on negative density\t";
                        std::cout << "Position =\t" << X << "\t" << Y << std::endl;
                        //MASS_DENSITY = 0.000001;
                        exit(0);
                }
                if (PRESSURE < 0.0){
                        std::cout << "B WARNING: Exiting on negative pressure\t";
                        std::cout << "Position =\t" << X << "\t" << Y << std::endl;
                        
                        // PRESSURE = 0.000001;
                        exit(0);
                }
                if (MASS_DENSITY_HALF < 0.0){
                        std::cout << "B WARNING: Exiting on negative half state density\t";
                        std::cout << "Position =\t" << X << "\t" << Y << std::endl;
                        //MASS_DENSITY_HALF = 0.000001;
                        exit(0);
                }
                if (PRESSURE_HALF < 0.0){
                        std::cout << "B WARNING: Exiting on negative half state pressure\t";
                        std::cout << "Position =\t" << X << "\t" << Y << std::endl;
                        // PRESSURE_HALF = 0.000001;
                        exit(0);
                }

        }

        // Calculate min timestep this cell requires
        void calc_next_dt(double DX, double CFL, double &NEXT_DT){
                double C_SOUND = sqrt(GAMMA*PRESSURE/MASS_DENSITY);
                double V_MAX = max_val(abs(X_VELOCITY)+C_SOUND,abs(Y_VELOCITY)+C_SOUND);
                NEXT_DT = 2.0*CFL*DUAL/(6.0*DX*V_MAX);
                //std::cout << NEXT_DT << std::endl;
        }

        double max_val(double A, double B){
                if(A>B){
                        return A;
                }else{
                        return B;
                }
        }

};