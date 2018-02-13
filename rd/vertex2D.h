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
                U_VARIABLES = array conatining values of vector U (see README)
*/

using namespace std;

class VERTEX{

private:

        double X,Y,DX,DY;
        double DUAL;
        double MASS_DENSITY,X_VELOCITY,Y_VELOCITY,PRESSURE,SPECIFIC_ENERGY;
        double U_VARIABLES[4],F_VARIABLES[4],DU[4];
        vector<int> ASSOC_TRIANG;

public:

        // setter functions preventing varaibles being changed accidentally
        // (no setter functions for U and F(U) as these are set by the other variables)
        void set_x(double NEW_X){
                X = NEW_X;
        }

        void set_y(double NEW_Y){
                Y = NEW_Y;
        }

        void set_dx(double NEW_DX){
                DX = NEW_DX;
        }

        void set_dy(double NEW_DY){
                DY = NEW_DY;
        }

        void set_mass_density(double NEW_MASS_DENSITY){
                MASS_DENSITY = NEW_MASS_DENSITY;
        }

        void set_x_velocity(double NEW_X_VELOCITY){
                X_VELOCITY = NEW_X_VELOCITY;
        }

        void set_y_velocity(double NEW_Y_VELOCITY){
                Y_VELOCITY = NEW_Y_VELOCITY;
        }

        void set_pressure(double NEW_PRESSURE){
                PRESSURE = NEW_PRESSURE;
        }

        void add_triang(int NEW_TRIANGLE){
                ASSOC_TRIANG.push_back(NEW_TRIANGLE);
        }

        //getter functions for eXtracting values of variables
        double get_x(){
                return X;
        }

        double get_y(){
                return Y;
        }

        double get_dx(){
                return DX;
        }

        double get_dy(){
                return DY;
        }

        double get_specific_energy(){
                return SPECIFIC_ENERGY;
        }

        double get_mass_density(){
                return MASS_DENSITY;
        }

        double get_x_velocity(){
                return X_VELOCITY;
        }

        double get_y_velocity(){
                return Y_VELOCITY;
        }


        double get_pressure(){
                return PRESSURE;
        }

        double get_u0(){
                return U_VARIABLES[0];
        }

        double get_u1(){
                return U_VARIABLES[1];
        }

        double get_u2(){
                return U_VARIABLES[2];
        }

        double get_u3(){
                return U_VARIABLES[3];
        }


        // functions to set up the specific energy varaible, as well as u and f arrays
        void setup_specific_energy(){
                double VEL_SQ_SUM = X_VELOCITY*X_VELOCITY + Y_VELOCITY*Y_VELOCITY;
                SPECIFIC_ENERGY = PRESSURE/((GAMMA-1.0)*MASS_DENSITY) + VEL_SQ_SUM/2.0;
        }

        void calculate_dual(){
                DUAL = 6.0/(2.0*3.0)*DX*DY;
        }

        void prim_to_con(){
                U_VARIABLES[0] = MASS_DENSITY;
                U_VARIABLES[1] = MASS_DENSITY * X_VELOCITY;
                U_VARIABLES[2] = MASS_DENSITY * Y_VELOCITY;
                U_VARIABLES[3] = MASS_DENSITY * SPECIFIC_ENERGY;
        }

        // convert conserved variables to primitive variables
        void con_to_prim(){
                MASS_DENSITY = U_VARIABLES[0];
                X_VELOCITY = U_VARIABLES[1]/MASS_DENSITY;
                Y_VELOCITY = U_VARIABLES[2]/MASS_DENSITY;
                SPECIFIC_ENERGY = U_VARIABLES[3]/MASS_DENSITY;
        }

        // recacluate PRESSURE based on updated primitive varaibles
        void recalculate_pressure(){
                double VEL_SQ_SUM = X_VELOCITY*X_VELOCITY + Y_VELOCITY*Y_VELOCITY;
                PRESSURE = (GAMMA-1.0) * MASS_DENSITY * (SPECIFIC_ENERGY - VEL_SQ_SUM/2.0);
        }

        // reset the changes in primative variables
        void reset_du(){
                DU[0] = 0.0;
                DU[1] = 0.0;
                DU[2] = 0.0;
                DU[3] = 0.0;
        }

        // update DU with value from face
        void update_du(double NEW_DU[4]){
                DU[0] = DU[0] + NEW_DU[0];
                DU[1] = DU[1] + NEW_DU[1];
                DU[2] = DU[2] + NEW_DU[2];
                DU[3] = DU[3] + NEW_DU[3];
        }

        // change primaitive vriables based on accumulated DU from all faces of cell
        void update_u_variables(){
                U_VARIABLES[0] = U_VARIABLES[0] + DU[0];
                U_VARIABLES[1] = U_VARIABLES[1] + DU[1];
                U_VARIABLES[2] = U_VARIABLES[2] + DU[2];
                U_VARIABLES[3] = U_VARIABLES[3] + DU[3];
        }

        // calculate min timestep this cell requires
        void calc_next_dt(double DX, double CFL, double &NEXT_DT){
                double C_SOUND;
                C_SOUND = sqrt(GAMMA*PRESSURE/MASS_DENSITY);
                NEXT_DT = CFL*(DX/(C_SOUND+abs(X_VELOCITY)+abs(Y_VELOCITY)));
        }

};