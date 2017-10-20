/*      class containing conserved and primative variables of fluid at the position of the centre
                x = position
                mass_density = mass_density of material in cell
                velocity = velocity of material in cell
                pressure = pressure in cell
                specific_energy = specific_energy of cell (?) 
                u_variables = array conatining values of vector U (see README)
                f_variables = array of values for vector F(U) (see README)
*/

using namespace std;

class centre{

private:

        double x;
        double mass_density,velocity,pressure,specific_energy;
        double u_variables[3],f_variables[3],du[3];
        vector<int> assoc_cells;

public:

        // setter functions preventing varaibles being changed accidentally
        // (no setter functions for U and F(U) as these are set by the other variables)
        void set_x(double new_x){
                x = new_x;
        }

        void set_mass_density(double new_mass_density){
                mass_density = new_mass_density;
        }

        void set_velocity(double new_velocity){
                velocity = new_velocity;
        }

        void set_pressure(double new_pressure){
                pressure = new_pressure;
        }

        void add_cell(int new_centre){
                assoc_cells.push_back(new_centre);
        }


        //getter functions for extracting values of variables
        double get_x(){
                return x;
        }

        double get_specific_energy(){
                return specific_energy;
        }

        double get_mass_density(){
                return mass_density;
        }

        double get_velocity(){
                return velocity;
        }

        double get_pressure(){
                return pressure;
        }

        double get_u0(){
                return u_variables[0];
        }

        double get_u1(){
                return u_variables[1];
        }

        double get_u2(){
                return u_variables[2];
        }

        double get_f0(){
                return f_variables[0];
        }

        double get_f1(){
                return f_variables[1];
        }

        double get_f2(){
                return f_variables[2];
        }
        

        // functions to set up the specific energy varaible, as well as u and f arrays
        void setup_specific_energy(){
                specific_energy = pressure/((GAMMA-1.0)*mass_density)+velocity*velocity/2.0;
        }

        void con_to_prim(){
                u_variables[0] = mass_density;
                u_variables[1] = mass_density*velocity;
                u_variables[2] = mass_density*specific_energy;
        }

        void setup_f_variables(){
                f_variables[0] = mass_density * velocity;
                f_variables[1] = mass_density * velocity * velocity + pressure;
                f_variables[2] = (mass_density * specific_energy + pressure) * velocity;
        }

        // convert conserved variables to primitive variables
        void prim_to_con(){
                mass_density = u_variables[0];
                velocity = u_variables[1]/mass_density;
                specific_energy = u_variables[2]/mass_density;
        }

        // recacluate pressure based on updated primitive varaibles
        void recalculate_pressure(){
                pressure = (GAMMA-1.0) * mass_density * (specific_energy - (velocity*velocity)/2.0);
        }

        // reset the changes in primative variables
        void reset_du(){
                du[0] = 0.0;
                du[1] = 0.0;
                du[2] = 0.0;
        }

        // update du with value from face
        void update_du(double new_du[3]){
                du[0] = du[0] + new_du[0];
                du[1] = du[1] + new_du[1];
                du[2] = du[2] + new_du[2];
        }

        // change primaitive vriables based on accumulated du from all faces of cell
        void update_u_variables(){
                u_variables[0] = u_variables[0] + du[0];
                u_variables[1] = u_variables[1] + du[1];
                u_variables[2] = u_variables[2] + du[2];
        }

        // calculate min timestep this cell requires
        void calc_next_dt(double dx, double cfl, double &next_dt){
                double c_sound;
                c_sound = sqrt(GAMMA*pressure/mass_density);
                next_dt = cfl*(dx/(c_sound+abs(velocity)));
        }

};