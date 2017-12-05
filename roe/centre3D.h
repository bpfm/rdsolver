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

        double x,y,z;
        double mass_density,pressure,specific_energy;
        double x_velocity,y_velocity,z_velocity;
        double u_variables[5],f_variables[5],du[5];
        vector<int> assoc_cells;

public:

        // setter functions preventing varaibles being changed accidentally
        // (no setter functions for U and F(U) as these are set by the other variables)
        void set_x(double new_x){
                x = new_x;
        }

        void set_y(double new_y){
                y = new_y;
        }

        void set_z(double new_z){
                z = new_z;
        }

        void set_mass_density(double new_mass_density){
                mass_density = new_mass_density;
        }

        void set_x_velocity(double new_velocity){
                x_velocity = new_velocity;
        }

        void set_y_velocity(double new_velocity){
                y_velocity = new_velocity;
        }

        void set_z_velocity(double new_velocity){
                z_velocity = new_velocity;
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

        double get_y(){
                return y;
        }

        double get_z(){
                return z;
        }

        double get_specific_energy(){
                return specific_energy;
        }

        double get_mass_density(){
                return mass_density;
        }

        double get_x_velocity(){
                return x_velocity;
        }

        double get_y_velocity(){
                return y_velocity;
        }

        double get_z_velocity(){
                return z_velocity;
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

        double get_u3(){
                return u_variables[3];
        }

        double get_u4(){
                return u_variables[4];
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

        double get_f3(){
                return f_variables[3];
        }

        double get_f4(){
                return f_variables[4];
        }


        // functions to set up the specific energy varaible, as well as u and f arrays
        void setup_specific_energy(){
                double vel_sq_sum = x_velocity*x_velocity + y_velocity*y_velocity + z_velocity*z_velocity;
                specific_energy = pressure/((GAMMA-1.0)*mass_density) + vel_sq_sum/2.0;
        }

        void prim_to_con(){
                u_variables[0] = mass_density;
                u_variables[1] = mass_density * x_velocity;
                u_variables[2] = mass_density * y_velocity;
                u_variables[3] = mass_density * z_velocity;
                u_variables[4] = mass_density * specific_energy;
        }

        void setup_f_variables(){
                f_variables[0] = mass_density * x_velocity;
                f_variables[1] = mass_density * x_velocity * x_velocity + pressure;
                f_variables[2] = mass_density * x_velocity * y_velocity;
                f_variables[3] = mass_density * x_velocity * z_velocity;
                f_variables[4] = (mass_density * specific_energy + pressure) * x_velocity;
        }

        // convert conserved variables to primitive variables
        void con_to_prim(){
                mass_density = u_variables[0];
                x_velocity = u_variables[1]/mass_density;
                y_velocity = u_variables[2]/mass_density;
                z_velocity = u_variables[3]/mass_density;
                specific_energy = u_variables[4]/mass_density;
        }

        // recacluate pressure based on updated primitive varaibles
        void recalculate_pressure(){
                double vel_sq_sum = x_velocity*x_velocity + y_velocity*y_velocity + z_velocity*z_velocity;
                pressure = (GAMMA-1.0) * mass_density * (specific_energy - (vel_sq_sum)/2.0);
        }

        // reset the changes in primative variables
        void reset_du(){
                du[0] = 0.0;
                du[1] = 0.0;
                du[2] = 0.0;
                du[3] = 0.0;
                du[4] = 0.0;
        }

        // update du with value from face
        void update_du(double new_du[3]){
                du[0] = du[0] + new_du[0];
                du[1] = du[1] + new_du[1];
                du[2] = du[2] + new_du[2];
                du[3] = du[3] + new_du[3];
                du[4] = du[4] + new_du[4];
        }

        // change primaitive vriables based on accumulated du from all faces of cell
        void update_u_variables(){
                u_variables[0] = u_variables[0] + du[0];
                u_variables[1] = u_variables[1] + du[1];
                u_variables[2] = u_variables[2] + du[2];
                u_variables[3] = u_variables[3] + du[3];
                u_variables[4] = u_variables[4] + du[4];
        }

        // calculate min timestep this cell requires (pick smallest value for all faces)
        void calc_next_dt(double dx, double cfl, double &next_dt){
                int i;
                double c_sound;
                double next_dt_pick[3];
                c_sound = sqrt(GAMMA*pressure/mass_density);
                next_dt_pick[0] = cfl*(dx/(c_sound+abs(x_velocity)));
                next_dt_pick[1] = cfl*(dx/(c_sound+abs(y_velocity)));
                next_dt_pick[2] = cfl*(dx/(c_sound+abs(z_velocity)));
                
                next_dt=next_dt_pick[0];

                for(i=1;i<3;i++){if(next_dt_pick[i] < next_dt){next_dt = next_dt_pick[i];}}
        }

};