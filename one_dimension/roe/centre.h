/* 	class containing conserved and primative variables of fluid at the position of the centre
		x = position
		mass_density = mass_density of material in cell
		velocity = velocity of material in cell
		pressure = pressure in cell
		energy_density = energy_density of cell (?) 
		u_variables = array conatining values of vector U (see README)
		f_variables = array of values for vector F(U) (see README)
*/

using namespace std;

class centre{

private:

	double x;
	double mass_density,velocity,pressure,energy_density;
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

	double get_energy_density(){
		return energy_density;
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
	

	//functions to set up energy density varaible, and u and f arrays
	void setup_energy_density(){
		double gamma=1.4;
		energy_density = pressure/((gamma-1.0)*mass_density)+velocity*velocity/2.0;
	}

	void con_to_prim(){
		u_variables[0] = mass_density;
		u_variables[1] = mass_density*velocity;
		u_variables[2] = mass_density*energy_density;
	}

	void prim_to_con(){
		mass_density = u_variables[0];
		velocity = u_variables[1]/mass_density;
		energy_density = u_variables[2]/mass_density;
	}

	void recalculate_pressure(){
		double gamma=1.4;
//		cout << x << " " << pressure << endl;
		pressure = (energy_density*(gamma-1.0)*mass_density) - ((velocity*velocity)/2.0);
		if(x==13.5){cout << x << " " << pressure << " " << energy_density << " " << velocity << " " << mass_density << endl;}
	}

	void setup_f_variables(){
		f_variables[0] = mass_density * velocity;
		f_variables[1] = mass_density * velocity * velocity + pressure;
		f_variables[2] = (mass_density * energy_density + pressure) * velocity;
	}

	void reset_du(){
		du[0] = 0.0;
		du[1] = 0.0;
		du[2] = 0.0;
	}

	void update_du(double new_du[3]){
		du[0] = du[0] + new_du[0];
		du[1] = du[1] + new_du[1];
		du[2] = du[2] + new_du[2];
	}

	void update_u_variables(){
		u_variables[0] = u_variables[0] + du[0];
		u_variables[1] = u_variables[1] + du[1];
		u_variables[2] = u_variables[2] + du[2];
	}

	void calc_next_dt(double dx, double cfl, double &next_dt){
		double c_sound,gamma=1.4;
		c_sound = sqrt(gamma*pressure/mass_density);
		next_dt = cfl*(dx/(c_sound+abs(velocity)));
	}

};