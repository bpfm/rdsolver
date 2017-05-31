/* 	class containing conserved and primative variables of fluid at the position of the vertex
		x = position
		mass = mass contained in cell
		momentum = momentum of material in cell
		energy_density = energy_density mass_density (?) of cell
		mass_density = mass_density of material in cell
		velocity = velocity of material in cell
		pressure = pressure in cell
		u_variables = array conatining values of vector U (see README)
		f_variables = array of values for vector F(U) (see README)
*/

class vertex{

private:

	float x;
//	float mass,momentum,energy;
	float mass_density,velocity,pressure,energy_density;
	float u_variables[3],f_variables[3],du[3];
//	float residual;
//	float residual_cont_0,residual_cont_1;

public:

	// setter functions preventing varaibles being changed accidentally
	// (no setter functions for U and F(U) as these are set by the other variables)
	void set_x(float new_x){
		x = new_x;
	}

//	void set_mass(float new_mass){
//		mass = new_mass;
//	}

//	void set_momentum(float new_momentum){
//		momentum = new_momentum;
//	}

//	void set_energy(float new_energy){
//		energy = new_energy;
//	}

	void set_mass_density(float new_mass_density){
		mass_density = new_mass_density;
	}

	void set_velocity(float new_velocity){
		velocity = new_velocity;
	}

	void set_pressure(float new_pressure){
		pressure = new_pressure;
	}

//	void set_residual_cont_0(float new_res){
//		residual_cont_0 = new_res;
//	}

//	void set_residual_cont_1(float new_res){
//		residual_cont_1 = new_res;
//	}

	//getter functions for extracting values of variables
	float get_x(){
		return x;
	}

//	float get_mass(){
//		return mass;
//	}

//	float get_momentum(){
//		return momentum;
//	}

//	float get_energy(){
//		return energy;
//	}

	float get_energy_density(){
		return energy_density;
	}

	float get_mass_density(){
		return mass_density;
	}

	float get_velocity(){
		return velocity;
	}

	float get_pressure(){
		return pressure;
	}

	float get_u0(){
		return u_variables[0];
	}

	float get_u1(){
		return u_variables[1];
	}

	float get_u2(){
		return u_variables[2];
	}

	float get_f0(){
		return f_variables[0];
	}

	float get_f1(){
		return f_variables[1];
	}

	float get_f2(){
		return f_variables[2];
	}

//	float get_residual_cont_0(){
//		return residual_cont_0;
//	}

//	float get_residual_cont_1(){
//		return residual_cont_1;
//	}

	//functions to set up energy density varaible, and u and f arrays
	void setup_energy_density(){
		energy_density = pressure/(0.4*mass_density);
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

	void setup_f_variables(){
		f_variables[0] = mass_density*velocity;
		f_variables[1] = mass_density*velocity*velocity + pressure;
		f_variables[2] = (mass_density*energy_density+pressure) * velocity;
	}

	void reset_du(){
		du[0] = 0.0;
		du[1] = 0.0;
		du[2] = 0.0;
	}

	void update_du(float new_du[3]){
		du[0] = du[0] + new_du[0];
		du[1] = du[1] + new_du[1];
		du[2] = du[2] + new_du[2];
	}

	void update_u_variables(){
		u_variables[0] = u_variables[0] + du[0];
		u_variables[1] = u_variables[1] + du[1];
		u_variables[2] = u_variables[2] + du[2];
	}

};