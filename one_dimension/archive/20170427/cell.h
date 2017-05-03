/* 	class containing conserved and primative variables of fluid
		x = position
		mass = mass contained in cell
		momentum = momentum of material in cell
		energy = energy density (?) of cell
		density = density of material in cell
		velocity = velocity of material in cell
		pressure = pressure in cell
*/

extern float conserve(float q);

class cell{

private:

	float x,vertex_lower,vertex_upper;
	float mass,momentum,energy;
	float density,velocity,pressure;
	float u_variables[3],fluxes[3];
	float residual;

public:

	//setter functions preventing varaibles being changes accidentally
	void set_x(float new_x){
		x = new_x;
	}

	void set_vertex_lower(float new_vertex_lower){
		vertex_lower = new_vertex_lower;
	}

	void set_vertex_upper(float new_vertex_upper){
		vertex_upper = new_vertex_upper;
	}

	void set_mass(float new_mass){
		mass = new_mass;
	}

	void set_momentum(float new_momentum){
		momentum = new_momentum;
	}

	void set_energy(float new_energy){
		energy = new_energy;
	}

	void set_density(float new_density){
		density = new_density;
	}

	void set_velocity(float new_velocity){
		velocity = new_velocity;
	}

	void set_pressure(float new_pressure){
		pressure = new_pressure;
	}


	//getter functions for extracting values of variables
	float get_x(){
		return x;
	}

	float get_vertex_lower(){
		return vertex_lower;
	}

	float get_vertex_upper(){
		return vertex_upper;
	}

	float get_mass(){
		return mass;
	}

	float get_momentum(){
		return momentum;
	}

	float get_energy(){
		return energy;
	}

	float get_density(){
		return density;
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
		return fluxes[0];
	}

	float get_f1(){
		return fluxes[1];
	}

	float get_f2(){
		return fluxes[2];
	}


	//functions to convert between conserved and primative variables
	void con_to_prim(){
		float gam = 1.8;
		density = mass/fabs(vertex_upper - vertex_lower);
		velocity = momentum/mass;
		pressure = density*(gam-1.0)*energy;
	}

	void prim_to_con(){
		float gam = 1.8;
		mass = density*(fabs(vertex_upper - vertex_lower));
		momentum = density*velocity;
		energy = pressure/(density*(gam-1.0));
	}


	//function to set up variable and flux arrays
	void setup_u_variables(){
		u_variables[0] = density;
		u_variables[1] = density*velocity;
		u_variables[2] = density*energy;
	}

	void setup_fluxes(){
		fluxes[0] = density*velocity;
		fluxes[1] = density*velocity*velocity + pressure;
		fluxes[2] = (density*energy+pressure) * velocity;
	}

};