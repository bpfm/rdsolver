/* 	class containing conserved and primative variables of fluid at the position of the vertex
		x = position
		mass = mass contained in cell
		momentum = momentum of material in cell
		energy = energy density (?) of cell
		density = density of material in cell
		velocity = velocity of material in cell
		pressure = pressure in cell
		u_variables = array conatining values of vector U (see README)
		f_variables = array of values for vector F(U) (see README)
*/

class vertex{

private:

	float x;
	float mass,momentum,energy;
	float density,velocity,pressure;
	float u_variables[3],f_variables[3];
	float residual;

public:


	// setter functions preventing varaibles being changed accidentally
	// (no setter functions for U and F(U) as these are set by the other variables)
	void set_x(float new_x){
		x = new_x;
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
		return f_variables[0];
	}

	float get_f1(){
		return f_variables[1];
	}

	float get_f2(){
		return f_variables[2];
	}

	//functions to set up variable and flux arrays
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