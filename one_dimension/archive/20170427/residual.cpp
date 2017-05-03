extern float simpson(float n_points,float lower, float upper, std::function<float (float)> func);

//equivaent to Phi
float element_residual(float lower,float upper, float q0, float q1){
	int n_points = 100;
	float residual;

	residual = simpson(n_points, lower, upper,conserve);

	return residual;
}

//equivalent to Phi_i
float nodal_residual(){
	float residual;

	return residual;
}

//equivalent to r(Q)
float conserve(float q){
	float r,delta_q,delta_t,delta_f,delta_x;

	r = delta_q/delta_t + delta_f/delta_x;
	return r;
}

//equivalent to N_i(x)
float basis_function(float x, float x_point){
	float out; 

	out = 1.0 + (x-x_point)/(2.0);

	return out;
}

//equivalent to Q_h
float approximation(float x,float lower,float upper,float lower_value,float upper_value){
	int i;
	float q_app,x_point,value;

	for (i=0;i<2;i++){
		if(i==0){
			x_point = lower;
			value = lower_value;
		}else{
			x_point = upper;
			value = upper_value;
		}
		q_app += basis_function(x,x_point)*value;
	}

	return q_app;
}