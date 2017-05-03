float simpson(int n_points,float lower, float upper, std::function<float (float)> func){
	float integral;
	float lower_arr,upper_arr;
	int i;

	for(i=0;i<n_points;i++){
		lower_arr = lower+(upper-lower)*float(i)/float(n_points);
		upper_arr = lower+(upper-lower)*float(i+1)/float(n_points);
		integral += (upper_arr-lower_arr)/6.0*(func(lower_arr)+4.0*func((lower_arr+upper_arr)/2)+func(upper_arr));
	}

	return integral;
}

float trapezium(int n_points, float lower, float upper, float lower_value, float upper_value){
	float integral;
	float lower_arr,upper_arr;
	int i;

	for(i=0;i<n_points;i++){
		lower_arr = lower+(upper-lower)*float(i)/float(n_points);
		upper_arr = lower+(upper-lower)*float(i+1)/float(n_points);
		integral += (upper_arr-lower_arr)*lower_value+0.5*(lower_value-upper_value);
	}

	return integral;
}