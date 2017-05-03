float compact_euler(float u,float delta_t,float delta_f,float delta_x,float function(float)){
	float delta_u;

	delta_u = (-1.0)*(delta_f/delta_x)*delta_t;

	return delta_u;
}

