/* class containing values associated with the cell
	*vertex_0 = pointer to vertex ID of lower vertex
	*vertex_1 = pointer to vertex ID of upper vertex
	element_residual = value of residual of this element
	nodal_residual = value of residual distributed from this element to that node
*/

class cell{

private:

	int *vertex_0,*vertex_1;
	float element_residual,weight_0,weight_1;

public:

	void set_vertex_0(int *new_vertex_0){
		vertex_0 = new_vertex_0;
	}

	void set_vertex_1(int *new_vertex_1){
		vertex_1 = new_vertex_1;
	}

	int get_vertex_0(){
		return vertex_0;
	}

	int get_vertex_1(){
		return vertex_1;
	}

	/* 	calculate element residual (based on equation 7 of Deconinick et al.)
			res = element residual 
	*/
	void calc_element_residual(vertex vertex_lower, vertex vertex_upper){
		float res,u,u_plus,q,q_plus,q_half,f,f_plus,f_half,a_half,h_half;

		q = 0.0;		// ignore source terms for now 
		q_plus = 0.0;

		u = vertex_lower.get_density();
		u_plus = vertex.get_density();

		f = vertex_lower.get_f0();
		f_plus = vertex_upper.get_f0();

		h_half = (vertex_lower.get_x()+vertex_upper.get_x())/2.0;

		q_half = (q+q_plus)/2.0;
		a_half = (f_plus-f)/(u_plus-u);
		f_half = (f+f_plus)/2.0 - sign(a_half)/2.0*(f_plus-f);

		res = f_plus-f-q_half*h_half;

		element_residual = res;
	}


	// returns 1.0 for positive numbers and -1.0 for negative numbers
	float sign(float x){

		if(x>0.0){
			return 1.0;
		}elseif(x<0.0){
			return -1.0;
		}elseif(x==0.0){
			cout << "ERROR: SIGN FUNCTION UNDEFINED AT 0" << endl;
			return 0.0;
		}
		
	}

}