/* class containing values associated with the cell
	*vertex_0 = pointer to vertex ID of lower vertex
	*vertex_1 = pointer to vertex ID of upper vertex
	element_residual = value of residual of this element
	nodal_residual = value of residual distributed from this element to that node
*/

using namespace std;

class cell{

private:

	vertex *vertex_0,*vertex_1;
	float element_residual,weight_0,weight_1;

public:

	void set_vertex_0(vertex* new_vertex_0){
		vertex_0 = new_vertex_0;
	}

	void set_vertex_1(vertex* new_vertex_1){
		vertex_1 = new_vertex_1;
	}

	vertex* get_vertex_0(){
		return vertex_0;
	}

	vertex* get_vertex_1(){
		return vertex_1;
	}

	float get_element_residual(){
		return element_residual;
	}

	/* 	calculate element residual (based on equation 7 of Deconinick et al.)
			res = element residual 
	*/
	void calc_element_residual(){
		float res,u,u_plus,q,q_plus,q_half,f,f_plus,f_half,a_half,h_half;

		q = 0.0;		// ignore source terms for now 
		q_plus = 0.0;

		u = vertex_0->get_mass_density();
		u_plus = vertex_1->get_mass_density();

		f = vertex_0->get_f0();
		f_plus = vertex_1->get_f0();

		h_half = (vertex_0->get_x()+vertex_1->get_x())/2.0;

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
		}else if(x<0.0){
			return -1.0;
		}else{
			//cout << "ERROR: SIGN FUNCTION UNDEFINED AT 0" << endl;
			return 0.0;
		}
		
	}

	float calcute_weights(){
		float beta_0,beta_1;

		beta_0 = 1.0/2.0
	}

};