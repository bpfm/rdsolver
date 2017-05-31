/* class containing values associated with the cell
	*vertex_0 = pointer to vertex ID of lower vertex
	*vertex_1 = pointer to vertex ID of upper vertex
*/

using namespace std;

class cell{

private:

	vertex *vertex_0,*vertex_1;

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

	//float get_element_residual(){
	//	return element_residual;
	//}

	/* 	calculate element residual (based on equation 7 of Deconinick et al.)
			res = element residual 
	*/
	void distribute_residual(float dx, float dt){
		int i;
		float res[3],u[3],u_plus[3],q[3],q_plus[3],q_half[3],f[3],f_plus[3],f_half[3],a_half[3],beta_0[3],beta_1[3],du_0[3],du_1[3];
		float h_half;

		u[0] = vertex_0->get_u0();
		u_plus[0] = vertex_1->get_u0();

		u[1] = vertex_0->get_u1();
		u_plus[1] = vertex_1->get_u1();

		u[2] = vertex_0->get_u2();
		u_plus[2] = vertex_1->get_u2();

		f[0] = vertex_0->get_f0();
		f_plus[0] = vertex_1->get_f0();

		f[1] = vertex_0->get_f1();
		f_plus[1] = vertex_1->get_f1();

		f[2] = vertex_0->get_f2();
		f_plus[2] = vertex_1->get_f2();

		h_half = (vertex_0->get_x()+vertex_1->get_x())/2.0;

		for(i=0;i<3;i++){
			q[i] = 0.0;		// ignore source terms for now 
			q_plus[i] = 0.0;

			q_half[i] = (q[i]+q_plus[i])/2.0;
			a_half[i] = (f_plus[i]-f[i])/(u_plus[i]-u[i]);
			f_half[i] = (f[i]+f_plus[i])/2.0 - sign(a_half[i])/2.0*(f_plus[i]-f[i]);

			res[i] = f_plus[i]-f[i]-q_half[i]*h_half;

			beta_0[i] = 0.5*(1.0-(a_half[i])/abs(a_half[i]));
			beta_1[i] = 0.5*(1.0+(a_half[i])/abs(a_half[i]));

			du_0[i] = beta_0[i]*res[i]*dt/dx;
			du_1[i] = beta_1[i]*res[i]*dt/dx;
		}

		cout << "weights" << " " << beta_0[0] << " " << beta_1[0] << " " << a_half[0] << " " << u[0] << " " << u_plus[0] << endl;

		vertex_0->update_du(du_0);
		vertex_1->update_du(du_1);

	}

	// returns 1.0 for positive numbers and -1.0 for negative numbers
	float sign(float x){
		if(x>0.0){
			return 1.0;
		}else if(x<0.0){
			return -1.0;
		}else{
			return 0.0;
		}
		
	}

};