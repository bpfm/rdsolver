/* class containing values associated with the cell
	*vertex_0 = pointer to vertex ID of lower vertex
	*vertex_1 = pointer to vertex ID of upper vertex
	element_residual = value of residual of this element
	nodal_residual = value of residual distributed from this element to that node
*/

class cell{

private:

	int *vertex_0,*vertex_1;
	float element_residual,nodal_residual_1,nodal_residual_2;

public:

	void set_vertex_0(int *new_vertex_0){
		vertex_0 = new_vertex_0;
	}

	void set_vertex_1(int *new_vertex_1){
		vertex_1 = new_vertex_1;
	}

	float calc_element_residual(){

	}

	float calc_nodal_residual(int vertex){

	}

}