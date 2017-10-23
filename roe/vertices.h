#include vertex.h

class vertices{
private:

	std::vector<vertex> points;

public:

	void add_member(vertex new_vertex){
		points.pushback(new_vertex);
	}

	vector get_member(int i){
		return points[i];
	}

};