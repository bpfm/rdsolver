/* class containing values associated with the face
        *VERTEX_0 = pointer to left VERTEX
        *VERTEX_1 = pointer to right VERTEX
*/

using namespace std;

class TRIANGLE{

private:

        int ID;
        VERTEX *VERTEX_0,*VERTEX_1,*VERTEX_2;

public:

        void set_id(int NEW_ID){
                ID = NEW_ID;
        }

        void set_vertex_0(VERTEX* NEW_VERTEX){
                VERTEX_0 = NEW_VERTEX;
        }

        void set_vertex_1(VERTEX* NEW_VERTEX){
                VERTEX_1 = NEW_VERTEX;
        }

        void set_vertex_2(VERTEX* NEW_VERTEX){
                VERTEX_2 = NEW_VERTEX;
        }

        VERTEX* get_vertex_0(){
                return VERTEX_0;
        }

        VERTEX* get_vertex_1(){
                return VERTEX_1;
        }

        VERTEX* get_vertex_2(){
                return VERTEX_2;
        }

        void calculate_change(double DX, double DT, double T){
                double U_N[3][4],U_HALF[3][4];
                double DU0[4],DU1[4],DU2[4];

                //cout << VERTEX_0->get_x() << "\t" << VERTEX_0->get_y() << endl;

                // U_N[0][0] = VERTEX_0->get_u0();
                // U_N[0][1] = VERTEX_0->get_u1();
                // U_N[0][2] = VERTEX_0->get_u2();
                // U_N[0][3] = VERTEX_0->get_u3();

                // U_N[1][0] = VERTEX_1->get_u0();
                // U_N[1][1] = VERTEX_1->get_u1();
                // U_N[1][2] = VERTEX_1->get_u2();
                // U_N[1][3] = VERTEX_1->get_u3();

                // U_N[2][0] = VERTEX_2->get_u0();
                // U_N[2][1] = VERTEX_2->get_u1();
                // U_N[2][2] = VERTEX_2->get_u2();
                // U_N[2][3] = VERTEX_2->get_u3();

                //cout << VERTEX_0->get_x() << endl;

                // VERTEX_0->update_du(DU0);
                // VERTEX_1->update_du(DU1);
                // VERTEX_2->update_du(DU2);

                return ;
        }

        // returns Roe average of left and right states
        double roe_avg(double L1, double L2, double R1, double R2){
                double AVG;
                AVG = (sqrt(L1)*L2+sqrt(R1)*R2)/(sqrt(L1)+sqrt(R1));
                return AVG;
        }

};