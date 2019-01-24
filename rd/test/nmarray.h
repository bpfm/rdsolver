using namespace std;

class NMARRAY{

private:
        double *ARRAY;
public:

        int N,M;

        // Define constructor
        NMARRAY(int i, int j){
                N = i;
                M = j;
                ARRAY = new double[N * M];
                cout << "Finished initializing array" << endl;
        }

        // Set array element [i][j]
        void set_element(int i, int j,double VAL){ARRAY[M*i+j] = VAL;}

        // Get array element [i][j]
        double get_element(int i, int j){return ARRAY[M*i+j];}

        // Define destructor
        ~NMARRAY(){
                delete[] ARRAY;
                cout << "Freeing memory" << endl;
        }

};

