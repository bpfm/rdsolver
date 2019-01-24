#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main(){
        int n_points = 10;
        double x_side = 10.0;
        double x,y,dx,dy;

        dx = x_side/double(n_points);
        dy = sqrt(3.0)*dx/2.0;

        ofstream positions;

        positions.open("positions.txt");

        for(int i=0;i<n_points;i++){
                for(int j=0;j<n_points;j++){
                        if((j % 2)== 0){
                                x=double(i)*dx;
                                y=double(j)*dy;
                                cout << i << "\t" << j % 2 << "\t" << x << "\t" << y << endl;
                        }else{
                                x=(double(i)+0.5)*dx;
                                y=double(j)*dy;
                                cout << i << "\t" << j % 2 << "\t" << x << "\t" << y << endl;
                        }
                        
                        positions << x << "\t" << y << endl;
                }
        }

        positions.close();

        return 0;
}