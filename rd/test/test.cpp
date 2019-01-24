#include "nmarray.h"
#include <iostream>

using namespace std;

int main(){
  nmarray trialarray(4,3);
  int ncount=4,mcount=3;

  cout << "1" << endl;

  for (int i = 0; i < ncount; ++i){
    for (int j = 0; j < mcount; ++j)
    {
      trialarray.set_element(i,j,i+j);
    }
  }



  cout << "2" << endl;

  for (int i = 0; i < ncount; ++i){
    for (int j = 0; j < mcount; ++j)
    {
      std::cout << i << "\t" << j << "\t" << trialarray.get_element(i,j) << std::endl;
    }
  }

  return 0;
}