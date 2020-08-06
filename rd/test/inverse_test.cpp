#include <iostream>

#include "cblas.h"
#include "lapacke.h"
#include "inverse.cpp"

int main(){
  double INFLOW_MINUS_SUM[5][5];

  INFLOW_MINUS_SUM[0][0] = 1;
  INFLOW_MINUS_SUM[0][1] = 2;
  INFLOW_MINUS_SUM[0][2] = 0;
  INFLOW_MINUS_SUM[0][3] = 6;
  INFLOW_MINUS_SUM[0][4] = 3;
  
  INFLOW_MINUS_SUM[1][0] = 5;
  INFLOW_MINUS_SUM[1][1] = 0;
  INFLOW_MINUS_SUM[1][2] = 2;
  INFLOW_MINUS_SUM[1][3] = 6;
  INFLOW_MINUS_SUM[1][4] = 9;

  INFLOW_MINUS_SUM[2][0] = 8;
  INFLOW_MINUS_SUM[2][1] = 3;
  INFLOW_MINUS_SUM[2][2] = 5;
  INFLOW_MINUS_SUM[2][3] = 0;
  INFLOW_MINUS_SUM[2][4] = 1;

  INFLOW_MINUS_SUM[3][0] = 0;
  INFLOW_MINUS_SUM[3][1] = 8;
  INFLOW_MINUS_SUM[3][2] = 2;
  INFLOW_MINUS_SUM[3][3] = 4;
  INFLOW_MINUS_SUM[3][4] = 3;

  INFLOW_MINUS_SUM[4][0] = 9;
  INFLOW_MINUS_SUM[4][1] = 2;
  INFLOW_MINUS_SUM[4][2] = 1;
  INFLOW_MINUS_SUM[4][3] = 0;
  INFLOW_MINUS_SUM[4][4] = 7;

  mat_inv(&INFLOW_MINUS_SUM[0][0],5,0,0,0);
}
