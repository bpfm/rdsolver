#define IC 0

#define TWO_D

#define N_SNAP 10

#define CLOSED

// #define DEBUG

// #define FIXED_DT

#define LDA_SCHEME

// #define N_SCHEME

int N_POINTS_X = 400;
int N_POINTS_Y = 50;

// Sod Shcok Tube
#if IC == 0
double CFL = 0.1;
double T_TOT = 0.1;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 30.0;
double SIDE_LENGTH_Y = 3.0;
#endif

// Sine Wave Tube
#if IC == 1
double CFL = 0.1;
double T_TOT = 500.0;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 50.0;
double SIDE_LENGTH_Y = 10.0;
#endif

// Sedov Blast Wave
#if IC == 2
double CFL = 0.01;
double T_TOT = 0.1;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 50.0;
double SIDE_LENGTH_Y = 10.0;
#endif

// Gaussian pulse advection
#if IC == 3
double CFL = 0.1;
double T_TOT = 10.0;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 50.0;
double SIDE_LENGTH_Y = 1.0;
#endif

// Sod Shcok Tube (Varied in Y)
#if IC == 4
double CFL = 0.1;
double T_TOT = 0.1;
double GAMMA = 1.4;
double SIDE_LENGTH = 30.0;
#endif
