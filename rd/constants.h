#define IC 0

#define TWO_D

#define N_SNAP 100

//#define DEBUG

int N_POINTS = 100;

// Sod Shcok Tube
#if IC == 0
double CFL = 0.1;
double T_TOT = 0.01;
double GAMMA = 1.4;
double SIDE_LENGTH = 50.0;
#endif

// Sine Wave Tube
#if IC == 1
double CFL = 0.1;
double T_TOT = 100.0;
double GAMMA = 1.4;
double SIDE_LENGTH = 50.0;
#endif

// Sedov Blast Wave
#if IC == 2
double CFL = 0.1;
double T_TOT = 0.0001;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH = 50.0;
#endif