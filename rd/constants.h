#define IC 2

#define TWO_D

#define N_SNAP 10

// #define DEBUG
// #define ASTRIX_COMP

#define GENERATE_IC
// #define READ_IC           // doesn't work yet

// #define CLOSED
// #define REFLECTIVE        // doesn't work yet

// #define FIXED_DT

#define OFFSET_GRID
// #define EQUILATERAL_GRID

// #define LDA_SCHEME
#define N_SCHEME

int N_POINTS_X = 64;
int N_POINTS_Y = 64;

double RANDOM_LVL = 0.4;	 // doesn't work for non zero values

// Sod Shock Tube
#if IC == 0
double CFL = 0.01;
double T_TOT = 0.1;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 0.05;
#endif

// Sine Wave Tube
#if IC == 1
double CFL = 0.1;
double T_TOT = 5.0;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 0.1;
#endif

// Sedov Blast Wave
#if IC == 2
double CFL = 0.1;
double T_TOT = 5.0;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 10.0; // if altered, change setup.cpp as well
double SIDE_LENGTH_Y = 10.0;
#endif

// Gaussian pulse advection (x-direction)
#if IC == 3
double CFL = 0.2;
double T_TOT = 1.0;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 0.05;
#endif

// Sod Shock Tube (Varied in Y)
#if IC == 4
double CFL = 0.01;
double T_TOT = 0.1;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 1.0;
#endif

// Uniform flow
#if IC == 5
double CFL = 0.01;
double T_TOT = 1.0;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 30.0;
double SIDE_LENGTH_Y = 30.0;
#endif

// 2D Noh problem
#if IC == 6
double CFL = 0.01;
double T_TOT = 2.0;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 2.0;
double SIDE_LENGTH_Y = 2.0;
#endif

// Gaussian pulse advection (y-direction)
#if IC == 7
double CFL = 0.01;
double T_TOT = 0.1;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 0.05;
double SIDE_LENGTH_Y = 1.0;
#endif

// KH instability (x flow)
#if IC == 8
double CFL = 0.1;
double T_TOT = 1.0;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 1.0;
#endif

// KH instability (y flow)
#if IC == 9
double CFL = 1.0;
double T_TOT = 2.0;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 1.0;
#endif


double GAMMA_1 = GAMMA - 1.0;
double GAMMA_2 = GAMMA - 2.0;

double BLAST_VERTICES = 14.0;
int POINT_CHECK = 0;
