#define IC 0

#define TWO_D

#define N_SNAP 10

// #define DEBUG
// #define FORCE_ASTRIX
// #define DEBUG_MOMENTUM

#define CLOSED
// #define REFLECTIVE        # doesn't work yet

#define FIXED_DT

// #define LDA_SCHEME
#define N_SCHEME

int N_POINTS_X = 200;
int N_POINTS_Y = 10;


// Sod Shock Tube
#if IC == 0
double CFL = 0.01;
double T_TOT = 0.1;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 1.0;
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
double SIDE_LENGTH_X = 10.0;
double SIDE_LENGTH_Y = 10.0;
#endif

// Gaussian pulse advection
#if IC == 3
double CFL = 0.01;
double T_TOT = 1.0;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 5.0;
double SIDE_LENGTH_Y = 1.0;
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

double GAMMA_1 = GAMMA - 1.0;
double GAMMA_2 = GAMMA - 2.0;
