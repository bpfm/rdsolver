#define IC 3

//#define ONE_D
#define THREE_D

//#define BOUNDARY "PERIODIC"
#define BOUNDARY "CLOSED"

#define FLUX_LIMITER 0
#define FLUX_LIMITER_TYPE "MINMOD"

#define FORCE_ANALYTIC_PULSE 1

#define N_SNAP 5

int N_POINTS=50;

// Sod Shcok Tube
#if IC == 0
double cfl = 0.1;
double t_tot = 0.5;
double GAMMA = 1.4;
double SIDE_LENGTH = 50.0;
#endif

// Sine Wave
#if IC == 1
double cfl = 0.1;
double t_tot = 100.0;
double GAMMA = 1.4;
double SIDE_LENGTH = 50.0;
#endif

// Sedov Blast Wave
#if IC == 2
double cfl = 0.1;
double t_tot = 0.2;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH = 50.0;
#endif

// Gaussian Pulse
#if IC == 3
double cfl = 0.1;
double t_tot = 1.0;
double GAMMA = 1.4;
double SIDE_LENGTH = 50.0;
#endif

// Square wave
#if IC == 4
double cfl = 0.1;
double t_tot = 5.0;
double GAMMA = 1.4;
double SIDE_LENGTH = 50.0;
#endif

// Steady flow
#if IC == 5
double cfl = 0.1;
double t_tot = 0.1;
double GAMMA = 1.4;
double SIDE_LENGTH = 50.0;
#endif

// Linear Waves Test (ATHENA)
#if IC == 6
double cfl = 0.1;
double t_tot = 0.1;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH = 1.0;
#endif

// Linear Waves Test (ATHENA)
#if IC == 7
double cfl = 0.1;
double t_tot = 0.038;
double GAMMA = 1.4;
double SIDE_LENGTH = 1.0;
#endif
