#define IC 10

//-----------------------------------------
/* set dimensionality */
//-----------------------------------------
#define TWO_D

//-----------------------------------------
/* set umber of snapshots */
//-----------------------------------------
#define N_SNAP 10

//-----------------------------------------
/* debug flag for debug output */
//-----------------------------------------
// #define DEBUG

//-----------------------------------------
/* define flag to read positions and triangles from file (only option currently) */
//-----------------------------------------
#define READ_IC           // doesn't work yet


//-----------------------------------------
/* define flag for using either QHULL or CGAL triangulation */
//-----------------------------------------
// #define QHULL_IC
#define CGAL_IC

//-----------------------------------------
/* define boundary conditions (none for periodic) */
//-----------------------------------------
// #define CLOSED
#define PERIODIC
// #define REFLECTIVE        // doesn't work yet

//-----------------------------------------
/* define flag for fixed timestep */
//-----------------------------------------
// #define FIXED_DT

//-----------------------------------------
/* define type of grid (none for square grid of vertices) */
//-----------------------------------------
// #define OFFSET_GRID
// #define EQUILATERAL_GRID

//-----------------------------------------
/* define distribution scheme */
//-----------------------------------------
#define LDA_SCHEME
// #define N_SCHEME
// #define BLENDED

//-----------------------------------------
/* set order of scheme (none for 2nd order) */
//-----------------------------------------
#define FIRST_ORDER


// Sod Shock Tube (Varied in X)
#if IC == 0
double CFL = 0.1;
double T_TOT = 0.1;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 2.0;
double SIDE_LENGTH_Y = 2.0;
#endif

// Sod Shock Tube (Varied in Y)
#if IC == 1
double CFL = 0.1;
double T_TOT = 0.1;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 2.0;
double SIDE_LENGTH_Y = 2.0;
#endif

// Sine Wave Tube
#if IC == 2
double CFL = 0.1;
double T_TOT = 5.0;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 2.0;
double SIDE_LENGTH_Y = 2.0;
#endif

// Sedov Blast Wave
#if IC == 3
double CFL = 1.0;
double T_TOT = 1.0;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 10.0; // if altered, change setup.cpp as well
double SIDE_LENGTH_Y = 10.0;
#endif

// Gaussian pulse advection (x-direction)
#if IC == 4
double CFL = 0.1;
double T_TOT = 0.1;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 1.0;
#endif

// Gaussian pulse advection (y-direction)
#if IC == 5
double CFL = 0.1;
double T_TOT = 0.1;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 0.05;
double SIDE_LENGTH_Y = 1.0;
#endif

// Uniform flow
#if IC == 6
double CFL = 0.1;
double T_TOT = 1.0;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 2.0;
double SIDE_LENGTH_Y = 2.0;
#endif

// 2D Noh problem
#if IC == 7
double CFL = 0.5;
double T_TOT = 0.1;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 2.0;
double SIDE_LENGTH_Y = 2.0;
#endif

// KH instability (x flow)
#if IC == 8
double CFL = 0.5;
double T_TOT = 2.0;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 1.0;
#endif

// KH instability (y flow)
#if IC == 9
double CFL = 0.5;
double T_TOT = 2.0;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 1.0;
#endif

// KH instability - smoothed (x flow)
#if IC == 10
double CFL = 0.5;
double T_TOT = 1.0;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 1.0;
#endif

// KH instability - smoothed (y flow)
#if IC == 11
double CFL = 0.5;
double T_TOT = 1.0;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 1.0;
#endif

// Blob test !!! NOT WORKING !!!
#if IC == 12
double CFL = 0.01;
double T_TOT = 0.01;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 1.0;
#endif

// Grav Test !!! NOT WORKING !!!
#if IC == 13
double CFL = 0.5;
double T_TOT = 100.0;
double GAMMA = 1.4;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 1.0;
#endif

double GAMMA_1 = GAMMA - 1.0;
double GAMMA_2 = GAMMA - 2.0;

double BLAST_VERTICES = 46.0;
int POINT_CHECK = 1;

double GRAV   = 6.63e-11;
double MSOLAR = 1.989e30;

#ifdef FIXED_DT
double DT_FIX = 0.00001
#endif
