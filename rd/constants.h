//---------------OA--------------------------
/* choose hydro test */
//-----------------------------------------

// #define SODX
// #define SODY
// #define SINEX
#define SEDOV
// #define GAUSSX
// #define GAUSSY
// #define UNIFORM
// #define NOH
// #define KHX
// #define KHY
// #define KHXSMOOTH
// #define KHYSMOOTH
// #define BLOB
// #define GRAVITY

//-----------------------------------------
/* set dimensionality */
//-----------------------------------------
#define TWO_D
#define THREE_D

//-----------------------------------------
/* set umber of snapshots */
//-----------------------------------------
#define N_SNAP 100

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
// #define LDA_SCHEME
#define N_SCHEME
// #define BLENDED

//-----------------------------------------
/* set order of scheme (none for 2nd order) */
//-----------------------------------------
#define FIRST_ORDER

// #define SELF_GRAVITY // !!! NOT PERIODIC !!!
// #define PARA_RES
// #define PARA_UP


// Sod Shock Tube (Varied in X)
#ifdef SODX
double CFL = 0.1;
double T_TOT = 0.1;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 2.0;
double SIDE_LENGTH_Y = 2.0;
#endif

// Sod Shock Tube (Varied in Y)
#ifdef SODY
double CFL = 0.1;
double T_TOT = 0.1;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 2.0;
double SIDE_LENGTH_Y = 2.0;
#endif

// Sine Wave Tube
#ifdef SINX
double CFL = 0.1;
double T_TOT = 5.0;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 2.0;
double SIDE_LENGTH_Y = 2.0;
#endif

// Sedov Blast Wave
#ifdef SEDOV
double CFL = 0.5;
double T_TOT = 0.01;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 10.0; // if altered, change setup.cpp as well
double SIDE_LENGTH_Y = 10.0;
#endif

double BLAST_E_TOT = 0.0;
double R_BLAST = 0.1;
int POINT_CHECK = 0;

// Gaussian pulse advection (x-direction)
#ifdef GAUSSX
double CFL = 0.1;
double T_TOT = 0.1;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 1.0;
#endif

// Gaussian pulse advection (y-direction)
#ifdef GAUSSY
double CFL = 0.1;
double T_TOT = 0.1;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 1.0;
#endif

// Uniform flow
#ifdef UNIFORM
double CFL = 0.1;
double T_TOT = 1.0;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 2.0;
double SIDE_LENGTH_Y = 2.0;
#endif

// 2D Noh problem
#ifdef NOH
double CFL = 0.4;
double T_TOT = 0.1;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 2.0;
double SIDE_LENGTH_Y = 2.0;
#endif

// KH instability (x flow)
#ifdef KHX
double CFL = 0.4;
double T_TOT = 2.0;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 1.0;
#endif

// KH instability (y flow)
#ifdef KHY
double CFL = 0.4;
double T_TOT = 2.0;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 1.0;
#endif

// KH instability - smoothed (x flow)
#ifdef KHXSMOOTH
double CFL = 0.4;
double T_TOT = 0.5;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 1.0;
#endif

// KH instability - smoothed (y flow)
#ifdef KHYSMOOTH
double CFL = 0.4;
double T_TOT = 2.0;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 1.0;
double SIDE_LENGTH_Y = 1.0;
#endif

// Blob test
#ifdef BLOB
double CFL = 0.4;
double T_TOT = 10.0;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 20.0;
double SIDE_LENGTH_Y = 20.0;
#endif

// Grav Test !!! NOT WORKING !!!
#ifdef GRAVITY
double CFL = 0.4;
double T_TOT = 10.0;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 10.0;
double SIDE_LENGTH_Y = 10.0;
#endif

double GAMMA_1 = GAMMA - 1.0;
double GAMMA_2 = GAMMA - 2.0;

double GRAV = 6.67e-11;
double MSOLAR = 1.989e+30;

int N_TBINS = 8;

#ifdef FIXED_DT
double DT_FIX = 0.00005;
#endif
