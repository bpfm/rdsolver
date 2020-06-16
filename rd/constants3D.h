//-----------------------------------------
/* set dimensionality */
//-----------------------------------------
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
#define LDA_SCHEME
// #define N_SCHEME
// #define BLENDED

//-----------------------------------------
/* set order of scheme (none for 2nd order) */
//-----------------------------------------
#define FIRST_ORDER

double CFL = 0.1;
double T_TOT = 0.1;
double GAMMA = 5.0/3.0;
double SIDE_LENGTH_X = 2.0;
double SIDE_LENGTH_Y = 2.0;
double SIDE_LENGTH_Z = 2.0;

double GAMMA_1 = GAMMA - 1.0;
double GAMMA_2 = GAMMA - 2.0;

double GRAV = 6.67e-11;
double MSOLAR = 1.989e+30;
