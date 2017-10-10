#define IC 3

double cfl = 0.5;

// 1D Sod Shcok Tube
#if IC == 0
double t_tot = 0.5;
double GAMMA = 1.4;
#endif

// 1D Sine Wave
#if IC == 1
double t_tot = 20000.0;
double GAMMA = 1.4;
#endif

// 1D Sedov Blast Wave
#if IC == 2
double t_tot = 1.0;
double GAMMA = 5.0/3.0;
#endif

// 1D Gaussian Pulse
#if IC == 3
double t_tot = 0.5;
double GAMMA = 1.4;
#endif