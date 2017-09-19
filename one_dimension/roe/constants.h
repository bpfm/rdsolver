#define IC 3

double cfl = 0.1;


#if IC == 0
double t_tot = 1.0;
double GAMMA = 1.4;
#endif

#if IC == 1
double t_tot = 20000.0;
double GAMMA = 1.4;
#endif

#if IC == 2
double t_tot = 1.0;
double GAMMA = 5.0/3.0;
#endif

#if IC == 3
double t_tot = 1.0;
double GAMMA = 1.4;
#endif