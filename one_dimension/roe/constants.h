#define IC 2

double t_tot = 0.001;

#if IC == 0
double GAMMA = 1.4;
#endif

#if IC == 1
double GAMMA = 1.4;
#endif

#if IC == 2
double GAMMA = 5.0/3.0;
#endif