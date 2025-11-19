/* wave1d.c
   Explicit finite-difference solver for 1D wave equation (d'Alembert)
   - options: standing mode, pluck (initial displacement), strike (initial velocity), two pulses
   - boundary: both fixed (Dirichlet) or left-fixed right-free (one fixed)
   - outputs CSV file "output.csv" with columns: t, x, y
   - compile: gcc -O2 -o wave1d wave1d.c -lm
*/

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

/* Utility for safe malloc */
void *xmalloc(size_t n) {
    void *p = malloc(n);
    if (!p) { fprintf(stderr,"Out of memory\n"); exit(1); }
    return p;
}

/* initial displacement functions */
double ic_standing(double x, double L, int mode, double A) {
    return A * sin(mode * M_PI * x / L);
}

double ic_pluck_tri(double x, double L, double x0, double width, double A) {
    // triangular pulse centered at x0 with half-width 'width'
    double d = fabs(x - x0);
    if (d > width) return 0.0;
    return A * (1.0 - d/width);
}

double ic_gaussian(double x, double x0, double sigma, double A) {
    double r = (x-x0)/sigma;
    return A * exp(-0.5 * r * r);
}

/* initial velocity */
double vel_gaussian(double x, double x0, double sigma, double A) {
    return ic_gaussian(x,x0,sigma,A);
}

/* write triple t,x,y to file */
void write_triple(FILE *f, double t, double x, double y) {
    fprintf(f,"%.8g,%.8g,%.8g\n", t, x, y);
}

int main() {
    /* ----- User parameters (you can change or ask input) ----- */
    double L = 1.0;            // string length (m)
    int Nx = 201;              // number of spatial points (odd recommended)
    double T0 = 30.0;          // tension (N)
    double mu0 = 0.1;          // linear mass density (kg/m)
    double v = sqrt(T0 / mu0); // wave speed (m/s)
    double dx = L / (Nx - 1);
    double r = 0.9;            // target CFL v*dt/dx (choose <=1)
    double dt = r * dx / v;
    double tmax = 0.6;         // simulation duration (s)
    int save_every = 5;        // write output every N time steps

    /* Menu choices (set here or read from stdin) */
    printf("1) Standing mode (harmonic)\n");
    printf("2) Pluck: initial displacement (triangular or gaussian)\n");
    printf("3) Strike: initial velocity (localized)\n");
    printf("4) Two opposite pulses\n");
    printf("Choose case (1-4): ");
    int case_choice=0;
    if (scanf("%d", &case_choice)!=1) case_choice = 2;

    printf("Boundary condition: 1=both fixed, 2=left fixed, right free: ");
    int bc_choice=1;
    if (scanf("%d", &bc_choice)!=1) bc_choice = 1;

    /* allocate arrays */
    int N = Nx;
    double *y_prev = (double*) xmalloc(sizeof(double)*N); // y^{n-1}
    double *y_curr = (double*) xmalloc(sizeof(double)*N); // y^{n}
    double *y_next = (double*) xmalloc(sizeof(double)*N); // y^{n+1}
    double *u0 = (double*) xmalloc(sizeof(double)*N);     // initial velocity u(x,0)

    /* initialize to zero */
    for (int i=0;i<N;i++){ y_prev[i]=y_curr[i]=y_next[i]=0.0; u0[i]=0.0; }

    /* Choose initial conditions based on case */
    if (case_choice == 1) {
        int mode = 1;
        double A = 0.01;
        printf("Mode number (n>=1): ");
        if (scanf("%d",&mode)!=1) mode = 1;
        printf("Amplitude (m): ");
        if (scanf("%lf",&A)!=1) A=0.01;
        for (int i=0;i<N;i++) {
            double x = i * dx;
            y_curr[i] = ic_standing(x,L,mode,A);
            u0[i] = 0.0;
        }
    } else if (case_choice == 2) {
        // pluck
        int shape = 1;
        double A = 0.02, x0=0.25, width=0.1, sigma=0.05;
        printf("Pluck shape 1=triangular 2=gaussian: ");
        if (scanf("%d",&shape)!=1) shape=1;
        printf("Amplitude (m): ");
        if (scanf("%lf",&A)!=1) A=0.02;
        printf("Center x0 (m): ");
        if (scanf("%lf",&x0)!=1) x0 = 0.25;
        if (shape==1) {
            printf("Half-width (m): ");
            if (scanf("%lf",&width)!=1) width=0.1;
            for (int i=0;i<N;i++){
                double x=i*dx;
                y_curr[i] = ic_pluck_tri(x,L,x0,width,A);
                u0[i]=0.0;
            }
        } else {
            printf("Gaussian sigma (m): ");
            if (scanf("%lf",&sigma)!=1) sigma=0.05;
            for (int i=0;i<N;i++){
                double x=i*dx;
                y_curr[i] = ic_gaussian(x,x0,sigma,A);
                u0[i]=0.0;
            }
        }
    } else if (case_choice == 3) {
        // strike (initial velocity)
        double A=1.0, x0=0.25, sigma=0.05;
        printf("Strike amplitude (m/s): ");
        if (scanf("%lf",&A)!=1) A=1.0;
        printf("Center x0 (m): ");
        if (scanf("%lf",&x0)!=1) x0=0.25;
        printf("Sigma (m): ");
        if (scanf("%lf",&sigma)!=1) sigma=0.05;
        for (int i=0;i<N;i++){
            double x=i*dx;
            y_curr[i] = 0.0;
            u0[i] = vel_gaussian(x,x0,sigma,A);
        }
    } else {
        // two pulses
        double A=0.02, sigma=0.03;
        double xL = 0.2, xR = 0.8;
        printf("Two pulses amplitude (m): ");
        if (scanf("%lf",&A)!=1) A=0.02;
        printf("sigma (m): ");
        if (scanf("%lf",&sigma)!=1) sigma=0.03;
        for (int i=0;i<N;i++){
            double x=i*dx;
            y_curr[i] = ic_gaussian(x,xL,sigma,A) + ic_gaussian(x,xR,sigma,A);
            u0[i]=0.0;
        }
    }

    /* enforce boundary displacement for initial condition (if fixed ends) */
    if (bc_choice == 1) {
        y_curr[0] = 0.0;
        y_curr[N-1] = 0.0;
    } else {
        // left fixed:
        y_curr[0] = 0.0;
    }

    /* compute dx, dt, and check CFL */
    dx = L/(N-1);
    /* keep dt as chosen above: dt = r * dx / v; recompute r from dt */
    r = v * dt / dx;
    if (r >= 1.0) {
        fprintf(stderr,"CFL condition violated: r = %g >= 1. Reduce dt or increase Nx.\n", r);
        fprintf(stderr,"You can recompile or edit dt/r in the source.\n");
        // exit or clamp? exit to be safe
        exit(1);
    }
    printf("Using N=%d dx=%.6g dt=%.6g v=%.6g CFL r=%.6g\n", N, dx, dt, v, r);

    int Nt = (int) ceil(tmax / dt);
    printf("Number of time steps Nt = %d (tmax=%.6g)\n", Nt, tmax);

    /* Build first time step y_prev -> y_curr -> y_next
       y_prev must represent y^{n-1}. For n=0 (initial), we need y^{-1} implicitly.
       We'll use the standard first-step formula:
       y^1 = y^0 + dt*u0 + 0.5*r^2*(y^0_{i+1}-2y^0_i+y^0_{i-1})
       So set y_prev = y^0 (for convenience), compute y_next as y^1,
       then advance: y_prev <- y_curr; y_curr <- y_next, etc.
    */
    for (int i=0;i<N;i++) y_prev[i] = y_curr[i]; // store y^0 in y_prev

    /* compute first step y_next = y^1 */
    for (int i=1;i<N-1;i++) {
        double lap = y_curr[i+1] - 2.0*y_curr[i] + y_curr[i-1];
        y_next[i] = y_curr[i] + dt * u0[i] + 0.5 * r * r * lap;
    }
    /* boundaries */
    if (bc_choice == 1) {
        y_next[0] = 0.0;
        y_next[N-1] = 0.0;
    } else {
        y_next[0] = 0.0;
        // right free: Neumann zero slope -> y[N-1] = y[N-2]
        y_next[N-1] = y_next[N-2];
    }

    /* open output file */
    FILE *fout = fopen("output.csv","w");
    if (!fout) { perror("output.csv"); exit(1); }
    fprintf(fout,"t,x,y\n");

    /* save initial snapshot (t=0) */
    double t = 0.0;
    for (int i=0;i<N;i++) write_triple(fout, t, i*dx, y_curr[i]);

    /* save y^1 (first step) at t=dt */
    t = dt;
    for (int i=0;i<N;i++) write_triple(fout, t, i*dx, y_next[i]);

    /* advance: now set y_prev = y_curr; y_curr = y_next; */
    for (int i=0;i<N;i++) { y_prev[i]=y_curr[i]; y_curr[i]=y_next[i]; }

    /* time loop */
    for (int n=1; n < Nt; n++) {
        t = n * dt;
        /* interior points */
        for (int i=1;i<N-1;i++) {
            y_next[i] = 2.0*y_curr[i] - y_prev[i] + r*r*(y_curr[i+1] - 2.0*y_curr[i] + y_curr[i-1]);
        }
        /* boundaries */
        if (bc_choice == 1) {
            y_next[0]=0.0; y_next[N-1]=0.0;
        } else {
            y_next[0]=0.0;
            y_next[N-1] = y_next[N-2]; // zero slope
        }
        /* every save_every steps write snapshot at time t+dt */
        if ((n+1) % save_every == 0) {
            double tout = (n+1)*dt;
            for (int i=0;i<N;i++) write_triple(fout, tout, i*dx, y_next[i]);
        }
        /* roll arrays: y_prev <- y_curr; y_curr <- y_next */
        for (int i=0;i<N;i++) { y_prev[i]=y_curr[i]; y_curr[i]=y_next[i]; }
    }

    fclose(fout);
    printf("Simulation done. Output in output.csv (columns t,x,y). Use Octave/Matlab to plot.\n");
    printf("Example (Octave):\n");
    printf(" >> D = csvread('output.csv',1,0); # skip header\n");
    printf(" >> t = unique(D(:,1)); x = unique(D(:,2));\n");
    printf(" >> %% To plot snapshot for a particular time idx:\n");
    printf(" >> tidx = round(length(t)/3); Ym = D(D(:,1)==t(tidx), 3); plot(x, Ym, '-o');\n");
    printf(" >> %% Or animate by looping on times and plotting.\n");

    /* free memory */
    free(y_prev); free(y_curr); free(y_next); free(u0);

    return 0;
}