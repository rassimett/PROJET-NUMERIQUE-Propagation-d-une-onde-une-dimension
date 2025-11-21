/* wave1d_project.c
   Projet : Propagation d'une onde 1D - Paul Sabatier
   - Menu complet conforme au sujet (I, II, III)
   - Génère des fichiers CSV pour visualisation Octave
   - Compilation : gcc -O2 -o wave1d_project wave1d_project.c -lm
   NOTE: spécification fournie dans le PDF : /mnt/data/B_Projet_Propagation_Onde_1D.pdf
*/

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

/* utilitaire malloc sûr */
void *xmalloc(size_t n) {
    void *p = malloc(n);
    if (!p) { fprintf(stderr, "Erreur : memoire insuffisante.\n"); exit(1); }
    return p;
}

/* ---------- conditions initiales ---------- */

/* mode propre : selon BC
   bc==1 -> fixe-fixe : sin(n*pi*x/L)
   bc==2 -> fixe-libre : sin((2*n-1)*pi*x/(2L))
*/
double mode_propre(double x, double L, int n, double A, int bc) {
    if (bc == 1) return A * sin(n * M_PI * x / L);
    else return A * sin((2 * n - 1) * M_PI * x / (2.0 * L));
}

double ic_pince_tri(double x, double L, double x0, double largeur, double A) {
    double d = fabs(x - x0);
    if (d > largeur) return 0.0;
    return A * (1.0 - d / largeur);
}

double ic_gauss(double x, double x0, double sigma, double A) {
    double r = (x - x0) / sigma;
    return A * exp(-0.5 * r * r);
}

double vitesse_gauss(double x, double x0, double sigma, double A) {
    return ic_gauss(x, x0, sigma, A);
}

/* écriture CSV t,x,y */
void ecrire(FILE *f, double t, double x, double y) {
    fprintf(f, "%.8g,%.8g,%.8g\n", t, x, y);
}

/* ---------- noyau de simulation (schéma explicite) ----------
   Paramètres d'entrée :
     L, Nx, T0, mu, dt, tmax, save_every, bc
     init_choice : 1=mode propre (only used for I.1), 2=pluck, 3=strike, 4=two pulses
   Pour les cas spéciaux on passe des paramètres supplémentaires via structure simple
*/
typedef struct {
    int mode;      // pour mode propre
    double A;      // amplitude (déplacement ou vitesse)
    double x0;     // position
    double width;  // largeur (triangulaire)
    double sigma;  // sigma pour gauss
    int shape;     // 1=triangulaire,2=gaussienne
} InitParams;

void simulate_case(const char *outname,
                   double L, int Nx, double T0, double mu,
                   double dt, double tmax, int save_every,
                   int bc, int init_choice, InitParams p)
{
    int N = Nx;
    double dx = L / (N - 1);
    double v = sqrt(T0 / mu);
    double r = v * dt / dx;
    if (r >= 1.0) {
        fprintf(stderr, "CFL non respecte r=%g >=1 pour %s\n", r, outname);
        return;
    }

    int Nt = (int) ceil(tmax / dt);

    double *y_prev = xmalloc(sizeof(double)*N);
    double *y = xmalloc(sizeof(double)*N);
    double *y_next = xmalloc(sizeof(double)*N);
    double *u0 = xmalloc(sizeof(double)*N);

    for (int i=0;i<N;i++){ y_prev[i]=y[i]=y_next[i]=u0[i]=0.0; }

    /* initial conditions */
    if (init_choice == 1) {
        // mode propre
        for (int i=0;i<N;i++){
            double x = i*dx;
            y[i] = mode_propre(x, L, p.mode, p.A, bc);
            u0[i] = 0.0;
        }
    } else if (init_choice == 2) {
        // pluck (deformation initiale, vitesse nulle)
        for (int i=0;i<N;i++){
            double x = i*dx;
            if (p.shape == 1) y[i] = ic_pince_tri(x, L, p.x0, p.width, p.A);
            else y[i] = ic_gauss(x, p.x0, p.sigma, p.A);
            u0[i] = 0.0;
        }
    } else if (init_choice == 3) {
        // strike (deplacement initial nul, vitesse non nulle)
        for (int i=0;i<N;i++){
            double x = i*dx;
            y[i] = 0.0;
            u0[i] = vitesse_gauss(x, p.x0, p.sigma, p.A);
        }
    } else if (init_choice == 4) {
        // two pulses (two gaussians)
        double xL = L*0.2, xR = L*0.8;
        for (int i=0;i<N;i++){
            double x = i*dx;
            y[i] = ic_gauss(x, xL, p.sigma, p.A) + ic_gauss(x, xR, p.sigma, p.A);
            u0[i] = 0.0;
        }
    }

    /* apply BC at t=0 */
    y[0]=0.0;
    if (bc==1) y[N-1]=0.0;

    /* prepare first time step */
    for (int i=0;i<N;i++) y_prev[i] = y[i];
    for (int i=1;i<N-1;i++){
        double lap = y[i+1] - 2.0*y[i] + y[i-1];
        y_next[i] = y[i] + dt*u0[i] + 0.5 * r * r * lap;
    }
    if (bc==1) {
        y_next[0]=0.0; y_next[N-1]=0.0;
    } else {
        y_next[0]=0.0;
        y_next[N-1] = y_next[N-2]; // Neumann approx for free end
    }

    /* open output file */
    FILE *f = fopen(outname, "w");
    if (!f) { perror("fopen"); goto cleanup; }
    fprintf(f,"t,x,y\n");

    double t = 0.0;
    for (int i=0;i<N;i++) ecrire(f, t, i*dx, y[i]);
    t = dt;
    for (int i=0;i<N;i++) ecrire(f, t, i*dx, y_next[i]);

    /* roll */
    for (int i=0;i<N;i++){ y_prev[i]=y[i]; y[i]=y_next[i]; }

    /* main time loop */
    for (int n=1; n<Nt; n++){
        double tout = (n+1)*dt;
        for (int i=1;i<N-1;i++){
            y_next[i] = 2.0*y[i] - y_prev[i] + r*r*(y[i+1] - 2.0*y[i] + y[i-1]);
        }
        if (bc==1) { y_next[0]=0.0; y_next[N-1]=0.0; }
        else { y_next[0]=0.0; y_next[N-1] = y_next[N-2]; }

        if ((n+1) % save_every == 0) {
            for (int i=0;i<N;i++) ecrire(f, tout, i*dx, y_next[i]);
        }

        for (int i=0;i<N;i++){ y_prev[i]=y[i]; y[i]=y_next[i]; }
    }

    fclose(f);

cleanup:
    free(y_prev); free(y); free(y_next); free(u0);
    return;
}

/* ---------- menus helpers ---------- */

int ask_int(const char *prompt, int def) {
    int v;
    printf("%s", prompt);
    if (scanf("%d", &v) != 1) { v = def; }
    return v;
}
double ask_double(const char *prompt, double def) {
    double v;
    printf("%s", prompt);
    if (scanf("%lf", &v) != 1) { v = def; }
    return v;
}

/* ---------- programme principal : menu complet ---------- */

int main() {
    printf("=== PROJET : Propagation d'une onde 1D (Equation de d'Alembert) ===\n");
    printf("Sujet (pdf) : /mnt/data/B_Projet_Propagation_Onde_1D.pdf\n\n");

    /* paramètres par défaut phys. et numériques */
    double L = 1.0;
    int Nx = 201;
    double T0 = 30.0;
    double mu = 0.1;
    double v = sqrt(T0 / mu);
    double dx = L / (Nx - 1);
    double default_r = 0.9;
    double dt = default_r * dx / v;
    double tmax = 2.0;
    int save_every = 1;

    while (1) {
        printf("\nMenu principal :\n");
        printf("1) Corde fixee aux deux extremites (I)\n");
        printf("2) Corde fixee a gauche, libre a droite (II)\n");
        printf("3) Influence du pas en temps (III)\n");
        printf("4) Modifier parametres (L, Nx, T0, mu, tmax)\n");
        printf("0) Quitter\n");
        int mainc = ask_int("Votre choix : ", 0);

        if (mainc == 0) { printf("Au revoir.\n"); break; }

        if (mainc == 4) {
            L = ask_double("L (m) ? [1.0] : ", L);
            Nx = ask_int("Nx (points) ? [201] : ", Nx);
            T0 = ask_double("Tension T0 (N) ? [30] : ", T0);
            mu = ask_double("masse lineique mu (kg/m) ? [0.1] : ", mu);
            tmax = ask_double("tmax (s) ? [2.0] : ", tmax);
            dx = L / (Nx - 1);
            v = sqrt(T0 / mu);
            dt = default_r * dx / v;
            printf("Parametres mis a jour. dt calcule = %g (r=%g)\n", dt, default_r);
            continue;
        }

        if (mainc == 1 || mainc == 2) {
            int bc = (mainc==1) ? 1 : 2;
            printf("\n--- Vous avez choisi : %s ---\n",
                   (bc==1) ? "Corde fixe-fixe (I)" : "Corde fixe-libre (II)");
            printf("1) Ondes stationnaires (I.1) [visualiser modes propres]\n");
            printf("2) Propagation (I.2 / II) [deformation initiale]\n");
            int sub = ask_int("Choix : ", 1);

            if (sub == 1) {
                // modes propres
                int mode = ask_int("Mode n (>=1) ? [1]: ", 1);
                double A = ask_double("Amplitude (m) ? [0.01]: ", 0.01);
                InitParams p = {mode, A, 0.0, 0.1, 0.05, 1};
                char fname[256];
                snprintf(fname, sizeof(fname), "output_mode_bc%d_n%d.csv", bc, mode);
                double dx_local = L/(Nx-1);
                double dt_local = default_r * dx_local / sqrt(T0/mu);
                simulate_case(fname, L, Nx, T0, mu, dt_local, tmax, save_every, bc, 1, p);
                printf("Fichier genere : %s (mode %d, BC=%d)\n", fname, mode, bc);
                printf("-> Visualiser avec visualize_wave.m en ajustant amplitude.\n");
            } else if (sub == 2) {
                printf("Propagation : choix du type initial :\n");
                printf("1) Deformation initiale, vitesse nulle: pincement/impulsion/train/2 impulsions\n");
                printf("2) Deformation nulle, vitesse initiale (choc)\n");
                int typ = ask_int("Choix : ", 1);
                if (typ == 1) {
                    printf("a) Pincement triangulaire\nb) Pincement gaussien\nc) Impulsion gaussienne\n");
                    printf("d) Deux impulsions symetriques\n");
                    char choice[8];
                    printf("Choix (a/b/c/d) [a]: ");
                    if (scanf("%s", choice) != 1) strcpy(choice,"a");
                    InitParams p; memset(&p,0,sizeof(p));
                    p.A = ask_double("Amplitude (m) ? [0.02]: ", 0.02);
                    p.x0 = ask_double("Position x0 (m) ? [0.25]: ", 0.25);
                    p.width = ask_double("Largeur (pour triangulaire) ? [0.1]: ", 0.1);
                    p.sigma = ask_double("sigma (pour gaussienne) ? [0.05]: ", 0.05);
                    p.shape = (choice[0]=='b') ? 2 : 1;
                    char fname[256];
                    if (choice[0]=='d') {
                        p.sigma = ask_double("sigma pour les deux impulsions ? [0.03]: ", 0.03);
                        snprintf(fname, sizeof(fname),"output_two_pulses_bc%d.csv", bc);
                        simulate_case(fname, L, Nx, T0, mu, dt, tmax, save_every, bc, 4, p);
                    } else if (choice[0]=='c') {
                        snprintf(fname, sizeof(fname),"output_pulse_bc%d.csv", bc);
                        // use single gaussian pulse at x0
                        simulate_case(fname, L, Nx, T0, mu, dt, tmax, save_every, bc, 2, p);
                    } else {
                        // a or b = pinch
                        snprintf(fname, sizeof(fname),"output_pluck_bc%d.csv", bc);
                        simulate_case(fname, L, Nx, T0, mu, dt, tmax, save_every, bc, 2, p);
                    }
                    printf("Fichier genere. Visualisez avec Octave.\n");
                } else {
                    // strike
                    InitParams p; memset(&p,0,sizeof(p));
                    p.A = ask_double("Amplitude vitesse (m/s) ? [1.0]: ", 1.0);
                    p.x0 = ask_double("Position x0 (m) ? [0.25]: ", 0.25);
                    p.sigma = ask_double("sigma (m) ? [0.05]: ", 0.05);
                    char fname[256];
                    snprintf(fname, sizeof(fname),"output_strike_bc%d.csv", bc);
                    simulate_case(fname, L, Nx, T0, mu, dt, tmax, save_every, bc, 3, p);
                    printf("Fichier genere : %s\n", fname);
                }
            } else {
                printf("Choix inconnu.\n");
            }
        } else if (mainc == 3) {
            /* Influence du pas en temps : on fait varier r et on enregistre */
            printf("\n--- Etude influence du pas en temps (III) ---\n");
            double rmin = 0.3, rmax = 0.99;
            int Nr = ask_int("Nombre de valeurs r a tester ? [4]: ", 4);
            rmin = ask_double("r min (>=0) ? [0.3]: ", 0.3);
            rmax = ask_double("r max (<1) ? [0.99]: ", 0.99);
            if (rmin < 0) rmin = 0.3;
            if (rmax >= 1.0) rmax = 0.99;
            printf("Quel cas utiliser pour l'etude ?: 1=plucker (deformation), 2=strike (vitesse) [1]: ");
            int cas = ask_int("", 1);
            InitParams p; memset(&p,0,sizeof(p));
            if (cas==1) {
                p.A = ask_double("Amplitude deformation (m) ? [0.02]: ", 0.02);
                p.x0 = ask_double("x0 (m) ? [0.25]: ", 0.25);
                p.sigma = ask_double("sigma (m) ? [0.05]: ", 0.05);
            } else {
                p.A = ask_double("Amplitude vitesse (m/s) ? [1.0]: ", 1.0);
                p.x0 = ask_double("x0 (m) ? [0.25]: ", 0.25);
                p.sigma = ask_double("sigma (m) ? [0.05]: ", 0.05);
            }
            for (int k=0;k<Nr;k++){
                double rcur = rmin + (rmax-rmin) * k / (double) (Nr-1);
                double dx_local = L/(Nx-1);
                double dt_local = rcur * dx_local / sqrt(T0/mu);
                char fname[256];
                snprintf(fname, sizeof(fname), "output_r_%g_case%d.csv", rcur, cas);
                int init_choice = (cas==1)?2:3;
                simulate_case(fname, L, Nx, T0, mu, dt_local, tmax, save_every, 1, init_choice, p);
                printf("generé r=%g -> %s\n", rcur, fname);
            }
            printf("Etude terminee. Comparez fichiers (par ex. amplitude, stability).\n");
        } else {
            printf("Choix invalide.\n");
        }
    }

    return 0;
}
