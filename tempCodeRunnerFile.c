/* wave1d.c
   Simulation numérique d’une onde 1D (équation de d’Alembert)
   Méthode : différences finies explicites (schéma classique)
   - Cas possibles : mode propre, corde pincée, choc (vitesse initiale), deux impulsions
   - Conditions limites : deux extrémités fixes ou une fixe + une libre
   - Le programme génère un fichier "output.csv" contenant t, x, y
   - Compilation : gcc -O2 -o wave1d wave1d.c -lm
*/

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

/* petite fonction pour éviter les problèmes d’allocation */
void *xmalloc(size_t n) {
    void *p = malloc(n);
    if (!p) {
        fprintf(stderr, "Erreur : memoire insuffisante.\n");
        exit(1);
    }
    return p;
}

/* ---------------- Fonctions pour les conditions initiales ---------------- */

double ic_mode(double x, double L, int n, double A) {
    return A * sin(n * M_PI * x / L);
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

/* écriture dans le fichier CSV */
void ecrire(FILE *f, double t, double x, double y) {
    fprintf(f, "%.8g,%.8g,%.8g\n", t, x, y);
}

/* ============================= PROGRAMME PRINCIPAL ============================= */

int main() {

    /* ---------------- Paramètres principaux ---------------- */

    double L = 1.0;      // longueur de la corde (m)
    int Nx = 201;        // nombre de points spatiaux
    double T0 = 30.0;    // tension (N)
    double mu = 0.1;     // masse linéique (kg/m)
    double v = sqrt(T0 / mu); // célérité
    double dx = L / (Nx - 1);

    double r = 0.9;      // nombre CFL (<=1)
    double dt = r * dx / v;

    double tmax = 0.6;   // durée de simulation
    int save_every = 5;  // on enregistre un état tous les N pas de temps

    /* ---------------- Menu utilisateur ---------------- */

    printf("=== Simulation d'une onde 1D (Equation de d'Alembert) ===\n\n");

    printf("1) Mode propre (sinus)\n");
    printf("2) Corde pincee (triangulaire ou gaussienne)\n");
    printf("3) Choc : vitesse initiale\n");
    printf("4) Deux impulsions\n");
    printf("\nVotre choix : ");

    int choix = 0;
    if (scanf("%d", &choix) != 1) choix = 2;

    printf("\nConditions limites :\n");
    printf("1) Les deux extremites fixes\n");
    printf("2) Gauche fixe, droite libre\n");
    printf("Votre choix : ");

    int bc = 1;
    if (scanf("%d", &bc) != 1) bc = 1;

    /* ---------------- Allocations mémoire ---------------- */

    int N = Nx;
    double *y_prec = xmalloc(sizeof(double)*N);
    double *y     = xmalloc(sizeof(double)*N);
    double *y_suiv = xmalloc(sizeof(double)*N);
    double *u0    = xmalloc(sizeof(double)*N);

    for (int i = 0; i < N; i++)
        y_prec[i] = y[i] = y_suiv[i] = u0[i] = 0.0;

    /* ---------------- Conditions initiales ---------------- */

    if (choix == 1) {

        /* Mode propre */
        int mode = 1;
        double A = 0.01;

        printf("Numero du mode (n>=1) : ");
        if (scanf("%d", &mode) != 1) mode = 1;

        printf("Amplitude (m) : ");
        if (scanf("%lf", &A) != 1) A = 0.01;

        for (int i=0; i<N; i++) {
            double x = i * dx;
            y[i] = ic_mode(x, L, mode, A);
        }

    } else if (choix == 2) {

        /* Corde pincée */
        int forme = 1;
        double A = 0.02, x0 = 0.25, largeur = 0.1, sigma = 0.05;

        printf("Forme : 1=triangulaire, 2=gaussienne : ");
        if (scanf("%d", &forme) != 1) forme = 1;

        printf("Amplitude (m) : ");
        if (scanf("%lf", &A) != 1) A = 0.02;

        printf("Position du pic x0 (m) : ");
        if (scanf("%lf", &x0) != 1) x0 = 0.25;

        if (forme == 1) {
            printf("Largeur (m) : ");
            if (scanf("%lf", &largeur) != 1) largeur = 0.1;

            for (int i=0; i<N; i++) {
                double x = i * dx;
                y[i] = ic_pince_tri(x, L, x0, largeur, A);
            }

        } else {
            printf("Sigma (m) : ");
            if (scanf("%lf", &sigma) != 1) sigma = 0.05;

            for (int i=0; i<N; i++) {
                double x = i * dx;
                y[i] = ic_gauss(x, x0, sigma, A);
            }
        }

    } else if (choix == 3) {

        /* Choc = vitesse initiale */
        double A = 1.0, x0 = 0.25, sigma = 0.05;

        printf("Amplitude vitesse (m/s) : ");
        if (scanf("%lf", &A) != 1) A = 1.0;

        printf("Position x0 (m) : ");
        if (scanf("%lf", &x0) != 1) x0 = 0.25;

        printf("Sigma (m) : ");
        if (scanf("%lf", &sigma) != 1) sigma = 0.05;

        for (int i=0; i<N; i++) {
            double x = i * dx;
            y[i] = 0.0;
            u0[i] = vitesse_gauss(x, x0, sigma, A);
        }

    } else {

        /* Deux impulsions */
        double A=0.02, sigma=0.03;
        double xL=0.2, xR=0.8;

        printf("Amplitude (m) : ");
        if (scanf("%lf", &A) != 1) A = 0.02;

        printf("Sigma (m) : ");
        if (scanf("%lf", &sigma) != 1) sigma = 0.03;

        for (int i=0; i<N; i++) {
            double x = i * dx;
            y[i] = ic_gauss(x, xL, sigma, A) + ic_gauss(x, xR, sigma, A);
        }
    }

    /* ---------------- Conditions limites au temps t=0 ---------------- */

    if (bc == 1) {
        y[0] = 0.0;
        y[N-1] = 0.0;
    } else {
        y[0] = 0.0;
        /* extrémité droite libre traitée plus loin */
    }

    /* ---------------- Vérification du CFL ---------------- */

    r = v * dt / dx;
    if (r >= 1.0) {
        fprintf(stderr, "Erreur : condition CFL non respectee (r>=1).");
        exit(1);
    }

    int Nt = (int) ceil(tmax / dt);
    printf("\nParametres : N=%d, dt=%.6g, dx=%.6g, r=%.6g, Nt=%d\n",
           N, dt, dx, r, Nt);

    /* ---------------- Premier pas de temps ---------------- */

    for (int i=0; i<N; i++)
        y_prec[i] = y[i];

    for (int i=1; i<N-1; i++) {
        double lap = y[i+1] - 2*y[i] + y[i-1];
        y_suiv[i] = y[i] + dt*u0[i] + 0.5*r*r*lap;
    }

    if (bc == 1) {
        y_suiv[0] = 0.0;
        y_suiv[N-1] = 0.0;
    } else {
        y_suiv[0] = 0.0;
        y_suiv[N-1] = y_suiv[N-2];
    }

    /* ---------------- Ouverture du fichier CSV ---------------- */

    FILE *f = fopen("output.csv", "w");
    if (!f) {
        perror("Impossible d'ouvrir output.csv");
        exit(1);
    }
    fprintf(f, "t,x,y\n");

    double t = 0.0;
    for (int i=0; i<N; i++)
        ecrire(f, t, i*dx, y[i]);

    t = dt;
    for (int i=0; i<N; i++)
        ecrire(f, t, i*dx, y_suiv[i]);

    /* Mise à jour : y_prec ← y, y ← y_suiv */
    for (int i=0; i<N; i++) {
        y_prec[i] = y[i];
        y[i] = y_suiv[i];
    }

    /* ---------------- Boucle en temps ---------------- */

    for (int n=1; n<Nt; n++) {

        for (int i=1; i<N-1; i++) {
            y_suiv[i] = 2*y[i] - y_prec[i] + r*r*(y[i+1] - 2*y[i] + y[i-1]);
        }

        if (bc == 1) {
            y_suiv[0] = 0.0;
            y_suiv[N-1] = 0.0;
        } else {
            y_suiv[0] = 0.0;
            y_suiv[N-1] = y_suiv[N-2];
        }

        if ((n+1) % save_every == 0) {
            double t_enr = (n+1)*dt;
            for (int i=0; i<N; i++)
                ecrire(f, t_enr, i*dx, y_suiv[i]);
        }

        for (int i=0; i<N; i++) {
            y_prec[i] = y[i];
            y[i] = y_suiv[i];
        }
    }

    fclose(f);

    printf("\nSimulation terminee. Resultats -> output.csv\n");
    printf("OCTAVE Pour Visualiser\n\n");

    free(y_prec);
    free(y);
    free(y_suiv);
    free(u0);

    return 0;
}
