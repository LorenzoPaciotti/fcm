#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#define c 2 //numero di centri di cluster
#define n 4 //numero di punti totale in input
double m = 2.0; //fuzzification
#define d 2 //dimensioni spaziali
double epsilon = 0.001; //minima distanza per arrestare

double x[n][d];
double v[c][d];
double u[c][n];

double distanza(double a[d], double b[d]) {
    double somma = 0;
    int i;
    for (i = 0; i < d; i++)
        somma += pow(a[i] - b[i], 2);
    return sqrt(somma);
}

double aggiorna() {
    int i, j, k, l;
    double esp = 2.0 / (m - 1);

    for (i = 0; i < c; i++)
        for (j = 0; j < n; j++) { /* calcolo u[i][j] */
            double somma = 0;
            double num = distanza(x[j], v[i]);
            for (k = 0; k < c; k++) {
                double den = distanza(x[j], v[k]);
                double fraz = num / den;
                somma += pow(fraz, esp);
            }
            u[i][j] = 1.0 / somma;
        }

    double max_diff = 0;
    for (i = 0; i < c; i++) {
        double den = 0;
        for (j = 0; j < n; j++)
            den += pow(u[i][j], m);
        double diff = 0;
        for (k = 0; k < d; k++) { /* calcolo v[i][k] */
            double num = 0;
            for (j = 0; j < n; j++)
                num += pow(u[i][j], m) * x[j][k];
            double old_v = v[i][k];
            v[i][k] = num / den;
            diff += pow(old_v - v[i][k], 2);
        }
        if (diff > max_diff) max_diff = diff;
    }
    return max_diff;
}

void scriviSituazione() {
    int i, j;

    printf("matrice U\n");
    for (i = 0; i < c; i++) {
        for (j = 0; j < n; j++)
            printf("%lf ", u[i][j]);
        printf("\n");
    }
    printf("\n");
    printf("matrice V\n");
    for (i = 0; i < c; i++) {
        for (j = 0; j < d; j++)
            printf("%lf ", v[i][j]);
        printf("\n");
    }
    printf("\n");
    /* scrivere un numero per avanzare */
    scanf("%*d");
}

void test() {
    int i, j;
    double max_diff;

    srand48(3);

    x[0][0] = 1.0;
    x[0][1] = 2;
    x[1][0] = 1.1;
    x[1][1] = 2;
    x[2][0] = 8.0;
    x[2][1] = 4;
    x[3][0] = 7.9;
    x[3][1] = 4.1;

    for (i = 0; i < c; i++)
        for (j = 0; j < d; j++)
            v[i][j] = 10 * drand48() - 5;
    scriviSituazione();

    do {
        max_diff = aggiorna();
        scriviSituazione();
    } while (max_diff > epsilon);
}

int main() {
    test();
}

