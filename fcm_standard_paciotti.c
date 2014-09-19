#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#define c 2 //numero di centri di cluster
#define n 4 //numero di punti totale in input
double m = 2.0; //fuzzification
#define d 2 //dimensioni spaziali
double epsilon = 0.001; //minima distanza per arrestare
double distanze[c]; //vettore con le dist fra centroidi dopo l'aggiornamento


double X[n][d]; //dati input
double U[c][n]; //partition matrix
double V[c][d]; //matr centroidi
double max;

void stampaMatrice(int righe, int col, double mat[][col]) {
    int i, j;
    printf("\n\n");
    for (i = 0; i < righe; i++) {
        for (j = 0; j < col; j++) {
            printf("%lf", mat[i][j]);
            printf(" ");
        }
        puts("");
    }
}

void stampaVett(double x[]) {
    puts("____");
    printf("|%lf|\n", x[0]);
    printf("|%lf|\n", x[1]);
    puts("____");
}

double calcDistanza(double a[d], double b[d]) {
    double ris = 0;
    int i;
    for (i = 0; i < d; i++)
        ris += pow(a[i] - b[i], 2);
    return sqrt(ris);
}

void prodottoScalareVettore(double scal, double vett_in[], double vett_out[]) {

    int i;
    for (i = 0; i < d; i++) {
        vett_out[i] = vett_in[i] * scal;
    }
}

double maxDistCentroidi() {
    double max = 0.0;
    int i;
    for (i = 0; i < c; i++) {
        if (distanze[i] > max)
            max = distanze[i];
    }
    return max;
}

int main(int argc, char** argv) {
    /*
     * INPUT: X, c, m
     * OUTPUT: U,V
     */
    int i;
    int j;

    srand48(3);

    X[0][0] = 1.0;
    X[0][1] = 2;
    X[1][0] = 1.1;
    X[1][1] = 2;
    X[2][0] = 8.0;
    X[2][1] = 4;
    X[3][0] = 7.9;
    X[3][1] = 4.1;

    for (i = 0; i < c; i++)
        for (j = 0; j < d; j++)
            V[i][j] = 10 * drand48() - 5;

    puts("\ninizializzazione matrice V:");
    stampaMatrice(c, d, V);
    puts("fine init V");

    int contPassi = 0;
    max = 0.0;
    do {
        //CALCOLO PARTITION MATRIX
        printf("PASSO: %d\n\n", ++contPassi);
        for (i = 0; i < c; i++) {
            for (j = 0; j < n; j++) {
                double esponente = 2.0 / (m - 1.0);
                double denom = 0.0;
                double dist_x_j__v_i = calcDistanza(X[j], V[i]);

                int k; //SOMMATORIA 1
                for (k = 0; k < c; k++) {
                    double dist_xj_vk = calcDistanza(X[j], V[k]);
                    //if (dist_xj_vk <= 0) dist_xj_vk = DBL_MIN * 2;//NaN
                    denom += pow((dist_x_j__v_i / dist_xj_vk), esponente);
                }
                U[i][j] = 1.0 / denom;
            }
        }


        //RICALCOLO POSIZIONE CENTROIDI
        int i, j, z, k;
        double old[d];
        double denom;
        for (i = 0; i < c; i++) {
            for (z = 0; z < d; z++)
                old[z] = V[i][z]; //per confronto diff
            double denom = 0.0;
            for(j=0; j<n;j++)//sommatoria denom
                denom += pow(U[i][j], m);
            for (k = 0; k < d; k++) {
                double num = 0.0;
                for (j = 0; j < n; j++) {//SOMMATORIA 2 e 3
                    num += X[j][k] * pow(U[i][j], m);
                }
                V[i][k] = num / denom;
            }
            distanze[i] = pow(calcDistanza(V[i], old), 2.0);
        }



        puts("matrice X:");
        stampaMatrice(n, d, X);
        puts("");
        puts("matrice U:");
        stampaMatrice(c, n, U);
        puts("");
        puts("matrice V:");
        stampaMatrice(c, d, V);
        puts("");
        stampaVett(distanze);
        puts("");
        printf("MAX: %lf", maxDistCentroidi());
        printf("\n\n\n\n\n");

        sleep(2);
    } while (maxDistCentroidi() > epsilon); //
    return (0);
}