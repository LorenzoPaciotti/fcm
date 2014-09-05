#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#define c 2 //numero di centri di cluster
#define n 4 //numero di punti totale in input
double m = 2.0; //fuzzification
#define d 2 //dimensioni spaziali
double epsilon = 0.001; //minima distanza per arrestare
double miaDist[2]; //vettore delle distanze fra centroidi


double *X[n]; //dati input  (d x n)
double *U[n]; //partition matrix   (c x n)
double *V[c]; //matr centroidi  (d x c)

void stampaMatrice(int righe, int col, double *mat[]) {
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

double normaVett(double x[]) {
    double somma = 0.0;
    int i;
    for (i = 0; i < d; i++) {
        somma += pow(x[i], 2);
    }
    return sqrt(somma);
}

void diffVett(double x[], double y[], double ris[]) {
    int i;
    for (i = 0; i < d; i++) {
        ris[i] = x[i] - y[i];
    }
}

void stampaVett(double x[]) {
    puts("____");
    printf("|%lf|\n", x[0]);
    printf("|%lf|\n", x[1]);
    puts("____");
}

double calcDistanza2D(double x[], double y[]) {
    //d(x,y) = norm(x-y)
    double z[d];
    diffVett(x, y, z);
    return normaVett(z);
}

int maxSuVett(double vet[]) {

    return 0; //ritornare l'indice con il massimo

}

void prodottoScalareVettore(double scal, double vett_in[], double vett_out[]) {

    int i;
    for (i = 0; i < d; i++) {
        vett_out[i] = vett_in[i] * scal;
    }
}

void ricalcCentroide(int i_ext, double *U[], double *X[], double *V[]) {

    int j;

    for (j = 0; j < n; j++) {
        prodottoScalareVettore(pow(U[i_ext][j], m), X[j], V[i_ext]);
        prodottoScalareVettore((1 / (pow(U[i_ext][j], m))), V[i_ext], V[i_ext]); //
    }
}

double calc_u_ij(double x_j[], double v_i[]) {
    double denom;
    double dist_x_j__v_i = calcDistanza2D(x_j, v_i);

    int k;
    for (k = 0; k < c; k++) {
        denom += (dist_x_j__v_i / calcDistanza2D(x_j, V[k]));
    }
    denom = pow(denom, (2.0 / (m - 1.0)));

    return (1.0 / denom);
}

int main(int argc, char** argv) {
    /*
     * INPUT: X, c, m
     * OUTPUT: U,V
     */
    int i;
    int j;
    //alloc V ha dimensione dxc
    for (i = 0; i < c; i++) {
        V[i] = malloc(sizeof (double)*d);
    }
    //init V
    V[0][0] = 2.0;
    V[0][1] = 2.0;
    V[1][0] = 0.0;
    V[1][1] = 0.0;

    //alloc U ha dimensione cxn
    for (i = 0; i < n; i++) {
        U[i] = malloc(sizeof (double)*c);
    }

    //alloc X ha dimensione dxn
    for (i = 0; i < n; i++) {
        X[i] = malloc(sizeof (double)*d);
    }
    //init X
    X[0][0] = 0.1;
    X[0][1] = 0.9;
    X[1][0] = 4.5;
    X[1][1] = 7.4;
    X[2][0] = 0.3;
    X[2][1] = 200;
    X[3][0] = 2000;
    X[3][1] = 1111;




    while (1) {
        for (i = 0; i < c; i++) {
            for (j = 0; j < n; j++) {
                U[i][j] = calc_u_ij(X[j], V[i]);
            }
        }


        for (i = 0; i < c; i++) {
            ricalcCentroide(i, U, X, V);
        }

        puts("matrice X:");
        stampaMatrice(n, d, X);
        puts("");
        puts("matrice U:");
        stampaMatrice(n, c, U);
        puts("");
        puts("matrice V:");
        stampaMatrice(c, d, V);
        printf("\n\n\n\n\n");

        sleep(2);
    }

    return (0);
}
