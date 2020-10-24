#include <iostream>
#include <math.h>
using namespace std;

int N = 10;
double h = 1.0 / ((double)N);

FILE*outr = fopen("koef.txt", "w");
double F(double k);
double K(double x);
double P(double x);
double Y_true(double x);
void kf1(double *a, double *b, double *c);
void kf2(double *Alpha, double *Beta, double *a, double *b, double *c);



int main() {
    FILE*out = fopen("outprogonka.txt", "w");
    double Alpha[N + 1], Beta[N + 1], result[N + 1], a[N+1], b[N+1], c[N+1];
    double sum = 0.0;
    for (int i = 0; i < N + 1; i++) {
        result[i] = 0.0;
        Alpha[i] = 0.0;
        Beta[i] = 0.0;
        a[i] = 0.0;
        b[i] = 0.0;
        c[i] = 0.0;
    }
    c[0] = 1;
    b[0] = 0;
    //fprintf(outr, "c[%d] = %.5f\tb[%d] = %.5f\n",0, c[0], 0, b[0]);
    kf1(a, b, c);
    
    Alpha[1] = b[0] / c[0];
    Beta[1] = F(0) / (c[0]);
    //Alpha[2] = b[1] / c[1];
    //Beta[2] = F(h) / (c[1]);
    fprintf(outr, "Alfa[%d] = %.5f\tBeta[%d] = %.5f\n",2, Alpha[2], 2, Beta[2]);
    kf2(Alpha, Beta, a, b, c);

    result[N] = (F(N) + a[N] * Beta[N]) / (c[N] - a[N] * Alpha[N]);

    for (int i = N - 1; i > 0; i--) {
        result[i] = Alpha[i + 1] * result[i + 1] + Beta[i + 1];
    }

        sum += 0.5 * h * (result[N] - Y_true(N*h))*(result[N] - Y_true(N*h));
        for (int i = 1; i < N; i++) sum += h * (result[i] - Y_true(i*h))*(result[i] - Y_true(i*h));
        sum += 0.5 * h * (result[0] - Y_true(0*h))*(result[0] - Y_true(0*h));
         printf("||result - y_true||h = %0.10f\n ", sqrt(sum));
        for (int i = 0; i < N+1; i++)
               fprintf(out, "%0.10f\t%0.10f\t%0.10f\t%0.10f\n", i*h, result[i], Y_true(i*h), fabs(Y_true(i*h) - result[i]));

    return 0;
}



double F(double k) {
    if(k == 0)
        return 0.;
     if(k == N) {
         double y = sin(M_PI/2 *k*h);
         double y1 = 0;//M_PI/2 *cos(M_PI/2 *k*h);
    double y2 = -M_PI/2 * M_PI/2 *sin(M_PI/2 *k*h);
         double k1 = 1;//(1+k*h);
         double p = 1;//(1+k*h*k*h);
        return -h*(y2*k1  + y1 - p*y)/(2*K(N*h));}
    else {
        double y= sin(M_PI/2 *k*h);
        double y1 = 0;//M_PI/2 *cos(M_PI/2 *k*h);
        double y2 = -M_PI/2 * M_PI/2 *sin(M_PI/2 *k*h);
        double k1 = 1;//(1+k*h);
        double p = 1;//(1+k*h*k*h);
        return -(y2*k1  + y1 - p*y);
    }

}


double K(double x) {
   //return 1 + x;
   return 1;
}

double P(double x) {
    //return 1 + x*x;
    return 1;
}

double Y_true(double x) {
    return sin(M_PI/2 *x);
}

void kf1(double *a, double *b, double *c) {
    for (int k = 1; k < N; k++) b[k] = K(k*h + 0.5*h) / (h*h);
    for (int k = 1; k < N; k++) a[k] = K(k*h - 0.5*h) / (h*h);
    for (int k = 1; k < N ; k++) c[k] = K(k*h + 0.5*h) / (h*h) + K(k*h - 0.5*h) / (h*h) + P(k*h);
    a[N] = 1/(h);
    c[N] = 1 / (h) + h/2 *P(N*h)/K(N*h);
    fprintf(outr, "c[%d] = %.5f\t b[%d] = %.5f\n", 1, c[1], 1, b[1]);
    for (int j = 2; j < N; j++) {
    fprintf(outr, "a[%d] = %.5f\t c[%d] = %.5f \tb[%d] = %.5f\n", j, a[j], j, c[j], j, b[j]);
    }
    fprintf(outr, "a[%d] = %.5f \tc[%d] = %.5f\n",  N, a[N], N, c[N]);
}

void kf2(double *Alpha, double *Beta, double *a, double *b, double *c) {
    for (int i = 1; i < N; i++) {
        Alpha[i + 1] = b[i] / (c[i] - Alpha[i] * a[i]);
        Beta[i + 1] = (F(i) + a[i] * Beta[i]) / (c[i] - Alpha[i] * a[i]);
    }
    for (int j = 2; j < N; j++) {
        fprintf(outr, "Alfa[%d] = %.5f \tBeta[%d] = %.5f\n", j+1, Alpha[j+1], j+1, Beta[j+1]);
    }
}



