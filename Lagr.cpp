#include <iostream>
#include <math.h>
#include <locale.h>
using namespace std;
double func(double x);
double P(int n, double *p, double x);
void in_matr(int n, double *matr, double *x);
void print_matr(FILE*fmatr, int n, double *matr);
void swap(int n, int i, int k, double *a);
double gauss(int n, double *a, double *b, double *x);
double F(int i, int n, double y, double *x);
double L(int n, double y, double *x);
//void print_2(FILE*f, int n, int m, double *x, double *koef);
void print_3(FILE*f, int n, int m, double *x, double *koef);

int main() {
 //part 1
    FILE*fout = fopen("out1.txt", "wt");
    double a, b;
    int n;
    cin >> a >> b >> n;
    double *x = new double[n];
    double *f = new double[n];
    fprintf(fout, "%d\n", n);
    double k = (a+b)/2., p = (b-a)/2.;
    
    for(int i = 1; i < n+1; i++) {
        //x[i-1] = a + (b-a)*(i-1)/(n-1);
        x[i-1] = k + p*cos((2.*M_PI*i - 1.)/(2.*n));
        f[i-1] = func(x[i-1]);
        fprintf(fout, "%0.10f\t%0.10f\n", x[i-1], f[i-1]);
    }
    fclose(fout);
    
//part 2
    FILE*fmatr = fopen("matr.txt", "wt");
    //FILE*fout2 = fopen("out2.txt", "wt");
    double *matr = new double[n*n];
    double *koef = new double[n];
    in_matr(n, matr, x);
    print_matr(fmatr, n, matr);
    
    gauss(n, matr, f, koef);
    //проверка Гаусса и вывод коэфф тут в первой строке
    /*for(int i = 0; i < n; i++) printf("%0.10f\n", koef[i]);
        for(int i = 0; i < n; i++)
           for(int j = 0; j < n; j++)
               matr[i*n + j] = 0.;
    in_matr(n, matr, x);
    double tmp = 0.;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++)
            tmp += matr[i*n + j]*koef[j];
        printf("%0.10f   %0.10f\n", tmp, func(x[i]));
        tmp = 0.;
    }*/
    fclose(fmatr);
    delete []matr;
    delete []f;
    int m = 3; //разбиение маленького отрезка
    //print_2(fout2, n, m, x, koef);
    //fclose(fout2);
    
//part 3
    FILE*fout3 = fopen("out3.txt", "wt");
    print_3(fout3, n, m, x, koef);
    fclose(fout3);
    delete []x;
    delete []koef;
    return 0;
}

double func(double x) {
    return 1./(25*x*x + 1.); //1. + pow(x, 5);
}

double P(int n, double *p, double x) {
    double tmp = p[n-1];
    for(int i = 1; i < n; i++)
        tmp = p[n-1-i] + x*tmp;
    return tmp;
}

void in_matr(int n, double *matr, double *x) {
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            matr[n*i + j] = pow(x[i], j);
}
void print_matr(FILE*fmatr, int n, double *matr){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++)
               fprintf(fmatr, "%0.10f\t", matr[n*i + j]);
        fprintf(fmatr, "\n");
    }
}

void swap(int n, int i, int k, double *a) {
    double tmp;
    for (int j = 0; j < n; ++j) {
        tmp = a[i * n + j];
        a[i * n + j] = a[k * n + j];
        a[k * n + j] = tmp;
    }
}

double gauss(int n, double *a, double *b, double *x) {
    int i, j, k, indMax;
    double tmp, max;

    for (i = 0; i < n; ++i) {
        max = fabs(a[i * n + i]);
        indMax = i;
        for (j = i + 1; j < n; ++j)
            if (max < fabs(a[j * n + i])){
                max = fabs(a[j * n + i]);
                indMax = j;
            }
        if(indMax != i) {
            swap(n, i, indMax, a);
            tmp = b[i];
            b[i] = b[indMax];
            b[indMax] = tmp;
        }
        if (fabs(a[i * n + i]) < 1e-15)
            return -1;
        tmp = 1.0 / a[i * n + i];
        for (j = i; j < n; ++j)
            a[i * n + j] *= tmp;
        b[i] *= tmp;

        for (j = i + 1; j < n; ++j) {
            tmp = a[j * n + i];
            for (k = i; k < n; ++k)
                a[j * n + k] -= a[i * n + k] * tmp;
            b[j] -= b[i] * tmp;
        }
    }

    for (i = n - 1; i >= 0; --i) {
        tmp = b[i];
        for (j = i + 1; j < n; ++j)
            tmp -= a[i * n + j] * x[j];
        x[i] = tmp;
    }
    return 0;
}

double F(int i, int n, double y, double *x) {
    double tmp = 1.;
    for(int j = 0; j < n; j++)
        if(j != i)
            tmp = ((y - x[j])/(x[i] - x[j]))*tmp;
    return tmp;
}

double L(int n, double y, double *x) {
    double tmp = 0.;
    for(int i = 0; i < n; i++)
        tmp += func(x[i])*F(i, n, y, x);
    return tmp;
}

/*void print_2(FILE*f, int n, int m, double *x, double *koef){
    double k;
    for(int i = 0; i < n-1; i++) {
        k = (x[i+1] - x[i])/m;
        for(int j = 0; j < m; j++) {
            fprintf(f, "%0.10f\t%0.10f\n", x[i] + j*k, func(x[i] + j*k));
            //fprintf(f, "x = %0.10f\tf(x) = %0.10f\tP(x) = %0.10f\tf(x) - P(x) =%0.10e\n", x[i] + j*k, func(x[i] + j*k), P(n, koef, x[i] + j*k), fabs(func(x[i] + j*k)-P(n, koef, x[i] + j*k)));
        }
        if(i == n-2)
            fprintf(f, "%0.10f\t%0.10f\n", x[n-1], func(x[n-1]));
           // fprintf(f, "x = %0.10f\tf(x) = %0.10f\tP(x) = %0.10f\tf(x) - P(x) = %0.10e\n", x[n-1], func(x[n-1]), P(n, koef, x[n-1]), fabs(func(x[n-1])-P(n, koef, x[n-1])));
    }
}
*/

void print_3(FILE*f, int n, int m, double *x, double *koef){
    double k;
    for(int i = 0; i < n-1; i++) {
        k = (x[i+1] - x[i])/m;
        for(int j = 0; j < m; j++) {
            fprintf(f, "x = %0.10f\tf(x) = %0.10f\tP(x) = %0.10f\tf(x) - P(x) = %0.10e\tL(x) = %0.10f\tL(x) - P(x) = %0.10e\n", x[i] + j*k, func(x[i] + j*k), P(n, koef, x[i] + j*k), fabs(func(x[i] + j*k)-P(n, koef, x[i] + j*k)), L(n, x[i] + j*k, x), fabs(L(n, x[i] + j*k, x)-P(n, koef, x[i] + j*k)));
        }
        if(i == n-2)
            fprintf(f, "x = %0.10f\tf(x) = %0.10f\tP(x) = %0.10f\tf(x) - P(x) = %0.10e\tL(x) = %0.10f\tL(x) - P(x) = %0.10e\n", x[n-1], func(x[n-1]), P(n, koef, x[n-1]), fabs(func(x[n-1])-P(n, koef, x[n-1])), L(n, x[n-1], x), fabs(L(n, x[n-1], x)-P(n, koef, x[n-1])));
    }
}
