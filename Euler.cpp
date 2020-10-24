#include <stdio.h>
#include <iostream>
#include <math.h>
#include <locale.h>

#define eps 1e-4
#define eps1 1e-6
#define eps2 1e-4
#define h_min 1e-8
#define h_max 1e-1

void True_y(double *y, double x);
void Func(double *f, double *y_curr, double x_k);
void Step(double *f, double *y_new, double *y_curr, double x_k, double h);
void Step_2(double* f, double *y_new1, double *y_new2, double *y_curr, double xk, double &h_curr, double &h_rec);
double Scalar_product(double *y_1, double *y_2, double h);


 void True_y(double *y, double x) {
     y[0] = exp(-(x*x - 4*x + 4));
    //y[1] = ...
}

 
void Func(double *f, double *y, double x) {
    f[0] = -2*(x - 2)*y[0];
    //f[1] = ...
}


 void Step(double *f, double *y_new, double *y_curr, double x_k, double h) {
    Func(f, y_curr, x_k);
    for(int i = 0; i < 1; i++)
        y_new[i] = y_curr[i] + h*f[i]; //Явный Эйлер
}

 
double Scalar_product(double *y_1, double *y_2, double h) {
    double tmp = 0.0;
    for(int i = 0; i < 1; i++)
        tmp += (y_1[i] - y_2[i])*(y_1[i] - y_2[i])*h;
    return sqrt(tmp);
}


 void Step_2(double* f, double *y_new1, double *y_new2, double *y_curr, double xk, double &h_curr, double &h_rec) {
     double c[1];
     Step(f, y_new1, y_curr, xk, h_curr);
     Step(f, y_new2, y_curr, xk, h_curr / 2.0);
     for(int i = 0; i < 1; i++) c[i]=y_new2[i];
     Step(f, y_new2, c, xk + (h_curr / 2.0), h_curr / 2.0);
     while (Scalar_product(y_new1, y_new2, h_curr) > eps2){
        h_curr /= 2.0;
        Step(f, y_new1, y_curr, xk, h_curr);
        Step(f, y_new2, y_curr, xk, h_curr / 2.0);
        for(int i = 0; i < 1; i++) c[i]=y_new2[i];
        Step(f, y_new2, c, xk + (h_curr / 2.0), h_curr / 2.0);
    }
    h_rec = h_curr;
}
 
 
 int main () {
     FILE*file = fopen("out.txt", "wt");
     FILE*file1 = fopen("out1.txt", "wt");
     FILE*file2 = fopen("out2.txt", "wt");
     double h_curr = 0.1, h_rec = 0.1, x_k = 0.0;
     double y_curr[1], y[1], y_new1[1], y_new2[1], f[1]; //стат. массивы длины 1
     y_curr[0] = exp(-4);
     
     while ((20.0 - x_k) > eps) {
         Step_2(f, y_new1, y_new2, y_curr, x_k, h_curr, h_rec);
         True_y(y, x_k);
         fprintf(file, "%0.10f\t%0.10f\n", x_k, h_rec);
         fprintf(file1, "%0.10f\t%0.10f\n", x_k, y_new1[0]);
         fprintf(file2, "%0.10f\t%0.10f\t%0.10f\t%0.10f\n", x_k, h_rec, y_new1[0], y[0]);
         if (Scalar_product(y_new1,y_new2,h_curr) < eps1)
             h_curr *= 2;
         if(h_curr > h_max)
             h_curr = h_max;
         if(h_curr < h_min)
             return -1;
         for (int i = 0; i < 1; i++) {
             y_curr[i] = y_new1[i];
             y_new2[i] = y_new1[i];
         }
         x_k += h_rec;
     }
     return 0;
 }
