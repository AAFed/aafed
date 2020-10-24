#include <iostream>
#include <math.h>

using namespace std;

#define eps1 1e-12
#define eps2 1e-10
#define hmax 1e-2
#define hmin 1e-18
#define N 2
#define X  3.141592653589793238//7//10

double* F(double *f, double *y_curr, double xk);
double scalar(double *y1, double *y2);
void RK(double *y_curr, double *y_new, double *f, double h, double xk);
void Step(double* f, double *y_new1, double *y_new2, double *y_curr, double xk, double &h_curr, bool &flag);
void y_true(double *y, double xk);

int main() {
    FILE*f0 = fopen("xh.txt", "wt");
    FILE*f1 = fopen("outr1.txt", "wt");
    FILE*f2 = fopen("outr2.txt", "wt");
    bool flag = true;
    double h_curr = 0.01, xk = 0.0;
    double f[N], y[N], y_curr[N], y_new1[N], y_new2[N];
    //y_curr[0] = exp(-5);
    y_curr[0] = 0.;
    y_curr[1] = 1.0;

    //начальные условия для полиномов, сделать h_curr=0.9
   //y_curr[0] = 1.0;
   //y_curr[1] = 1.0;

    y_true(y, xk);
    fprintf(f0, "%0.10f\t%0.10f\n", xk, h_curr);
    fprintf(f1, "%0.10f\t%0.10f\t%0.10f\t%0.10f\n", xk, y_curr[0],  y[0], fabs(y[0]-y_curr[0]));
    fprintf(f2, "%0.10f\t%0.10f\t%0.10f\t%0.10f\n", xk, y_curr[1], y[1], fabs(y[1]-y_curr[1]));

    while (flag) {
        Step(f, y_new1, y_new2, y_curr, xk, h_curr, flag);
        xk += h_curr;
        y_true(y, xk);
        fprintf(f0, "%0.15f\t%0.15f\n", xk, h_curr);
        fprintf(f1, "%0.15f\t%0.15f\t%0.15f\t%0.15f\n", xk, y_new1[0],  y[0], fabs(y_new1[0]-y[0]));
        fprintf(f2, "%0.15f\t%0.15f\t%0.15f\t%0.15f\n", xk, y_new1[1], y[1], fabs(y_new1[1]-y[1]));
        for (int i = 0; i < N; i++) y_curr[i] = y_new1[i];
        if (scalar(y_new1, y_new2) < eps1)
            h_curr *= 2.0;

        if (xk + h_curr > X) {
            h_curr = X - xk;
            flag = false;
        }
    }

    RK(y_curr, y_new1, f, h_curr, xk);
    xk += h_curr;
    y_true(y, xk);

    fprintf(f0, "%0.15f\t%0.10f\n", xk, h_curr);
    fprintf(f1, "%0.15f\t%0.15f\t%0.15f\t%0.15f\n", xk, y_new1[0],  y[0], fabs(y_new1[0]-y[0]));
    fprintf(f2, "%0.15f\t%0.15f\t%0.15f\t%0.15f\n", xk, y_new1[1], y[1], fabs(y_new1[1]-y[1]));
    return 0;
}


double* F(double *f, double *y_curr, double xk) {

    f[0] = y_curr[1];
    f[1] = -y_curr[0];

    //f[0] = 1 + 4 * xk*xk*xk; //должно считаться точно так как порядок метода s=3 и решением будет полином четвертой степени
    //f[1] = 1 + 5 * xk*xk*xk*xk;
    //f[0] = (-10)*(xk-1)*exp(-5*(xk-1)*(xk-1));

    return f;
}

double scalar(double *y1, double *y2) {
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += (y1[i] - y2[i])*(y1[i] - y2[i]);
    }
    return sqrt(sum);
}

void RK(double *y_curr, double *y_new, double *f, double h, double xk) {

    double k1[N], k2[N], k3[N], k4[N], tmp[N];

    for (int i = 0; i < N; i++) k1[i] = h*F(f, y_curr, xk)[i];
    for (int i = 0; i < N; i++) tmp[i] = y_curr[i] + k1[i] / 2;
    for (int i = 0; i < N; i++) k2[i] = h*F(f, tmp, xk + h/2)[i];
    for (int i = 0; i < N; i++) tmp[i] = y_curr[i] + (k1[i]+k2[i])/4;
    for (int i = 0; i < N; i++) k3[i] = h*F(f, tmp, xk + h/2)[i];
    for (int i = 0; i < N; i++) tmp[i] = y_curr[i] - k2[i]+2*k3[i];
    for (int i = 0; i < N; i++) k4[i] =  h*F(f, tmp, xk + h)[i];
    for (int i = 0; i < N; i++) y_new[i] = y_curr[i] + (k1[i] + 4*k3[i] + k4[i]) / 6;

}


void Step(double* f, double *y_new1, double *y_new2, double *y_curr, double xk, double &h_curr, bool &flag) {

    //double tmp[N];

    //if (h_curr > hmax) h_curr = hmax;

    RK(y_curr, y_new1, f, h_curr, xk);

    //уменьшение шага, убрать для полиномов
    /*RK(y_curr, y_new2, f, h_curr*0.5, xk);
    for (int i = 0; i < N; i++) tmp[i] = y_new2[i];
    RK(tmp, y_new2, f, h_curr*0.5, xk+h_curr*0.5);

    while (scalar(y_new1, y_new2) > eps2) {
        h_curr *= 0.5;
        RK(y_curr, y_new1, f, h_curr, xk);
        RK(y_curr, y_new2, f, h_curr*0.5, xk);
        for (int i = 0; i < N; i++) tmp[i] = y_new2[i];
        RK(tmp, y_new2, f, h_curr*0.5, xk + h_curr * 0.5);
    }// до сюда

    if (h_curr < hmin) {
        cout << "аварйная остановка" << endl;
        flag = false;
    }*/
}

void y_true(double *y, double xk) {

    y[0] = sin(xk);
    y[1] = cos(xk);

    //y[0] = 1 + xk + pow(xk, 4);
    //y[1] = 1 + xk + pow(xk, 5);
    //y[0] = exp(-5*(xk-1)*(xk-1));
}
