#include <iostream>
#include <math.h>
#include <fstream>

#define EPS1 1e+10
#define EPS2 1e-15
#define acc 14


using namespace std;
ofstream out;
int N;

void func(double*x, double*f) {
    f[0] = x[0]*x[0]/8. - sin(x[1]/6.);
    f[1] = sin(x[0]*x[1]/4.) + sin(x[1])- exp(x[2]);
    f[2] = sin(x[0]*x[1]) + cos(x[1]) + exp(x[2]);
    //f[0] = exp(x[0])*sin(x[1] * M_PI) - 1.;
    //f[1] = cos(x[1] * M_PI) + sin(x[0]);
}

double norm(double*f) {
    double tmp = 0.;
    for(int i = 0; i < N; i++) tmp += f[i]*f[i];
    return sqrt(tmp);
}

void Newton(double*x, double*f, double*A1){
    for(int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            x[i] -= A1[i*N+j]*f[j];
}

int InvertA(double*a, double*x) {
    double tmp1, tmp2;
    
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            x[i * N + j] = (double)(i == j);

    for (int i = 0; i < N; ++i) {
        tmp1 = 0.0;
        for (int j = i + 1; j < N; ++j)
            tmp1 += a[j * N + i] * a[j * N + i];
        
        tmp2 = sqrt(tmp1 + a[i * N + i] * a[i * N + i]);
        if (tmp2 < 1e-30) return -1;
        a[i * N + i] -= tmp2;
        
        tmp1 = sqrt(tmp1 + a[i * N + i] * a[i * N + i]);
        if (tmp1 < 1e-30) {
            a[i * N + i] += tmp2;
            continue;
        }
        tmp1 = 1. / tmp1;
        for (int j = i; j < N; ++j) a[j * N + i] *= tmp1;

        for (int k = i + 1; k < N; ++k) {
            tmp1 = 0.0;
            for (int j = i; j < N; ++j)
                tmp1 += a[j * N + i] * a[j * N + k];

            tmp1 *= 2.0;
            for (int j = i; j < N; ++j)
                a[j * N + k] -= tmp1 * a[j * N + i];
        }

        for (int k = 0; k < N; ++k) {
            tmp1 = 0.0;
            for (int j = i; j < N; ++j)
                tmp1 += a[j * N + i] * x[j * N + k];

            tmp1 *= 2.0;
            for (int j = i; j < N; ++j)
                x[j * N + k] -= tmp1 * a[j * N + i];
        }

        a[i * N + i] = tmp2;
    }

    for (int k = 0; k < N; ++k)
        for (int i = N - 1; i >= 0; --i) {
            tmp1 = x[i * N + k];
            for (int j = i + 1; j < N; ++j)
                tmp1 -= a[i * N + j] * x[j * N + k];
            x[i * N + k] = tmp1 / a[i * N + i];
        }
    
    return 0;
}

int result(double*x) {
    int k = 0, flag = 1;
    double*A = new double[N*N];
    double*A1 = new double[N*N];
    double* x1 = new double[N];
    double* f = new double[N];
    double* f1 = new double[N];
    
    out.open("output.txt");
    while(true) {
        k += 1;
        func(x, f);
        out << "Итерация: " << k << endl;
        if(norm(f) < EPS2) break;
        if(k >= 100) {
            flag = -1;
            break;
         }
         for(int i = 0; i < N; i++) out << setprecision(acc) << fixed << x[i] << "    ";
         out << endl << endl;
         for(int i = 1; i < N; i++) x1[i] = x[i];
         
         for(int i = 0; i < N; i++)
             for(int j = 0; j < N; j++){
                 x1[j] = 1./EPS1 + x[j];
                 func(x1,f1);
                 A[i*N+j] =(f1[i]-f[i])*EPS1;
                 x1[j] = x[j];
             }
         
         InvertA(A,A1);
         Newton(x,f,A1);
     }
     for(int i = 0; i < N; i++) out << setprecision(acc) << fixed << x[i] << "    ";
    
     out << endl << endl;
     if(flag == 1) {
         out << "Решение = " << endl;
         for(int i = 0; i < N; i ++) out << setprecision(acc) << fixed << x[i] << "    ";
     }
     else out << "Метод не сходится к решению" << endl;
    
    delete[]x1;
    delete[]f1;
    delete[]f;
    delete[]A1;
    delete[]A;
    return 0;
    
}


int main(void) {
    N = 3;//cout << "Enter N:" << endl; cin >> N;
    double* x = new double[N];
    cout << "Enter x0 vector:" << endl; //0.4 0.4 N=2; 2.3 2.5 1 N=3
    for(int i = 0; i < N; i++) cin >> x[i];
    result(x);
    delete[]x;
}
