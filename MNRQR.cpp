#include <iostream>
#include <math.h>

using namespace std;

FILE* output;
FILE* check_matr;
FILE* check_vect;
int n, N, aa, bb, flag;

double func(double x);
void print_matr(int m, int k, double*matr);
void print_vect(int m, double*v);

void in_xb(double *x, double *b);
void in_A(double*a, double *x);
void transpose(double*a, double*at);
void in_AtA(double *x, double *a, double *at, double*ata);
void in_Atb(double *b, double*at, double*atb);
void in_Qt(double*qt, double*x);

void in_alfa(double*alfa);
void in_F(double *alfa, double*c, double *x, double*F);
void print_res(double*x, double*b, double*c, double*F, double*alfa);
double Gauss(int N, double *c, double *ata, double *atb);


int main(void) {
    check_matr = fopen("check_matr.txt", "wt");
    check_vect = fopen("check_vect.txt", "wt");

    cout << "Enter mode: 0 - MNR, 1 - QR" << endl;
    cin >> flag;
    cout << "Enter N, n, a, b:" << endl;
    cin >> N >> n >> aa >> bb; //100 10 -1 3 100 20 -1 1
    
    double* x = new double[N];
    double* b = new double[N];
    double* c = new double[n];
    double* F = new double[N];
    double* alfa = new double[N];
    double* A = new double[N*n];
    double* At = new double[n*N]; //Qt
    double* AtA = new double[n*n]; //R
    double* Atb = new double[n]; //Qtb
    
    in_xb(x, b);
    in_A(A, x);
    if(flag == 0) transpose(A, At);
    if(flag == 1) in_Qt(At,x);
    in_AtA(x, A, At, AtA);
    in_Atb(b, At, Atb);
    Gauss(n, c, AtA, Atb);
    print_res(x, b, c, F, alfa);
    
    delete[]x;
    delete[]b;
    delete[]c;
    delete[]alfa;
    delete[]F;
    delete[]A;
    delete[]At;
    delete[]AtA;
    delete[]Atb;
        
    fclose(check_matr);
    fclose(check_vect);
}




double func(double x){
    //return 1.0 / (25 * x * x + 1.0);
    //return fabs(x);
    return exp(x);
}

void print_matr(int m, int k, double*matr) {
    fprintf(check_matr, "\n\n\n");
    for(int i = 0; i < m; i ++) {
        for(int j = 0; j < k; j++)
            fprintf(check_matr, "%0.10f\t", matr[i*k + j]);
        fprintf(check_matr, "\n");
    }
}

void print_vect(int m, double*v) {
    fprintf(check_vect, "\n\n\n");
    for(int i = 0; i < m; i ++)
            fprintf(check_vect, "%0.10f\n", v[i]);
}

void in_xb(double *x, double *b){
    for(int i = 0; i < N; i++) {
        x[i] = (double)aa + ((double)bb - (double)aa)*i/((double)N - 1.);
        b[i] = func(x[i]);
    }
    //print_vect(N, x);
    //print_vect(N, b);
}

void in_A(double*a, double *x) {
    for(int i = 0; i < N; i++)
        for(int j = 0; j < n; j++)
            a[i*n + j] = pow(x[i], j);
    //print_matr(N, n, A);
}

void transpose(double*a, double*at) {
    for(int i = 0; i < N; i++)
        for(int j = 0; j < n; j++)
            at[j * N + i] = a[i * n + j];
    //print_matr(n, N, At);
}

void in_AtA(double *x, double *a, double *at, double*ata) {
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++) {
            ata[i * n + j] = 0.0;
            for(int k = 0; k < N; k++)
                ata[i * n + j] += at[i * N + k] * a[k * n + j];
        }
    
    //print_matr(n, n, AtA);
}

void in_Atb(double *b, double*at, double*atb) {
    for(int i = 0; i < n; i++) {
        atb[i] = 0.0;
        for(int j = 0; j < N; j++)
            atb[i] += at[i * N + j] * b[j];
    }
    //print_vect(n, Atb);
}

void in_alfa(double*alfa){
    for(int i = 0; i < N; i++) {
        alfa[i] = 1.;
    }
    //alfa[31] = ((double)bb - (double)aa) / 3.5;
    //alfa[32] = ((double)bb - (double)aa) / 5;
    //alfa[33] = ((double)bb - (double)aa) / (double)N;
    //alfa[34] = ((double)bb - (double)aa) / 5;
    //alfa[35] = ((double)bb - (double)aa) / 3.5;
    //print_vect(N, alfa);
}

void in_F(double *alfa, double*c, double *x, double*F) {
    for(int i = 0; i < N; i++) {
        F[i] = 0.0;
        for(int j = 0; j < n; j++)
            F[i] += alfa[i] * c[j] * pow(x[i], j);
    }
}

void print_res(double*x, double*b, double*c, double*F, double*alfa) {
    double inf = 0.0;
    
    in_alfa(alfa);
    in_F(alfa, c, x, F);
    
    if(flag == 0) output = fopen("outAtAx=Atb.txt", "wt");
    if(flag == 1) output = fopen("outQR.txt", "wt");
    
    fprintf(output, "x\t\t\tf\t\t\tF\t\t\t|f - F|\n");
    for (int i = 0; i < N; i++){
        inf += alfa[i] * (b[i] - F[i]) * (b[i] - F[i]);
        fprintf(output, "%0.15f\t%0.15f\t%0.15f\t%0.15e\n", x[i], b[i], F[i], (b[i] - F[i]));
    }
    fprintf(output, "inf = %0.15e\t", inf);
    fclose(output);
}

void in_Qt(double*qt, double*x) {
    double*q = new double[N*n];
    double tmp = 0.;
    in_A(q, x);
    for (int k = 0; k < n; k++){
        for(int j = 0; j < k; j++){
            for (int i = 0; i < N; i++)
                tmp += q[i * n + k] * q[i * n + j];
            for (int i = 0; i < N; i++)
                q[i * n + k] -= tmp * q[i * n + j];
            tmp = 0.0;
        }
        for(int i = 0; i < N; i++)
            tmp += q[i * n + k] * q[i * n + k];
        for (int i = 0; i < N; i++)
            q[i * n + k] /= sqrt(tmp);
        tmp = 0.0;
    }
    transpose(q, qt);
    //print_matr(n, N, qt);
    delete[]q;
}


double Gauss(int m, double *c, double *ata, double *atb){
    int i_max;
    double tmp, max;
    
    for (int i = 0; i < m; i++){
        max = fabs(ata[i * m + i]);
        i_max = i;
        for (int j = i + 1; j < m; j++){
            if (max < fabs(ata[j * m + i])){
                max = fabs(ata[j * m + i]);
                i_max = j;
            }
        }
        if (fabs(ata[i * m + i]) < 1e-15) return -1;
        if (i_max != i) {
            for (int j = i; j < m; j++) {
                tmp = ata[i * m + j];
                ata[i * m + j] = ata[i_max * m + j];
                ata[i_max * m + j] = tmp;
            }
            tmp = atb[i];
            atb[i] = atb[i_max];
            atb[i_max] = tmp;
        }
        
        tmp = 1.0 / ata[i * m + i];
        for (int j = i; j < m; ++j)
            ata[i * m + j] *= tmp;
        atb[i] *= tmp;

        for (int j = i + 1; j < m; j++) {
            tmp = ata[j * m + i];
            for (int k = i; k < m; k++)
                ata[j * m + k] -= ata[i * m + k] * tmp;
            atb[j] -= atb[i] * tmp;
        }
    }

    for (int i = m - 1; i >= 0; --i) {
        tmp = atb[i];
        for (int j = i + 1; j < m; j++)
            tmp -= ata[i * m + j] * c[j];
        c[i] = tmp;
    }
    //print_vect(n, c);
    return 0;
}

