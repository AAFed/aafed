#include <iostream>
#include <fstream>
#include <math.h>
#include <locale.h>
//#pragma warning(disable: 4996)
using namespace std;
FILE*out;
FILE *out1 = fopen("outgauss.txt", "wt");

int nx = 12; //OX
int ny = 12; //OY
double ax = 0., bx = 0.1, ay = 0., by = 0.1;
int p = 3; //для Паде

double func(double x, double y);
double Gauss(int n, double *a, double *b, double *x);
void Lagr_Interpolation(double *y, int *Index, double *F, double *L);
double Ph(int i, int n, double y, double *x);
double Lagr(int n, double y, double *Y, double *FY);
void Pade_Interpolation(double *x, int *Index, double *F, double *Pade);
double R(int p, int q, double *a, double *b, double x);
void PrintMatrix(int m, int n, double *Matrix, int flag);

int main(void) {
    out = fopen("result.txt", "wt");
    double *x = new double[nx];
    double *y = new double[ny];
    int *Index = new int[nx*ny];
    double *F = new double[nx*ny]; //all function values
    double *L = new double[nx*ny]; //Lagr Interpolation
    double *Pade = new double[nx*ny]; //Pade Interpolation
    double *Res = new double[nx*ny]; //Result = (L + Pade)/2
    double *Error = new double[nx*ny];
    
    for (int i = 0; i < nx*ny; i++) {
        Index[i] = 0;
        F[i] = 0.;
        L[i] = 0.;
        Pade[i] = 0.;
        Res[i] = 0.;
        Error[i] = 0.;
    }
    
    //x, y
    for (int i = 0; i < nx; i++) x[i] = ax + (bx - ax)*(double)i / ((double)nx - 1.);
    for (int i = 0; i < ny; i++) y[i] = ay + (by - ay)*(double)i / ((double)ny - 1.);
    
    //Index
    fprintf(out, "Index Matrix:\n");
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < nx; j++) {
            if (i == 0 || j == 0 || i == ny - 1 || j == nx - 1) Index[i*nx + j] = 1;
            else Index[i*nx + j] = (i + j) % 2;
            fprintf(out, "%d\t", Index[i*nx + j]);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "\n\n\n");
    
    //F
    for (int i = 0; i < ny; i++)
        for (int j = 0; j < nx; j++)
            F[i*nx + j] = func(x[j], y[i]);

    //Лагранж по столбцам
    Lagr_Interpolation(y, Index, F, L);
    
    //Паде по строкам
    Pade_Interpolation(x, Index, F, Pade);
    
    //Res, Error
    for (int i = 0; i < ny; i++)
        for (int j = 0; j < nx; j++) {
            Res[i*nx + j] = (L[i*nx + j] + Pade[i*nx + j]) / 2;
            Error[i*nx + j] = fabs(Res[i*nx + j] - F[i*nx + j]);
        }

    PrintMatrix(ny, nx, F, 1);
    PrintMatrix(ny, nx, L, 2);
    PrintMatrix(ny, nx, Pade, 3);
    PrintMatrix(ny, nx, Res, 4);
    PrintMatrix(ny, nx, Error, 5);
    delete[]F;
    delete[]L;
    delete[]Pade;
    delete[]Res;
    delete[]Error;
    fclose(out);
}

double func(double x, double y) {
    //return exp(x + y);
    //return pow(x, 4) + pow(y, 4);
    return exp(-pow(x,2) - pow(y,2));
}

void Lagr_Interpolation(double *y, int *Index, double *F, double *L) {
    int count = 0, tmp = 0;
    double *Y, *FY;
    
    for (int i = 0; i < nx; i++) {
        
        //известные значения из F записываем в L сразу
        count = 0;
        for (int j = 0; j < ny; j++) {
            if (Index[j*nx + i] != 0) {
                L[j*nx + i] = F[j*nx + i];
                count++;
            }
        }
        
        Y = new double[count];
        FY = new double[count];
        for (int s = 0; s < count; s++) {
            Y[s] = 0.;
            FY[s] = 0.;
        }
        
        //известные y и значения функции в них
        tmp = 0;
        for (int s = 0; s < ny; s++)
            if (Index[s*nx + i] != 0) {
                FY[tmp] = F[s*nx + i];
                Y[tmp] = y[s];
                tmp++;
            }
        
        //неизвестные восстанавливаем
        for (int s = 0; s < ny; s++)
            if (Index[s*nx + i] == 0)
                L[s*nx + i] = Lagr(count, y[s], Y, FY);
        
        delete[]Y;
        delete[]FY;
    }
}

double Ph(int i, int n, double y, double *x) {
    double tmp = 1.;
    for(int j = 0; j < n; j++)
        if(j != i)
            tmp = ((y - x[j])/(x[i] - x[j]))*tmp;
    return tmp;
}

double Lagr(int n, double y, double *Y, double *FY) {
    double tmp = 0.;
    for(int i = 0; i < n; i++)
        tmp += FY[i]*Ph(i, n, y, Y);
    return tmp;
}

double R(int p, int q, double *a, double *b, double x) {
    
    double tmp1 = 0.;
    tmp1 = b[q];
    for (int i = 1; i < q + 1; i++)
        tmp1 = b[q - i] + x * tmp1;

    double tmp2 = 0.;
    tmp2 = a[p];
    for (int i = 1; i < p + 1; i++)
        tmp2 = a[p - i] + x * tmp2;

    return tmp2 / tmp1;
}

void Pade_Interpolation(double *x, int *Index, double *F, double *Pade) {
    
    double *X, *FX, *BA, *B, *Solve, *a, *b;
    int count = 0, tmp = 0;
    
    for (int i = 0; i < ny; i++) {
        
        count = 0;
        for (int j = 0; j < nx; j++)
            if (Index[i*nx + j] != 0) {
                Pade[i*nx + j] = F[i*nx + j];
                count++;
            }

        X = new double[count];
        FX = new double[count];
        for (int k = 0; k < count; k++) {
            X[k] = 0.;
            FX[k] = 0.;
        }
        
        tmp = 0;
        for (int s = 0; s < nx; s++)
            if (Index[i*nx + s] != 0) {
                FX[tmp] = F[i*nx + s];
                X[tmp] = x[s];
                tmp++;
            }
        
        int q = count - p - 1;
        BA = new double[count * count];
        for (int k = 0; k < count*count; k++) BA[i] = 0.;

        for (int k = 0; k < count; k++)
            for (int l = 0; l < count; l++) {
                if(l < q) BA[k*count + l] = FX[k] * pow(X[k], l + 1);
                else BA[k*count + l] = -1.*pow(X[k], l - q);
            }
        
        /*cout << "BA Matrix:" << endl;
        for (int i = 0; i < count; i++) {
            for (int j = 0; j < count; j++)
                cout << BA[i * count + j] << " ";
            cout << endl;
        }*/

        B = new double[count];
        Solve = new double[count];
        for (int k = 0; k < count; k++) {
            B[i] = 0.;
            Solve[i] = 0.;
        }

        a = new double[p + 1];
        b = new double[q + 1];
        for (int k = 0; k < p + 1; k++) a[k] = 0.;
        for (int k = 0; k < q + 1; k++) b[k] = 0.;

        for (int k = 0; k < count; k++) B[k] = -1.*FX[k];
        Gauss(count, BA, B, Solve);
        /*for (int k = 0; k < count; k++) fprintf(out1, "%0.5f\t", Solve[k]);
        fprintf(out1, "\n\n");*/
        //for (int k = 0; k < count; k++) cout<<Solve[i]<<endl;

        b[0] = 1.;
        for (int w = 1; w < q; w++) {
            b[w] = Solve[w - 1];
            //cout << b[w] <<"  ";
            
        }
        for (int w = q; w < count; w++) a[w - q] = Solve[w];
        for (int u = 0; u < nx; u++) {
            if (Index[i*nx + u] == 0) Pade[i*nx + u] = R(p, q, a, b, x[u]);
        }

        /*delete[]X;
        delete[]FX;
        delete[]a;
        delete[]b;
        delete[]B;
        delete[]BA;
        delete[]Solve;*/
    }

    
}

double Gauss(int n, double *a, double *b, double *x) {
    int i, j, k, indMax = 0;
    double tmp = 0., max = 0.;

    for (i = 0; i < n; ++i) {
        max = fabs(a[i * n + i]);
        indMax = i;
        for (j = i + 1; j < n; ++j)
            if (max < fabs(a[j * n + i])) {
                max = fabs(a[j * n + i]);
                indMax = j;
            }
        if (indMax != i) {
            tmp = 0.;
            for (int j = 0; j < n; ++j) {
                tmp = a[i * n + j];
                a[i * n + j] = a[indMax*n + j];
                a[indMax*n + j] = tmp;
            }
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

void PrintMatrix(int m, int n, double *Matrix, int flag) {
    if (flag == 0) fprintf(out, "Matrix:\n");
    if (flag == 1) fprintf(out, "Function Values Matrix:\n");
    if (flag == 2) fprintf(out, "L Interpolation Matrix:\n");
    if (flag == 3) fprintf(out, "Pade Interpolation Matrix:\n");
    if (flag == 4) fprintf(out, "Result Matrix:\n");
    if(flag > -1 && flag < 5) {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++)
                fprintf(out, "%0.4e\t", Matrix[i*n + j]);
            fprintf(out, "\n");
            }
        fprintf(out, "\n\n\n");
    }
    if (flag == 5) { fprintf(out, "Error Matrix:\n");
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++)
                fprintf(out, "%0.4e\t", Matrix[i*n + j]);
            fprintf(out, "\n");
        }
        fprintf(out, "\n\n\n");
    }
}

