#include <iostream>
#include <math.h>

using namespace std;

FILE* out;

double ax = 0.5, bx = 1.5, ay = 0.5, by = 1.5, hx = 0., hy = 0.;
int nx = 0, ny = 0, basisfunc = 4;

double func(double x, double y);
void InputMatrices(double *RealValuesMatrix, int *KnownValuesIndices);
void CubedSplain(double *RealValuesMatrix, int *KnownValuesIndices, double *SplainMatrix);
double CubedSplainFormula(int tmp, int s, double *X, double *Y, double *M);
void InputL2InterpolationMatrices(int KnownValues, double*A, double *Q, double*Qt, double*QtA, double*QtY, double*X, double*Y);
void L2Interpolation(double *RealValuesMatrix, int *KnownValuesIndices, double *L2InterpolationMatrix);
void Result(double *RealValuesMatrix, double *SplainMatrix, double *L2InterpolationMatrix, double *ResultMatrix, double *ErrorMatrix);

void Progonka(int M, double*result, double*C, double*P);
double Gauss(int n, double *a, double *b, double *x);
void Null(double n, double*a);
void PrintMatrix(int m, int n, double *Matrix, int flag);

int main(void) {
    out = fopen("output.txt", "wt");
    //cout << "Enter [ax, bx] , [ay, by]:" << endl;
    //cin >> ax >> bx >> ay >> by;
    cout << "Enter nx, ny:" << endl;
    cin >> nx >> ny;
    
    int *KnownValuesIndices = new int[nx*ny];
    double *RealValuesMatrix = new double[nx*ny];
    double *SplainMatrix = new double[nx*ny];
    double *L2InterpolationMatrix = new double[nx*ny];
    double *ResultMatrix = new double[nx*ny];
    double *ErrorMatrix = new double[nx*ny];

    Null(nx*ny, RealValuesMatrix);
    Null(nx*ny, SplainMatrix);
    Null(nx*ny, L2InterpolationMatrix);
    Null(nx*ny, ResultMatrix);
    Null(nx*ny, ErrorMatrix);

    InputMatrices(RealValuesMatrix, KnownValuesIndices);
    CubedSplain(RealValuesMatrix, KnownValuesIndices, SplainMatrix);
    L2Interpolation(RealValuesMatrix, KnownValuesIndices, L2InterpolationMatrix);
    Result(RealValuesMatrix, SplainMatrix, L2InterpolationMatrix, ResultMatrix, ErrorMatrix);

    delete[]KnownValuesIndices;
    delete[]SplainMatrix;
    delete[]L2InterpolationMatrix;
    delete[]ResultMatrix;
    delete[]ErrorMatrix;
    fclose(out);
}

double func(double x, double y) {
    return x + x*x + y + y*y;
    //return exp(- x*x - y*y);
}

void InputMatrices(double *RealValuesMatrix, int *KnownValuesIndices) {

    hx = (bx - ax) / ((double)nx - 1.);
    hy = (by - ay) / ((double)ny - 1.);

    for (int i = 0; i < nx*ny; i++) KnownValuesIndices[i] = 0;
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++) {
            if (j == 0 || j == ny - 1 || i == 0 || i == nx - 1) KnownValuesIndices[j*nx + i] = 1;
            else KnownValuesIndices[j*nx + i] = (i + j) % 3;
        }

    fprintf(out, "Known Values Indices Matrix:\n");
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < nx; j++)
            fprintf(out, "%d    ", KnownValuesIndices[i*nx + j]);
        fprintf(out, "\n");
    }
    fprintf(out, "\n\n\n");

    for (int i = 0; i < ny; i++)
        for (int j = 0; j < nx; j++)
            RealValuesMatrix[i*nx + j] = func(ax + j * hx, ay + i * hy);
    //PrintMatrix(ny, nx, RealValuesMatrix, 1);
}

void CubedSplain(double *RealValuesMatrix, int *KnownValuesIndices, double *SplainMatrix) {

    int KnownValues, tmp;
    double *X, *Y, *C, *M, *d;
    double *ProgonkaSolution;
    //double*GaussSolution;

    for (int i = 0; i < nx; i++) {

        KnownValues = 0;
        for (int j = 0; j < ny; j++)
            if (KnownValuesIndices[j*nx + i] != 0) {
                SplainMatrix[j*nx + i] = RealValuesMatrix[j*nx + i];
                KnownValues++;
            }

        X = new double[KnownValues];
        Y = new double[KnownValues];
        Null(KnownValues, X);
        Null(KnownValues, Y);
       
        tmp = 0;
        for (int s = 0; s < ny; s++)
            if (KnownValuesIndices[s*nx + i] != 0) {
                Y[tmp] = RealValuesMatrix[s*nx + i];
                X[tmp] = ay + s * hy;
                tmp++;
            }

        C = new double[(KnownValues - 2)*(KnownValues - 2)];
        Null((KnownValues - 2)*(KnownValues - 2), C);
        for (int s = 0; s < KnownValues - 2; s++)
            for (int j = 0; j < KnownValues - 2; j++) {
                if (s == j) C[s*(KnownValues - 2) + j] = (X[s + 2] - X[s]) / 3.;
                else if (j == s + 1) C[s*(KnownValues - 2) + j] = (X[j + 1] - X[j]) / 6.;
                else if (s == j + 1) C[s*(KnownValues - 2) + j] = (X[s + 1] - X[s]) / 6.;
                else C[s*(KnownValues - 2) + j] = 0.0;
            }
        //PrintMatrix(KnownValues - 2, KnownValues - 2, C, 0);

        d = new double[KnownValues - 2];
        Null((KnownValues - 2), d);
        for (int s = 0; s < KnownValues - 2; s++) d[s] = (Y[s + 2] - Y[s + 1]) / (X[s + 2] - X[s + 1]) - (Y[s + 1] - Y[s]) / (X[s + 1] - X[s]);

        /*GaussSolution = new double[KnownValues - 2];
        Null((KnownValues - 2), GaussSolution);
        Gauss(KnownValues - 2, C, d, GaussSolution);*/
        ProgonkaSolution = new double[KnownValues - 2];
        Null((KnownValues - 2), ProgonkaSolution);
        Progonka(KnownValues - 2, ProgonkaSolution, C, d);

        M = new double[KnownValues];
        Null(KnownValues, M);
        M[0] = 0.;
        M[KnownValues - 1] = 0.;
        for (int s = 1; s < KnownValues - 1; ++s) {
            //M[s] = GaussSolution[s - 1];
            M[s] = ProgonkaSolution[s - 1];
        }

        tmp = 0;
        for (int s = 0; s < ny; s++) {
            if (KnownValuesIndices[s*nx + i] != 0) tmp++;
            else if (KnownValuesIndices[s*nx + i] == 0) {
                SplainMatrix[s*nx + i] = CubedSplainFormula(tmp, s, X, Y, M);
            }
        }

        delete[]X;
        delete[]Y;
        delete[]M;
        delete[]d;
        //delete[]GaussSolution;
        delete[]ProgonkaSolution;
        delete[]C;
    }
    //PrintMatrix(ny, nx, SplainMatrix, 2);
}

double CubedSplainFormula(int tmp, int s, double *X, double *Y, double *M) {
    return M[tmp - 1] * pow(X[tmp] - (ay + s * hy), 3) / (6. * (X[tmp] - X[tmp - 1])) +
        M[tmp] * pow((ay + s * hy) - X[tmp - 1], 3) / (6. * (X[tmp] - X[tmp - 1])) +
        (Y[tmp - 1] - (M[tmp - 1] * (X[tmp] - X[tmp - 1])*(X[tmp] - X[tmp - 1])) / 6.) * (X[tmp] - (ay + s * hy)) / (X[tmp] - X[tmp - 1]) +
        (Y[tmp] - (M[tmp] * (X[tmp] - X[tmp - 1])*(X[tmp] - X[tmp - 1])) / 6.0) * (ay + s * hy - X[tmp - 1]) / (X[tmp] - X[tmp - 1]);
}

void InputL2InterpolationMatrices(int KnownValues, double*A, double*Q, double*Qt, double*QtA, double*QtY, double*X, double*Y) {
    
    for (int i = 0; i < KnownValues; i++)
        for (int j = 0; j < basisfunc; j++) {
            A[i*basisfunc + j] = pow(X[i], j);
            Q[i*basisfunc + j] = pow(X[i], j);
        }

     double tmp = 0.; //Q
     for (int k = 0; k < basisfunc; k++){
         for(int j = 0; j < k; j++){
             for (int i = 0; i < KnownValues; i++)
                 tmp += Q[i * basisfunc + k] * Q[i * basisfunc + j];
             for (int i = 0; i < KnownValues; i++)
                 Q[i * basisfunc + k] -= tmp * Q[i * basisfunc + j];
             tmp = 0.0;
         }
         for(int i = 0; i < KnownValues; i++)
             tmp += Q[i * basisfunc + k] * Q[i * basisfunc + k];
         for (int i = 0; i < KnownValues; i++)
             Q[i * basisfunc + k] /= sqrt(tmp);
         tmp = 0.0;
    }

    for (int i = 0; i < basisfunc; i++) //Qt
        for (int j = 0; j < KnownValues; j++)
            Qt[i * KnownValues + j] = Q[j * basisfunc + i];

    for (int i = 0; i < basisfunc; i++) //QtA
        for (int j = 0; j < basisfunc; j++) {
            QtA[i * basisfunc + j] = 0.;
            for (int k = 0; k < KnownValues; k++)
                QtA[i * basisfunc + j] += Qt[i * KnownValues + k] * A[k * basisfunc + j];
        }

    for (int i = 0; i < basisfunc; i++) { //QtY
        QtY[i] = 0.0;
        for (int j = 0; j < KnownValues; j++)
            QtY[i] += Qt[i * KnownValues + j] * Y[j];
    }
}


void L2Interpolation(double *RealValuesMatrix, int *KnownValuesIndices, double *L2InterpolationMatrix) {
    
    int KnownValues, tmp1, tmp2;
    double *X, *Y, *X2, *A, *Q, *Qt, *QtA, *QtY, *GaussSolution;

    for (int i = 0; i < ny; i++) {

        KnownValues = 0;
        for (int j = 0; j < nx; j++)
            if (KnownValuesIndices[i*nx + j] != 0) {
                L2InterpolationMatrix[i*nx + j] = RealValuesMatrix[i*nx + j];
                KnownValues++;
            }

        X = new double[KnownValues];
        Y = new double[KnownValues];
        Null(KnownValues, X);
        Null(KnownValues, Y);
        X2 = new double[nx - KnownValues];
        Null(nx - KnownValues, X2);
        tmp1 = 0; tmp2 = 0;
        for (int j = 0; j < nx; j++) {
            if (KnownValuesIndices[i*nx + j] != 0) {
                Y[tmp1] = RealValuesMatrix[i*nx + j];
                X[tmp1] = ax + j * hx;
                tmp1++;
            }
            else if (KnownValuesIndices[i*nx + j] == 0) {
                X2[tmp2] = ax + j * hx;
                tmp2++;
            }
        }

        A = new double[KnownValues*basisfunc];
        Q = new double[KnownValues*basisfunc];
        Qt = new double[basisfunc*KnownValues];
        QtA = new double[basisfunc*basisfunc];
        QtY = new double[basisfunc];
        GaussSolution = new double[basisfunc];
        Null(KnownValues*basisfunc, A);
        Null(KnownValues*basisfunc, Q);
        Null(KnownValues*basisfunc, Qt);
        Null(basisfunc*basisfunc, QtA);
        Null(basisfunc, QtY);
        Null(basisfunc, GaussSolution);

        InputL2InterpolationMatrices(KnownValues, A, Q, Qt, QtA, QtY, X, Y);
        Gauss(basisfunc, QtA, QtY, GaussSolution);

        tmp2 = 0;
        for (int j = 0; j < nx; j++) {
            if (KnownValuesIndices[i*nx + j] == 0) {
                for (int k = 0; k < basisfunc; k++)
                    L2InterpolationMatrix[i*nx + j] += GaussSolution[k] * pow(X2[tmp2], k);

                tmp2++;
            }
        }

        delete[]X;
        delete[]Y;
        delete[]X2;
        delete[]A;
        delete[]Q;
        delete[]Qt;
        delete[]QtA;
        delete[]QtY;
        delete[]GaussSolution;
    }
    //PrintMatrix(ny, nx, L2InterpolationMatrix, 3);
}

void Null(double n, double*a) {
    for (int i = 0; i < n; i++) a[i] = 0.;
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


void Result(double *RealValuesMatrix, double *SplainMatrix, double *L2InterpolationMatrix, double *ResultMatrix, double *ErrorMatrix) {
    
    for (int i = 0; i < ny; i++)
        for (int j = 0; j < nx; j++) {
            ResultMatrix[i*nx + j] = (SplainMatrix[i*nx + j] + L2InterpolationMatrix[i*nx + j]) / 2;
            ErrorMatrix[i*nx + j] = fabs(ResultMatrix[i*nx + j] - RealValuesMatrix[i*nx + j]);
        }

    PrintMatrix(ny, nx, RealValuesMatrix, 1);
    PrintMatrix(ny, nx, SplainMatrix, 2);
    PrintMatrix(ny, nx, L2InterpolationMatrix, 3);
   // PrintMatrix(ny, nx, ResultMatrix, 4);
    PrintMatrix(ny, nx, ErrorMatrix, 5);
}

void PrintMatrix(int m, int n, double *Matrix, int flag) {
    if (flag == 0) fprintf(out, "Matrix:\n");
    if (flag == 1) fprintf(out, "Real Values Matrix:\n");
    if (flag == 2) fprintf(out, "Splain Matrix:\n");
    if (flag == 3) fprintf(out, "L2 Interpolation Matrix:\n");
    if (flag == 4) fprintf(out, "Result Matrix:\n");
    if (flag == 5) fprintf(out, "Error Matrix:\n");
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            fprintf(out, "%.5f    ", Matrix[i*n + j]);
        fprintf(out, "\n");
    }
    fprintf(out, "\n\n\n");
}

void Progonka(int M, double*result, double*C, double*P) {

    double *Alpha, *Beta, *a, *b, *c;
    int N;
    N = M - 1;
    Alpha = new double[N+1];
    Beta = new double[N+1];
    a = new double[N];
    b = new double[N];
    c = new double[N+1];
    Null(N+1, Alpha);
    Null(N+1, Beta);
    Null(N, a);
    Null(N, b);
    Null(N+1, c);
    Null(N+1, result);
    
    //заполнение a b c
    for (int k = 0; k < N+1; k++) {
        for (int j = 0; j < N+1; j++) {
            if (j == k + 1) b[k] = C[k*M + j];
            else if (j == k) c[k] = C[k*M + k];
            else if (j == k - 1) a[k] = C[k*M + j];
        }
    }
    
    //заполнение alpha beta и решения
    Alpha[1] = b[0] / c[0];
    Beta[1] = P[0] / c[0];
    for (int i = 1; i < N; i++) {
        Alpha[i + 1] = b[i] / (c[i] - Alpha[i] * a[i]);
        Beta[i + 1] = (P[i] + a[i] * Beta[i]) / (c[i] - Alpha[i] * a[i]);
    }

    result[N] = (P[N] + a[N] * Beta[N]) / (c[N] - a[N] * Alpha[N]);
    for (int i = N - 1; i >= 0; i--)
        result[i] = Alpha[i + 1] * result[i + 1] + Beta[i + 1];

    delete[]Alpha;
    delete[]Beta;
    delete[]a;
    delete[]b;
    delete[]c;
}
