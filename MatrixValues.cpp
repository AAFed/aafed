#include <iostream>
#include <math.h>
#include <locale.h>

using namespace std;
int N;
double h;

void matrix(double *matr);
double Y(int m, int k); //- собств. ф-ии; m = {1,..,N} - номер в-ра, k={0,...,N} - номер координаты
double L(int m); //-собственные значения
double Norm(double *matr, int m); //относительная погр-ть ||Ay+Ly||/||Ly||

int main() {
    setlocale(LC_ALL, "rus");
    int max_deviation_number = 0, max_scalar_number1 = 0, max_scalar_number2 = 0;
    double max_deviation = 0.0, tmp = 0.0, max_scalar = 0.0;
    printf("Enter N\n");
    cin>>N;
    
    h = 1/ ((double)N);
    double *matr = new double[N*N]; //минус 1ая строка и столбец, y0=0
    
    for (int i = 0; i < N * N; i++)
        matr[i] = 0;
    matrix(matr);
//вывод св сз и матрицы
  /*  printf("СОБСТВЕННЫЕ ЗНАЧЕНИЯ:\n");
    for (int m = 1; m < N+1; m++)
        printf("%0.15f\n", L(m));

    printf("СОБСТВЕННЫЕ ВЕКТОРЫ:\n");
    for (int m = 1; m < N; m++) {
            for (int k = 0; k < N + 1; k++) {
                printf("%0.15f\n", Y(m, k));
            }
        printf("\n");
    }

    printf("МАТРИЦА:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%f  ", matr[j + i * N]);
        }
        printf("\n\n");
    }
*/
// ск. произв.имеет вид
//(x,y)=<сумма по j={0,..,N-1}> x_j*y_j*h + 1/2*x_N*y_N*h
    
//ск.произв. разных в-ров + поиск макс.
    for (int m = 1; m < N + 1; m++) {
        for (int q = m+1; q < N + 1 ; q++) {
                for (int i = 1; i < N + 1; i++){
                    if(i == N)
                        tmp += (1.0/(2*N)) * Y(m, i) * Y(q, i);
                    else
                        tmp += (1.0/N) * Y(m, i) * Y(q, i);
                }

                if (fabs(tmp) > fabs(max_scalar)) {
                    max_scalar = tmp;
                    max_scalar_number1 = m;
                    max_scalar_number2 = q;
                }
                tmp = 0.0;
        }
    }

    //ск. произв. для одинак. в-ров
 /*   for (int m = 1; m < N + 1; m++) {
        for (int i = 1; i < N + 1; i++){
            if(i == N)
                tmp += (h/2)* Y(m, i) * Y(m, i);
            else
                tmp += h * Y(m, i) * Y(m, i);
            }
        printf("same %0.10f \n", tmp);
        tmp = 0.0;
    }
  */
   //||Ay+Ly||/||Ly|| - относит.погр-ть (отклонение)
    for (int m = 1; m < N+1; m++) {
        tmp = Norm(matr, m);
        if (tmp > max_deviation) {
            max_deviation = tmp;
            max_deviation_number = m;
        }
        tmp = 0.0;
    }
    printf("МАКСИМАЛЬНОЕ ОТКЛОНЕНИЕ: %e НА ВЕКТОРE: %d \n", max_deviation, max_deviation_number);
    printf("МАКСИМАЛЬНОЕ СКАЛЯРНОЕ ПРОИЗВЕДЕНИЕ: %e НА ВЕКТОРАХ: %d  %d \n", max_scalar, max_scalar_number1, max_scalar_number2);
    FILE*f = fopen("out.txt", "w");
    for(int i = 0; i < N; i++) {
        fprintf(f, "%f\t", i*h);
        fprintf(f, "%f\n", Y(1, i));
        
    }
    return 0;
}

void matrix(double *matr) {
    for (int i = 0; i < N-1; i++) {
        for (int j = 0; j < N-1; j++) {
            if (i == j)
                matr[j + i * N] = -2/(h*h);
            else if ((i - j == 1) || (j - i == 1))
                matr[j + i * N] = 1.0/(h*h);
        }
    }
    matr[N*(N-1) - 1] = 1.0/(h*h);
    matr[N*N - 1] = -2.0/(h*h);
    matr[N*N - 2] = 2.0/(h*h);
}

double Y(int m, int k) {
    return (sqrt(2))*sin(M_PI*(2 * m - 1)*k / (2 * N));
}

double L(int m) {
    return (4/(h*h))*sin(M_PI*(2 * m - 1) / (4 * N))*sin(M_PI*(2 * m - 1) / (4 * N));
}

double Norm(double* matr, int m) {
    double tmp = 0.0;
    double *a = new double[N]; //массив для вектора=A*y
    double norma1 = 0.0, norma2 = 0.0;
    for (int i = 0; i < N; ++i) {
        for (int j = 1; j < N+1; j++) {
            tmp += matr[j - 1 + i * N] * Y(m, j);
        }
        a[i] = tmp;
        tmp = 0.0;
    }
    for (int i = 1; i < N+1; i++) {
        norma1 += (a[i-1] + L(m)*Y(m, i))*(a[i-1] + L(m)*Y(m, i)); //||AY+LY||
        norma2 += (L(m)*Y(m, i))*(L(m)*Y(m, i)); //||LY||
    }
    delete[]a;
    return sqrt(norma1 / norma2);
}

