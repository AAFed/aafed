#include <iostream>
#include <math.h>
#include <locale.h>
using namespace std;

int p = 1;
int N = 10;
double h = 1/(double)N;
double F(double x);
double Y(int m, int k);
double L(int m);
double Y_true(double x);
void Fourier(double *d);


int main() {
    FILE*out = fopen("outfurie.txt", "w");
    double result[N+1], d[N+1];
    double sum = 0.0;
    double max=0.0, min=1.0;
    for (int i = 0; i < N+1; i++){
        result[i] = 0.;
        d[i] = 0.;
    }

    int p1 = 0, q1 = 0;

    Fourier(d);

    for (int i = 1; i < N+1; i++) {
        for (int j = 1; j < N+1; j++) {
            result[i] += d[j]*Y(j, i);//b = A y
        }
    }

    cout << endl;

    cout << "Scalar: " << endl;

    for (int i = 1; i < N+1; i++) {
        for (int j =1; j < N; j++) {
            sum += h * Y(i, j)*Y(i, j);
        }
        sum += 0.5 * h * Y(i, N) * Y(i, N);
        cout << "(" << i << " , " << i << ")" << " = " << sum << endl;

        if (sum > max) {
            max = sum;
            p1 = i;
        }
        if (sum < min) {
            min = sum;
            q1 = i;
        }
        sum = 0.0;
    }

    cout << "Max scalar: " << "(" << p1 << " , " << p1 << ")" << " = " << max << endl;
    cout << "Min scalar: " << "(" << q1 << " , " << q1 << ")" << " = " << min << endl;
    for(int i = 0; i < N+1; i++) fprintf(out, "%.10f\t%.10f\t%.10f\n", i*h, result[i], Y_true(i * h));
    for (int i = 1; i < N; i++) {
        sum += h * (Y_true(i*h) - result[i])*(Y_true(i*h) - result[i]);
    }
    sum += 0.5 * h * (Y_true(N*h) - result[N])*(Y_true(N*h) - result[N]);
    printf("||result - y_real||h = %0.10f\n", sqrt(sum));

    return 0;
}



double F(double x) {
    return (M_PI/2 * M_PI/2 )* sin(M_PI/2 *x) + sin(M_PI/2 *x);
}


double Y(int m, int k) {
    return sqrt(2)*sin(M_PI*(2 * m - 1) * k * h/2);
}

double L(int m) {
    return pow(2*sin(M_PI*(2 * m - 1) * h/ 4)/h,2);
}

double Y_true(double x) {
    return sin(M_PI/2 *x);
}

void Fourier(double *d) {
    for (int j = 1; j < N+1; j++) {
        for (int i = 1; i < N; i++) {
            d[j] += h * F(i*h)*Y(j, i);
        }
        d[j] += 0.5 * h * F(N*h)*Y(j, N);
        d[j] /=( L(j) + p);
    }
}

