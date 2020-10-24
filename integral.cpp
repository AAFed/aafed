#include <iostream>
#include <math.h>

using namespace std;

FILE*out;
double aa, bb;
int N, flag;

double func(double x);
double formula(double a, double b, double x1, double x2, double x3);
double max_derivative(double a, double b);
double derivative(double x);
double integral1(double x);
void result1();

double fxy(double x, double y);
double arg(int i, double c);
double I2(int i, int j);
double integral2(double x);
void result2();


int main(void){
    cout<< "Enter mode: 0 [dim=1], 1 [dim = 2]; Enter N:" << endl;
    cin >> flag >> N; //1000 10 100
    //aa = 0.5; bb = 1.5;
    //aa = - M_PI; bb = M_PI;
    aa = 0.; bb = 1.;
    //aa = -1.;
    //bb = 1.;
    if(flag == 0) result1();
    if(flag == 1) result2();
    
}


double func(double x){
    return pow(x,10);
    //return exp(x);
    //return sin(2. * x) * sin(5. * x);
    //return sin(2. * x) * sin(2. * x);
   //return pow(sin((pow(1.10, 3)) * x), 2);
    //return 1./(sqrt(1 - x*x));
}

double formula(double a, double b, double x1, double x2, double x3) {
    return (1./18.)*(b - a)*(5.* func(x1) + 8.*func(x2) + 5.*func(x3));
}

double integral1(double x){
    return (1.0 / 11.0) * pow(x, 11);
    //return exp(x);
    //return (1.0 / 42.0) * (7*sin(3.*x) - 3. * sin(7. * x));
   // return (1./8.)*(4*x - sin(4.*x));
    //return (1./ 2.) * x - 0.187829 * sin(2.662 * x);
    //return asin(x);
}

double derivative(double x){
    return 151200.0 * pow(x, 4);
    //return exp(x);
    //return 58460.0 * cos(2. * x) * cos(5. * x) - 59189.0 * sin(2. * x) * sin(5. * x);
    //return 2048*cos(4.*x);
    //return 177.917 * pow(cos(1.331 * x), 2) - 177.917 * pow(sin(1.331 * x), 2);
    //return 4725.*x*x/pow(sqrt(1. - x*x), 9) + 225./pow(sqrt(1. - x*x), 7) + 10395.*x*x*x*x*x*x/pow(sqrt(1. - x*x), 13) + 14175.*x*x*x*x/pow(sqrt(1. - x*x), 11);
}

double max_derivative(double a, double b){
    int n = 1000;
    double max = 0.0;
    for (int i = 0; i < n; i++){
        if (fabs(derivative(a + ((double)i / (double)n) * (b - a))) > max){
            max = fabs(derivative(a + (double)i / (double)n) * (b - a));
        }
    }
    return max;
}

double fxy(double x, double y){
    //return x + exp(y);
    return x * exp(y);
}

double arg(int i, double c) {
    double h = (bb-aa)/(double)N;
    return aa+((double)i + c)*h;
}

double I2(int i, int j) {
    double A = fxy(arg(j, 0.5),arg(i, 0.));
    double B = fxy(arg(j, 1.), arg(i, 0.5));
    double C = fxy(arg(j, 0.5), arg(i, 0.5));
    double tmp1 = A + B + C;
    double D = fxy(arg(j, 0.),arg(i, 0.5));
    double F = fxy(arg(j, 0.5), arg(i, 1.));
    double G = fxy(arg(j, 0.5), arg(i, 0.5));
    double tmp2 = D + F + G;
    return  (1./6.)*(tmp1 + tmp2);
}
 

double integral2(double x) {
    double I = 0.;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            I += I2(i, j)*pow(x, 2);
    return I;
}

void result1(){
    out = fopen("integr1.txt", "wt");
    double Sn = 0., Rn = 0., If = 0., R, lnC, C, p, accuracy;
    double x1, x2, x3, ak, bk;
    for (int k = 0; k < N; ++k){
        ak = aa + ((bb - aa)/(double)N)*(double)k;
        bk = aa + ((bb - aa)/(double)N)*((double)k+1.);
        x1 = 0.5*(ak + bk) - 0.5*(bk - ak)*sqrt(3./5.);
        x2 = 0.5*(ak + bk);
        x3 = 0.5*(ak + bk) + 0.5*(bk - ak)*sqrt(3./5.);
        If += (integral1(bk) - integral1(ak));
        Sn += formula(ak, bk, x1, x2, x3);
        Rn = (1./737280.)*pow(bk - ak, 7)*max_derivative(ak, bk);
    }
    R = log(1./fabs(Sn - If));
    lnC = double(int(R) % int(log(N)));
    C = 1.0 / exp((lnC));
    p = (R - lnC)/log(N);
    accuracy = C/exp(p*log(N));
    fprintf(out, "a = %0.15f, b = %0.15f, N = %d\n\n", aa, bb, N);
    fprintf(out, "Sn = %0.30f\nIf = %0.30f\n|Sn - If| = %.5e\nRn = %.5e\n", Sn, If, fabs(Sn - If), Rn);
    fprintf(out, "|Sn - If| = %.30f ~ %.10f / %d^%.10f = %.5e\naccuracy = %.5e\n", fabs(Sn - If), C, N, p, accuracy, fabs(fabs(Sn - If) - accuracy));
    fclose(out);
}

void result2() {
    out = fopen("integr2.txt", "wt");
    double h = (bb - aa) / (double)N;
    double p = 0., real_I = 0., I = 0.;
    real_I = 0.5 * (exp(1.) - 1.);
    //real_I = exp(1.) - 1. / 2.;
    I = integral2(h);
    p = log(fabs(real_I - I)) / log(h);
    fprintf(out, "a = %f, b = %f, N = %d\n\n", aa, bb, N);
    fprintf(out, "I = %.20f\ntrue I = %.20f\n|I - true I| = %.10e = %f ^ %f\n", I, real_I, fabs(real_I - I), h, p);
    fclose(out);
}



