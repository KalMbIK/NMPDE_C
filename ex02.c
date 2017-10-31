//#define TYPE int
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "myNumPy.h"

int N = 1001, M = 11;
double x0=0, x1 = 1;
double T = 0.2;
double a_diff = 1;

double u0(double x){
    return sin(M_PI*x);
}

double f2(double t){
    return exp(-M_PI*M_PI*t);
}

double exactSolution(double x){
    return u0(x)*f2(T);
}

void step(double r, double la, double a, double b, double cK,
          double cK1, double * f, double * d, double * u, double * v){
    f[0] = 0;
    d[0] = 0;
    //f - initialization
    for (int i = 1; i < N; i++){
        f[i] = r*u[i-1]+(1-la)*u[i]+r*u[i+1];
        if (i > 1)
            d[i] = (f[i]-a*d[i-1])/(b-a*cK);
        else
            d[i] = (f[i]-a*d[i-1])/(b-a*cK1);
    }
    f[N-1] = 0;
    d[N-1] = 0;
    v[N-1] = d[N-1];
    for (int i = N-2; i >= 0; i--)
        v[i] = d[i]-cK*u[i];

    double  * temp = u;
    u = v;
    v = temp;
}

double getError(double* a, double* b, int size){
    double * delta = subtr(a,b,size);
    double max = 0, temp = 0;
    for (int i = 0; i < size; i++){
        temp = fabs(delta[i]);
        if (temp > max)
            max = temp;
    }
    free(delta);
    return max;
}

int main() {
    double h = (x1-x0)/(N-1);
    double t = T/(M-1);
    double la = a_diff*t/h/h;
    double r = la/2;
    double * gr = linspace(x0,x1,N);
    double * u = vectorize(&u0,gr,N);
    double * exact = vectorize(&exactSolution,gr,N);
    double * v = newArray(N);
    double * f = newArray(N);
    double * d = newArray(N);
    //Consts initialization
    double a = -r;
    double b = 1+la;
    double c = -r;
    double cK1 = c/b;
    double cK = c/(b-a*c);

    for (int i = 0; i < M; i++){
        step(r,la,a,b,cK,cK1,f,d,u,v);
    }
//    printArray(u,N);
    printf("%lg\n",getError(u,exact,N));

    free(gr);
    free(u);
    free(v);
    free(f);
    free(d);
    free(exact);
    return 0;
}