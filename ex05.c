//
// Created by pusheen on 29.11.17.
//
//#define TYPE int
#include "myNumPyV2.h"
#include <string.h>
#define EPS 1E-6
#define EPS_ZERO(x) ((x) < EPS ? 0: x)

size_t N = 1000;
double x0=0, x1 = 1;
double T = 1;

double u0(double x){
    return 1+sin(2*M_PI*x)/2;
}
double uu0(double x){
    if ((x<=.1)||(x>=.3))
        return 0;
    else
        return 1;
}

double alpha1(double x){
    return fabs(x);
}
double alpha2(double x){
    return 2*fabs(x);
}

double f(double u){
    return u*u/2;
}
double fprime(double u){
    return u;
}

double getTau(double *abar, double *alpha, double *temp, size_t size, double h){
    double divisor = (abar[size-1]+alpha[size-1])+(alpha[0]-abar[0]);
    if (divisor!=0){
        temp[0]=h/divisor;
    } else {
        temp[0]=1;
    }
    for (int i=1; i<size; i++){
        divisor=(abar[i-1]+alpha[i-1])+(alpha[i]-abar[i]);
        if (divisor!=0){
            temp[i]=h/divisor;
        } else {
            temp[i]=1;
        }
    }
    double res = getMinElement(temp, size);
    return 2*res;
}

void solver(double step, double* u, double * v, size_t size, double (*alphaFun)(double)){
    double * flow = newArray(size);
    double * alpha = newArray(size);
    double * abar = newArray(size);
    double * numFlow = newArray(size);
    double * temp = newArray(size);
    double tau = 0;
    double t = 0;
    double h = step;

    while (t < T){
        vectorize(&f,u,flow,size);
        //abar and alpha
//        if (fabs(u[0]-u[size-1])>EPS)
//            abar[0] = (flow[0]-flow[size-1])/(u[0]-u[size-1]);
//        else
//            abar[0] = fprime(u[size-1]);
        for (int i = 0; i < size-1; i++){
            if (fabs(u[i+1]-u[i])>EPS)
                abar[i] = (flow[i+1]-flow[i])/(u[i+1]-u[i]);
            else
                abar[i] = fprime(u[i]);
        }
        if (fabs(u[0]-u[size-1])>EPS)
            abar[size-1] = (flow[0]-flow[size-1])/(u[0]-u[size-1]);
        else
            abar[size-1] = fprime(u[size-1]);
        vectorize(alphaFun,abar,alpha,size);

        //numflow
        for (int i = 0; i < size-1; i++){
            numFlow[i]=(flow[i]+flow[i+1])/2-alpha[i]/2*(u[i+1]-u[i]);
        }
        numFlow[size-1]=(flow[size-1]+flow[0])/2-alpha[size-1]/2*(u[0]-u[size-1]);

        tau = getTau(abar, alpha, temp, size, step);
        t+=tau;

        //nextIteration
        v[0]=u[0]-tau/h*(numFlow[0]-numFlow[size-1]);
        for (int i = 1; i < size; i++){
            v[i]=u[i]-tau/h*(numFlow[i]-numFlow[i-1]);
        }

        double * tt = v;
        v = u;
        u = tt;
    }

    free(flow);
    free(numFlow);
    free(alpha);
    free(abar);
    free(temp);
}

int main(void){
    double h = (x1-x0)/(N);
    double * gr = newArray(N+1);
    double * u = newArray(N+1);
    double * v1 = newArray(N+1);
    double * v2 = newArray(N+1);
    linspace(x0,x1,gr,N+1);

    vectorize(&uu0,gr,u,N+1);
    solver(h,u,v1,N,&alpha1);
    vectorize(&uu0,gr,u,N+1);
    solver(h,u,v2,N,&alpha2);

    //filename, size of grid, number of vectors to print, ... - pointers to vectors
    vectorsToCsv("test"".csv",N,3,gr,v1,v2);

    free(gr);
    free(u);
    free(v1);
    free(v2);
    return 0;
}