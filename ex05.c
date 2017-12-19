//
// Created by pusheen on 17.12.17.
//
//#define TYPE int
#include "myNumPyV2.h"
#define EPS 1E-6
//#define EPS_ZERO(x) ((x) < EPS ? 0: x)

size_t N = 100;
double x0=0, x1 = 1;
double T = 1;

// InitCond
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

// Calculate the value of the new time step
double getTau(const double *abar, const double *alpha, double *temp, size_t size, double h){
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
    return 2*getMinElement(temp, size);
}

void updateValues(double* flux, double* abar, double* alpha, double* numFlux, double* u, size_t size, double (*alphaFun)(double)){
    //flux
    vectorize(&f,u,flux,size);
    //abar and alpha
    for (int i = 0; i < size-1; i++){
        if (fabs(u[i+1]-u[i])>EPS)
            abar[i] = (flux[i+1]-flux[i])/(u[i+1]-u[i]);
        else
            abar[i] = fprime(u[i]);
    }
    if (fabs(u[0]-u[size-1])>EPS)
        abar[size-1] = (flux[0]-flux[size-1])/(u[0]-u[size-1]);
    else
        abar[size-1] = fprime(u[size-1]);
    vectorize(alphaFun,abar,alpha,size);

    //numflux
    for (int i = 0; i < size-1; i++){
        numFlux[i]=(flux[i]+flux[i+1])/2-alpha[i]/2*(u[i+1]-u[i]);
    }
    numFlux[size-1]=(flux[size-1]+flux[0])/2-alpha[size-1]/2*(u[0]-u[size-1]);
}

void solver(double step, double* u, double * v, size_t size, double (*alphaFun)(double)){
    double * flux = newArray(size);
    double * numFlux = newArray(size);
    double * abar = newArray(size);
    double * alpha = newArray(size);
    double * temp = newArray(size);
    double tau = 0;
    double t = 0;

    while (t < T){
        updateValues(flux,abar,alpha,numFlux,u,size,alphaFun);
        //time
        tau = getTau(abar, alpha, temp, size, step);
        t+=tau;

        //nextIteration
        v[0]=u[0]-tau/step*(numFlux[0]-numFlux[size-1]);
        for (int i = 1; i < size; i++){
            v[i]=u[i]-tau/step*(numFlux[i]-numFlux[i-1]);
        }

        double * tt = v;
        v = u;
        u = tt;

        printf("%lg ",t);
    }

    free(flux);
    free(numFlux);
    free(abar);
    free(alpha);
    free(temp);
}

int main(void){
    double h = (x1-x0)/(N);
    double * gr = newArray(N+1);
    double * u = newArray(N+1);
    double * v1 = newArray(N+1);
    double * v2 = newArray(N+1);
    linspace(x0,x1,gr,N+1);

    // To test the first problem - use the u0
    // To perform the second test - use the uu0
    vectorize(&u0,gr,u,N+1);
    solver(h,u,v1,N,&alpha1);
    vectorize(&u0,gr,u,N+1);
    solver(h,u,v2,N,&alpha2);

    //filename, size of grid, number of vectors to print, ... - pointers to vectors
    vectorsToCsv("test.csv",N,3,gr,v1,v2);

    free(gr);
    free(u);
    free(v1);
    free(v2);
    return 0;
}