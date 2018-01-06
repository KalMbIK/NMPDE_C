//
// Created by pusheen on 17.12.17.
//
//#define TYPE int
#include "myMatrices.h"
#define EPS 1E-6
#define CONST_C 1.
#define CONST_CFL 0.2
//#define EPS_ZERO(x) ((x) < EPS ? 0: x)

size_t N = 100;
size_t DIM = 3;
double x0=0, x1 = 1;
double T = .2;
double gFactor = 1.4;

double initFieldDensity(double x){
    if (x < .5)
        return 1;
    else
        return .125;
}

double initFieldVelocity(double x){
    return 0.;
}

double initFieldEnergy(double x){
    if (x < .5)
        return 2.5;
    else
        return .25;
}

double momentum(double ro, double v){
    return ro*v;
}

double pressure(double* state, size_t size){
    return (gFactor-1)*(state[2]-state[1]*state[1]/state[0]/2);
}

void flux(double* u, double p, double* res){
    res[0]=u[1];
    double uu = u[1]/u[0];
    res[1]=p+uu*u[1];
    res[2]=uu*(u[2]+p);
}

double** initCondition(double* grid, size_t size){
    double **tmp = newMatrix(size, DIM);
    double *v = newArray(size);
    vectorize(&initFieldVelocity, grid ,v, size);
    vectorize(&initFieldDensity, grid,tmp[0], size);
    vectorize2D(&momentum, tmp[0], v, tmp[1], size);
    vectorize(&initFieldEnergy, grid,tmp[2], size);
    double **res = newMatrixTransposed(tmp,size,DIM);
    free(v);
    freeMatrix(tmp);
    return res;
}

double speedOfSound(double ro, double p){
    return sqrt(gFactor*p/ro);
}

void solverF1(double step, double **u, double **v, size_t size){
    double lambda = 0;
    double tau = 0, t = 0;
    double *press = newArray(size);
    double *velocity = newArray(size);
    double *spOs = newArray(size);
    double *temp = newArray(size);
    double *alphas = newArray(size-1);
    double **fluxes = newMatrix(DIM,size);
    double **numFluxes = newMatrix(DIM,size-1);

    while (t < T){
        matrixReduce(&pressure,u,3,press,size);

        for (int i = 0; i < size; i++){
            velocity[i]=u[i][1]/u[i][0];
            spOs[i]=speedOfSound(u[i][0],press[i]);
        }
        vectorize(&fabs,velocity,temp,size);
        add(temp,spOs,temp,size);
        for(int i = 0; i < size-1; i++)
            alphas[i]=CONST_C*fmax(temp[i],temp[i+1]);

        lambda = getMaxElement(alphas,size-1)/CONST_C;
        tau=CONST_CFL*step/lambda;
        t += tau;

        for (int i = 0; i < size; i++)
            flux(u[i],press[i],fluxes[i]);

        for (int i = 0; i < size-1; i++){
            subtr(u[i],u[i+1],numFluxes[i],DIM);
            constMult(numFluxes[i],alphas[i],numFluxes[i],DIM);
            add(fluxes[i],numFluxes[i],numFluxes[i],DIM);
            add(fluxes[i+1],numFluxes[i],numFluxes[i],DIM);
            constMult(numFluxes[i],.5,numFluxes[i],DIM);
        }

        subtr(numFluxes[0],fluxes[0],v[0],DIM);
        for (int i = 1; i < size-1; i++){
            subtr(numFluxes[i],numFluxes[i-1],v[i],DIM);
        }
        subtr(fluxes[size-1],numFluxes[size-2],v[size-1],DIM);

        for (int i = 0; i < size; i++){
            constMult(v[i],tau/step,v[i],DIM);
            subtr(u[i],v[i],u[i],DIM);
        }
    }

    free(press);
    free(velocity);
    free(spOs);
    free(temp);
    free(alphas);
    freeMatrix(fluxes);
    freeMatrix(numFluxes);
}

int main(void){
    double * gr = newArray(N+1);
    linspace(x0,x1,gr,N+1);
    double h = gr[1]-gr[0];
    double ** u = initCondition(gr,N+1);
    double ** v = newMatrix(DIM, N+1);
    solverF1(h, u, v, N + 1);

    double ** toFile = newMatrixTransposed(u,DIM,N+1);
    vectorsToCsv("test.csv",N+1,DIM+1,gr,toFile[0],toFile[1],toFile[2]);
    freeMatrix(toFile);
    freeMatrix(u);
    free(gr);
    return 0;
}