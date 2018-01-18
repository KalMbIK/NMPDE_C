//
// Created by pusheen on 17.12.17.
//
//#define TYPE int
#include "myMatrices.h"
#define EPS 1E-6
#define CONST_C 1.
#define CONST_C_ROE 1.
#define CONST_Ala 1.
#define CONST_CFL 0.2
//#define EPS_ZERO(x) ((x) < EPS ? 0: x)

size_t N = 1000;
size_t DIM = 3;
double x0=0, x1 = 1;
double T = .2;
double gFactor = 1.4;

// The following functions define the initial field distribution
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

double speedOfSound(double ro, double p){
    return sqrt(gFactor*p/ro);
}

double velocity(double* state, size_t size){
    return state[1]/state[0];
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
    vectorize(&initFieldVelocity, grid, v, size);
    vectorize(&initFieldDensity, grid, tmp[0], size);
    vectorize2D(&momentum, tmp[0], v, tmp[1], size);
    vectorize(&initFieldEnergy, grid, tmp[2], size);
    double **res = newMatrixTransposed(tmp, size, DIM);
    free(v);
    freeMatrix(tmp);
    return res;
}

void numFlux1(double **u, double **fluxes, double **numFluxes, double *alphas, double *temp, size_t size){
    for (int i = 0; i < size-1; i++){
        alphas[i]=CONST_C*fmax(temp[i],temp[i+1]);
        subtr(u[i],u[i+1],numFluxes[i],DIM);
        constMult(numFluxes[i],alphas[i],numFluxes[i],DIM);
        add(fluxes[i],numFluxes[i],numFluxes[i],DIM);
        add(fluxes[i+1],numFluxes[i],numFluxes[i],DIM);
        constMult(numFluxes[i],.5,numFluxes[i],DIM);
    }
}
void myAvg(double *ro, double *q, double *res, size_t size){
    for (int i = 0; i < size-1; i++){
        res[i] = (ro[i]*q[i]+ro[i+1]*q[i+1])/(ro[i]+ro[i+1]);
    }
}
void fillEig(double* eig, double** R, double c){
    eig[0] = fabs(R[1][0]);
    eig[1] = fabs(R[1][1]);
    eig[2] = fabs(R[1][2]);
    //ENTROPY FIX
    double cc = CONST_Ala*c;
    for (int i = 0; i < DIM; i++)
        if (eig[i]<cc)
            eig[i]=(cc+eig[i]*eig[i]/cc)/2;
}
void fillR(double** R, double u, double c, double h){
    double t = u*c;
    R[0][0]=1;
    R[0][1]=1;
    R[0][2]=1;
    R[1][0]=u-c;
    R[1][1]=u;
    R[1][2]=u+c;
    R[2][0]=h-t;
    R[2][1]=u*u/2;
    R[2][2]=h+t;
}
void fillR_1(double** R, double u, double c, double h){
    double G = (gFactor-1)/c/c;
    double tt = u*u-h;
    double c_1 = 1/c;
    double uc_1 = u*c_1;
    double Gu = G*u;
    double R00 = 1+G*tt;
    R[0][0]=(R00+uc_1)/2.;
    R[2][0]=(R00-uc_1)/2.;
    R[0][1]=-(Gu+c_1)/2.;
    R[2][1]=-(Gu-c_1)/2.;
    R[0][2]=G/2;
    R[2][2]=R[0][2];
    R[1][0]=-G*tt;
    R[1][1]=Gu;
    R[1][2]=-G;
}

void numFlux2(double **u, double **fluxes, double **numFluxes, double *vel, double *press,
              double *sqRo, double* enthalpy,
              double *uHat, double* hHat, double* cHat, double **R, double **R_1, double *eig, double *temp, size_t size){
    for (int i = 0; i < size; i++){
        enthalpy[i] = (u[i][2]+press[i])/u[i][0];
        sqRo[i] = sqrt(u[i][0]);
    }

    myAvg(sqRo,vel,uHat,size);
    myAvg(sqRo,enthalpy,hHat,size);
    for(int i = 0; i < size-1; i++){
        cHat[i]=sqrt((gFactor-1)*(hHat[i]-1/2*uHat[i]*uHat[i]));
    }

    for (int i = 0; i < size-1; i++){
        //fill eig, R_1, R
        fillR(R,uHat[i],cHat[i],hHat[i]);
        fillR_1(R_1,uHat[i],cHat[i],hHat[i]);
        fillEig(eig,R, cHat[i]);
        subtr(u[i+1],u[i],temp,DIM);
        matrixVector(R_1,temp,numFluxes[i],DIM,DIM);
        diagMatrixVector(eig,numFluxes[i],temp,DIM);
        matrixVector(R,temp,numFluxes[i],DIM,DIM);

        constMult(numFluxes[i],-CONST_C_ROE,numFluxes[i],DIM);
        add(fluxes[i],numFluxes[i],numFluxes[i],DIM);
        add(fluxes[i+1],numFluxes[i],numFluxes[i],DIM);
        constMult(numFluxes[i],.5,numFluxes[i],DIM);
    }

}

void solver(double step, double **u, double **v, size_t size){
    double lambda = 0;
    double tau = 0, t = 0;
    double *press = newArray(size);
    double *vel = newArray(size);
    double *spOs = newArray(size);
    double *temp = newArray(size);
    double *alphas = newArray(size-1);
    double **fluxes = newMatrix(DIM,size);
    double **numFluxes = newMatrix(DIM,size-1);


    double *sqRo = newArray(size);
    double *enthalpy = newArray(size);
    double *uHat = newArray(size-1);
    double *hHat = newArray(size-1);
    double *cHat = newArray(size-1);

    double **R = newMatrix(DIM,DIM);
    double **R_1 = newMatrix(DIM,DIM);
    double *eig = newArray(DIM);


    while (t < T){
        matrixReduce(&pressure,u,DIM,press,size);
        matrixReduce(&velocity,u,DIM,vel,size);

        for (int i = 0; i < size; i++){
            spOs[i]=speedOfSound(u[i][0],press[i]);
            flux(u[i],press[i],fluxes[i]);
        }

        vectorize(&fabs,vel,temp,size);
        add(temp,spOs,temp,size);

        lambda = getMaxElement(temp,size);
        tau=CONST_CFL*step/lambda;
        t += tau;

        // COMPUTATION OF THE NUM_FLUXES
        // TO CHOOSE WHICH ONE OF THE FOLLOWING FLUXES TO USE, SIMPLY HIDE ANOTHER ONE IN THE COMMENTS
//        numFlux1(u,fluxes,numFluxes,alphas,temp,size);
        numFlux2(u,fluxes,numFluxes,vel,press,sqRo,enthalpy,uHat,hHat,cHat,R,R_1,eig,temp,size);

        // SCHEME CORE
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
    free(vel);
    free(spOs);
    free(temp);
    free(alphas);
    freeMatrix(fluxes);
    freeMatrix(numFluxes);


    free(sqRo);
    free(enthalpy);
    free(cHat);
    free(hHat);
    free(uHat);

    freeMatrix(R);
    freeMatrix(R_1);
    free(eig);

}

int main(void){
    double * gr = newArray(N+1);
    linspace(x0,x1,gr,N+1);
    double h = gr[1]-gr[0];
    double ** u = initCondition(gr,N+1);
    double ** v = newMatrix(DIM, N+1);
    solver(h, u, v, N + 1);
    double * pr = newArray(N+1);
    double * vel = newArray(N+1);
    matrixReduce(&pressure,u,DIM,pr,N+1);
    matrixReduce(&velocity,u,DIM,vel,N+1);

    double ** toFile = newMatrixTransposed(u,DIM,N+1);
    vectorsToCsv("test.csv",N+1,DIM+1,gr,toFile[0],vel,pr);

    freeMatrix(toFile);
    freeMatrix(u);
    free(gr);
    free(pr);
    free(vel);
    return 0;
}