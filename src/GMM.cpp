/*************************************************************************
 *     Deconvolution using GMM algorithm
 *     Copyright (C) Esteban del Castillo Pérez
 *     esteban.delcastillo@urv.cat
 *     abril 2023
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **************************************************************************/
#include "GMM.h"

//Constructor
Cgmm::Cgmm(){ }


//Gauss Mixture Model: applies the EM algorithm with Gaussians
//Arguments:
// Data in X/Y (peak magnitude), number of Gaussians for decomposition
// initial parameters for each Gaussian, maximum number of iterations
// convergence value = relative change in likelihood between iterations of the algorithm
//return
// number of iterations performed
int Cgmm::gmm(GMM_STRUCT *gmm)
    {
    int    iterCount=0;         //iteration limit
    double sumaPxk, sumaY=0, sumaTotalPxk;    //accumulated
    double L=1.0, Lold;         //likelihood

    //memory to hold partial results
    double newMean[DECONV_MAX_GAUSSIAN], newSigma[DECONV_MAX_GAUSSIAN], newWeight[DECONV_MAX_GAUSSIAN];
    double  *pxk_p[DECONV_MAX_GAUSSIAN];

    //memory reservation and initialization
    double *sumaPxp_p=new double[gmm->size];
    for(int g=0; g<gmm->nGauss; g++)
        {
        pxk_p[g]=new double[gmm->size];
        gmm->status[g]=0;
        }

    for(int iy=0; iy<gmm->size; iy++) //accumulated in magnitude
        sumaY+=gmm->y[iy];

    while(true) //implementation of the EM algorithm (expectation-maximization)
    {
    //Espectation: probability that the variable x belongs to the Gaussian 'g'
    Lold=L;
    L=1.0;
    sumaTotalPxk=0.0;
    for(int ix=0; ix<gmm->size; ix++)
        {
        sumaPxp_p[ix]=0;
        sumaPxk=0.0;
        for(int g=0; g<gmm->nGauss; g++)
            {
            double A=(gmm->x[ix]-gmm->params[g].mean)/gmm->params[g].sigma;
            pxk_p[g][ix]=gmm->params[g].weight*exp(-A*A*0.5)/(2.506628*gmm->params[g].sigma); //likelihod
            sumaPxk+=pxk_p[g][ix]; //accumulated likelihood for each Gaussian over each mz
            }
        sumaPxp_p[ix]=sumaPxk;  //note for later use
        L*=sumaPxk; //likelihood
        sumaTotalPxk+=sumaPxk; //total accumulated
        for(int g=0; g<gmm->nGauss; g++)
            if(sumaPxk==0) pxk_p[g][ix]=1e-30; else pxk_p[g][ix]/=sumaPxk;   //likelihood->probability normalization
        }

    //Maximization in EM algorithm
    //new weights for each Gaussian
    for(int g=0; g<gmm->nGauss; g++)
        {
        newWeight[g]=0;
        for(int ix=0; ix<gmm->size; ix++)
            {
            newWeight[g]+=pxk_p[g][ix]*gmm->y[ix];
            }
        }

    //new mean for each Gaussian
    for(int g=0; g<gmm->nGauss; g++)
        {
        newMean[g]=0;
        for(int ix=0; ix<gmm->size; ix++)
            {
            newMean[g]+=pxk_p[g][ix]*gmm->y[ix]*gmm->x[ix];
            }
        newMean[g]/=newWeight[g];
        }

    //new sigma for each Gaussian
    for(int g=0; g<gmm->nGauss; g++)
        {
        newSigma[g]=0;
        for(int ix=0; ix<gmm->size; ix++)
            {
            newSigma[g]+=pxk_p[g][ix]*gmm->y[ix]*(gmm->x[ix]-newMean[g])*(gmm->x[ix]-newMean[g]);
            }
        if(newWeight[g]==0) newSigma[g]=gmm->limits[g].maxSigma;
        else {newSigma[g]/=newWeight[g]; newSigma[g]=sqrt(newSigma[g]);}
        }

    //normalization of weights
    for(int g=0; g<gmm->nGauss; g++)
        if(sumaY>0) newWeight[g]/=sumaY;
        else newWeight[g]=0; //gmm->limits[g].maxWeight;

    //new values become old ones
    for(int g=0; g<gmm->nGauss; g++)
        {
        gmm->params[g].weight=newWeight[g];
        gmm->params[g].mean  =newMean  [g];
        gmm->params[g].sigma =newSigma [g];
        }

    //limit control
    //it is noted if any limit is reached
    for(int g=0; g<gmm->nGauss; g++) //para cada gaussiana
        {
        if(gmm->params[g].weight<gmm->limits[g].minWeight)
          {gmm->params[g].weight=gmm->limits[g].minWeight; gmm->status[g]|=0x01;} //límite inferior en peso alcanzado
        if(gmm->params[g].weight>gmm->limits[g].maxWeight)
          {gmm->params[g].weight=gmm->limits[g].maxWeight; gmm->status[g]|=0x02;} //límite superior en peso alcanzado
        if(gmm->params[g].mean<gmm->limits[g].minMean)
          {gmm->params[g].mean=gmm->limits[g].minMean;     gmm->status[g]|=0x04;} //límite inferior en media alcanzado
        else if(gmm->params[g].mean>gmm->limits[g].maxMean)
          {gmm->params[g].mean=gmm->limits[g].maxMean;     gmm->status[g]|=0x08;} //límite superior en media alcanzado
        if(gmm->params[g].sigma<gmm->limits[g].minSigma)
          {gmm->params[g].sigma=gmm->limits[g].minSigma;   gmm->status[g]|=0x10;} //límite inferior en sigma alcanzado
        else if(gmm->params[g].sigma>gmm->limits[g].maxSigma)
          {gmm->params[g].sigma=gmm->limits[g].maxSigma;   gmm->status[g]|=0x20;} //límite superior en sigma alcanzado
        }

     //likelihod control for convergence of the algorithm
    double rel=(L-Lold)/L;
    if(iterCount>1 && fabs(rel)<gmm->relChange)
        break;

    if(++iterCount>=gmm->maxIter) break;
    }

    //factor to magnify the Gaussians obtained
    if(sumaTotalPxk==0)  sumaTotalPxk=1e-30; //casi cero
    gmm->yFactor=sumaY/sumaTotalPxk; //factor buscado

    //error promedio
    double err, qErr=0;
    for(int ix=0; ix<gmm->size; ix++)
        {
        err=gmm->y[ix]-sumaPxp_p[ix]*gmm->yFactor;
        qErr+=err*err; //quadratic error
        }
    qErr=sqrt(qErr/gmm->size);  //root mean square error
    gmm->quality=1.0-qErr/sumaY;  //relative fit quality in magnitude

    //free reserved memory
    if(sumaPxp_p) delete [] sumaPxp_p;
    for(int g=0; g<gmm->nGauss; g++)
        if(pxk_p[g]) delete []pxk_p[g];

  return iterCount; //returns the number of iterations performed
  }


