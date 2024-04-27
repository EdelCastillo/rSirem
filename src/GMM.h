/***********************************************
 *   Copyright (C) 2018 by Esteban del Castillo
 *   esteban.delcastillo@urv.cat
 *
 *   Project
 *      December 2018
************************************************/
//GMM Gaussian Mixture Model
//Decomposes a peak into sum of Gaussians

#ifndef MALDI_GMM_H
#define MALDI_GMM_H

#include <stdio.h>
#include <math.h>

#define DECONV_MAX_GAUSSIAN 100

typedef struct{                 //deconvolution overview for a compound magnitude peak
    float
                mean,           //mean of this Gaussian
                sigma,          //standard deviation of this Gaussian
                weight,         //relative weight of this Gaussian
                yFactor;        //factor to determine the final magnitude of the Gaussian
}GAUSSIAN;

typedef struct      //GMM parameters
    {
    int maxIterations,  //maximum iterations of the algorithm.
        maxGaussians;   //maximum supported Gaussians.
    double
        relChange;      //relative change in likelihood between iterations of the algorithm
    }GMM_PARAMS;

typedef struct      //gaussian parameters
{
    double  weight,
            mean,
            sigma;
}GAUSS_PARAMS;

typedef struct  //gaussian parameters array
    {
    int size;
    GAUSS_PARAMS params[100];
    }GAUSS_GROUP;

typedef struct  //limits on variations
{
    double  minWeight,
            maxWeight,
            minMean,
            maxMean,
            minSigma,
            maxSigma;
}GAUSS_LIMITS;

typedef struct
{
    float   *x,                         //X coordinates
            *y;                         //Y coordinates
    int     size,                       //array size
            maxIter,                    //maximum number of iterations
            nGauss;                     //resulting Gaussians
    char    status[DECONV_MAX_GAUSSIAN];//resulting status
    double  relChange;                  //relative change in likelihood between iterations of the algorithm
    double  yFactor,                    //factor to recover the original magnitude
            quality;                    //fit quality

    GAUSS_PARAMS params[DECONV_MAX_GAUSSIAN]; //initial parameters
    GAUSS_LIMITS limits[DECONV_MAX_GAUSSIAN]; //limits on variations
}GMM_STRUCT;

//class for the implementation of the Gauss Mixture (GM) algorithm.
class Cgmm
  {
  public:
    //constructor
    Cgmm();

     //Establishes the decomposition of a peak magnitude into sum of Gaussians using the EM algorithm
     //(EM -> Expectation-Maximization)
     //Arguments:
     // Data in X/Y (peak magnitude), number of Gaussians for decomposition
     // initial parameters for each Gaussian, maximum number of iterations
     // convergence value = relative change in likelihood between iterations of the algorithm
     //return
     // the results in the input gmm structure itself.
     // number of iterations performed
    int gmm(GMM_STRUCT *gmm);

  };

#endif

