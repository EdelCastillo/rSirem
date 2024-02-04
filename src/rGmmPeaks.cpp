/*************************************************************************
 *     rSirem - R package for MSI data deconvolution
 *     Copyright (C) november 2023, Esteban del Castillo Pérez
 *     esteban.delcastillo@urv.cat
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

#include <Rcpp.h>
#include "siremPeaks.h"
#include "gmmPeaks.h"

using namespace Rcpp;

//From the information on concentration peaks and sirem
//concentration information is deconvolved.
//Arguments:
//data: (adapted for the output of siremPeaksR()).
// "magnitudes"->concentration for each scan.
// "magnitudePeaks" -> concentration peaks.
// "siremPeaks" -> sirem peaks.
// "unitedMagnitudePeaks" -> contiguous magnitude peaks.
// "unitedSiremPeaks" -> siren peaks inside each concentration peak.
//magnitudes:
// magnitudes to adjust in the deconvolution.
//minMeanPxMag:
// minimum averaged magnitude to consider a peak valid.
//mzAxis:
// mass vector if you want the result in Daltons. If not given, the result is in scans.
//Return:
// Matrix with info of the Gaussians resulting from deconvolution. columns: mean, sigma and value.

// [[Rcpp::export]]
NumericMatrix rGmmPeaks(List data, NumericVector magnitudes, float minMeanPxMag, NumericVector mzAxis)
{
  NumericMatrix retFail(2,2); //used for return with error.
  retFail(1,1)=-1;
//  List retFail=List::create(-1);

  //The information passed from R is adapted to C++.
  NumericVector meanMagVec=data["scanMean"];
  NumericMatrix magPeaksMx=data["magnitudePeaks"];
  NumericMatrix siremPeaksMx=data["siremPeaks"];
  NumericMatrix uMagPeaksMx=data["unitedMagnitudePeaks"];
  NumericMatrix uSiremPeaksMx=data["unitedSiremPeaks"];

  float*magVec_p=0; //hosts the information on magnitudes (concentrations).
  float*meanMagVec_p=0; //hosts the averaged information of magnitudes (concentrations).
  int nMagPeaks=magPeaksMx.nrow(); //number of magnitude peaks.
  CONVOLVED_PEAKS  *cnvPeaks_p=0; //hosts all the peak information.
  int nMeanMagVec=meanMagVec.size(); //vector size.
  int nMagVec=magnitudes.size(); //concentration data size.
   
  try
  {

    if(nMeanMagVec>0 && nMagVec==0)
    {
    meanMagVec_p=new float[nMeanMagVec]; //memory reservation.
    for(int i=0; i<nMeanMagVec; i++)
      meanMagVec_p[i]=(float)meanMagVec[i]; //copy data (mean of concentrations).
       
    }

    if(nMagVec>0)
    {
    magVec_p=new float[nMagVec]; //memory reservation.
    for(int i=0; i<nMagVec; i++)
      magVec_p[i]=(float)magnitudes[i]; //copy data (concentration).
    }
    
    if(nMagPeaks>0)
    {
    cnvPeaks_p=new CONVOLVED_PEAKS[nMagPeaks]; //memory reservation.
    for(int mp=0; mp<nMagPeaks; mp++) //initialization.
      cnvPeaks_p[mp].etpPeaks_p=0;
    }
  }
  catch(const std::bad_alloc& e)
  {
    printf("Error reserving memory: %s\n",e.what());
    return retFail;
  }

  //  int nSiremPeaks=uSiremPeaksMx.nrow(); //number of sirem peaks.
  if(uSiremPeaksMx.nrow() != nMagPeaks) //agreement analysis.
    {printf("Data size consistency error.\n"); return retFail;}

  //movement of data from the passed list to the CONVOLVED_PEAKS structure.
  for(int mp=0; mp<nMagPeaks; mp++)
  {
    cnvPeaks_p[mp].magPeaks.low =magPeaksMx(mp,0); //magnitude peak.
    cnvPeaks_p[mp].magPeaks.max =magPeaksMx(mp,1);
    cnvPeaks_p[mp].magPeaks.high=magPeaksMx(mp,2);

    int A=uSiremPeaksMx(mp, 0); //low  sirem peak index.
    int B=uSiremPeaksMx(mp, 1); //high sirem peak index.
    int nSiremP=B-A+1; //number of sirem peaks inside the magnitude peak.
    try
      {
      cnvPeaks_p[mp].etpPeaks_p=new Peaks::ION_INDEX[nSiremP]; //memory reservation.
      }
    catch(const std::bad_alloc& e)
      {
        printf("Error reserving memory: %s\n",e.what());
        return retFail;
      }

    //sirem peaks interior to the magnitude peak.
    for(int sirP=A, k=0; sirP<=B; sirP++, k++)
    {
      cnvPeaks_p[mp].etpPeaks_p[k].low =siremPeaksMx(sirP, 0); //low peak index.
      cnvPeaks_p[mp].etpPeaks_p[k].max =siremPeaksMx(sirP, 1); //maximum peak index.
      cnvPeaks_p[mp].etpPeaks_p[k].high=siremPeaksMx(sirP, 2); //high peak index.
    }
    cnvPeaks_p[mp].nEtpPeaks=nSiremP; //number of sirem peaks interior to the magnitude peak.
  }
  //class for deconvolution.
  GmmPeaks gmmPeaks(minMeanPxMag);

  if(magnitudes.length()==0) //if no magnitude info is passed, the average magnitude is used.
  {
    GROUP_F pxMag;
    pxMag.set=meanMagVec_p;
    pxMag.size=nMeanMagVec;
    gmmPeaks.setMagnitudes(&pxMag);
  }

  GAUSSIAN deconvIn, deconvOut, *deconv_p; //if it is required to adapt the mass axis.
  deconv_p=&deconvIn;
  
  NumericMatrix ret(siremPeaksMx.nrow(), 3); //the results matrix is created and sized.
  colnames(ret)=CharacterVector({"mean", "sigma", "value"}); //column names.
  int gaussianIndex=0; //indices for each deconvolved Gaussian.
  float *mzAxis_p=0; //mass axis
  for(int uPeak=0; uPeak<uMagPeaksMx.nrow(); uPeak++)
    {
    int mPeakLow   =uMagPeaksMx(uPeak, 0);    //lower simple peak.
    int mPeakHigh  =uMagPeaksMx(uPeak, 1);    //upper simple peak.
    int lowMzIndex =magPeaksMx(mPeakLow, 0);  //lower mass index.
    int highMzIndex=magPeaksMx(mPeakHigh, 2); //higher mass index.
    int mzSize=highMzIndex-lowMzIndex+1;

    int nPeaks=mPeakHigh-mPeakLow+1; //number of simple peaks.
    
    //The information of simple peaks of magnitude and entropy that make up the segment is established.
    //The class that generates the peaks is connected to the class that uses them for deconvolution.
    gmmPeaks.setPeaks(&cnvPeaks_p[mPeakLow], nPeaks); //pointer to the beginning of peaks to be treated.

    //The magnitude information that must be adjusted is established.
    if(nMagVec>0) //new magnitude (px or ROI)
      {
      GROUP_F pxMag;
      pxMag.set=magVec_p+lowMzIndex;
      pxMag.size=mzSize;
      gmmPeaks.setMagnitudes(&pxMag); //magnitude peaks.
      }
    else   //averaged magnitude.
      {
      GROUP_F pxMag;
      pxMag.set=meanMagVec_p+lowMzIndex;
      pxMag.size=mzSize;
      gmmPeaks.setMagnitudes(&pxMag); //magnitude peaks.
      }

    //The Gaussians are formed whose sum reproduces the magnitude peaks.
   if(gmmPeaks.gmmDeconvolution()<0) //deconvolution
      {
      printf("ERROR: limits exceeded. The aim is to deconvolve a peak composed of more than %d Gaussians.\n", DECONV_MAX_GAUSSIAN);
      return retFail;
      }
    bool unitConver=false; //true If a change of units should be made(scans to Daltons)

    if(mzAxis.size()>0) //if the vector for unit conversion exists.
    {
      unitConver=true;
      mzAxis_p=new float[mzAxis.size()];
      for(int i=0; i<mzAxis.size(); i++)
        mzAxis_p[i]=(float)mzAxis(i);
    }
    if(unitConver && (mzAxis.size()!=nMeanMagVec)) //if the dimensions are not correct.
      {printf("warning: the dimensions of mzAxis and data[magnitudes] must match.\n"); return retFail;}

    //info of Gaussians, adapted in magnitude, that make up the segment of simple magnitude peaks.
    for(int i=0; i<gmmPeaks.getDeconvNumber(); i++)
    {
      deconvIn=gmmPeaks.getDeconv(i); //get a gaussian.
      if(unitConver)
      {
        gmmPeaks.gaussConversion(&deconvIn, &deconvOut, mzAxis_p+lowMzIndex, nMeanMagVec); //unit conversion (scans to Daltos).
        deconv_p=&deconvOut;
      }
      else
      {
        deconv_p=&deconvIn;
      }
      //the Gaussians are returned.
      ret(gaussianIndex, 0)=(double)deconv_p->mean;
      ret(gaussianIndex, 1)=(double)deconv_p->sigma;
      ret(gaussianIndex, 2)=(double)deconv_p->yFactor*deconv_p->weight;
      gaussianIndex++;
  
    }
  }
//  return ret;

  //Reserved memory is released.
  if(meanMagVec_p) delete [] meanMagVec_p;
  if(magVec_p) delete [] magVec_p;
  if(cnvPeaks_p) 
  {
    for(int mp=0; mp<nMagPeaks; mp++)
      if(cnvPeaks_p[mp].etpPeaks_p) delete [] cnvPeaks_p[mp].etpPeaks_p;
    delete [] cnvPeaks_p;
  }
  if(mzAxis_p) delete [] mzAxis_p;

  return ret;
}


