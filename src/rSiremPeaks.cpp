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

using namespace Rcpp;

//The concentration information and positions of each pixel determine an image.
//From the image, and certain parameters, the measurement of sirem (derived from entropy) is determined.
//Then the peaks associated with the concentration and sirem information are determined.
//Arguments:
//data:
//  Matrix with as many columns (scans) as images to be analyzed and as many rows as pixels in the images.
//    info of concentration of each pixel on each scan
//params:
// algorithm         -> 0 = sirem; 1=entropy
// cutLeves          -> cut levels to apply to each image to generate binary images (percentiles).
// magSensitivity    -> variable sensitivity depending on the concentration level[0:maxConcentration] (bounded between limits).
// siremSensitivity  -> variable sensitivity depending on the sirem level[0:1] (bounded between limits).
// minMeanPxMag      -> minimum averaged concentration value to be considered.
// minSectionDensity -> minimum active pixel density value in the image to be considered.
// noiseLevel        -> absolute noise Level[0:maxConcentration]: lower values are considered null.
// pxCoord_X         -> X coordinates of each pixel.
// pxCoord_Y         -> Y coordinates of each pixel.
// referenceType     -> reference image type:
//                      0=obtained from concentration info (data);
//                      1=given in argument(scanReference);
//                      2=all pixels have significant concentration.
// scanReference     -> reference image to determine the area to consider.
// tileSide          -> number of side pixels of the square tile: 1(1x1), 2(2x2), 3(3x3), 4(4x4)
//return:
// "magnitudes"           -> concentration averaged over the data of each column (image).
// "sirem"                -> sirem measurement associated with each image.
// "magnitudePeaks"       -> concentration peaks: low, maximum, high (they are indices to the scans).
// "unitedMagnitudePeaks" -> united concentration peaks: low and high.
// "siremPeaks"           -> sirem peaks: low, maximum, high (they are indices to the scans).
// "unitedSiremPeaks"     -> sirem peaks inside each concentration peak.

// [[Rcpp::export]]
List rSiremPeaks(NumericMatrix data, NumericMatrix pxCoord, List params, NumericVector scanReference)
{
  List retFalse=List::create(Named("ret")=-1);
  int nScan=data.ncol();
  int nPixels=data.nrow();

  float *magnitudes_p=0;
  if(nScan<=0 || nPixels<=0)
    {printf("Warning: no data detected\n"); return retFalse;}
  
  //The passed information is copied to the SIREM_PEAKS_INFO structure.
  SIREM_PEAKS_INFO siremPeaksInfo;
  
  NumericVector v=params["algorithm"];
  siremPeaksInfo.algorithm=(char) v[0];
  
  v=params["minMeanPxMag"];
  siremPeaksInfo.minMeanPxMag=v[0];
  
  v=params["minSectionDensity"];
  siremPeaksInfo.minSectionDensity=v[0];
  
  v=params["noiseLevel"];
  siremPeaksInfo.noiseLevel=v[0];
  if(siremPeaksInfo.noiseLevel==0) siremPeaksInfo.noiseLevel=1e-6;
  if(siremPeaksInfo.noiseLevel<0) 
  {printf("warning: noiseLevel must be within range [0:maxConcentration]\n"); return retFalse;}
  
  siremPeaksInfo.nScans=nScan;
  
  v=params["tileSide"];
  siremPeaksInfo.tileSide=(char)v[0];
  if(siremPeaksInfo.tileSide>4 || siremPeaksInfo.tileSide<1)
  {printf("warning: tileSide must be within range [1:4]\n");}
  
  v=params["referenceType"];
  siremPeaksInfo.referenceType=(char)v[0];
  if(siremPeaksInfo.referenceType>2 || siremPeaksInfo.referenceType<0)
  {printf("warning: referenceType must be within range [0:2]\n"); return retFalse;}
  
  //cutting levels.
  NumericVector cutLevels = as<NumericVector>(params["cutLevels"]);
  if(cutLevels.length()<=0)
  {printf("warning: cutLevels must be greater than zero\n"); return retFalse;}
  
  siremPeaksInfo.cutLevels.set=0;
  siremPeaksInfo.cutLevels.set=new float[cutLevels.length()];
  for(int i=0; i<cutLevels.length(); i++)
    siremPeaksInfo.cutLevels.set[i]=cutLevels[i];
  siremPeaksInfo.cutLevels.size= cutLevels.length();
  
  NumericVector etpSen = as<NumericVector>(params["siremSensitivity"]);
  if(etpSen[1]>1 || etpSen[0]<0)
  {printf("warning: siremSensitivity must be within range [0:1]\n"); return retFalse;}
  siremPeaksInfo.etpSenMax   = etpSen[1];
  siremPeaksInfo.etpSenMin   = etpSen[0];
  
  NumericVector magSen = as<NumericVector>(params["magSensitivity"]);
  if(magSen[0]<0)
  {printf("warning: magSensitivity must be within range [0:maxConcentration]\n"); return retFalse;}
  siremPeaksInfo.magSenMax   = magSen[1];
  siremPeaksInfo.magSenMin   = magSen[0];
  
  try
  {
    //The reference concentration to be used is analyzed to determine the useful area.
    if(siremPeaksInfo.referenceType==0) //reference is maximun in data
    {
      //An array is generated with the maximum value of each row (pixel).
     if(nPixels>0)
        siremPeaksInfo.mzInfo.set=new float[nPixels];
      double maxValue;
      for(int i=0; i<nPixels; i++)
      {
        maxValue=0;
        for(int scan=0; scan<nScan; scan++)
        {
          if(data(i,scan)>maxValue) maxValue=data(i,scan); //value maximun
        }
        if(maxValue>(double)siremPeaksInfo.noiseLevel) //if it overcomes the noise.
          siremPeaksInfo.mzInfo.set[i]=(float)maxValue;
        else
          siremPeaksInfo.mzInfo.set[i]=0.0; //It is canceled if it does not exceed the noise.
      }
      siremPeaksInfo.mzInfo.size=nPixels;
    }
    else if(siremPeaksInfo.referenceType==1) //the reference is given as an argument.
    {
    siremPeaksInfo.mzInfo.set=0;
    if(scanReference.length()>0)
      {
      siremPeaksInfo.mzInfo.set=new float[scanReference.length()];
      for(int i=0; i<scanReference.length(); i++)
          siremPeaksInfo.mzInfo.set[i]=scanReference[i];
      }
      siremPeaksInfo.mzInfo.size=scanReference.length();
    }
    else //reference is maximun possible
    {
      siremPeaksInfo.mzInfo.set=0;
      siremPeaksInfo.mzInfo.size=0;
    }

    //copy of pixel coordinate info.
    if(pxCoord.ncol()<=0 || pxCoord.nrow()<=0)
      {printf("warning: pxCoord must be greater than zero\n"); return retFalse;}
    siremPeaksInfo.pxCoord.set=0;
    if(nPixels != pxCoord.nrow())
        {
        printf("ERROR: lack of consistency in size between data elements and params.\n");
        return (retFalse);
        }
    if(nPixels>0)
      {
        siremPeaksInfo.pxCoord.set=new PIXEL_XY[nPixels];
        for(int i=0; i<nPixels; i++)
            {
            siremPeaksInfo.pxCoord.set[i].x=pxCoord(i, 0);
            siremPeaksInfo.pxCoord.set[i].y=pxCoord(i, 1);
            }
      }
      siremPeaksInfo.pxCoord.size=nPixels;
      
      magnitudes_p=new float[nPixels]; //maintains the concentrations of each scan.
    }

    catch(const std::bad_alloc& e)
      {
        printf("Error reserving memory: %s\n",e.what());
        return retFalse;
      }

      SiremPeaks siremPeaks(&siremPeaksInfo);//SIREM/Entropy & peaks class
    GROUP_F nextScan;
//  printf("nScan:%d px:%d\n", nScan, nPixels);
    
    //The images to be analyzed are launched (concentrations of one scan) and sirem is obtained.
    for(int scan=0; scan<nScan; scan++)
        {
        for(int i=0; i<nPixels; i++)
          {
            magnitudes_p[i]=(float)data(i,scan);
          }
        nextScan.set=magnitudes_p;
        nextScan.size=nPixels;
        siremPeaks.newMz(&nextScan); //a new image for analysis. get sirem measures
        }

    if(siremPeaks.getPeaksList()<0)  //get peaks info
    {printf("Warning: there are no peaks in data.\n"); return  retFalse;}
    
    int ncPeaks=siremPeaks.getCompoundPeakNumber(); //if there are no simple peak segments.
    if(ncPeaks<=0)
    {printf("Warning: there are no peaks in data.\n"); return retFalse;}

    //info on magnitudes averaged over each scan.
    NumericVector mag(nScan);
    for(int i=0; i<nScan; i++) //return average magnitudes
      mag[i]=siremPeaks.getMeanMagnitude(i);

    //sirem matrix: rows=scan; columns=cut level.
    NumericMatrix siremMx(nScan, siremPeaksInfo.cutLevels.size); //dimensioned.
    for(int cut=0; cut<siremPeaksInfo.cutLevels.size; cut++)
    {
      float *cut_p=siremPeaks.getSiremCut(cut);
      for(int i=0; i<nScan; i++)
        siremMx(i,cut)=cut_p[i]; //rows=scan; columns=cut level.
    }

    //info on composite magnitude peaks (united by their valleys).
    UNITED_PEAKS uPeaks; //peaks composed of simple concentration peaks.
    NumericMatrix cMagPeaks(ncPeaks, 2);
    for(int cp=0; cp<ncPeaks; cp++) //for each segment composed of simple magnitude peaks.
    {
      uPeaks=siremPeaks.getCompoundPeak(cp);
      cMagPeaks(cp,0)=uPeaks.peakLow; //index to the first magnitude peak.
      cMagPeaks(cp,1)=uPeaks.peakHigh;//index to the last magnitude peak.
    }

    //size info for variable sizing.
    Peaks::ION_INDEX peak;
    int nMagPeaks=siremPeaks.getSinglePeaksNumber(); //number of magnitude peaks.
    int nEtpPeaks=0;
    int nEtp=0;
    for(int i=0; i<nMagPeaks; i++)
      nEtpPeaks+=siremPeaks.getEtpPeaksNumber(i); //number of sirem peaks.
    
    //info of magnitude peaks and sirem peaks associated with each magnitude peak.
    NumericMatrix sPeaks(nEtpPeaks, 3), cSiremPeaks(nMagPeaks, 2);
    NumericMatrix magPeaks(nMagPeaks, 3);
    for(int i=0, k=0; i<nMagPeaks; i++)
      {
      peak=siremPeaks.getSingleMagPeak(i);
      magPeaks(i,0)=peak.low;  //low to peak magnitude ratio.
      magPeaks(i,1)=peak.max;  //index to maximum peak magnitude.
      magPeaks(i,2)=peak.high; //high ratio to peak magnitude.
      nEtp=siremPeaks.getEtpPeaksNumber(i);
      cSiremPeaks(i, 0)=k;        //sirem peaks low  associated with this magnitude peak
      cSiremPeaks(i, 1)=k+nEtp-1; //sirem peaks high associated with this magnitude peak
      for(int etp=0; etp<nEtp; etp++)
        {
        peak=siremPeaks.getEtpSinglePeak(i, etp);
        sPeaks(k+etp, 0)=peak.low;  //index to sirem low to peak.
        sPeaks(k+etp, 1)=peak.max;  //index to maximum sirem magnitude.
        sPeaks(k+etp, 2)=peak.high; //index to sirem high to peak.
        }
      k+=nEtp;
      }

    //returned list
    List ret=List::create(Named("scanMean")=mag, Named("sirem")=siremMx, Named("magnitudePeaks")=magPeaks,
                          Named("unitedMagnitudePeaks")=cMagPeaks, Named("siremPeaks")=sPeaks, Named("unitedSiremPeaks")=cSiremPeaks);

    //Reserved memory is released.
    if(siremPeaksInfo.cutLevels.set) delete [] siremPeaksInfo.cutLevels.set;
    if(siremPeaksInfo.pxCoord.set) delete [] siremPeaksInfo.pxCoord.set;
    if(magnitudes_p) delete [] magnitudes_p;
    return ret;
}

