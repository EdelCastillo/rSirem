/*************************************************************************
 *     Deconvolution from sirem peaks
 *     Copyright (C) octubre 2023 Esteban del Castillo Pérez
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

#include "gmmPeaks.h"

//Constructor:
//Argument:
//minMeanPeakMagnitude -> minimum averaged magnitude for a peak to be deconvolved
GmmPeaks::GmmPeaks(float minMeanPeakMagnitude)
{
    m_minMeanPeakMagnitude=minMeanPeakMagnitude;
    m_deconv_p=0;
    m_nDeconv=0;
    m_sGmm.x=0;
    m_sGmm.y=0;
    m_mag_p=0;
}

//Destructor: frees reserved memory.
GmmPeaks::~GmmPeaks()
{
//printf("...ini GmmPeaks destructor\n");
    if(m_deconv_p)  delete [] m_deconv_p;
//printf("...end GmmPeaks destructor\n");

}

//magnitudes associated with all the peaks to be treated
//Arguments:
//mag_p -> structure with magnitude information associated with each scan
void GmmPeaks::setMagnitudes(GROUP_F *mag_p)
{
    m_mag_p=mag_p->set;
    m_nMag=mag_p->size;
}

//peaks to consider
//Arguments:
//cnvPeaks_p -> pointer to array of the set of simple peaks that constitute the composite peak
//nUnitedPeaks -> #simple peaks
void GmmPeaks::setPeaks(CONVOLVED_PEAKS *cnvPeaks_p, int nUnitedPeaks)
{
    m_cnvPeaks_p=cnvPeaks_p;
    m_nUnitedPeaks=nUnitedPeaks;

    if(m_sGmm.x) {delete []m_sGmm.x; m_sGmm.x=0;}
    if(m_sGmm.y) {delete []m_sGmm.y; m_sGmm.y=0;}
    int range=0;
    range=cnvPeaks_p[nUnitedPeaks-1].magPeaks.high-cnvPeaks_p[0].magPeaks.low+1; //#composite peak scans
    range+=2;   //to accommodate the extreme zeros
    m_sGmm.x=0;
    m_sGmm.y=0;
    try
        {
        m_sGmm.x=new float [range]; //mz
        m_sGmm.y=new float [range]; //magnitude
        }
    catch(const std::bad_alloc& e)
        {
        printf("Error reserving memory: %s\n",e.what());
        return;
        }

    m_sGmm.size=range;
    m_sGmm.maxIter=100;     //maximum iterations of the GMM algorithm
    m_sGmm.relChange=0.001; //minimum change between iterations to stop the algorithm
    m_sGmm.quality=0.9;     //fit quality
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
//prepare the necessary information before calling the GMM algorithm: initialize the m_sGmm structure
//Requires the reference to the composite peak to be deconvolved and the associated magnitudes
//Return:
//  -1 if the average magnitude within the peak does not reach a minimum (m_minMeanPeakMagnitude)
//  -2 if gaussian number>DECONV_MAX_GAUSSIAN
int GmmPeaks::iniGMM()
{
    int nGauss=0;
    m_etpHits=0;
    int etpLow, etpHigh, mzSize;
    double sumaMagEpeak, sumaMagMpeak;
    Peaks::ION_INDEX *magPeak_p, *etpPeak_p;

    int lowPeakIndex, highPeakIndex, lowMzIndex, highMzIndex;
    lowPeakIndex =0;  //reference to the lower simple peak
    highPeakIndex=m_nUnitedPeaks-1;//reference to the upper simple peak

    if(lowPeakIndex>highPeakIndex) return -1; //limit control

    lowMzIndex =m_cnvPeaks_p[lowPeakIndex] .magPeaks.low;  //lower mass index
    highMzIndex=m_cnvPeaks_p[highPeakIndex].magPeaks.high; //upper mass index
    mzSize=highMzIndex-lowMzIndex+1;

    if(lowMzIndex>=highMzIndex) return -1; //limit control

    //sum of the magnitudes within the magnitude peak
    sumaMagMpeak=0;
    float maxValue=0;
    for(int i=0; i< mzSize; i++)
        {
        if(m_mag_p[i] > maxValue) maxValue=m_mag_p[i];
        sumaMagMpeak+=m_mag_p[i];
        }
    if(maxValue < m_minMeanPeakMagnitude)
        {m_sGmm.quality=-1; nGauss=0; return -1;}

    //The array of magnitudes and the array of masses are created
    //A zero is artificially added to the data on each side of the peak: it improves the fit.
     m_sGmm.x[0]=0.0;//To the left of the peak
    m_sGmm.y[0]=0.0;
    int j=1;
    for(int i=0; i< mzSize; i++, j++)
        {
        m_sGmm.x[j]=i; //magnitude peak
        m_sGmm.y[j]=m_mag_p[i];
        }
    //A la derecha del pico
    m_sGmm.x[j]=highMzIndex+1;
    m_sGmm.y[j]=0.0;
    m_sGmm.size=mzSize+2; //data size increased by two
    nGauss=0; //num of gausianas
    int iCount=-1;
    m_nDeconv=0;
    
    //in each single magnitude peak there may be several entropy peaks, or none.
    //are made to fit within the limits of the peak magnitude
    for(int peak=lowPeakIndex; peak<=highPeakIndex; peak++)
        {
        magPeak_p=&m_cnvPeaks_p[peak].magPeaks; //pointer to simple peak magnitude info
        if(m_cnvPeaks_p[peak].nEtpPeaks>0) //if there is any entropy peak in the magnitude peak
            {
            //for each entropy peak the information is prepared for the GMM algorithm
            for(int i=0; i<m_cnvPeaks_p[peak].nEtpPeaks; i++)
                {
                iCount++; m_nDeconv++;
                etpPeak_p=&m_cnvPeaks_p[peak].etpPeaks_p[i]; //pointer to simple peak entropy info

                //the low and high limits of the entropy peak are adjusted if it is not completely contained in the magnitude peak
                etpLow=etpPeak_p->low;
                if(etpLow < magPeak_p->low) etpLow=magPeak_p->low;

                etpHigh=etpPeak_p->high;
                if(etpHigh > magPeak_p->high) etpHigh=magPeak_p->high;

                //sum of magnitudes at peak entropy
                sumaMagEpeak=0;
                for(int j=etpLow; j<=etpHigh; j++)
                    sumaMagEpeak+=m_mag_p[j-lowMzIndex];
                if(sumaMagEpeak<m_minMeanPeakMagnitude)  continue; //??
                m_etpHits|=(unsigned long)1<<iCount;

                //the parameters of the Gaussian are initialized:
                //The mean is made to coincide with the maximum value of the entropy peak
                m_sGmm.params[nGauss].mean=etpPeak_p->max-lowMzIndex;
                //Let sigma be half the width of the entropy peak.
                m_sGmm.params[nGauss].sigma=((float)etpHigh-etpLow)/2.0;
                //The weight is made to be the quotient between sum_of_each_peak_entropy and sum_of_peak_magnitude
                if(sumaMagMpeak==0) m_sGmm.params[nGauss].weight=0;
                else m_sGmm.params[nGauss].weight=sumaMagEpeak/sumaMagMpeak;

                //the limits of the Gaussian parameters are established:
                //The minimum value for the mean is the minimum value of the peak ntropy
                m_sGmm.limits[nGauss].minMean=etpLow-lowMzIndex-2;
                //The maximum value for the mean is the maximum value of the entropy peak
                m_sGmm.limits[nGauss].maxMean=etpHigh-lowMzIndex+2;
                //The minimum value for sigma is made to be half a scan
                float minSigma=0.5;
                m_sGmm.limits[nGauss].minSigma=minSigma;
                //The maximum value for sigma is the width of the entropy peak
                m_sGmm.limits[nGauss].maxSigma=(etpHigh-etpLow)*1.0;
                //Control of the minimum value for sigma
                if(m_sGmm.params[nGauss].sigma<=minSigma)  m_sGmm.params[nGauss].sigma=minSigma;
                //The minimum weight is zero
                m_sGmm.limits[nGauss].minWeight=0.0;
                //The maximum weight is one
                m_sGmm.limits[nGauss].maxWeight=1.0;
                nGauss++;
                if(nGauss>=DECONV_MAX_GAUSSIAN) //limit control
                  return -2;
                }
            }
        //If there are no entropy peaks within the magnitude peak, the magnitude info is used
        else //if(nGauss==0)
            {
            //the parameters of the Gaussian are initialized. Same as before only with magnitude values
            magPeak_p=&m_cnvPeaks_p[peak].magPeaks;
            m_sGmm.params[nGauss].mean=magPeak_p->max;
            m_sGmm.params[nGauss].sigma=((float)magPeak_p->high-magPeak_p->low)/2.0;
            m_sGmm.params[nGauss].weight=1; //sum_of_each_ion/sum_of_peak

            //the limits of the parameters of the Gaussian are established
            m_sGmm.limits[nGauss].minMean=magPeak_p->low;
            m_sGmm.limits[nGauss].maxMean=magPeak_p->high;
            float minSigma=1.0;
            m_sGmm.limits[nGauss].minSigma=minSigma;
            m_sGmm.limits[nGauss].maxSigma=(m_sGmm.limits[nGauss].maxMean-m_sGmm.limits[nGauss].minMean)*1.0;
            if(m_sGmm.params[nGauss].sigma<=minSigma)  m_sGmm.params[nGauss].sigma=minSigma;
            m_sGmm.limits[nGauss].minWeight=0;
            m_sGmm.limits[nGauss].maxWeight=1;
            m_etpHits=1;
            nGauss=1;
            m_nDeconv++;
            }
        }
    m_sGmm.maxIter=nGauss*40; //maximum iterations for adjustment in GMM
    return nGauss; //number of Gaussians that will make up the peak magnitude
}

//Deconvolution
//From a compound magnitude peak and its entropy peaks,
//obtains the Gaussians that compose it: EM Algorithm
//establishes the array m_deconv_p[] with m_nDeconv elements with info of each Gaussian ordered from lowest to highest mass
//if an element of m_deconv_p[] has its values null -> non-deconvolved entropy peak.
//If a peak could not be deconvolved, only one structure appears, its field nGauss=0 and its field quality=-1.
//The mean and sigma information of the Gaussians is in units of scans
//return:
//Number of gaussians
int GmmPeaks::gmmDeconvolution()
{
    Cgmm cGmm; //Class for GMM

    int nGauss=0;

    nGauss=iniGMM(); //structure initialization for GMM

    if(nGauss>0) //if it can be decomposed
        {
        //the Gaussians that make up the magnitude peak are constructed
        m_sGmm.nGauss=nGauss;
        cGmm.gmm(&m_sGmm);
        }
    else if(nGauss==-1)//if the average magnitude within the peak does not reach a minimum (m_minMeanPeakMagnitude)
        {
        m_sGmm.quality=-1;
        nGauss=0; //a structure appears in the file without deconvolution info
        m_nDeconv=0;
        }
    else   //if gaussian number>DECONV_MAX_GAUSSIAN;
        return -1;

    for(int g=0; g<nGauss; g++) //magnitude error control
        {
        if(isnan(m_sGmm.params[g].mean))      {m_sGmm.params[g].mean=-1;   m_sGmm.quality=-1;}
        if(isnan(m_sGmm.params[g].sigma))     {m_sGmm.params[g].sigma=-1;  m_sGmm.quality=-1;}
        if(isnan(m_sGmm.params[g].weight))    {m_sGmm.params[g].weight=-1; m_sGmm.quality=-1;}
        }

    if(isnan(m_sGmm.quality)) m_sGmm.quality=-1; //limit control

    if(m_deconv_p) {delete [] m_deconv_p; m_deconv_p=0;} //old memory is released
    try
        {
        m_deconv_p=new GAUSSIAN[m_nDeconv];   //new memory is reserved
        //m_deconv_p=new GAUSSIAN[nGauss];   //new memory is reserved
        }
    catch(const std::bad_alloc& e)
        {
        printf("Error reserving memory: %s\n",e.what());
        return -1;
        }
    //m_nDeconv=nGauss;

    int lowMzIndex, highMzIndex;
    lowMzIndex =m_cnvPeaks_p[0].magPeaks.low;       //index at lower mz
    highMzIndex=m_cnvPeaks_p[m_nUnitedPeaks-1].magPeaks.high;  //index at upper mz
    if(highMzIndex<=lowMzIndex)                     //limit control
            {lowMzIndex=0; highMzIndex =0;}

    //saves deconvolution info.
    //m_sGmm returns the constructed Gaussians but this function returns all of them; 
    //that is to say not all entropy peaks may generate Gaussians. Those who don't
    //they set all their values to zero
//    for(int g=0; g<nGauss; g++)
    for(int g=0, k=0; g<m_nDeconv; g++)
      {
      if((m_etpHits|=((unsigned long)1<<g))!=0)          //if the gaussian exist
        {
        m_deconv_p[g].mean=m_sGmm.params[k].mean;       //mean value
        m_deconv_p[g].sigma=m_sGmm.params[k].sigma;     //standard deviation
        m_deconv_p[g].weight=m_sGmm.params[k].weight;   //weighted weight
        m_deconv_p[g].yFactor=m_sGmm.yFactor;           //factor for reconstruction
        k++;
        }
      else
        {
        m_deconv_p[g].mean=0;       
        m_deconv_p[g].sigma=0;     
        m_deconv_p[g].weight=0;   
        m_deconv_p[g].yFactor=0;           
        }
      }
      
    m_quality=m_sGmm.quality;           //fit quality
    
//return nGauss; //all good
return m_nDeconv; //all good
}

//unit conversion.
//mean and sigma have units of scans and must be converted to Daltons.
//yFactor must adapt to the new units.
//proceeds to linear interpolation.
//Arguments:
// deconvIn_p -> pointer to structure with original units (scans).
// deconOut_p -> pointer to structure with units of Daltons.
// mzAxis_p   -> pointer to array with Dalton units.
//               Each new value is a new scan. Association between mz axes and magnitudes.
// mzAxisSize -> mz array size. It must coincide with the size of the magnitudes to be analyzed.
void GmmPeaks::gaussConversion(GAUSSIAN *deconvIn_p, GAUSSIAN *deconvOut_p, float *mzAxis_p, int mzAxisSize)
{
  double delta, offset, mean, sigma;
  int tmp;

  //mean adjustment
  mean=deconvIn_p->mean;
  tmp=(int)mean; //lower integer value.

  //interval in Daltons between two consecutive neighboring scans of mean.
  if(tmp+1 < mzAxisSize) //if it is within range.
    delta=mzAxis_p[tmp+1]-mzAxis_p[tmp];//posterior delta.
  else
    delta=mzAxis_p[tmp]-mzAxis_p[tmp-1];//previous delta.

  offset=delta*(mean-(double)tmp); //displacement with respect to the origin of the compound peak.
  deconvOut_p->mean=mzAxis_p[tmp]+offset; //converted value.

  //sigma adjustment, as in mean.
  sigma=deconvIn_p->sigma;
  tmp=(int)sigma;
  offset=delta*(sigma-(double)tmp);
  deconvOut_p->sigma=tmp*delta+offset;

  //Adjustment in Y. yFactor is increased by the quotient factor between the sigma in Daltons and in scans.
  deconvOut_p->yFactor = deconvIn_p->yFactor*deconvOut_p->sigma/deconvIn_p->sigma;
  deconvOut_p->yFactor/=(deconvOut_p->sigma*2.50662827); //normalizing factor (area=1.0).
  deconvOut_p->weight=deconvIn_p->weight; //The weight remains the same.
}


//Returns information about the sirem peaks that are candidates for generating Gaussians
unsigned long GmmPeaks::getEtpHits()
    {return m_etpHits;}

//Returns the number of Gaussians that integrate a compound magnitude peak
int GmmPeaks::getDeconvNumber()
    {return m_nDeconv;}

//Returns the Gaussian indexed by 'index' housed in a compound magnitude peak
GAUSSIAN GmmPeaks::getDeconv(int index)
    {return m_deconv_p[index];}

//Returns the direction of a Gaussian indexed by 'index' housed in a compound magnitude peak
GAUSSIAN *GmmPeaks::getDeconv_p(int index)
    {return &m_deconv_p[index];}

//Returns the quality of the Gaussian fit with the composite magnitude info
float GmmPeaks::getQuality()
    {return m_quality;}

//Returns a pointer to the magnitude information handled by the composite peak
float *GmmPeaks::getMagnitude()
    {return m_mag_p;}

//Returns an element of the magnitude information handled by the composite peak
float GmmPeaks::getMagnitude(int index)
    {return m_mag_p[index];}

//Returns the number of magnitude elements in the composite peak
int GmmPeaks::getMagnitudeNumber()
    {return m_nMag;}
