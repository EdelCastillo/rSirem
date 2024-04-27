/*************************************************************************
 *     Sirem and peaks
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

#include "siremPeaks.h"

//Constructor: initialization
SiremPeaks::SiremPeaks(SIREM_PEAKS_INFO *siremContex)
{
    m_siremContex_p=siremContex;
     m_siremCuts_p=0;

    //mz to determine the  number of tiles to use
    m_siremInfo.goodMz.set =siremContex->mzInfo.set;
    m_siremInfo.goodMz.size=siremContex->mzInfo.size;

    //coordinates of each pixel
    m_siremInfo.pxCoord.set=siremContex->pxCoord.set;
    m_siremInfo.pxCoord.size=siremContex->pxCoord.size;
    m_siremInfo.tileSide=4;//tesela de 4x4
    if(siremContex->algorithm==0)
        m_siremInfo.SIREM=true;
    else m_siremInfo.SIREM=false;
    
    m_cnvPeaks1_p=0;
    m_cnvPeaks2_p=0;
    m_magPeaks_p=0;
    m_siremIndex=0;
    if(siremContex->cutLevels.size<=0)
      {printf("Warning: cut level size must be greater than zero\n"); return;}
    if(siremContex->nScans<=0)
      {printf("Warning: no data detected\n"); return;}
    
    try
    {
        //memory to house sirem information about each cut level and each mz
        m_siremCuts_p=new float*[siremContex->cutLevels.size];
        for(int i=0; i<siremContex->cutLevels.size; i++)
            {m_siremCuts_p[i]=0; m_siremCuts_p[i]=new float[siremContex->nScans];}

        m_sirem_p=0;
        m_sirem_p=new Sirem(&m_siremInfo);

        //memory to hold average magnitude values for each mz
        m_meanMag_p=0;
        m_meanMag_p=new float[siremContex->nScans];
        m_unitedPeaks_p=0;
    }
    catch(const std::bad_alloc& e)
        {
        printf("Error reserving memory: %s\n",e.what());
        return;
        }
}

//destructor: free reserved memory
SiremPeaks::~SiremPeaks()
{
//printf("...ini SiremPeaks destructor\n");
    if( m_siremCuts_p)
        {
        for(int i=0; i<m_siremContex_p->cutLevels.size; i++)
            if(m_siremCuts_p[i]) delete[] m_siremCuts_p[i];
        delete m_siremCuts_p;
        }
    if(m_sirem_p)       delete m_sirem_p;
    if(m_meanMag_p)     delete m_meanMag_p;
    if(m_cnvPeaks1_p)   delete m_cnvPeaks1_p;
    if(m_cnvPeaks2_p)   delete m_cnvPeaks2_p;
    if(m_magPeaks_p)    delete m_magPeaks_p;
    if(m_unitedPeaks_p) delete [] m_unitedPeaks_p;
//printf("...end SiremPeaks destructor\n");
}

//Obtain sirem for each cut level of a mz, the results remain in the m_siremCuts_p matrix
//if the average value does not reach minMeanPxMag, sirem is aborted
//sirem is organized in a matrix of as many rows as cut levels have been applied and as many columns as mz are analyzed
//Arguments:
//  mzInfo_p -> vector with magnitude info for each pixel
//return:
//  The first free column in the matrix m_siremCuts_p
int SiremPeaks::newMz(GROUP_F *mzInfo_p)
{
    double meanValue=0;
    //average value of the magnitudes for this scan
    for(int i=0; i<m_siremContex_p->pxCoord.size; i++)
        meanValue+=mzInfo_p->set[i];
    meanValue/=m_siremContex_p->pxCoord.size;
    m_meanMag_p[m_siremIndex]=meanValue;

    if(meanValue<m_siremContex_p->minMeanPxMag) //If it does not exceed a minimum, the entropy is nullified
        {
        for(int cut=0; cut<m_siremContex_p->cutLevels.size; cut++)
            m_siremCuts_p[cut][m_siremIndex]=0;
        }

    else //the entropy is obtained for each cut
      {
        for(int cut=0; cut<m_siremContex_p->cutLevels.size; cut++)
            {
            m_siremCuts_p[cut][m_siremIndex]=m_sirem_p->getSirem(mzInfo_p, m_siremContex_p->cutLevels.set[cut], m_siremContex_p->noiseLevel, m_siremContex_p->minSectionDensity);
            //mzRelPxBins.set[k]=m_sirem_p->getBinaryImageDensity();
            }
      }
    m_siremIndex++; //updates to the first free index in the matrix.
    return m_siremIndex;
}

//Returns the SIREM value corresponding to the passed image.
//Arguments:
//magnitudes_p -> pointer to the magnitude of each pixel.
//magSize -> Number of pixels.
//cutLevel -> cut level (percentile).
//return
// 0 if the mean magnitude is less than minMeanPxMag.
// 0 if the image density is less than minSectionDensity.
// SIREM otherwise.
float SiremPeaks::getSirem(float *magnitudes_p, int magSize, float cutLevel)
{
    double meanValue=0;
    float sirem;
    GROUP_F magInfo;
    magInfo.set=magnitudes_p;
    magInfo.size=magSize;
    //average value of the magnitudes for this scan.
    for(int i=0; i<magSize; i++)
        meanValue+=magnitudes_p[i];
    meanValue/=m_siremContex_p->pxCoord.size;
    m_meanMag_p[m_siremIndex]=meanValue;

    if(meanValue<m_siremContex_p->minMeanPxMag) //If it does not exceed a minimum, the entropy is nullified.
        return 0;
    else //the entropy is obtained for each cut.
        sirem=m_sirem_p->getSirem(&magInfo, cutLevel, m_siremContex_p->noiseLevel, m_siremContex_p->minSectionDensity);
    return sirem;
}

//Get a list of sirem simple peaks associated with each magnitude simple peak
//Get the list of compound peaks: simple peaks joined by their valleys
//makes use of the magnitude information and the m_siremCuts_p structure
//Return -1 if unable to do so
int SiremPeaks::getPeaksList()
{
    if(m_siremIndex!=m_siremContex_p->nScans) return -1;
    Peaks *etpPeaks1_p=0, *etpPeaks2_p=0;
    int nEtpBlocks;
  
    try
        {
        //get the list of peaks from the magnitude
        //each block in the list contains three indices: mzLow, mzMax and mzHigh
        m_magPeaks_p=new Peaks(m_meanMag_p, m_siremContex_p->nScans, m_siremContex_p->magSenMin,
            m_siremContex_p->magSenMax, m_siremContex_p->noiseLevel);
        m_nMagBlocks=m_magPeaks_p->get(0, m_siremContex_p->nScans-1); //#magnitude peaks
        if(m_nMagBlocks<=0) return -1; //no peaks in data

        //get list of peaks from entropy
        etpPeaks1_p=new Peaks(m_siremCuts_p[0], m_siremContex_p->nScans, m_siremContex_p->etpSenMin, m_siremContex_p->etpSenMax, m_siremContex_p->etpSenMin);
        nEtpBlocks=etpPeaks1_p->get(0, m_siremContex_p->nScans-1); //#peaks from entropy

        //as there may be entropy information for several cuts, it is necessary to create a single array of peaks from all of them.
        //two lists are created: one maintains the history and the other the current list
        m_cnvPeaks1_p=new CONVOLVED_PEAKS[m_nMagBlocks]; //array de estructuras con info de picos simples
        m_cnvPeaks2_p=new CONVOLVED_PEAKS[m_nMagBlocks];
        for(int i=0; i<m_nMagBlocks; i++)
            {m_cnvPeaks1_p[i].etpPeaks_p=0;  m_cnvPeaks2_p[i].etpPeaks_p=0;}
           
        
        //for each magnitude peak, associate the contained entropy peaks. The result is m_cnvPeaks1_p
        splitEtpPeaks(m_cnvPeaks1_p, m_magPeaks_p->m_mzIndex_p, m_nMagBlocks, etpPeaks1_p->m_mzIndex_p, nEtpBlocks);
//      return(0);      
        
        //get the list of peaks from the entropy.
        //Won't be useful if there is only one cut
        etpPeaks2_p=new Peaks(m_siremCuts_p[0], m_siremContex_p->nScans, m_siremContex_p->etpSenMin, m_siremContex_p->etpSenMax,
            m_siremContex_p->etpSenMin);
        nEtpBlocks=etpPeaks2_p->get(0, m_siremContex_p->nScans-1); //#picos de entropía
        }
    catch(const std::bad_alloc& e)
        {
        printf("Error reserving memory: %s\n",e.what());
        return -1;
        }

    for(int cut=1; cut<m_siremContex_p->cutLevels.size; cut++)
        {
      //sirem peaks are arranged with magnitude peaks
        etpPeaks2_p->setMagnitude(m_siremCuts_p[cut], m_siremContex_p->nScans); //sets the magnitudes of sirem to extract its peaks
        nEtpBlocks=etpPeaks2_p->get(0, m_siremContex_p->nScans-1); //peaks
        splitEtpPeaks(m_cnvPeaks2_p, m_magPeaks_p->m_mzIndex_p, m_nMagBlocks, etpPeaks2_p->m_mzIndex_p, nEtpBlocks); //sorting on m_cnvPeaks2_p

        //Possible new peaks (m_cnvPeaks2_p) are added to the previous unified list (m_cnvPeaks1_p)
        peaksMix(m_cnvPeaks1_p, m_cnvPeaks2_p, m_nMagBlocks);

        //exchanging pointers to iterate
        Peaks *tmpPeaks_p=etpPeaks1_p;
        etpPeaks1_p=etpPeaks2_p;
        etpPeaks2_p=tmpPeaks_p;

        //exchanging pointers to iterate
        CONVOLVED_PEAKS *tmpCnv_p=m_cnvPeaks1_p;
        m_cnvPeaks1_p=m_cnvPeaks2_p;
        m_cnvPeaks2_p=tmpCnv_p;
        }

    //The simple peaks that are connected by their valleys are determined to form compound peaks.
    if(m_siremContex_p->cutLevels.size==1) //si solo hay un corte
        {
        unitedPeaks(m_cnvPeaks1_p, m_nMagBlocks);
        m_cnvPeaks_p=m_cnvPeaks1_p;
        }
    else
        {
        unitedPeaks(m_cnvPeaks2_p, m_nMagBlocks);
        m_cnvPeaks_p=m_cnvPeaks2_p;
        }
    if(etpPeaks1_p) delete etpPeaks1_p;
    if(etpPeaks2_p) delete etpPeaks2_p;
    return 0;
}

//The lists of entropy peaks from the Peaks class are analyzed and organized
//An array of CONVOLVED_PEAKS structures is initialized for all simple peaks of magnitude
//for each magnitude peak, the entropy peaks contained are noted
//Arguments:
// cnvPeaks_p -> target entropy peak sorting under a magnitude peak
// magPeaks -> pointer to magnitude peaks
// nMagPeaks -> #magnitude peaks
// etpPeak_p -> pointer to entropy peaks
// nEtpPeaks -> #entropy peaks
//Returns the total number of final entropy peaks
int SiremPeaks::splitEtpPeaks(CONVOLVED_PEAKS *cnvPeaks_p, Peaks::ION_INDEX *magPeaks_p, int nMagPeaks, Peaks::ION_INDEX *etpPeaks_p, int nEtpPeaks)
{
    int lastEtpPeak, ePeak, nEtpPeaksConv, totalPeaks=0;

    //we look for the first entropy index located within the magnitude peak
    lastEtpPeak=-1;
    for(int i=0; i<nEtpPeaks; i++) //tour of entropy peaks
      if(etpPeaks_p[i].max >= magPeaks_p[0].low && etpPeaks_p[i].max <= magPeaks_p[0].high)
        {lastEtpPeak=i; break; }
    if(lastEtpPeak==-1) //no entropy peak is within the magnitude peak
        {cnvPeaks_p[0].nEtpPeaks=0;}
    for(int mPeak=0; mPeak<nMagPeaks; mPeak++) //for all magnitude peaks
        {
        for(int i=lastEtpPeak; i<nEtpPeaks; i++) //tour of entropy peaks to find the first one into mPeak
            if(etpPeaks_p[i].max >= magPeaks_p[mPeak].low && etpPeaks_p[i].max <= magPeaks_p[mPeak].high)
            {lastEtpPeak=i; break;}

        cnvPeaks_p[mPeak].magPeaks.low =magPeaks_p[mPeak].low; //magnitude peak
        cnvPeaks_p[mPeak].magPeaks.max =magPeaks_p[mPeak].max;
        cnvPeaks_p[mPeak].magPeaks.high=magPeaks_p[mPeak].high;

//        if(!hit) //if there are no entropy peaks either within or above this mag peak
//        {cnvPeaks_p[mPeak].nEtpPeaks=0; continue;}
        
        //number of entropy peaks in each magnitude peak
        cnvPeaks_p[mPeak].nEtpPeaks=0; //default.
        nEtpPeaksConv=0;  //#entropy peaks
        totalPeaks++;
        
        //in each magnitude peak there can be several entropy peaks, or none.
        //we look for those magnitude peaks that include several entropy peaks and count them
        //si el pico de sirem coincide con el valle derecho   del pico de mag, no se incluye 
        //si el pico de sirem coincide con el valle izquierdo del pico de mag, sí se incluye 
        //se debe hacer un análisis más detallado
        for(ePeak=lastEtpPeak; ePeak<nEtpPeaks; ePeak++) //tour of entropy peaks
            {
            //if the maximum value of the entropy peak is contained in the magnitude peak
            if(etpPeaks_p[ePeak].max >= magPeaks_p[mPeak].low && etpPeaks_p[ePeak].max <= magPeaks_p[mPeak].high)
                nEtpPeaksConv++;
            else if(etpPeaks_p[ePeak].max > magPeaks_p[mPeak].high) //outside the peak magnitude
                break;
            }
        
        //Entropy peaks are associated with each magnitude peak to which they correspond.
        //if two maxima are very close, the second is discarded. The tube effect is taken into account (small continuous increments)
        try
            {
            if(nEtpPeaksConv>0) //if there are peaks for deconvolution
                {
                //memory to accommodate all peaks plus a possible inserted
                cnvPeaks_p[mPeak].etpPeaks_p=new Peaks::ION_INDEX [nEtpPeaksConv+1];

                for(int i=0; i<nEtpPeaksConv; i++)
                    {
                    cnvPeaks_p[mPeak].etpPeaks_p[i].low =etpPeaks_p[i+lastEtpPeak].low;
                    cnvPeaks_p[mPeak].etpPeaks_p[i].max =etpPeaks_p[i+lastEtpPeak].max;
                    cnvPeaks_p[mPeak].etpPeaks_p[i].high=etpPeaks_p[i+lastEtpPeak].high;
                    cnvPeaks_p[mPeak].nEtpPeaks++;
                    }
                }

            else //without entropy peaks: magnitude information is assigned
                {
                cnvPeaks_p[mPeak].etpPeaks_p=new Peaks::ION_INDEX [1]; //memoria para alojar un pico
                cnvPeaks_p[mPeak].etpPeaks_p[0].low =magPeaks_p[mPeak].low; //pico de magnitud
                cnvPeaks_p[mPeak].etpPeaks_p[0].max =magPeaks_p[mPeak].max;
                cnvPeaks_p[mPeak].etpPeaks_p[0].high=magPeaks_p[mPeak].high;
                cnvPeaks_p[mPeak].nEtpPeaks=1;
                }
            lastEtpPeak+=nEtpPeaksConv; //base for each magnitude peak
            totalPeaks +=nEtpPeaksConv; //total peaks
            }
        catch(const std::bad_alloc& e)
            {
            printf("Error reserving memory: %s\n",e.what());
            return -1;
            }
        }

    //ensures that a magnitude peak can be shaped by at least one entropy peak.
    //If the sigma of the entropy peaks do not cover the magnitude peak, a new peak is inserted
    analizeMagPeak(cnvPeaks_p, magPeaks_p, nMagPeaks);
    return totalPeaks;
}

//the possibility of a magnitude peak being left without an entropy peak near its maximum is analyzed.
//If so, the range of neighboring peaks is expanded to cover it
//near: that the width of the entropy peak does not contain the maximum magnitude.
//Arguments:
// cnvPeaks_p -> info on magnitude and entropy peaks
// magPeaks_p -> pointer to magnitude peaks
// nMagPeaks -> #magnitude peaks
void SiremPeaks::analizeMagPeak(CONVOLVED_PEAKS *cnvPeaks_p, Peaks::ION_INDEX *magPeaks_p, int nMagPeaks)
{
    //No spike is inserted; the range of neighboring peaks is expanded to cover it
    bool hit;
    for(int mPeak=0; mPeak<nMagPeaks; mPeak++) //for all magnitude peaks
        {
        hit=false;
        for(int p=0; p<cnvPeaks_p[mPeak].nEtpPeaks; p++) //for all entropy peaks
            {
          //if the maximum in magnitude is contained in the entropy peak
            if(cnvPeaks_p[mPeak].etpPeaks_p[p].low<magPeaks_p[mPeak].max && cnvPeaks_p[mPeak].etpPeaks_p[p].high>magPeaks_p[mPeak].max)
                {hit=true; break;}
            }
        if(!hit)
            {
            if(cnvPeaks_p[mPeak].nEtpPeaks==1) //just a peak of entropy: ranks are equalized
                {
                cnvPeaks_p[mPeak].etpPeaks_p[0].low =magPeaks_p[mPeak].low;
                cnvPeaks_p[mPeak].etpPeaks_p[0].high=magPeaks_p[mPeak].high;
                }
            else
                {
                //neighboring peaks are identified
                int pLeft=-1, pRight=-1;
                for(int p=0; p<cnvPeaks_p[mPeak].nEtpPeaks; p++) //for all entropy peaks into mPeak
                    if(cnvPeaks_p[mPeak].etpPeaks_p[p].high<=magPeaks_p[mPeak].max) pLeft=p;
                    else if(cnvPeaks_p[mPeak].etpPeaks_p[p].low>=magPeaks_p[mPeak].max) {pRight=p; break;}
                if(pLeft!=-1 && pRight==-1) //There are no entropy peaks to the right of the magnitude peak
                    cnvPeaks_p[mPeak].etpPeaks_p[pLeft].high =magPeaks_p[mPeak].high; //expands to the right
                if(pLeft==-1 && pRight!=-1) //There are no entropy peaks to the left of the magnitude peak
                    cnvPeaks_p[mPeak].etpPeaks_p[pRight].low =magPeaks_p[mPeak].low; //expands to the left
                if(pLeft!=-1 && pRight!=-1) //There are entropy peaks to the left and right of the magnitude peak.
                    {
                    cnvPeaks_p[mPeak].etpPeaks_p[pLeft].high =magPeaks_p[mPeak].max+1; //expands to the right
                    cnvPeaks_p[mPeak].etpPeaks_p[pRight].low =magPeaks_p[mPeak].max-1; //expands to the left
                    }
                  /*             
                  if(m_siremContex_p->vervose)
                    {
                    if(pLeft!=-1 && pRight==-1)
                        printf("...Ajustado [%d %d %d]\n", cnvPeaks_p[mPeak].etpPeaks_p[pLeft].low, cnvPeaks_p[mPeak].etpPeaks_p[pLeft].max, cnvPeaks_p[mPeak].etpPeaks_p[pLeft].high);
                    if(pLeft==-1 && pRight!=-1)
                        printf("...Ajustad: [%d %d %d]\n", cnvPeaks_p[mPeak].etpPeaks_p[pRight].low, cnvPeaks_p[mPeak].etpPeaks_p[pRight].max, cnvPeaks_p[mPeak].etpPeaks_p[pRight].high);
                    if(pLeft!=-1 && pRight!=-1)
                        printf("...Ajustados: [%d %d %d][%d %d %d]\n", cnvPeaks_p[mPeak].etpPeaks_p[pLeft].low, cnvPeaks_p[mPeak].etpPeaks_p[pLeft].max, cnvPeaks_p[mPeak].etpPeaks_p[pLeft].high,
                            cnvPeaks_p[mPeak].etpPeaks_p[pRight].low, cnvPeaks_p[mPeak].etpPeaks_p[pRight].max, cnvPeaks_p[mPeak].etpPeaks_p[pRight].high);
                    }
                  */
                
                }
            }
        }
}

//Integrates two peak lists into a single list
//If the two lists are unequal in size, the largest of them is considered.
//If they are equal, the old list (cnvPeaks2_p) is considered but the possibility that some peak of cnvPeaks1_p is isolated is taken into account: it is included.
//If two peaks coincide in their max value (within the minScansPeaksGap range) the first one is included.
//Arguments:
// cnvPeaks1_p -> list of source and destination peaks
// cnvPeaks2_p -> list of source peaks
//Return true if OK
bool SiremPeaks::peaksMix(CONVOLVED_PEAKS *cnvPeaks1_p, CONVOLVED_PEAKS *cnvPeaks2_p, int nMagPeaks)
{
    int size, maxSize=0;
    if(cnvPeaks1_p==0 || cnvPeaks2_p==0 || nMagPeaks >=0) return false;

    for(int mPeak=0; mPeak<nMagPeaks; mPeak++) //greater number of peaks for memory reservation.
        {
        size=cnvPeaks1_p[mPeak].nEtpPeaks+cnvPeaks2_p[mPeak].nEtpPeaks;
        if(size>maxSize) maxSize=size;
        }

    int index, count;
    for(int mPeak=0; mPeak<nMagPeaks; mPeak++) //for all simple peaks.
        {
        count=0;
        printf("   A:");
        for(index=0; index<cnvPeaks2_p[mPeak].nEtpPeaks; index++)
            printf("[%d %d %d] ", cnvPeaks2_p[mPeak].etpPeaks_p[index].low,
                   cnvPeaks2_p[mPeak].etpPeaks_p[index].max,
                   cnvPeaks2_p[mPeak].etpPeaks_p[index].high);
        printf("    B:");
        for(int i=0; i<cnvPeaks1_p[mPeak].nEtpPeaks; i++, index++)
            printf("[%d %d %d] ", cnvPeaks1_p[mPeak].etpPeaks_p[i].low,
                   cnvPeaks1_p[mPeak].etpPeaks_p[i].max,
                   cnvPeaks1_p[mPeak].etpPeaks_p[i].high);
        printf("\n");

        Peaks::ION_INDEX *peaks_p=0;
        try
            {
            peaks_p=new Peaks::ION_INDEX[maxSize+1]; //temporary memory
            }
        catch(const std::bad_alloc& e)
            {
            printf("Error reserving memory: %s\n",e.what());
            return -1;
            }

        //Se elige a  la mayor lista
        if(cnvPeaks1_p[mPeak].nEtpPeaks > cnvPeaks2_p[mPeak].nEtpPeaks)         //peaks before > peaks now
            for(index=0; index<cnvPeaks1_p[mPeak].nEtpPeaks; index++)           //the old peaks are considered
                {
                peaks_p[index].low =cnvPeaks1_p[mPeak].etpPeaks_p[index].low;
                peaks_p[index].max =cnvPeaks1_p[mPeak].etpPeaks_p[index].max;
                peaks_p[index].high=cnvPeaks1_p[mPeak].etpPeaks_p[index].high;
                }
        else if(cnvPeaks1_p[mPeak].nEtpPeaks < cnvPeaks2_p[mPeak].nEtpPeaks)    //peaks before < peaks now
            for(index=0; index<cnvPeaks2_p[mPeak].nEtpPeaks; index++)           //new peaks are considered
                {
                peaks_p[index].low =cnvPeaks2_p[mPeak].etpPeaks_p[index].low;
                peaks_p[index].max =cnvPeaks2_p[mPeak].etpPeaks_p[index].max;
                peaks_p[index].high=cnvPeaks2_p[mPeak].etpPeaks_p[index].high;
                }
        else //si son iguales en tamaño
            {
            for(index=0; index<cnvPeaks2_p[mPeak].nEtpPeaks; index++)           //the old peaks are considered
                {
                peaks_p[index].low =cnvPeaks2_p[mPeak].etpPeaks_p[index].low;
                peaks_p[index].max =cnvPeaks2_p[mPeak].etpPeaks_p[index].max;
                peaks_p[index].high=cnvPeaks2_p[mPeak].etpPeaks_p[index].high;
                }

            //Also included in the new list are those old peaks located beyond the minimum distance from any other.
            for(int i=0; i<cnvPeaks1_p[mPeak].nEtpPeaks; i++)
                {
                if(insertPeak(peaks_p, &cnvPeaks1_p[mPeak].etpPeaks_p[i], index)!=-1)
                    index++;
                }
            }
        //se libera memoria
        if(cnvPeaks1_p[mPeak].etpPeaks_p) delete [] cnvPeaks1_p[mPeak].etpPeaks_p; //frees temporary memory.
        if(peaks_p) delete peaks_p;

        //cnvPeaks1_p becomes the new peak list
        cnvPeaks1_p[mPeak].etpPeaks_p=peaks_p;
        cnvPeaks1_p[mPeak].nEtpPeaks=index;

printf("Out->");
for(int i=0; i<cnvPeaks1_p[mPeak].nEtpPeaks; i++)
    printf("[%d %d %d] ", cnvPeaks1_p[mPeak].etpPeaks_p[i].low,
            cnvPeaks1_p[mPeak].etpPeaks_p[i].max,
            cnvPeaks1_p[mPeak].etpPeaks_p[i].high);
printf("\n");

        }

   return true;
}

//insert a peak into an array of peaks ordered from smallest to largest if it respects the minimum distance from its neighbors.
//there must be an extra position at the end of the array!!!
//Arguments:
// peak_p -> peak to insert
// etpPeaks_p -> array of peaks for insertion
// etpSize -> array size
//Return the index to the insertion point
int SiremPeaks::insertPeak(Peaks::ION_INDEX *etpPeak_p, Peaks::ION_INDEX *peak_p, int etpSize)
{
    int insertIndex=-1;
    for(int p=0; p<etpSize; p++) //search for insertion point
        {
        if(p==0 && peak_p->max < etpPeak_p[p].max) //below the first
            {
            if(abs(etpPeak_p[p].max-peak_p->max) <= m_siremContex_p->minScansPeaksGap) {insertIndex=-1; break;} //not inserted
            else insertIndex=p;
            }
        else if(p+1>=etpSize && peak_p->max > etpPeak_p[p].max)
            {
            if(abs(peak_p->max-etpPeak_p[p].max) <= m_siremContex_p->minScansPeaksGap) {insertIndex=-1; break;} //not inserted
            else insertIndex=p+1;
            }
        else if(peak_p->max > etpPeak_p[p].max && peak_p->max < etpPeak_p[p+1].max)
            {
            if(abs(peak_p->max-etpPeak_p[p].max) > m_siremContex_p->minScansPeaksGap &&
             abs(peak_p->max-etpPeak_p[p+1].max) > m_siremContex_p->minScansPeaksGap) //if the position exists
                insertIndex=p+1;
            else {insertIndex=-1; break;}
            }

        if(insertIndex!=-1)
            {
            for(int i=etpSize-1; i>=insertIndex; i--) //if found, all peaks are moved above it one position.
                {
                etpPeak_p[i+1].low =etpPeak_p[i].low;
                etpPeak_p[i+1].max =etpPeak_p[i].max;
                etpPeak_p[i+1].high=etpPeak_p[i].high;
                }
            //se inserta el pico dado
            etpPeak_p[insertIndex].low =peak_p->low;
            etpPeak_p[insertIndex].max =peak_p->max;
            etpPeak_p[insertIndex].high=peak_p->high;
            break;
            }
        }
    return insertIndex;
}

//The array of CONVOLVED_PEAKS structures, of simple peaks, is analyzed in search of possible compound peaks.
//A composite peak is a grouping of consecutive simple peaks joined by a valley of magnitude greater than m_noiseLevel.
//The info is in the array m_unitedPeaks_p[] of size m_nUnitedPeaks. The size is returned.
//Arguments:
// cnvPeaks_p -> array of structures with simple peaks info
// nMagPeaks -> array size
//Returns the number of compound peaks
int SiremPeaks::unitedPeaks(CONVOLVED_PEAKS *cnvPeaks_p, int nMagPeaks)
{
    if(m_unitedPeaks_p) {delete [] m_unitedPeaks_p; m_unitedPeaks_p=0;}
    try
        {
        m_unitedPeaks_p=new UNITED_PEAKS[nMagPeaks];
        }
    catch(const std::bad_alloc& e)
        {
        printf("Error reserving memory: %s\n",e.what());
        return -1;
        }
    m_nUnitedPeaks=0;
    bool hit;
    int mPeak=0;
    while(true)
        {
        //if the index is less than the maximum and is linked to the next peak
        if(mPeak<(nMagPeaks-1) && (cnvPeaks_p[mPeak].magPeaks.high==cnvPeaks_p[mPeak+1].magPeaks.low) &&
            m_meanMag_p[cnvPeaks_p[mPeak].magPeaks.high]>m_siremContex_p->noiseLevel)
            {
            //se toma nota
            m_unitedPeaks_p[m_nUnitedPeaks].peakLow =mPeak;
            m_unitedPeaks_p[m_nUnitedPeaks].peakHigh=mPeak+1;
            mPeak++;
            //perhaps there are more peaks linked in a chain
            hit=false;
            while(true)
                {
                //if the index is less than the maximum and is linked to the next peak.
                if(mPeak<(nMagPeaks-1) && (cnvPeaks_p[mPeak].magPeaks.high==cnvPeaks_p[mPeak+1].magPeaks.low) &&
                    m_meanMag_p[cnvPeaks_p[mPeak].magPeaks.high]>m_siremContex_p->noiseLevel)
                    {
                    //the end of the joint peak is updated.
                    m_unitedPeaks_p[m_nUnitedPeaks].peakHigh=mPeak+1;
                    mPeak++;
                    }
                //If you have reached the last peak and it is linked to the previous one: end.
                else if(mPeak>=(nMagPeaks-1))
                    {
                    mPeak++; m_nUnitedPeaks++;
                    hit=true;
                    break;
                    }
                //if you have not reached the end and are not linked to the next one.
                else {mPeak++; m_nUnitedPeaks++; break;} //two separate peaks.
                }
            if(hit) break;
            }

        else  //single peak
            {
            m_unitedPeaks_p[m_nUnitedPeaks].peakLow =mPeak;
            m_unitedPeaks_p[m_nUnitedPeaks].peakHigh=mPeak;
            m_nUnitedPeaks++;
            if(mPeak >= nMagPeaks-1) break;
            else mPeak++;
            }

        }
  return m_nUnitedPeaks;
}

//conversión de valores en un sistema no lineal
//Argumentos:
//  dataIn_p  -> viene en unidades en el rango [0:size-1]
//  dataOut_p -> son los datos convertidos
//  newUnit_p -> es un array con las unidades de destino
//  size      -> es el tamaño de todos los arrays
void SiremPeaks::unitConversion(float *dataIn_p, float *dataOut_p, float *newUnit_p, int size)
{
  float data, delta, offset;
  int tmp;

  //ajuste de mean
  for(int i=0; i<size; i++)
  {
    data=dataIn_p[i];
    tmp=(int)data;
    //delta se adecúa a la diferencial de los datos finales
    if(tmp+1 < size) //si dentro de rango
      delta=newUnit_p[tmp+1]-newUnit_p[tmp];//delta posterior
    else
      delta=newUnit_p[tmp]-newUnit_p[tmp-1];//delta anterior

    offset=delta*(data-(float)tmp); //desplazamiento respecto al origen del pico compuesto
    dataOut_p[i]=newUnit_p[tmp]+offset; //valor convertido

  }
}

//Returns a pointer to the CONVOLVED_PEAKS structure indexed by mPeak.
CONVOLVED_PEAKS* SiremPeaks::getCnvPeak(int mPeak)
    { return m_cnvPeaks_p+mPeak;}

//Returns the information of a simple magnitude peak.
Peaks::ION_INDEX SiremPeaks::getSingleMagPeak(int mPeak)
    {return m_cnvPeaks_p[mPeak].magPeaks;}

//Returns a pointer to the entropy peaks corresponding to a single magnitude peak.
Peaks::ION_INDEX *SiremPeaks::getEtpPeaks(int mPeak)
    {return m_cnvPeaks_p[mPeak].etpPeaks_p;}

//Returns a peak entropy, indexed by 'index', corresponding to a single peak magnitude.
Peaks::ION_INDEX SiremPeaks::getEtpSinglePeak(int mPeak, int index)
    {return m_cnvPeaks_p[mPeak].etpPeaks_p[index];}

//Returns the size of the array of entropy peaks corresponding to a single magnitude peak.
int SiremPeaks::getEtpPeaksNumber(int mPeak)
    {return m_cnvPeaks_p[mPeak].nEtpPeaks;}

//Returns the number of simple peaks.
int SiremPeaks::getSinglePeaksNumber()
    {return m_nMagBlocks;}

//Returns the simple low and high peaks that make up a composite peak.
UNITED_PEAKS SiremPeaks::getCompoundPeak(int peak)
    {return m_unitedPeaks_p[peak];}

//Returns the number of compound peaks.
int SiremPeaks::getCompoundPeakNumber()
    {return m_nUnitedPeaks;}

//Returns SIREM for a given cut and scan.
float SiremPeaks::getSiremSingleCut(int cut, int index)
    {return m_siremCuts_p[cut][index];}

//Returns a pointer to all SIREM information for a given cut.
float* SiremPeaks::getSiremCut(int cut)
    {return m_siremCuts_p[cut];}

//Returns the number of rows in the cut matrix (cut levels).
int SiremPeaks::getSiremCutRows()
    {return m_siremContex_p->cutLevels.size;}

//Returns the number of columns in the slice matrix (number of scans-mz).
int SiremPeaks::getSiremCutCols()
    {return m_siremContex_p->nScans;}

float *SiremPeaks::getMeanMagnitude()
    {return m_meanMag_p;}

//set average magnitudes (concentrations)
//Arguments:
//  meanMag -> average magnitudes
//  size    -> data size;
//return:
//  -1 if size not match siremContex_p->nScans; 0 if OK
int SiremPeaks::setMeanMagnitudes(float *meanMag, int size)
    {
    if(size!=m_siremContex_p->nScans) return -1;
    for(int i=0; i<size; i++)
        m_meanMag_p[i]=meanMag[i];
    m_siremIndex=size; //sample size
    return 0;
    }
    
float SiremPeaks::getMeanMagnitude(int scan)
    {return m_meanMag_p[scan];}

