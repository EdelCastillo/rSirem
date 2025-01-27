#ifndef GMM_PEAKS_H
#define GMM_PEAKS_H

#include "stdio.h"
#include <stdexcept>
#include "peaks.h"
#include "GMM.h"
//#include "gmmPeaks.h"
#include "siremPeaks.h"


//Decomposes a single composite peak, possibly formed by several simple peaks, into Gaussians.
class GmmPeaks
{
public:
    //Constructor:
    //Argument:
    //minMeanPeakMagnitude -> minimum averaged magnitude for a peak to be deconvolved
    GmmPeaks(float minMeanPeakMagnitude);

    //Destructor: frees reserved memory.
    ~GmmPeaks();

    //magnitudes associated with all the peaks to be treated
    //Arguments:
    //mag_p -> structure with magnitude information associated with each scan
    void setMagnitudes(GROUP_F *mag_p);

    //peaks to consider
    //Arguments:
    //cnvPeaks_p -> pointer to array of the set of simple peaks that constitute the composite peak
    //nUnitedPeaks -> #simple peaks
    void setPeaks(CONVOLVED_PEAKS *cnvPeaks_p, int nUnitedPeaks);

    //Deconvolution
    //From a compound magnitude peak and its entropy peaks,
    //obtains the Gaussians that compose it: EM Algorithm
    //establishes the array m_deconv_p[] with m_nDeconv elements with info of each Gaussian ordered from lowest to highest mass
    //If a peak could not be deconvolved, only one structure appears, its field nGauss=0 and its field quality=-1.
    //The mean and sigma information of the Gaussians is in units of scans
    //return:
    //Number of gaussians
    int gmmDeconvolution();

    //conversión de unidades
    //mean y sigma tienen unidades de scans y deben convertirse a Daltons
    //yFactor debe adaptarse a las nuevas unidades
    //se procede a una interpolación lineal
    //Argumentos:
    //  deconvIn_p -> puntero a estructura con unidades originales (scans)
    //  deconOut_p -> puntero a estructura con unidades de Daltons
    //  mzAxis_p   -> puntero a array con unidades de Dalton.
    //                Cada nuevo valor es un nuevo scan. Asociación entre ejes de mz y magnitudes
    //  mzAxisSize -> tamaño del array de mz. Debe coincidir con el tamaño de las magnitudes a analizar
    void gaussConversion(GAUSSIAN *deconvIn_p, GAUSSIAN *deconvOut_p, float *mzAxis_p, int mzAxisSize);

    //Returns information about the sirem peaks that are candidates for generating Gaussians
    unsigned long getEtpHits();

    //Returns the number of Gaussians that integrate a compound magnitude peak
    int getDeconvNumber();

    //Returns the Gaussian indexed by 'index' housed in a compound magnitude peak
    GAUSSIAN getDeconv(int index);

    //Retorna la dirección de una gausiana de un pico de magnitud compuesto
    GAUSSIAN *getDeconv_p(int index);

    //Returns the quality of the Gaussian fit with the composite magnitude info
    float getQuality();

    //Returns a pointer to the magnitude information handled by the composite peak
    float *getMagnitude();

    //Returns an element of the magnitude information handled by the composite peak
    float getMagnitude(int index);

    //Returns the number of magnitude elements in the composite peak
    int getMagnitudeNumber();

private:
    //prepare the necessary information before calling the GMM algorithm: initialize the m_sGmm structure
    //Requires the reference to the composite peak to be deconvolved and the associated magnitudes
    //Return -1 if the average magnitude within the peak does not reach a minimum (m_minMeanPeakMagnitude)
    int iniGMM();

    float           *m_mag_p,       //magnitude array
                    m_minMeanPeakMagnitude; //peaks of lower magnitude are discarded
    int             m_nMag;         //magnitude array size
    GMM_STRUCT      m_sGmm;         //structure for deconvolution of a peak
    CONVOLVED_PEAKS *m_cnvPeaks_p;  //all peaks: magnitude and entropy
    int             m_nUnitedPeaks; //number of entries in m_unitedPeaks_p
    GAUSSIAN        *m_deconv_p,    //pointer to gausians
                    m_meanDeconv[32]; //array de gausianas con valores promediados
    int             m_nDeconv;      //deconvolution number of gaussians
    unsigned long   m_etpHits;      //sirem peaks array
    float           m_quality;      //quality of the Gaussian fit
};

#endif
