#ifndef SIREM_PEAKS_H
#define SIREM_PEAKS_H

#include "stdio.h"
#include "relEntropy.h"
#include <stdexcept>
#include "peaks.h"
#include "sirem.h"
//#include "../GMM/GMM.hpp"

    enum CUT_LEVEL_TYPE             //ways to determine the cut level.
        {
        PERCENTIL,                  //according to a percentile.
        MAX_VALUE_ABSOLUTE,         //according to a percentage of the absolute maximum value.
        MAX_VALUE_PERCENTIL,        //according to a percentage of the maximum level obtained by a percentile.
        MEAN_VALUE                  //average value
        };

    typedef struct
        {
        CUT_LEVEL_TYPE  type;       //type of procedure to determine the cut level.
        bool        up;             //pixels that exceed or equal the cutoff level are extracted.
        float       percentil_1,    //percentile for the PERCENTILE type.
                    percentil_2,    //percentile for type MAX_VALUE_PERCENTIL.
                    maxAbsolute,    //percentage of the absolute maximum value.
                    meanFactor;     //factor applicable to the average value.
        int         nCuts,
                    cutIndex;       //indicates the cut to be made.
        float       cutLevels[20];
        }CUT_LEVEL;                 //info about the cut level.

    typedef struct                      //información general
        {
        GROUP_F     mzInfo,             //magnitude info for a representative scan. Used to determine the useful tiles in the image.
                    cutLevels;          //info about the cut to 3D info (pixel, masses, magnitude).
        GROUP_XY    pxCoord;            //coordinates of each pixel.
        int         nScans,             //#of mz to analyze. For memory reservation.
                    minScansPeaksGap;   //minimum separation in scans to consider separate peaks.
        char        tileSide;           //tile size for entropy.
        float       noiseLevel,         //noise level in magnitude info.
                    minMeanPxMag,       //minimum average value of the magnitude information of a peak for consideration.
                    minSectionDensity,  //Minimum density of a sectioned image for consideration.
                    magSenMin,          //minimum sensitivity in magnitude.
                    magSenMax,          //maximum sensitivity in magnitude.
                    etpSenMin,          //minimum sensitivity in entropy.
                    etpSenMax;          //maximum sensitivity in entropy.
        char        algorithm,          //SIREM/entropy
                    referenceType;      //reference image
        }SIREM_PEAKS_INFO;


    typedef struct                  //info on a simple peak.
        {
        Peaks::ION_INDEX
                    magPeaks,       //simple magnitude peak.
                    *etpPeaks_p;    //array of entropy and/or correlation peaks.
        int         nEtpPeaks;      //number of entries in etpPeaks_p.
        }CONVOLVED_PEAKS;

    typedef struct                  //info about a compound peak.
        {                           //If it is a peak composed of a single simple peak, both entries coincide.
        int         peakLow,        //index to the lower peak of the compound peak.
                    peakHigh;       //index to the upper peak of the composite peak.
        }UNITED_PEAKS;


class SiremPeaks
{
    public:
        //Constructor: initialization
        SiremPeaks(SIREM_PEAKS_INFO *siremContex_p);

        //destructor: frees reserved memory.
        ~SiremPeaks();

        //Obtains sirem for each cut level of a mz, the results remain in the m_siremCuts_p matrix.
        //if the average value does not reach minMeanPxMag, sirem is aborted.
        //sirem is organized in a matrix of as many rows as cut levels have been applied and as many columns as mz are analyzed.
        //Arguments:
        //mzInfo_p -> vector with magnitude info for each pixel
        //return:
        //The first free column in the matrix m_siremCuts_p
        int newMz(GROUP_F *mzInfo_p);

        //Returns the SIREM value corresponding to the passed image.
        //Arguments:
        //magnitudes_p -> pointer to the magnitude of each pixel.
        //magSize -> Number of pixels.
        //cutLevel -> cut level (percentile).
        //return
        // 0 if the mean magnitude is less than minMeanPxMag.
        // 0 if the image density is less than minSectionDensity.
        // SIREM otherwise.
        float getSirem(float *magnitudes_p, int magSize, float cutLevel);

        //Get a list of simple sirem peaks associated with each magnitude simple peak.
        //Gets the list of compound peaks: simple peaks joined by their valleys.
        //makes use of the magnitude information and the m_siremCuts_p structure.
        //Return -1 if unable to do so
        int getPeaksList();

        //conversión de valores en un sistema no lineal
        //Argumentos:
        //  dataIn_p  -> viene en unidades en el rango [0:size-1]
        //  dataOut_p -> son los datos convertidos
        //  newUnit_p -> es un array con las unidades de destino
        //  size      -> es el tamaño de todos los arrays
        void unitConversion(float *dataIn_p, float *dataOut_p, float *newUnit_p, int size);

        //Returns a pointer to the CONVOLVED_PEAKS structure indexed by mPeak.
        CONVOLVED_PEAKS* getCnvPeak(int mPeak);

        //Returns the information of a simple magnitude peak.
        Peaks::ION_INDEX getSingleMagPeak(int mPeak);

        //Returns a pointer to the entropy peaks corresponding to a single magnitude peak.
        Peaks::ION_INDEX *getEtpPeaks(int mPeak);

        //Returns a peak entropy, indexed by 'index', corresponding to a single peak magnitude.
        Peaks::ION_INDEX getEtpSinglePeak(int mPeak, int index);

        //Returns the size of the array of entropy peaks corresponding to a single magnitude peak.
        int getEtpPeaksNumber(int mPeak);

        //Returns the number of simple peaks.
        int getSinglePeaksNumber();

        //Returns the simple low and high peaks that make up a composite peak.
        UNITED_PEAKS getCompoundPeak(int peak);

        //Returns the number of compound peaks.
        int getCompoundPeakNumber();

        //Returns SIREM for a given cut and scan.
        float getSiremSingleCut(int cut=0, int index=0);

        //Returns a pointer to all SIREM information for a given cut.
        float* getSiremCut(int cut);

        //Returns the number of rows in the cut matrix (cut levels).
        int getSiremCutRows();

        //Returns the number of columns in the slice matrix (number of scans-mz).
        int getSiremCutCols();

        //Returns a pointer to the averaged values of all scans over all pixels.
        float *getMeanMagnitude();

        //Returns a value averaged from one scan over all pixels.
        float getMeanMagnitude(int scan);
        
        //set average magnitudes (concentrations)
        //Arguments:
        //  meanMag -> average magnitudes
        //  size    -> data size;
        //return:
        //  -1 if size not match siremContex_p->nScans; 0 if OK
        int setMeanMagnitudes(float *meanMag, int size);
          
    private:
        //The lists of entropy peaks from the Peaks class are analyzed and organized.
        //An array of CONVOLVED_PEAKS structures is initialized for all single magnitude peaks.
        //for each magnitude peak, the contained entropy peaks are noted.
        //Arguments:
        // cnvPeaks_p -> target entropy peak sorting under a magnitude peak
        // magPeaks -> pointer to magnitude peaks
        // nMagPeaks -> #magnitude peaks
        // etpPeak_p -> pointer to entropy peaks
        // nEtpPeaks -> #entropy peaks
        //Returns the total number of final entropy peaks
        int splitEtpPeaks(CONVOLVED_PEAKS *cnvPeaks_p, Peaks::ION_INDEX *magPeaks_p, int nMagPeaks, Peaks::ION_INDEX *etpPeaks_p, int
nEtpPeaks);

        //the possibility of a magnitude peak being left without an entropy peak near its maximum is analyzed.
        //If so, the range of neighboring peaks is expanded to cover it.
        //close: the width of the entropy peak does not contain the maximum magnitude
        //Arguments:
        // cnvPeaks_p -> info on magnitude and entropy peaks
        // magPeaks_p -> pointer to magnitude peaks
        // nMagPeaks -> #magnitude peaks
        void analizeMagPeak(CONVOLVED_PEAKS *cnvPeaks_p, Peaks::ION_INDEX *magPeaks_p, int nMagPeaks);

        //Integrates two peak lists into a single list.
        //If the two lists are unequal in size, the largest of them is considered.
        //If they are equal, the old list (cnvPeaks2_p) is considered but the possibility that some peak of cnvPeaks1_p is isolated is aken into account: it is included.
        //If two peaks coincide in their max value (within the minScansPeaksGap range) the first one is included.
        //Arguments:
        // cnvPeaks1_p -> list of source and destination peaks
        // cnvPeaks2_p -> list of source peaks
        //Return true if OK
        bool peaksMix(CONVOLVED_PEAKS *cnvPeaks1_p, CONVOLVED_PEAKS *cnvPeaks2_p, int nMagPeaks);

        //insert a peak into an array of peaks ordered from smallest to largest if it respects the minimum distance from its neighbors.
        //there must be an extra position at the end of the array!!!
        //Arguments:
        // peak_p -> peak to insert
        // etpPeaks_p -> array of peaks for insertion
        // etpSize -> array size
        //Return the index to the insertion point
        int insertPeak(Peaks::ION_INDEX *etpPeak_p, Peaks::ION_INDEX *peak_p, int etpSize);

        //The array of CONVOLVED_PEAKS structures, of simple peaks, is analyzed in search of possible compound peaks.
        //A composite peak is a grouping of consecutive simple peaks joined by a valley of magnitude greater than m_noiseLevel.
        //The info is in the array m_unitedPeaks_p[] of size m_nUnitedPeaks. The size is returned.
        //Arguments:
        // cnvPeaks_p -> array of structures with simple peaks info
        // nMagPeaks -> array size
        //Returns the number of compound peaks
        int unitedPeaks(CONVOLVED_PEAKS *cnvPeaks_p, int nMagPeaks);


        SIREM_PEAKS_INFO *m_siremContex_p; //general info
        int             m_siremIndex;   //index to the last free element of the cuts array.
        SIREM_INFO      m_siremInfo;    //structure with info for SIREM.
        Sirem           *m_sirem_p;     //pointer to RelEntropy class.
        Peaks           *m_magPeaks_p;  //pointer to Peaks class for magnitude.
        CONVOLVED_PEAKS *m_cnvPeaks1_p, //support pointers.
                        *m_cnvPeaks2_p,
                        *m_cnvPeaks_p;  //pointer to array of structures for information on all simple peaks.
        UNITED_PEAKS    *m_unitedPeaks_p; //array with information of simple peaks that should be considered united.
        int             m_nUnitedPeaks; //number of entries in m_unitedPeaks_p.
        float           **m_siremCuts_p,//matrix with sirem info for each cut and each mz.
                        *m_meanMag_p;   //average magnitude values over each mz.
        int             m_nMagBlocks;   //#magnitude peaks.
};

#endif
