#ifndef SIREM_H
#define SIREM_H

//#include "types.h"
#include "stdio.h"
#include "siremTypes.h"
#include "relEntropy.h"
#include <stdexcept>

  typedef struct
    {
    GROUP_F     goodMz;             //representative ion of the sample
    GROUP_XY    pxCoord;            //pixel coordinates
    char        tileSide;           //tile size for entropy
    bool        SIREM;              //true if SIREM algorithm is used; false if using Shannon entropy
    }SIREM_INFO;

//Class derived from class RelEntropy.
//From the magnitude information about each pixel of a scan and a cut level, generate a binary image.
//On it, the SIREM or relative entropy algorithm is applied and the measurement is returned.
//Attends to the level of noise and concentration of active pixels in the binary image
    
//Adapt the RelEntropy class to deal with the information of an ion
class Sirem: public RelEntropy
{
public:     
    //Constructor
    //Initialize the class
    Sirem(SIREM_INFO *sirem_p);
    
    //destructor:
    //free reserved memory
    ~Sirem();
    
    //Returns the SIREM value corresponding to the image passed in mzInfo_p
    //Arguments:
    //  mzInfo_p -> pointer to the magnitude of each pixel
    //  cutLevel -> cut level (percentile)
    //  noiseLevel -> Lower values are discarded for determining the cutoff level
    //  minDensity -> the binary image whose relative density(%) is lower, are discarded
    float getSirem(GROUP_F *mzInfo_p, float cutLevel, float noiseLevel, float minDensity);
    
    //returns the density of the binary image. Use after getSirem()
    float getBinaryImageDensity();
    
private:
    //free reserved memory.
    void outApp();
    
    //Initialize certain parameters.
    //The tiles with information to be used by the base class are established:
    //To do this, an ion is sectioned by its base.
    // It is suggested to use an ion with a large average magnitude that is present in the entire sample
    //Arguments:
    //  tile -> tile to use
    //Return:
    //  -1 if base class initialization failed
    int initClass(Tile tile);
    
    //leave in the m_imageSection array the pixels whose magnitude exceeds or equals the passed cut level (percentile)
    //Arguments:
    //  mzInfo_p -> pointer to the magnitude of each pixel
    //  cutLevel -> cut level (percentile)
    //  noiseLevel -> Lower values are discarded for determining the cutoff level
    //  minDensity -> the binary image whose relative density(%) is lower, are discarded
    //Return:
    //  the relative density
    float getSection(GROUP_F *mzInfo_p, float cutLevel, float noiseLevel);
    
    //returns the data value that makes the total of lower values match 'percentile', in percentage
    //the returned value is interpolated between the two neighboring quantities
    //Arguments:
    //  data -> array of data to evaluate
    //  size -> array size
    //  percentil -> desired percentile
    //  minvalue -> values less than minValue are not considered
    //Return:
    //  the percentile or -1 if memory allocation failed
    float getValueFromPercentil(float *data, int size, float percentil, float minValue);
    
    //returns in sort a chain of indexes to bufferIn such that the sequence is ordered from smallest to largest magnitude
    //Arguments:
    //  bufferIn -> input array
    //  sort     ->array with indexes to ordered data
    //  size     -> array size
    //Return:
    //  -1 if memory allocation failed
    int sortUp(float *bufferIn, int *sort, int size);
    
    //converts the array of pixels from a section 'm_imageSection' into a binary image in XY coordinates that is passed to the base class
    //return -1 if the info used is not good (m_imageEntropy.pxMapImage_p==0 || m_imageSection.size<0)
    int imageSection2binary();
    
    //overloaded: avoids the use of the argument
    int setRecBaseBinaryImage();
    
    IMAGE_ENTROPY       m_imageEntropy; //structure that houses the binary image from which entropy is extracted
    GROUP               m_imageSection; //set of pixels that make up a section
    SIREM_INFO          *m_sirem_p;     //pointer to structure passed to constructor
    float               m_density;      //binary image density
};

#endif
