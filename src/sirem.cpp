/*************************************************************************
 *     Sirem: Structure Inverse Relative Entropy Measure
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
#include "sirem.h"

//Constructor
//Initialize the class
Sirem::Sirem(SIREM_INFO *sirem_p)
{
    if(sirem_p==0) return;
    m_sirem_p=sirem_p;
    m_density=0;

    //the corners of the rectangle that houses the image are determined
    int maxX=0, maxY=0, minX=0x7FFFFFFF, minY=0x7FFFFFFF;
    for(int i=0; i<sirem_p->pxCoord.size; i++)
        {
        if(sirem_p->pxCoord.set[i].x>maxX) maxX=sirem_p->pxCoord.set[i].x;
        if(sirem_p->pxCoord.set[i].x<minX) minX=sirem_p->pxCoord.set[i].x;
        if(sirem_p->pxCoord.set[i].y>maxY) maxY=sirem_p->pxCoord.set[i].y;
        if(sirem_p->pxCoord.set[i].y<minY) minY=sirem_p->pxCoord.set[i].y;
        }

    //rectangle that houses the image
    m_imageEntropy.nCols=maxX-minX+1;
    m_imageEntropy.nRows=maxY-minY+1;
    m_imageSection.set=0;
    m_imageEntropy.pxMapImage_p=0;

    //memory reserve
    try
        {
        //memory to accommodate pixels greater than or equal to the cutoff level
        m_imageSection.set=new int[sirem_p->pxCoord.size];

        //memory to accommodate the binary image resulting from the section
        m_imageEntropy.pxMapImage_p=new bool*[m_imageEntropy.nRows];
        for(int y=0; y<m_imageEntropy.nRows; y++)
            {
            m_imageEntropy.pxMapImage_p[y]=0;
            m_imageEntropy.pxMapImage_p[y]=new bool[m_imageEntropy.nCols];
            }
        }
    catch(const std::bad_alloc& e)
        {
        printf("Error reserving memory: %s\n",e.what());
        m_error=true;
        return;
        }
    //tile size
    if(sirem_p->tileSide>0 && sirem_p->tileSide<=4)
        initClass((Tile)(sirem_p->tileSide-1)); //default initialization
    else
        initClass(T_4x4); //default initialization
}

//destructor:
//free reserved memory
Sirem::~Sirem()
{
    outApp();
}

//free reserved memory
void Sirem::outApp()
{
//printf("...ini Sirem destructor\n");
    if(m_imageSection.set) delete [] m_imageSection.set;
    if(m_imageEntropy.pxMapImage_p)
        {
        for(int y=0; y<m_imageEntropy.nRows; y++)
            if(m_imageEntropy.pxMapImage_p[y])      delete[] m_imageEntropy.pxMapImage_p[y];
        delete [] m_imageEntropy.pxMapImage_p;
        }
//printf("...end Sirem destructor\n");
}

//Initialize certain parameters.
//The tiles with information to be used by the base class are established:
//To do this, an ion is sectioned by its base.
// It is suggested to use an ion with a large average magnitude that is present in the entire sample
int Sirem::initClass(Tile tile)
{
    if(tile<0 || tile>3)
        {
        printf("Warning: unrecognized tile in Sirem class. Set as %dx%d\n", tile+1, tile+1);
        tile=T_4x4;
        }

    setTile(tile);

    m_imageEntropy.frame=MAX_AREA;
    m_imageEntropy.tile=tile;
    m_imageEntropy.tilesOffset_x=0; //tessellation starts in the corner
    m_imageEntropy.tilesOffset_y=0;
    m_imageEntropy.SIREM=m_sirem_p->SIREM;
    m_imageEntropy.noiseLevel=0; //There are several levels of image denoising. Not used
    setNoiseLevel(0);   //base class: no noise removed before parsing

    //if a reference image is not provided, all pixels are considered active.
    if(m_sirem_p->goodMz.set==0 || m_sirem_p->goodMz.size==0)
        {
        for(int j=0; j<m_sirem_p->pxCoord.size; j++) //for all pixels
            m_imageSection.set[j]=j;
        m_imageSection.size=m_sirem_p->pxCoord.size;
        }
    //if a reference image is provided.
    //the reference image is sectioned at the base to obtain the active pixels (with non-zero magnitude).
    //empty tiles are avoided by including a non-zero noise level.
    else
        getSection(&m_sirem_p->goodMz, 1.0, 1e-6); //base section (1% for magnitude percentile) => m_imageSection

    //The set of pixels that make up the section are converted to a binary image
    imageSection2binary();

    //the base class is initialized
    if(setRecBaseBinaryImage()) {m_error=true, outApp(); return -1;}

    return 0;
}

//Returns the SIREM value corresponding to the image passed in mzInfo_p
//Arguments:
//  mzInfo_p   -> pointer to the magnitude of each pixel
//  cutLevel   -> cut level (percentile)
//  noiseLevel -> Lower values are discarded for determining the cutoff level
//  minDensity -> the binary image whose relative density(%) is lower, are discarded
float Sirem::getSirem(GROUP_F *mzInfo_p, float cutLevel, float noiseLevel, float minDensity)
{
   m_density=getSection(mzInfo_p, cutLevel, noiseLevel);
   if(isnan(m_density) || m_density<0) m_density=0; //limit control
   if(m_density<minDensity/100.0) return 0.0;

    //The set of pixels that make up the section are converted to a binary image
    imageSection2binary();
    
    //get SIREM
    float relEntropy=RelEntropy::getRelEntropy();
    if(isnan(relEntropy)) //entropy cannot be determined
            relEntropy=1.0; //so that its inverse is 0

    //the inverse is done: 0=maximum entropy, 1=zero entropy
    return 1.0-relEntropy;
}

//Leave in the m_imageSection array the pixels whose magnitude exceeds or equals the passed cut level (percentile)
//Arguments:
//  mzInfo_p -> pointer to the magnitude of each pixel
//  cutLevel -> cut level (percentile)
//  noiseLevel -> Lower values are discarded for determining the cutoff level
//  minDensity -> the binary image whose relative density(%) is lower, are discarded
//Return:
//  the relative density
float Sirem::getSection(GROUP_F *mzInfo_p, float cutLevel, float noiseLevel)
{
    //The cutoff level to be applied is determined according to the percentile
    //only pixels with magnitude greater than noiseLevel are considered
    float level=getValueFromPercentil(mzInfo_p->set, mzInfo_p->size, cutLevel, noiseLevel); //exclude zeros
    if(level==0.0) level=1e-6; //avoid null value
    else if(level==-1) {m_imageSection.size=0; return 0;} //evita imágenes con valor máximo menor que el ruido
    
    //The list of pixels that exceed or equal in magnitude to the cutoff level is generated.
    m_imageSection.size=0; //start
    for(int j=0; j<mzInfo_p->size; j++) //for all pixels
        {
        if(mzInfo_p->set[j]>=level)
            m_imageSection.set[m_imageSection.size++]=j;
        }
    return (float)m_imageSection.size/mzInfo_p->size; //returns the relative density
}

//returns the data value that makes the total of lower values match 'percentile', in percentage
//the returned value is interpolated between the two neighboring quantities
//Arguments:
//  data -> array of data to evaluate
//  size -> array size
//  percentil -> desired percentile
//  minvalue -> values less than minValue are not considered
//Return:
//  the percentile or -1 if memory allocation failed
float Sirem::getValueFromPercentil(float *data, int size, float percentil, float minValue)
 {
     int *sortDataIndex_p=0;
     float a, b, x;
     int ia, ib, indexLow=0;

     if(size<1) return 0; //limit control
     try
        {
        sortDataIndex_p=new int[size];
        }
     catch(const std::bad_alloc& e)
        {
        printf("Error reserving memory: %s\n",e.what());
        return -1;
        }

     sortUp(data, sortDataIndex_p, size);
     if(minValue>0.0)
        {
        for(indexLow=0; indexLow<size; indexLow++)
            {
            if(data[sortDataIndex_p[indexLow]]>minValue) break;
            }
        if(indexLow==size-1) return -1.0;
        else
            size-=indexLow;
        }
     float indexSort=size*percentil/100.0 + indexLow;

     //interpolation and limit control
     ia=trunc(indexSort);  //previous index
     ib=ceil(indexSort);   //previous index
     if(ib>=size)  //limit control
        {ia-=1; ib-=1;}
     if(ia<0 || indexSort<0)
        return data[0];

     a=data[sortDataIndex_p[ia]];
     b=data[sortDataIndex_p[ib]];
     x=a+(b-a)*(indexSort-ia);

    if(sortDataIndex_p) delete []sortDataIndex_p;
    return x;
 }

//returns in sort a chain of indexes to bufferIn such that the sequence is ordered from smallest to largest magnitude
//Arguments:
//  bufferIn -> input array
//  sort     ->array with indexes to ordered data
//  size     -> array size
//Return:
//  -1 if memory allocation failed
int Sirem::sortUp(float *bufferIn, int *sort, int size)
{
    if(size<=0) return 0;
    int k=0, j;
    float  minor;
    float *list=0;
    try
        {
        list=new float[size];
        }
    catch(const std::bad_alloc& e)
        {
        printf("Error reserving memory: %s\n",e.what());
        return -1;
        }

    for(j=0; j<size; j++)
        list[j]=bufferIn[j]; //copy as changes occur in the content

    if(size==1) //If there is only one, it is already ordered. If there are none -> do nothing
        *sort=0;
    else
    do
        {
        minor=(float)0x7FFFFFFF; //minor>maxValue. The value greater than the maximum
        for(j=0; j<size; j++)
            {
            if((list[j])<minor)
                {
                minor=list[j];
                sort[k]=j;
                }
            }
        list[sort[k]]=(float)0x7FFFFFFF;
        k++;
        }
    while(k<size);
    if(list) delete [] list;
    return 0;
}

//converts the array of pixels from a section 'm_imageSection' into a binary image in XY coordinates
//return -1 if the info used is not good (m_imageEntropy.pxMapImage_p==0 || m_imageSection.size<0)
int Sirem::imageSection2binary()
{
    //verification that the necessary resources exist
    if(m_imageEntropy.pxMapImage_p==0 || m_imageSection.size<0) return -1;

    //image zeroing. Avoid influences from previous images
    for(int y=0; y<m_imageEntropy.nRows; y++)
        for(int x=0; x<m_imageEntropy.nCols; x++)
            m_imageEntropy.pxMapImage_p[y][x]=false;

    //Converting pixel list to image in XY coordinates
    for(int pxi=0; pxi<m_imageSection.size; pxi++)
        {
        int px=m_imageSection.set[pxi]; //pixel
        int x=m_sirem_p->pxCoord.set[px].x; //pixel coordinates
        int y=m_sirem_p->pxCoord.set[px].y;
        m_imageEntropy.pxMapImage_p[y][x]=true; //binary pixel image
        }
   RelEntropy::setRecBinaryImage(&m_imageEntropy); //image is set in base class
   return 0;
}

//overloaded: avoids the use of the argument
int Sirem::setRecBaseBinaryImage()
{
    return RelEntropy::initClass(&m_imageEntropy);
}

//returns the density of the binary image. Use after getSirem()
float Sirem::getBinaryImageDensity() {return m_density;}
