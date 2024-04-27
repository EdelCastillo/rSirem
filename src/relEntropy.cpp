/*************************************************************************
 *     Sirem base: Structure Inverse Relative Entropy Measure
 *     Copyright (C) agosto 2022 Esteban del Castillo Pérez
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


#include "relEntropy.h"

//Constructor
//initializes memory pointers and certain variables
RelEntropy::RelEntropy()
{
    m_tilesList_p=0;//the array of tiles does not exist
    m_evCount_p=0;

    m_SIREM=true; //the sirem algorithm is used and not the Shannon algorithm
    m_xOffset=m_yOffset=0; //without displacement of the tiles with respect to the origin of coordinates
    m_noiseLevel=0;  //no noise filter is run prior to the entropy analysis
    m_tile=T_4x4; //default tile
    for(int i=0; i<17; i++)
        m_split_p[i]=0;
    m_repeated.set=0;
}

//destructor
//memory is freed
RelEntropy::~RelEntropy()
{
//    printf("...relEntropy in\n");
    if(m_tilesList_p)    delete [] m_tilesList_p;
    if(m_evCount_p)     delete [] m_evCount_p;
    for(int i=0; i<17; i++)
            if(m_split_p[i]) delete [] m_split_p[i];
    if(m_repeated.set) delete [] m_repeated.set;

//    printf("...relEntropy out\n");
}

//receives a binary image corresponding to a section at the base of the primary image (maximum area image)
//sets the list of valid tiles (inside the image and filled with '1') and their size
//initialize the class
//Arguments:
//image_p => pointer to structure with information of the image to analyze
//Return -1 if I fail
int RelEntropy::initClass(IMAGE_ENTROPY *image_p)
{
    if(image_p==0)
        {printf("Error: no valid image\n"); return -1;}
    else
        m_image_p=image_p;

    //teselas
    setTile(image_p->tile);
    setNoiseLevel(image_p->noiseLevel);
    setTilesOffset(image_p->tilesOffset_x, image_p->tilesOffset_y);

    m_nXTiles=(m_image_p->nCols-m_xOffset)/m_tileSide; //tiles in X
    m_nYTiles=(m_image_p->nRows-m_yOffset)/m_tileSide; //tiles in Y
    m_totalTiles=m_nXTiles*m_nYTiles; //maximum size
    if(m_totalTiles<=0)
        {printf("warning: invalid image dimensions (%dx%d)\n", m_image_p->nRows, m_image_p->nCols); return -1;}
    try
    {
    //memoria para mantener las teselas completas interiores a la imagen
    if(m_tilesList_p) {delete []m_tilesList_p; m_tilesList_p=0;}
    m_tilesList_p=new bool[m_totalTiles];

    m_intoSampleTiles=0; //tiles inside the image

    //If a maximum area of the image is not provided, it is assumed equal to the entire given rectangle
    if(image_p->frame==REC_AREA)
        {
        for(int tile=0; tile<m_totalTiles; tile++)
            {m_tilesList_p[tile]=true; m_intoSampleTiles++; continue;}
        }
    else if(image_p->frame==MAX_AREA)//if a maximum area is provided
        {
        unsigned short mask_1=0x0001, mask_2=0x000F, mask_3=0x01FF, mask_4=0xFFFF, code;
        for(int tile=0; tile<m_totalTiles; tile++)
            {
            m_tilesList_p[tile]=false;
            code=getCodeTile(tile);      //code associated with the tile
            switch(m_tile)
                {
                case T_1x1:
                    {if(code==mask_1) {m_tilesList_p[tile]=true; m_intoSampleTiles++;}
                    break;}
                case T_2x2:
                    {if(code==mask_2) {m_tilesList_p[tile]=true; m_intoSampleTiles++;}
                    break;}
                case T_3x3:
                    {if(code==mask_3) {m_tilesList_p[tile]=true; m_intoSampleTiles++;}
                    break;}
                case T_4x4:
                    {if(code==mask_4) {m_tilesList_p[tile]=true; m_intoSampleTiles++;}
                    break;}
                }
            }
        }

    m_nPxTile=m_tileSide*m_tileSide;  //pixels on the tile

    m_nEvents=(int)ceil(pow(2.0, (double)m_nPxTile)); //Possible events (combinations of bits in the tile)
    if(m_evCount_p) {delete[]m_evCount_p; m_evCount_p=0;}
    m_evCount_p=new unsigned int[m_nEvents];

    //memory for code division structure according to its pixels
    //is reserved according to the combinations of 16 elements taken from n in n => 16!/[n!*(16-n)!]
    int sizes[]={16, 16, 120, 560, 1820, 4368, 8008, 11440, 12870, 11440, 8008, 4368, 1820, 560,
120, 16, 16};
    if(m_split_p[0]==0)   //if it has not been initialized before
        {
        for(int i=0; i<17; i++)
            m_split_p[i]=new SPLIT[sizes[i]];
        m_repeated.set=new int [m_intoSampleTiles+1]; //groupings of elements
        }
    m_repeated.size=0;

    }

    catch(const std::bad_alloc& e)
        {
        printf("Error reserving memory: %s\n",e.what());
        return -1;
        }
    return 0;
}

//Set the binary image for further analysis
//return -1 if the image is invalid (does not exist)
int RelEntropy::setRecBinaryImage(IMAGE_ENTROPY *image_p)
{
    if(image_p==0)
        {printf("Error: no valid image\n"); return -1;}
    else
        m_image_p=image_p;
    return 0;
}

//Sets the tile size used to analyze the image
void RelEntropy::setTile(Tile tile)
{
    switch(tile)
        {
        case T_1x1: m_tileSide=1; break;
        case T_2x2: m_tileSide=2; break;
        case T_3x3: m_tileSide=3; break;
        case T_4x4: m_tileSide=4; break;
        default:
            {
            printf("Warning: no valid tesela, so a 4x4-side_pixels tesela is assumed\n");
            m_tileSide=4;
            break;
            }
        }
    m_tile=tile;
}


//Determines the relative entropy of the image. Relative to the maximum possible.
//return the relative entropy in the range [0..1]:
//differential treatment is done according to the number of pixels that the tiles contain
//The lack of monotony at low concentrations of pixels in the image is corrected and the
//insensitivity of the measurement to noise is increased
float RelEntropy::getRelEntropy()
{
    m_intoSampleTiles=0;

    int canceledPxs=cancelImageNoise(); //removes isolated agglomerated pixels from the image

    //obtain the code of all the tiles in the image
    //the code is the sequence of bits that make up a tile
    m_intoSampleTiles=0;

    //determination of repeated events (tiles)
    for(unsigned int i=0; i<m_nEvents; i++)
        m_evCount_p[i]=0;  //resetting the repeat event counter

    //If no maximum area is provided, only tiles with pixels in the current image are considered.
    if(m_image_p->frame==ACTUAL_AREA)
        {
        unsigned short code=0;
        for(int tile=0; tile<m_totalTiles; tile++) //for each tile
            {
            code=getCodeTile(tile, 0);
            if(code==0) continue;
            m_evCount_p[code]++; //each position m_evCount_p[i] contains the times that event i has been found in the image
            m_intoSampleTiles++; //interior tiles.
            }
        }
    else //All tiles in the rectangle that includes the image or only those inside it are considered.
        {
        unsigned short code=0;
        for(int tile=0; tile<m_totalTiles; tile++) //for each tile
            {
            if(!m_tilesList_p[tile])//only the tiles completely within the image are of interest
                continue;
            code=getCodeTile(tile, 0);
             m_evCount_p[code]++; //each position m_evCount_p[i] contains the times that event i has been found in the image
            m_intoSampleTiles++; //interior tiles
            }
        }

    if(!m_image_p->SIREM) //if Shannon entropy is desired, unchanged
        {
        double H=0, prob;
        for(unsigned int j=0; j<m_nEvents; j++)
            {
            if(m_evCount_p[j]==0) continue;
            prob=(double)m_evCount_p[j]/m_intoSampleTiles;
            H+=prob*log2(prob); //entropy associated with the set of tiles with the same number of pixels
            }
        H/=log2(m_intoSampleTiles); //normalized entropy
        return -H; //the absolute value is returned
        }

    pxSplit(); //Repeated (non-exclusive) tiles are segregated according to their number of pixels

//    viewInfo();  //to view details of the processed data

    double relEntropy=0, prob;
    int newSize;
    m_totalExclusiveTiles=0;

    m_H[0]=0;  m_H[1]=0;
    for(int i=0; i<=m_nPxTile; i++)
        m_totalExclusiveTiles+=m_exclusiveSplit[i]; //total exclusive events (not repeated in the image)

    //elements of the repeated event sets that are noise are passed to the exclusive event set
    //shannon entropy non-monotonicity is eliminated and noise insensitivity is increased
    int tiles_p[17];
//    for(int i=0; i<17; i++) tiles_p[i]=0;

    //Noisy tiles in the image are identified and segregated according to the active pixels they contain.
    m_imageNoisyPxs=getNoisyPxs(tiles_p); //returns the total number of pxs considered noise
//m_imageNoisyPxs=0;

    //Noisy tiles join the list of exclusive tiles. Evaluated according to active pxs in tile
    for(int tilePxs=1; tilePxs<=m_nPxTile; tilePxs++)
        {
        int size=m_split_p[0][tilePxs].size; //SPLIT structure array size for tilePxs
        if(size==0) continue; //if there is no info
        int acu=0;
        int exclusiveTiles;
        newSize=0;

        exclusiveTiles=tiles_p[tilePxs]; //number of noisy tiles for tilePxs
        //Noisy ones are excluded from the list of tiles with the same number of pixels.
        for(int j=size-1; j>=0; j--)
            {
            int tmp=m_split_p[tilePxs][j].size; //tiles with the same code

            if(acu+tmp>exclusiveTiles) //If it cannot be excluded in its entirety, it must be divided
                {
                m_split_p[tilePxs][j].size-=exclusiveTiles-acu; //new size for number of repetitions of the event
                acu=exclusiveTiles; //tiles to exclude
                newSize=j+1; //index to list
                break;
                }
            else if(acu+tmp==exclusiveTiles) //the entire group is excluded
                {
                acu+=tmp;
                newSize=j; //index to list
                break;
                }
            else
                acu+=tmp; //acu maintains the number of tiles that are going to be separated
            }
        m_split_p[0][tilePxs].size=newSize;
        m_totalExclusiveTiles+=acu; //N0->tile that becomes part of the group of exclusive tiles
        }

    //Entropy of non-exclusive tiles
    getRepeatedGroups();  //speeds up the calculation
      for(int j=0; j<m_repeated.size; j++)
        {
        if(m_repeated.set[j]==0) continue;
        prob=(double)j/m_intoSampleTiles;
        if(prob>0)
            m_H[1]+=m_repeated.set[j]*prob*log2(prob); //entropy associated with the set of tiles with the same number of pixels
        else
            printf("ERROR in relEntropy class (probability<=0)\n");
        }

    //Entropy of exclusive tiles (empty+noise)
    prob=(double)m_intoSampleTiles;
    double log2prob=log2(prob);
    m_H[0]=(double)m_totalExclusiveTiles/prob*log2prob; //entropy associated with the set of exclusive tiles

    //entropía total
    relEntropy=(m_H[0]-m_H[1])/log2prob; //normalized total entropy
    return relEntropy; //returns the absolute value;
}

//establishes the maximum grouping of isolated pixels considered as noise:
//noiseLevel =0 -> no noise
//noiseLevel =1 -> groups of 1 pixel
//noiseLevel =2 -> groups of 3x3 pixels
//noiseLevel =3 -> groups of 5x5 pixels
void RelEntropy::setNoiseLevel(int noiseLevel) {m_noiseLevel=noiseLevel;}

//sets an offset on all image tiles
void RelEntropy::setTilesOffset(int x, int y)
{m_xOffset=x; m_yOffset=y;}

//returns the maximum number of active pixels included in the image
int RelEntropy::getMaxImagePixels() {return m_intoSampleTiles*m_nPxTile;}

//returns the number of tiles completely included in the image
int RelEntropy::getImageTiles() {return m_intoSampleTiles;}

//returns the active tile
Tile RelEntropy::getTile(){return m_tile;}

//Returns the part of the relative entropy corresponding to the exclusive tiles (noisiest empty ones)
//requires that before calling getRelEntropy()
float RelEntropy::getExclusiveRelEntropy() {return m_H[0]/log2(m_intoSampleTiles);}

//Returns the part of the relative entropy corresponding to non-exclusive tiles
//requires that before calling getRelEntropy()
float RelEntropy::getNonExclusiveRelEntropy() {return m_H[1]/log2(m_intoSampleTiles);}

//Returns the number of pixels considered noisy in the image
int RelEntropy::getImageNoisyPixels()  {return m_imageNoisyPxs;}

//returns the number of active and considered pixels in the image
int RelEntropy::getImageActivePixels()
{
    int totalActivePxs=0;
    //If no maximum area is provided, only tiles with pixels in the current image are considered.
    if(m_image_p->frame==ACTUAL_AREA)
        {
         for(int tile=0; tile<m_totalTiles; tile++) //for each tile
            totalActivePxs+=getPixelsTile(tile);
        }
    else //All tiles in the rectangle that includes the image or only those inside it are considered.
        {
        for(int tile=0; tile<m_totalTiles; tile++) //for each tile
            if(m_tilesList_p[tile])//only the tiles completely within the image are of interest
                totalActivePxs+=getPixelsTile(tile);
        }
return totalActivePxs;
}

//returns a code associated with the tile derived in a particular scan of the pixels that compose it.
//Code composition:
//horizontal axial sweep:
     //each pixel (bit) is extracted from the tile and placed in line starting from the upper left part of the tile (least weight bit).
     // read by rows from left to right. The least significant bit corresponds to the lower right corner of the tile.
//vertical axial sweep:
     //each pixel (bit) is extracted from the tile and placed in line starting from the upper left part of the tile (least weight bit).
     //it is read by columns from left to right. The least significant bit corresponds to the lower right corner of the tile.
//left-right diagonal sweep
     //each pixel (bit) is extracted from the tile and placed in line starting from the upper left part of the tile (least weight bit).
     //it is read following diagonal lines from left to right and from top to bottom. The least significant bit corresponds to the lower left corner of the tile.
//right-left diagonal sweep
     //each pixel (bit) is extracted from the tile and placed in line starting from the upper right part of the tile (least weight bit).
     //it is read following diagonal lines from right to left and from bottom to top. The least significant bit corresponds to the lower right corner of the tile.
//If the pixel is active, the corresponding bit is set to '1'
//Arguments:
//  tile -> tile of interest
//  type -> 0->axial horizontal; 1->vertical axial; 2->diagonal left-right; 3->diagonal right-left
unsigned short RelEntropy::getCodeTile(int tile, char type)
{
     int xIni, yIni, xTile, yTile;
   //coordinates of the upper left corner of the tile
    yTile=tile/m_nXTiles;
    xTile=tile%m_nXTiles;
    xIni=xTile*m_tileSide+m_xOffset;
    yIni=yTile*m_tileSide+m_yOffset;

    unsigned short code=0;
    int bitIndex=0;
    bool px;
    unsigned short byteMask;
    if(type==0)  //horizontal sweep
        {
        for(int y=0; y<m_tileSide; y++) //for each row
            for(int x=0; x<m_tileSide; x++)//for each column
                {
                px=m_image_p->pxMapImage_p[yIni+y][xIni+x]; //pixel
                byteMask=0x0001<<bitIndex; //the bit is positioned in the mask
                if(px==true)
                    code|= byteMask; //bit to '1' if active pixel
                bitIndex++;
                }
        }
    else if(type==1)  //barrido vertical
        {
        for(int x=0; x<m_tileSide; x++) //for each row
            for(int y=0; y<m_tileSide; y++)//for each column
                {
                px=m_image_p->pxMapImage_p[yIni+y][xIni+x]; //pixel
                byteMask=0x0001<<bitIndex; //the bit is positioned in the mask
                if(px==true)
                    code|= byteMask; //bit to '1' if active pixel
                bitIndex++;
                }
        }
    else if(type==2)  //barrido diagonal up-down
        {
        for(int f=0; f<m_tileSide; f++) //for each row
            for(int y=f, x=0; y>=0; y--, x++)//for each column
                {
                if(y>=m_tileSide) continue;
                px=m_image_p->pxMapImage_p[yIni+y][xIni+x]; //pixel
                byteMask=0x0001<<bitIndex; //the bit is positioned in the mask
                if(px==true)
                    code|= byteMask; //bit to '1' if active pixel
                bitIndex++;
                }
        for(int f=1; f<m_tileSide; f++) //for each row
            for(int y=f, x=m_tileSide-1; y<m_tileSide; y++, x--)//para cada columna
                {
                if(y>=m_tileSide) continue;
                px=m_image_p->pxMapImage_p[yIni+y][xIni+x]; //pixel
                byteMask=0x0001<<bitIndex; //the bit is positioned in the mask
                if(px==true)
                    code|= byteMask; //bit to '1' if active pixel
                bitIndex++;
                }
        }
    else if(type==3)  //barrido diagonal down-up
        {
        for(int f=m_tileSide-1; f>=0; f--) //for each row
            for(int y=0, x=f; x<m_tileSide; y++, x++)//for each column
                {
                if(y<0) continue;
                px=m_image_p->pxMapImage_p[yIni+y][xIni+x]; //pixel
                byteMask=0x0001<<bitIndex; //the bit is positioned in the mask
                if(px==true)
                    code|= byteMask; //bit to '1' if active pixel
                bitIndex++;
                }
        for(int f=1; f<m_tileSide; f++) //for each row
            for(int y=f, x=0; y<m_tileSide; y++, x++)//for each column
                {
                if(y<0) continue;
                px=m_image_p->pxMapImage_p[yIni+y][xIni+x]; //pixel
                byteMask=0x0001<<bitIndex; //the bit is positioned in the mask
                if(px==true)
                    code|= byteMask; //bit to '1' if active pixel
                bitIndex++;
                }

        }
    return code;
}

//Returns the number of bits to '1' in the tile
int RelEntropy::getPixelsTile(int tile)
{
     int xIni, yIni, xTile, yTile;
   //coordinates of the upper left corner of the tile
    yTile=tile/m_nXTiles;
    xTile=tile%m_nXTiles;
    xIni=xTile*m_tileSide;
    yIni=yTile*m_tileSide;

    int nPxs=0;
    bool px;
    for(int y=0; y<m_tileSide; y++) //for each row
        for(int x=0; x<m_tileSide; x++)//for each column
            {
            px=m_image_p->pxMapImage_p[yIni+y][xIni+x]; //pixel
            if(px==true)
                nPxs++; //bit to '1' if active pixel
            }
    return nPxs;
}

//Return the number of bits to '1' in the event
//the event must be of type unsigned short
int RelEntropy::getPixelsEvent(unsigned short event)
{
    unsigned short mask=0x0001;
    int index=0;
    int nPx=0;index=0;
    for(int j=0; j<m_nPxTile; j++)
           if(event & (mask<<index++)) nPx++;

    return nPx;
}

//returns true if the given tile is interior to the image as a whole
bool RelEntropy::isTileIntoImage(int tile)
{
    if(tile>=0 && tile<m_totalTiles)
        return m_tilesList_p[tile];
    return false;
}


//establishes a structure with information segregated from codes and sizes according to the number of pixels in the code
//m_split_p is an array of pointers to arrays of SPLIT structures of various sizes
//each SPLIT array maintains info associated with events with the same number of pixels.
//each SPLIT contains the code of a tile and the times it has been repeated in the image.
//The array[0] is different: it contains the indexes to the first free element of arrays 2 to 16 (repeated codes)
void RelEntropy::pxSplit()
{
    int nPxEvent, index;
    for(int i=0; i<17; i++) //The maximum tile considered is 4x4 = 16 possibilities in the number of pixels + the empty set
        {
        m_split_p[0][i].size=0; //zeroing index list to first free element
        m_split_p[1][i].size=0;
        m_exclusiveSplit[i]=0;
        }
    for(unsigned int event=0; event<m_nEvents; event++)
        {
        if(m_evCount_p[event]==0) continue;
        if(event==0) {m_exclusiveSplit[0]=m_evCount_p[event]; continue;}

        nPxEvent=getPixelsEvent(event); //pixels in code
        if(m_evCount_p[event]==1)
            m_exclusiveSplit[nPxEvent]++; //tiles of one px have special treatment: they become exclusive
        else
            {
            index=m_split_p[0][nPxEvent].size; //indexes to the first free element of each list
            m_split_p[nPxEvent][index].code=event; //code is updated
            m_split_p[nPxEvent][index].size=m_evCount_p[event]; //replays are updated
            m_split_p[0][nPxEvent].size++;
            }
        }

}

//establishes the array m_eventsByTilePxs[] indicating the number of repeated tiles in the image according to the number of active pixels
//used to present results
//return: the size of the array;
int RelEntropy::getEventsByTilePxs()
{
    int nPxEvent;
    for(int i=0; i<17; i++) //The maximum tile considered is 4x4 = 16 possibilities in the number of pixels + the empty set
        m_eventsByTilePxs[i]=0;

    for(unsigned int event=0; event<m_nEvents; event++)
        {
        if(m_evCount_p[event]==0) continue; //event not occurred => without interest
        nPxEvent=getPixelsEvent(event); //pixels in code
        m_eventsByTilePxs[nPxEvent]+=m_evCount_p[event]; //tiles in the image that contain the same number of active pixels
        }
   return 17;
}

//Handle m_split_p structures and set a global repeat count.
//The position in the array indicates the number of tiles with the same code, its content indicates
//the total number of tiles that repeat that pattern
//return: the size of the array
//Useful to speed up the calculation of entropy
int RelEntropy::getRepeatedGroups()
{
    int max=0;;
    for(int tilePxs=1; tilePxs<=m_nPxTile; tilePxs++)
        {
        for(int i=0; i<m_split_p[0][tilePxs].size; i++)
            if(m_split_p[tilePxs][i].size>max)
                max=m_split_p[tilePxs][i].size;
        }
    max++;

    for(int i=0; i<max; i++)
        m_repeated.set[i]=0;  //repeat counter initialization

    for(int tilePxs=1; tilePxs<=m_nPxTile; tilePxs++)
        for(int i=0; i<m_split_p[0][tilePxs].size; i++)
            m_repeated.set[m_split_p[tilePxs][i].size]++; //count the repetitions of each entry
    m_repeated.size=max;
    return max;

}

//returns the number of tiles in the image with the given number of active pixels
//arguments:
//  pixels=number of active pixels in the tile
//return
//  -1 if the argument is invalid
//  the number of tiles requested
int RelEntropy::tilesWithPixels(int pixels)
{
    if(pixels>m_nPxTile || pixels<0)
        return -1;
    int acu=0;
    for(unsigned int i=0; i<m_nEvents; i++)
        {
        if(getPixelsEvent(i)==pixels) acu+=m_evCount_p[i]; //number of tiles with the same number of pixels
        }
    return acu;
}


//presents detailed information, segregated according to the number of active pixels of the tiles
void RelEntropy::viewInfo()
{
    //number of non-repeated tiles (exclusive)
    printf("format: [nPx]number:\n");
    printf("exclusive:");
    for(int j=0; j<=m_nPxTile; j++) //based on the active pixels of the tiles
        printf("[%2d]%4d ",j, m_exclusiveSplit[j]);
    printf("\n");

    getEventsByTilePxs();
    printf(" repeated:");
    for(int j=0; j<=m_nPxTile; j++) //based on the active pixels of the tiles
        printf("[%2d]%4d ",j, m_eventsByTilePxs[j]-m_exclusiveSplit[j]);
    printf("\n");

    //number of total tiles (repeated + exclusive)
    printf("    total:");
    for(int j=0; j<=m_nPxTile; j++) //based on the active pixels of the tiles
        printf("[%2d]%4d ",j, m_eventsByTilePxs[j]);
    printf("\n");

    //for each tile size, number of tiles repeated on the indicated event
    printf(" repeated list format: [nPx,elements]->[repeated]0x_code\n");

    for(int i=1; i<=m_nPxTile; i++) //based on the active pixels of the tiles
    {
        int size=m_split_p[0][i].size; //elements in the array
        if(size==0) continue;

        int suma=0;
        printf("[%d,%d]->", i, size);
        for(int j=0; j<size; j++)
            {
            suma+=m_split_p[i][j].size; //total tiles repeated for that number of active pixels
            printf("[%d]%0x ",m_split_p[i][j].size, m_split_p[i][j].code);  //repeated tiles according to the code
            }
        printf("\n\trepeated tiles vs. total tiles:%d/%d\n", suma, m_intoSampleTiles);
    }
}

//sets the number of noisy tiles in the image ordered by number of active pixels (tiles_p)
//returns the number of pixels in the image considered noisy
//Clusters of pixels are detected based on four different scans of the tiles: two axial and two diagonal
//If these clusters do not exist, it is considered a noise tile
//Return:
//  Total number of noisy tiles
int RelEntropy::getNoisyPxs(int *tiles_p)
{
    unsigned short A;
    int pixelsTile, noisyPixelsTile, minNoisyPxs, totalNoisyPxs=0, direction;
    bool isNoisyTile, noisyTile=false;

    for(int i=0; i<=m_nPxTile; i++) tiles_p[i]=0;

    for(int tile=0; tile<m_totalTiles; tile++)
        {
        if(m_image_p->frame==ACTUAL_AREA || m_tilesList_p[tile])//if a maximum area is not provided or it is an interior tile
            {
            minNoisyPxs=0xFFFF;
            for(int i=0; i<4; i++) //codes for each of the 4 orientations
                {
                A=getCodeTile(tile, i); //The code associated with the tile is extracted (horizontal, vertical, etc.)
                isNoisyTile=getTileNoisyPxs(A, &pixelsTile, &noisyPixelsTile); //Is it noisy tile? returns the total pixels and noise pixels
                if(noisyPixelsTile<minNoisyPxs)
                    {
                    minNoisyPxs=noisyPixelsTile; //The minimum is of interest: they may be crowded in one orientation.
                    if(isNoisyTile) noisyTile=true; //in one orientation can be non-noisy tile
                    else noisyTile=false;
                    direction=i; //maximum crowding orientation
                    }
                }
            if(noisyTile) //noisy tile
                {
                tiles_p[pixelsTile]++; //noisy tiles are noted based on the number of pixels
                totalNoisyPxs+=pixelsTile; //Noisy pixels in the image are noted
                }
            }
        }
    return totalNoisyPxs;
}


//Returns the number of pixels considered as noise from the code associated with a tile.
//For the given code, count ones and zeros flanked by ones: ex: 001011101000-> 5 '1' and 2 '0'
//Those tiles that have one or two px and those that have
//  less than half of px and have any zeros in the sequence flanked by ones
//Arguments:
//  code          -> tile code
//  totalPixels   -> return the number of total active pixels on the tile
//  isolatedPixels-> return the number of isolated pixels in tile (with zeros before and after on the line)
//Return:
//  true if is noisy tile
bool RelEntropy::getTileNoisyPxs(unsigned short code, int *totalPixels, int *isolatedPixels)
{

    int ini, end, count1, count0;
    unsigned short A=0x0001;
    ini=-1; end=-1; count1=0; count0=0;
    *totalPixels=0;

    for(int i=0; i<m_nPxTile; i++) //search from the right of the first bit
        if(code & (A<<i))
            {ini=i; break;}

    for(int i=m_nPxTile-1; i>=0; i--) //search from the left of the first bit
        if(code & (A<<i))
            {end=i; break;}

    if(ini==-1 || end==-1)
        return 0; //no pixels on the tile

    for(int i=ini; i<=end; i++)
        {
        if(code & (A<<i)) count1++; //active pixels ('1')
        else count0++;              //zeros
        }

    //if the tile is less than 9 pixels, it is not considered noisy
    if(m_nPxTile<9)
        {*totalPixels=count1; isolatedPixels=0; return false;}

    //analysis of isolated pixels (with zeros before and after on the line)
    *totalPixels=count1;
    A=0xE000;
    unsigned short B=4000;
    int isolated=0;
    if((code & (B>>ini))==0) isolated++; //the first isolated pixel is analyzed
    for(int i=ini; i<end-3; i++)
        {
        if((A&code)==B) isolated++;
        A>>=i; B>>=i;
        }
    if((code & (B>>(end-2)))==0) isolated++; //the last isolated px is analyzed
    *isolatedPixels=isolated;

    //If there are only one or two pixels in the tile, they are considered noise
    //also, if the number of '1' agglomerates is less than half of the tile and contains some zero
    if(count1==1 || count1==2 || (count1<ceil(m_nPxTile/2.0) && count0>0)) //noise tile
        return true;

    return false;
}

//Groups of isolated pixels are canceled
// noiseLevel indicates the size of the grid considered as the maximum group of isolated pixels
// i.e. all pixels located in that grid with no pixels in the neighborhood are removed
// noiseLevel=0 ->does not exist
// noiseLevel=1 ->a single isolated pixel
// noiseLevel=2 ->isolated 3x3 grid
// noiseLevel=3 ->isolated 5x5 grid
//Returns the number of pixels removed from the image
//Return -1 if level<1
int RelEntropy::cancelImageNoise()
{
    PIXEL_2D *pxA_p=0, *pxB_p=0, *tmp_p;
    int nPxsIn, nPxsOut, count=0, size;
    if(m_noiseLevel<1) return -1;
    size=((m_noiseLevel-1)*2+3)*((m_noiseLevel-1)*2+3); //size of memory needed
    try
        {
        pxA_p=new PIXEL_2D[size]; pxB_p=new PIXEL_2D[size];
        }
    catch(const std::bad_alloc& e)
        {
        printf("Error reserving memory: %s\n",e.what());
        return -1;
        }

    for(int row=0; row<m_image_p->nRows; row++) //for all pixels in the image
        for(int col=0; col<m_image_p->nCols; col++)
            {
            if(m_image_p->pxMapImage_p[row][col]==false) continue; //if it is inactive pixel
            pxA_p[0].x=col; pxA_p[0].y=row; //start of sequence, base pixel
            nPxsIn=1;
            //chained search for new pixels. It ends when it detects an isolated block of size <= m_noiseLevel
            for(int i=1; i<=m_noiseLevel; i++)
                {
                nPxsOut=searchPxs(pxA_p, nPxsIn, pxB_p);
                if(nPxsOut==nPxsIn) //isolated package -> is removed
                    {
                    for(int i=0; i<nPxsIn; i++)
                        m_image_p->pxMapImage_p[pxA_p[i].y][pxA_p[i].x]=false; //canceled pixel
                    count+=nPxsIn;
                    break;
                    }
                nPxsIn=nPxsOut; //setting variables to iterate
                tmp_p=pxA_p;
                pxA_p=pxB_p;
                pxB_p=tmp_p;
                }
            }
    if(pxA_p) delete [] pxA_p;
    if(pxB_p) delete [] pxB_p;
    return count; //canceled pixels in the image
}


//Set in pxOut_p[] the coordinates of all the pxs included in a grid higher than the input grid
//The coordinates of the input px include the coordinates of all the px located on its periphery.
//Arguments:
//  pxIn_p[] -> array with the coordinates of the input pixels
//  pxOut_p[] -> array with the coordinates of the output pixels
//  nPxIn -> elements in pxIn_p[]
//return the size of pxOut_p[]
//  pxOut_p[] must have space reserved to accommodate (level-1)*2+3)^2 elements
int RelEntropy::searchPxs(PIXEL_2D *pxIn_p, int nPxIn, PIXEL_2D *pxOut_p)
{
    int X, Y, nNeigh, indexOut;
    PIXEL_2D neighbours_p[9];
    bool hit;
    for(int px=0; px<nPxIn; px++)
        {
        pxOut_p[px].x=pxIn_p[px].x; //we copy everything to the destination, necessary to detect new pixels
        pxOut_p[px].y=pxIn_p[px].y;
        }
    indexOut=nPxIn;

    for(int px=0; px<nPxIn; px++)
        {
        X=pxIn_p[px].x; Y=pxIn_p[px].y;
        if(m_image_p->pxMapImage_p[Y][X]==false) continue; //inactive pixel

        nNeigh=getNeighbours(Y, X, neighbours_p); //neighborhood of X/Y

        //Neighbor pixels that are not already in the output list are added
        for(int i=0; i<9; i++)
            {
            if(neighbours_p[i].x==-1 || i==4) continue; //inactive or central pixel
            hit=false;
            int tmpSize=indexOut;
            for(int j=0; j<tmpSize; j++) //it is analyzed if it is px not included in the list
                {
                if(neighbours_p[i].y==pxOut_p[j].y && neighbours_p[i].x==pxOut_p[j].x)
                    {hit=true; break;}
                }
            if(!hit)
                {
                pxOut_p[indexOut].x=neighbours_p[i].x; pxOut_p[indexOut].y=neighbours_p[i].y; //nuevo px
                indexOut++;
                }
            }
        }
    return indexOut; //output list size
}

//Set neighbor_p[] to list of active px on a 3x3 grid whose center is row/col
//Each element of neighbor_p[] is the coordinates of the active pixel.
//If the pixel is not activated, the coordinates are -1/-1
//Arguments:
// row -> row of the pixel image to analyze
// col -> column of the pixel image to analyze
// neighbor_p[] -> array of pixels in the 3x3 grid ordered from top to bottom and from left to right
//return the number of active pixels in the grid, including px
int RelEntropy::getNeighbours(int row, int col, PIXEL_2D *neighbour_p)
{
    int count=0;
    int maxCol=m_image_p->nCols-1;
    int maxRow=m_image_p->nRows-1;

    if(row>0 && col>0 && m_image_p->pxMapImage_p[row-1][col-1]==true)
        {neighbour_p[0].x=col-1; neighbour_p[0].y=row-1;}
    else
        {neighbour_p[0].x=-1; neighbour_p[0].y=-1;}

    if(row>0 && col>=0 && m_image_p->pxMapImage_p[row-1][col+0]==true)
        {neighbour_p[1].x=col+0; neighbour_p[1].y=row-1;}
    else
        {neighbour_p[1].x=-1; neighbour_p[1].y=-1;}

    if(row>0 && col<maxCol && m_image_p->pxMapImage_p[row-1][col+1]==true)
        {neighbour_p[2].x=col+1; neighbour_p[2].y=row-1;}
    else
        {neighbour_p[2].x=-1; neighbour_p[2].y=-1;}

    if(row>=0 && col>0 && m_image_p->pxMapImage_p[row+0][col-1]==true)
        {neighbour_p[3].x=col-1; neighbour_p[3].y=row+0;}
    else
        {neighbour_p[3].x=-1; neighbour_p[3].y=-1;}

    if(row>=0 && col>=0 && m_image_p->pxMapImage_p[row+0][col+0]==true)
        {neighbour_p[4].x=col+0; neighbour_p[4].y=row+0;}
    else
        {neighbour_p[4].x=-1; neighbour_p[4].y=-1;}

    if(row>=0 && col<maxCol && m_image_p->pxMapImage_p[row+0][col+1]==true)
        {neighbour_p[5].x=col+1; neighbour_p[5].y=row+0;}
    else
        {neighbour_p[5].x=-1; neighbour_p[5].y=-1;}

    if(row<maxRow && col>0 && m_image_p->pxMapImage_p[row+1][col-1]==true)
        {neighbour_p[6].x=col-1; neighbour_p[6].y=row+1;}
    else
        {neighbour_p[6].x=-1; neighbour_p[6].y=-1;}

    if(row<maxRow && col>=0 && m_image_p->pxMapImage_p[row+1][col+0]==true)
        {neighbour_p[7].x=col+0; neighbour_p[7].y=row+1;}
    else
        {neighbour_p[7].x=-1; neighbour_p[7].y=-1;}

    if(row<maxRow && col<maxCol && m_image_p->pxMapImage_p[row+1][col+1]==true)
        {neighbour_p[8].x=col+1; neighbour_p[8].y=row+1;}
    else
        {neighbour_p[8].x=-1; neighbour_p[8].y=-1;}

    for(int i=0; i<9; i++)
        if(neighbour_p[i].x!=-1) count++;
    return count;
}

//deactivates the active pixels in the tile considered noisy
//return:
// the number of tiles considered noisy
int RelEntropy::removeImageNoise()
{
    unsigned short A;
    int pixelsTile, noisyPixelsTile, minNoisyPxs, totalNoisyTiles=0, direction;
    bool isNoisyTile, noisyTile=false;
    int xTile, yTile, xIni, yIni;

    for(int tile=0; tile<m_totalTiles; tile++)
        {
        if(m_image_p->frame==ACTUAL_AREA || m_tilesList_p[tile])//if a maximum area is not provided or it is an interior tile
            {
            minNoisyPxs=0x7FFF;
            for(int i=0; i<4; i++) //codes for each of the 4 orientations
                {
                A=getCodeTile(tile, i); //The code associated with the tile is extracted (horizontal, vertical, etc.)
                isNoisyTile=getTileNoisyPxs(A, &pixelsTile, &noisyPixelsTile); //Is it noisy tile? returns the total pixels and the noise
                if(noisyPixelsTile<minNoisyPxs)
                    {
                    minNoisyPxs=noisyPixelsTile; //The minimum is of interest: they may be crowded in one orientation.
                    if(isNoisyTile) noisyTile=true; //in one orientation can be non-noisy tile
                    else noisyTile=false;
                    direction=i; //maximum crowding orientation
                    }
                }
            if(noisyTile) //noisy tile, its pixels are disabled
                {
                yTile=tile/m_nXTiles;
                xTile=tile%m_nXTiles;
                xIni=xTile*m_tileSide+m_xOffset;
                yIni=yTile*m_tileSide+m_yOffset;
                for(int y=0; y<m_tileSide; y++)
                    for(int x=0; x<m_tileSide; x++)
                       m_image_p->pxMapImage_p[yIni+y][xIni+x]=false;
                totalNoisyTiles++;
                }
            }
        }
    return totalNoisyTiles;
}




