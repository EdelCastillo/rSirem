/*************************************************************************
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

#ifndef IMAGE_REL_ENTROPY_H
#define IMAGE_REL_ENTROPY_H

#include <math.h>
#include <stdio.h>
#include <stdexcept>
#include "siremTypes.h"

    //accepted tile types
    enum Tile{
            T_1x1,      //tile of 1x1
            T_2x2,      //tile of 2x2
            T_3x3,      //tile of 3x3
            T_4x4       //tile of 4x4
    };
    
    enum  FRAME2{
            REC_AREA,   //all tiles inside the given rectangle are considered (nRows x nCols)
            MAX_AREA,   //The tiles obtained from the past image are considered as a reference.
                        //If no reference image is passed, all pixels are considered active.
            ACTUAL_AREA //tiles with activated pixels in the current image are considered
    };
    
    //basic information for the class
    typedef struct{
        FRAME2   frame;  //frame type
        Tile    tile;   //frame type
        char    noiseLevel; //establishes the maximum grouping of isolated pixels considered as noise:
                        //used with cancelImageNoise() function called before getting the entropies
                        //noiseLevel =0 -> no noise
                        //noiseLevel =1 -> groups of 1x1 pixel
                        //noiseLevel =2 -> groups of 3x3 pixels
                        //noiseLevel =3 -> groups of 5x5 pixels

        int     tilesOffset_x,//tiles are built from the given offset (default (0,0))
                tilesOffset_y; 
        bool    SIREM;  //the entropy of SIREM (true) or Shannon (false) is used
        int     nRows,  //rows in the image
                nCols;  //columns in the image
        bool    **pxMapImage_p; //image pixel map
    }IMAGE_ENTROPY;

    
    typedef struct{
        unsigned short code; //code associated with a tile
        int size;            //number of tiles with that code
    }SPLIT;         

    typedef struct{
        int x, y;       //coordinates of a pixel
    }PIXEL_2D;  
  
//Receives a binary image and returns its relative entropy[0..1] or SIREM, as indicated
//It does this by decomposing the image into tiles of a given size (1x1, 2x2, 3x3 or 4x4). Only tiles internal to the image are used
//A basic algorithm evaluates the number of times each tile is repeated to determine the probability of a certain event.
//An event is a particular distribution of pixels in the tile. Subsequently, the entropy equation according to Shannon is applied.
//In this algorithm, a couple of differentiated actions are proposed:
// 1) Tiles without pixels are now considered non-repeated tiles (exclusive or with a unique code in the image)
// 2) The tiles considered noisy, although they may be repeated, also add to the list of exclusive tiles
//With these actions you achieve:
// a) eliminate the non-monotonicity presented by Shannon entropy,
// b) the possibility of analyzing images with low pixel concentrations
// c) considerable insensitivity to noise in the image
//The entropy value is relativized to the maximum possible value.
    
//NOTE: before using the class, the setRecBaseBinaryImage(IMAGE_ENTROPY *image_p) function must be called with a particular image
//is the image that will serve as a basis to know the area of interest on which the tiles will be displayed
class RelEntropy
{
public:
    //Constructor
    //initializes memory pointers and certain variables
    RelEntropy();
    
    //destructor
    //image memory is freed
    ~RelEntropy();
    
    //receives a binary image corresponding to a section based on the primary image
    //sets the list of valid tiles (inside the image and filled with '1') and their size
    //initialize the class
    //Arguments:
    //  image_p => pointer to structure with image info
    //Return -1 if I fail    
    int initClass(IMAGE_ENTROPY *image_p);
    
    //Set the binary image for further analysis
    //return -1 if the image is invalid (does not exist)
    int setRecBinaryImage(IMAGE_ENTROPY *image_p);
    
    //Sets the tile size used to decompose the image
    void setTile(Tile tile);
    
    //Determines the relative entropy of the image. Relative to the maximum possible.
    //return the relative entropy in the range [0..1]:
    //differential treatment is done according to the pixels that house the tiles
    //This is to eliminate a lack of monotony at low concentrations of pixels in the image
    float getRelEntropy();

    //retorna cierto si el tile dado es interior a la imagen en su totalidad
    bool isTileIntoImage(int tile);
    
    //returns the number of tiles in the image with the given number of active pixels
    //arguments:
    //  pixels=number of active pixels in the tile
    //return
    //  -1 if the argument is invalid
    //the number of tiles requested   
    int tilesWithPixels(int pixels);
    
    //establishes the array m_repeatedSplit[] indicating the number of repeated tiles in the image according to the number of active pixels
    //return the size of the array;    
    int getEventsByTilePxs();
    
    //presents detailed information, segregated according to the number of active pixels of the tiles
    void viewInfo();
    
    //sets an offset on all image tiles
    void setTilesOffset(int x, int y);
    
    //returns the maximum number of active pixels in the image
    int getMaxImagePixels();
    
    //returns the number of active pixels in the image
    int getImageActivePixels();
    
    //Returns the number of pixels considered noisy in the image
    int getImageNoisyPixels();

    //returns the number of tiles completely included in the image
    int getImageTiles();
    
    //returns the active tile
    Tile getTile();
    
    //Returns the part of the relative entropy corresponding to the exclusive tiles (noisier empty ones)
    float getExclusiveRelEntropy();
    
    //Returns the part of the relative entropy corresponding to non-exclusive tiles
    float getNonExclusiveRelEntropy();
    
    //establishes the maximum grouping of isolated pixels considered as noise:
     //noiseLevel =0 -> no noise
     //noiseLevel =1 -> groups of 1 pixel
     //noiseLevel =2 -> groups of 3x3 pixels
     //noiseLevel =3 -> groups of 5x5 pixels
    void setNoiseLevel(int noiseLevel);

    //deactivates the active pixels in the tile considered noisy
    //return:
    // the number of tiles considered noisy
    int removeImageNoise();
    
    //variables
    bool    m_error;                //true if the constructor has failed
    int     m_totalTiles,           //total tiles of the rectangle that contains the image
            m_intoSampleTiles;      //interior tiles to the image. Those that are full of pxs in section at level 1
    Tile    m_tile;                 //tile used
    int     m_eventsByTilePxs[17],  //events based on the number of pixels in the tile
            m_exclusiveSplit[17];   //Exclusive tiles: no louder active pixels
    double  m_H[2];                 //entropy associated with each set from which information is segregated
    bool    m_SIREM;                //Structural Inverse Relativ Entropy Measure: true if Shannon's algorithm is not used
    
protected:
    
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
    //Arguments
    //tile -> tile of interest
    //type -> 0->axial horizontal; 1->vertical axial; 2->diagonal left-right; 3->diagonal right-left
    unsigned short getCodeTile(int tile, char type=0);
    
    //Retorna la cantidad de bits a '1' en el tile
    int getPixelsTile(int tile);
    
    //Return the number of bits to '1' in the event
    //the event must be of type unsigned short
    int getPixelsEvent(unsigned short event);
    
    //establishes a structure with information segregated from codes and sizes according to the number of pixels in the code
    //m_split_p is an array of pointers to arrays of SPLIT structures of various sizes
    //each SPLIT array maintains info associated with events with the same number of pixels.
    //each SPLIT contains the code of a tile and the times it has been repeated in the image.
    //The array[0] is different: it contains the indexes to the first free element of arrays 2 to 15 (repeated codes)
    void pxSplit();
    
    //Receives an array of m_split_p structures and establishes an ordering of its elements.
    //The result is in the m_repeatedGroups_p[] array with the following interpretation:
    //m_repeatedGroups_p[].code=number of repetitions of a code (repeated tiles)
    //m_repeatedGroups_p[].size=times a set of this same size is repeated (with other codes)
    //return the size of the array m_repeatedGroups_p[]
    int getRepeatedGroups();
    
    //sets the number of noisy tiles in the image ordered by number of pixels (tiles_p)
    //returns the number of pixels in the image considered noisy
    // Clusters of pixels are detected based on four different scans of the tiles: two axial and two diagonal
    //If those clusters do not exist, it is noise
    int getNoisyPxs(int *tiles_p);
    
    //Returns the number of pixels considered as noise from the code associated with a tile.
    //For the given code, count ones and zeros flanked by ones: ex: 001011101000-> 5 '1' and 2 '0'
    //Those tiles that have one or two px and those that have
    //less than half of px have any zeros in the sequence flanked by ones
    //int getTileNoisyPxs(unsigned short code);
    bool getTileNoisyPxs(unsigned short code, int *TotalPixels, int *isolatedPixels);
        
    //Groups of isolated pixels are cancelled.
    //Makes use of m_noiseLevel: indicates the size of the grid considered as the maximum group of isolated pixels
    // i.e. all pixels located in that grid with no pixels in the neighborhood are removed
    //Returns the number of pixels removed from the image
    //Return -1 if level<1
    int cancelImageNoise();
    
    //Set neighbor_p[] to list of active px on a 3x3 grid whose center is row/col
    //Each element of neighbor_p[] is the coordinates of the active pixel.
    //If the pixel is not activated, the coordinates are -1/-1
    //Arguments:
    // row -> row of the pixel image to analyze
    // col -> column of the pixel image to analyze
    // neighbor_p[] -> array of pixels in the 3x3 grid ordered from top to bottom and from left to right
    //return the number of active pixels in the grid, including px
    int getNeighbours(int row, int col, PIXEL_2D *neighbour_p);
    
    //Set in pxOut_p[] the coordinates of all the pxs included in a grid higher than the input grid
    //The coordinates of the input px include the coordinates of all the px located on its periphery.
    //Arguments:
    //pxIn_p[] -> array with the coordinates of the input pixels
    //pxOut_p[] -> array with the coordinates of the output pixels
    //nPxIn -> elements in pxIn_p[]
    //return the size of pxOut_p[]
    //pxOut_p[] must have space reserved to accommodate (level-1)*2+3)^2 elements
    int searchPxs(PIXEL_2D *pxIn_p, int nPxIn, PIXEL_2D *pxOut_p);
    
    //variables protegidas
    IMAGE_ENTROPY *m_image_p;
    int     m_tileSide,          //tile lateral length
            m_nXTiles,           //tiles in X
            m_nYTiles;           //tiles in Y            
    bool    *m_tilesList_p;      //houses the list of tiles inside the image
    int     m_nPxTile;           //number of pixels in the tile
    unsigned int *m_evCount_p,  //pointer to repeated event counter
            m_nEvents;          //total number of possible events
    SPLIT   *m_split_p[17];     //maintains lists of repeating tiles based on their pixel numbers
    int     m_xOffset,          //displacement in pixels of the start of the tessellation with respect to the reference frame
            m_yOffset;
    int     m_noiseLevel,       //indicates the maximum size of the clusters of isolated pixels that will be removed from the image
            m_imageNoisyPxs;    //total number of noisy pixels in the image
    LGROUP   m_repeated;        //count the repetitions of events
    int     m_totalExclusiveTiles; //total exclusive events (not repeated in the image)
};

#endif
