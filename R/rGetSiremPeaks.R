#########################################################################
#     rSirem - R package for MSI data deconvolution
#     Copyright (C) november 2023, Esteban del Castillo Pérez
#     esteban.delcastillo@urv.cat
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
############################################################################

#' rGetSiremPeaks:
#' get information about Sirem and its associated spikes.
#'
#' @param rMSIData   Sample data obtained from the file with rMSI2::LoadMsiData().
#' @param params     Parameters for sirem and for peaks.
#'    algorithm         -> 0 = sirem; 1=entropy.
#'    cutLeves          -> cut levels to apply to each image to generate binary images (percentiles).
#'    magSensitivity    -> variable sensitivity depending on the concentration level[0:maxConcentration] (bounded between limits).
#'    siremSensitivity  -> variable sensitivity depending on the sirem level[0:1] (bounded between limits).
#'    minMeanPxMag      -> minimum averaged concentration value to be considered.
#'    minSectionDensity -> minimum active pixel density value in the image to be considered.
#'    noiseLevel        -> absolute noise Level[0:maxConcentration]: lower values are considered null.
#'    pxCoord_X         -> X coordinates of each pixel.
#'    pxCoord_Y         -> Y coordinates of each pixel.
#'    referenceType     -> reference image type:
#'                      0=obtained from concentration info (data)
#'                      1=the image with the highest average value is adopted
#'                      2=all pixels have significant concentration.
#'    scanReference     -> reference image to determine the area to consider.
#'    tileSide          -> number of side pixels of the square tile: 1(1x1), 2(2x2), 3(3x3), 4(4x4)
#'    normalization     -> 0=none
#'                      -> 1=TIC
#'                      -> 2=RMS
#'                      -> 3=MAX
#'    relativeMag       -> !=0: magnitudes in %
#' @param initMass   Mass, in Daltons, corresponding to the first image to load.
#' @param finalMass  Mass, in Daltons, corresponding to the last  image to load.
#'
#' @return A list: 
#'    siremPeaks -> sirem and peak information; 
#'    massAxis   -> mass axis corresponding to the images.
#' @export
#' 
rGetSiremPeaks<-function(rMSIData, params, initMass, finalMass)
{
  normalizationType=params$normalization;
  
  if(normalizationType>3 || normalizationType<0)
  {cat("warning: normalization must be within range [0:3]\n"); return -1;}
  
  if(initMass<rMSIData$mass[1]) 
    {cat("Warning: the initial mass must be greather than", rMSIData$mass[1]); return(-1);}  
  if(finalMass>rMSIData$mass[length(rMSIData$mass)]) 
    {cat("Warning: the final mass must be less than or equal to", rMSIData$mass[length(rMSIData$mass)]); return(-1);}  
  initMassIndex <-rGetIndexFromMass(initMass, rMSIData$mass);   #index to the initial mass.
  finalMassIndex<-rGetIndexFromMass(finalMass, rMSIData$mass);  #index to the final mass.

  size=finalMassIndex-initMassIndex+1;
  if(size<=0) {print("Warning: the final mass must be greater than the initial mass."); return(-1);}  
  #the images are obtained (in columns).
  
  print("loading images..."); 
#  images<-rMSI2::load_imzMLImages(rMSIData, initMass, size);
#  images<-rImzML::getImages(rMSIData, initMass, size);
  
 images<-rLoadImages(rMSIData, rMSIData$mass[initMassIndex+1], size);

#...............................................
# loading from preloaded files. It is very fast
#...............................................
#    load("/home/esteban/MALDI/rSirem/images_C30_700_900.RData");
#    load("/home/esteban/MALDI/rSirem/images_C60_300_500.RData");
#    load("/home/esteban/MALDI/rSirem/images_C120_300_500.RData");
#            init <-rGetIndexFromMass(700, rMSIData$mass);   #index to the initial mass.
#            newInit=initMassIndex-init;
#            newEnd=size+newInit;
#            if(newInit==0)  {newInit=1;}
#            images<-images[, newInit:newEnd-1];
#            size=ncol(images);
# .....................................
  
  msg<-sprintf("OK -> pixels:%d scans:%d", nrow(images), ncol(images));
  print(msg);
  
  #images normalization as a percentage.
  maxValue=0;
  for(scan in 1:size)
    {
      if(normalizationType!=0) #if normalization should be used.
        images[, scan]= images[, scan]/as.matrix(rMSIData$normalizations[, normalizationType]); 
      if(mean(images[,scan])>maxValue) maxValue=mean(images[,scan]);
    }

  if(params$relativeMag!=0)
    {
    for(scan in 1:size)
      images[,scan]=images[,scan]*100.0/maxValue;
    }

  maxMean=0;
  maxMeanIndex=0;
  maxImage=0;
  #reference image:
  if(params$referenceType==1) #the image with the highest average value is adopted.
  {
    for(i in 1:length(rMSIData$mean))
      if(rMSIData$mean[i]>maxMean) {maxMean=rMSIData$mean[i]; maxMeanIndex=i;}
    
    if(maxMeanIndex==0){print("warning: no data");}
#    maxImage<-rMSI2::load_imzMLImages(rMSIData, rMSIData$mass[maxMeanIndex], 1); #load max image
    maxImage<-rLoadImages(rMSIData, rMSIData$mass[maxMeanIndex], 1); #load max image
#    maxImage<-rImzML::getImages(rMSIData, rMSIData$mass[maxMeanIndex], 1);
    
    #image normalization as a percentage.
    if(normalizationType!=0) #if normalization should be used.
        maxImage=maxImage/as.matrix(rMSIData$normalizations[, normalizationType]); 
    if(params$relativeMag!=0)
      maxImage=maxImage*100.0/mean(maxImage); 
    
    #sirem info + peaks are obtained.
    #1 is subtracted since coordinates starting at 0.0 are required.
    siremPeaksInfo<-rSiremPeaks(images, rMSIData$pos-1, params, maxImage);
  }
  #if a reference image is not provided, you must obtain:
  # params$referenceType=0 -> obtained from concentration info (data);
  # params$referenceType=2 -> all pixels have significant concentration.
  else #if a reference image is not provided.
  {
    #sirem info + peaks are obtained.
    #1 is subtracted since coordinates starting at 0.0 are required.
    #an empty vector is passed to use the data itself as a reference image.
    siremPeaksInfo<-rSiremPeaks(images, rMSIData$pos-1, params, numeric());
  }
  
  #initial index to the total mass axis.
  initMassIndex<-rGetIndexFromMass(initMass, rMSIData$mass);
  mzAxis=0;
  #Build the partial mass axis.
  mzAxis[1]=rMSIData$mass[initMassIndex];
  for(i in 2:size)
    mzAxis[i]=rMSIData$mass[initMassIndex+i];
  
  #Return list.
  siremInfo <- list(siremPeaks=siremPeaksInfo, massAxis=mzAxis);
  return (siremInfo);
}

#' rGetIndexFromMass
#' Returns the index to the full mass axis whose mass is closest to the given mass.
#' The algorithm of successive approximations is used.
#' @param rMSIData -> sample data obtained from the file with rMSI2::LoadMsiData().
#' @param mass -> mass in Daltons
#' 
#' @return The index to mass if OK; -1 if KO
#'  
#' 
rGetIndexFromMass<-function(mass, massAxis)
{
  nMassPoints=length(massAxis); #total length of the mass axis.
  if(nMassPoints==1){return (massAxis);}
  
  #index search to the matching mass or higher.
  initMassIndex=0;
  for(i in 1:nMassPoints) 
    if(massAxis[i]>=mass) {initMassIndex=i; break;}
  if(initMassIndex==0) {print("warning initMass out of range"); return -1;}
  
  #the index is updated to the nearest value.
  if(massAxis[initMassIndex]!=mass) 
  {
#    if(initMassIndex+1 <= nMassPoints) 
#    {delta=massAxis[initMassIndex+1]-massAxis[initMassIndex];} #mass increase.
#    else 
    {delta=massAxis[initMassIndex]-massAxis[initMassIndex-1];}
    
    if(massAxis[initMassIndex]-mass>delta/2.0) {initMassIndex=initMassIndex-1;}
  }
  return (initMassIndex);
}

#' rLoadImages
#' Returns an array with a set of images.
#'
#' @param rMSIData -> sample data obtained from the file with rMSI2::LoadMsiData().
#' @param initMass -> mass in Daltons corresponding to the first image.
#' @param size -> number of images.
#'
#' @return a matrix: column is the image corresponding to consecutive masses; rows contains the concentration data for each pixel.
#' 
rLoadImages<-function(rMSIData, initMass, size)
{
  initMassIndex<-rGetIndexFromMass(initMass, rMSIData$mass); #index to the initial mass.
  
  slices<-seq(initMassIndex, initMassIndex+size-1, by=1); #images.
  #  imgMx <- rMSI2::load_imzMLImage(rMSIData, initMass, size); #load the matrix from file.
  imgMx <- rMSI2::loadImageSliceFromCols(rMSIData, slices); #load the matrix from file.
  #  imgMx <- rImzML::getImages(rMSIData, )
  return (imgMx);
}
