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

#Returns a piece of the spectrum of a pixel.
#Arguments:
# rMSIData -> sample data obtained from the file with rMSI2::LoadMsiData().
# pixel -> pixel from which the partial spectrum is extracted.
# initMass -> mass in Daltons corresponding to the first image.
# size -> number of images.
#Return:
# A vector with the requested spectrum part.
rGetPixelChunkSpectra<-function(rMSIData, pixel, initMass, size)
{
  initMassIndex<-rGetIndexFromMass(rMSIData, initMass); #index to the initial mass.
  
  mySpectra <- rMSI2::loadImgChunkFromIds(rMSIData, pixel); #loads the full spectrum from file.

  chunkSpectra=0;
  for(i in 1:size)
    chunkSpectra[i]=mySpectra[1, initMassIndex+i];

  #partial spectrum.
  return (chunkSpectra);  
}