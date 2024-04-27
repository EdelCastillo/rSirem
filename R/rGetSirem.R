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

#' rGetSirem
#' Get information about sirem and its associated spikes.
#'
#' @param rMSIData   -> sample data obtained from the file with rMSI2::LoadMsiData().
#' @param params     -> parameters for sirem and for peaks.
#' @param initMass   -> mass, in Daltons, corresponding to the first image.
#' @param finalMass  -> mass, in Daltons, corresponding to the last image.
#'
#' @return a list: sirem->sirem information; scanMean->averaged values of each image; sirem->sirem value of each image.
#' @export
#' 
rGetSirem<-function(rMSIData, params, initMass, finalMass)
{
  if(initMass<rMSIData$mass[1]) 
  {cat("Warning: the initial mass must be greather than", rMSIData$mass[1]); return(-1);}  
  if(finalMass>rMSIData$mass[length(rMSIData$mass)]) 
  {cat("Warning: the final mass must be less than or equal to", rMSIData$mass[length(rMSIData$mass)]); return(-1);}  
  initMassIndex<-rSirem::rGetIndexFromMass(initMass, rMSIData$mass); #index to the initial mass.
  
  finalMassIndex<-rSirem::rGetIndexFromMass(finalMass, rMSIData$mass);  #index to the final mass.
  
  size=finalMassIndex-initMassIndex+1;
  if(size<=0) {print("Warning: the final mass must be greater than the initial mass."); return(-1);}  
  #the images are obtained (in columns).
  
  #se obtienen las imágenes (en columnas)
#  images<-rMSI2::load_imzMLImages(rMSIData, initMass, size);
  images<-rLoadImages(rMSIData, initMass, size);
#  images<-rImzML::getImages(rMSIData, initMass, size);
  
  #se obtiene la info de sirem + picos
  #se resta 1 dado que se requieren coordenadas iniciadas en 0,0
  #se pasa un vector vacío para usar los propios datos como imagen de referencia
  sirem<-rSirem::rSirem(images, rMSIData$pos-1, params, numeric());
  return(sirem);
}
