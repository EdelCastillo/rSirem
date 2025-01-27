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

#' rGetGaussians
#' Gets the Gaussians that decompose the given concentration data.
#' Makes use of the peak information returned by rGetSiremPeaks().
#'
#' @param rMSIData -> sample data obtained from the file with rMSI2::LoadMsiData().
#' @param siremInfo -> sirem data and peaks obtained from rGetSiremPeaks().
#' @param pixel -> pixel whose partial spectrum is deconvolved. if pixel<=0, averaged data is used
#' @param minMeanValue -> minimum averaged value of a peak to consider it valid.
#'
#' @return a list: 
#' gaussians -> matrix with the parameters of the Gaussians.
#' xAxis -> X axis (Daltons) (if error xAxis=0)
#' yAxis -> Y axis (concentration) (if error yAxis=0)
#' If pixel<=0 the averaged values are returned.
#' If pixel > 0, the part of the spectrum associated with the X axis is returned
#' @export
#' 
rGetGaussians<-function(rMSIData, siremInfo, pixel, minMeanValue)
{
  retFail<-list(gaussians=0, xAxis=0, yAxis=0);
  if(length(siremInfo$siremPeaks)==1 && siremInfo$siremPeaks[1]==-1)
     {print("Warning: there is no sirem data"); return(retFail);}
     
  initMass<-siremInfo$massAxis[1]; # mz inicial
  size=length(siremInfo$massAxis); # tamaño del eje X
  
  if(pixel>nrow(rMSIData$pos)) {print("Warning: pixel is out of range"); return(retFail);}
  
  if(pixel>0) # si interesa parte del espectro de un píxel específico
    {pixelMagnitudes<-rGetPixelChunkSpectra(rMSIData, pixel, initMass, size);}
  else # Se usan valores promediados
    {pixelMagnitudes<-vector();} #vector nulo
  
  #Descomposición en gausianas
  gaussiansInfo<-rGmmPeaks(siremInfo$siremPeaks, pixelMagnitudes, minMeanValue, siremInfo$massAxis);
  logic=gaussiansInfo[,1]!=0;
  gaussiansInfo=gaussiansInfo[logic,]; #removes possible null values from the end of the matrix.
  
#  logic=gaussiansInfo[,1]>0;
#  gaussiansInfo=gaussiansInfo[logic,];
  if(pixel>0) #se retorna las gausianas, el eje de mz y la parte del espectro asociado al eje mz 
    {gauss<-list(gaussians=gaussiansInfo, xAxis=siremInfo$massAxis, yAxis=pixelMagnitudes);}
  else #se retornan las gausianas, el eje de mz y los valores promediados
  {gauss<-list(gaussians=gaussiansInfo, xAxis=siremInfo$massAxis, yAxis=siremInfo$siremPeaks$scanMean);}
  
  return (gauss);
}
