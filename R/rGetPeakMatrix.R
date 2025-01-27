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


#' rGetPeaksMatrix
#' Build the peaks matrix.
#' Makes use of the gaussian information returned by rGetGaussians().
#'
#' @param rMSIData -> sample data obtained from the file with rMSI2::LoadMsiData().
#' @param gaussiansMatrix -> gaussians from rGetGaussians().
#'
#' @return a list: 
#' peakMatrix -> matrix[pixels, mass]
#' xAxis -> mz axis

#' @export
#' 
rGetPeaksMatrix<-function(rMSIData, gaussiansMatrix)
{
  retFail<-list(gaussians=0, xAxis=0, yAxis=0);
  nPeaks<-nrow(gaussiansMatrix); #gaussians number
  if(nPeaks==0)
    {print("Warning: there is no sirem data"); return(retFail);}
  
  nPixels<-nrow(rMSIData$pos); #pixels number
  peaksMatrix <-matrix(nrow = nPixels, ncol = nPeaks); #memory for peaks matrix
  
  for(index in 1:nPeaks)
  {
    #loading the magnitude info for the given mass over all pixels.
    imgMx<-rLoadImages(rMSIData, gaussiansMatrix[index, 1], 1);
    peaksMatrix[, index]<-imgMx[, 1];
  }
  ret<-list (peaksMatrix=peaksMatrix, xAxis=gaussiansMatrix[,1]);
  return (ret);  
}


