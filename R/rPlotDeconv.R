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

#' rPlotDeconv
#' Presents an image with the Gaussians passed as arguments,
#' along with the curve resulting from the sum of all of them (in red).
#' The blue circles indicate the original concentration data to be adjusted.
#'
#' @param drawInfo$gaussians:
#' column 1 -> mean values of each Gaussian.
#' column 2 -> standard deviation.
#' column 3 -> factor associated with each Gaussian.
#' @param drawInfo$axis:
#' column 1 -> Y axis (concentration data to adjust).
#' column 2 -> X axis (in Daltons).
#'
#' @return -1 if KO; 0 if OK
#' @export
#' 
rPlotDeconv<-function(drawInfo)
{
  if(length(drawInfo$xAxis)<2) {print("no data to draw"); return(-1);}
  
  nPoints<-length(drawInfo$xAxis); #number of original data on each axis.
  minX<-drawInfo$xAxis[1];         #minimum value on the X axis.
  maxX<-drawInfo$xAxis[nPoints];   #maximum value on the X axis.
  if(nPoints<2) {print("no data to draw."); return(-1);}
  
  nGauss<-nrow(drawInfo$gaussians); #number of Gaussians to represent.
  maxGaussValue=0;
  maxGaussIndex=0;
  minSigma<-min(drawInfo$gaussians[,2]); #required to determine the increments in X.
  
  Y=drawInfo$yAxis;
  X=drawInfo$xAxis;
  
  #Blue circles are drawn corresponding to the original data.
  plot(X, Y, type="p",col="blue",lwd=1, main="Deconvolution", xlab="mz(Daltons)", 
       ylab="Concentration", las=1, col.axis="black");  
  totalY=0; #used to represent the sum curve.
  
  deltaX<-minSigma/5; #5 points fit within the lowest sigma (resolution in X).
  if(deltaX<1e-6) deltaX=1e-6;
  
  X<-seq(minX, maxX, deltaX); # X axis.
  for(i in 1 : nGauss)
    {
    mean =drawInfo$gaussians[i, 1]; #gaussians parameters
    sigma=drawInfo$gaussians[i, 2];
    value=drawInfo$gaussians[i, 3];
      
    Y=value*exp(-((X-mean)^2)/(2*(sigma^2))); #magnitude of the current Gaussian.
    totalY=totalY+Y; #the sum is updated.
    color=((i-1)%%8)+1;
    lines(X, Y, type="l", col=color, lwd=1); #the current Gaussian is drawn.
    }
  lines(X, totalY, type="l", col="red", lwd=2); #the sum curve is drawn.
  return(0);
}
