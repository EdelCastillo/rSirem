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

# functions for graphical presentation.
#---------------------------------------------------------------------------

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


#' rPlotDeconv3
#' Presents an image with the Gaussians of two peaks,
#' along with the curve resulting from the sum of them (in red).
#' The blue circles indicate the original concentration data to be adjusted.
#' The drawInfo parameters come from rGetGaussians()
#' 
#' @param drawInfo1: gaussians info from low    resolution peak
#' @param drawInfo2: gaussians info from median resolution peak
#' @param drawInfo3: gaussians info from high   resolution peak
#' @param minMass: minimun mz
#' @param maxMass: maximun mz
#' @param rMSI2_peaks1: low    resolution peaks from rMSI peak matrix
#' @param rMSI2_peaks1: median resolution peaks from rMSI peak matrix
#' @param rMSI3_peaks1: high   resolution peaks from rMSI peak matrix
#' drawInfo$gaussians:
#' column 1 -> mean values of each Gaussian.
#' column 2 -> standard deviation.
#' column 3 -> factor associated with each Gaussian.
#' drawInfo$axis:
#' column 1 -> Y axis (concentration data to adjust).
#' column 2 -> X axis (in Daltons).
#'
#' @return -1 if KO; 0 if OK
#' @export
#' 
rPlotDeconv3<-function(drawInfo1, drawInfo2, drawInfo3, minMass=0, maxMass=0, rMSI2_peaks1=-1, rMSI2_peaks2=-1, rMSI2_peaks3=-1)
{
  #Los tres ejes X se unifican 
  minMass1=min(drawInfo1$xAxis);
  minMass2=min(drawInfo2$xAxis);
  minMass3=min(drawInfo3$xAxis);
  maxMass1=max(drawInfo1$xAxis);
  maxMass2=max(drawInfo2$xAxis);
  maxMass3=max(drawInfo3$xAxis);
  
  if(minMass==0) {minMass=min(c(minMass1, minMass2, minMass3));} #eje completo
  else 
  {
    minMass=minMass; 
    if(minMass<minMass1 & minMass<minMass2 & minMass<minMass3)
    {print("minMass is out of range"); return(-1);}
  }
  if(maxMass==0) {maxMass=max(c(maxMass1, maxMass2, maxMass3));}
  else 
  {
    maxMass=maxMass; 
    if(maxMass>maxMass1 & maxMass>maxMass2 & maxMass>maxMass3)
    {print("maxMass is out of range"); return(-1);}
  }
  
  if(maxMass<=minMass) {print("mass range is wrong"); return(-1);}
  
  #mz extremas de cada eje
  x1Range=c(minMass1, maxMass1);
  x2Range=c(minMass2, maxMass2);
  x3Range=c(minMass3, maxMass3);
  
  #se extrae un subconjunto de datos delimitado por minMz y maxMz
  #para arg1
  gauss1_logic=drawInfo1$gaussians[,1]>=minMass & drawInfo1$gaussians[,1]<=maxMass;
  mean1_sub=(drawInfo1$gaussians[,1])[gauss1_logic];
  sigma1_sub=(drawInfo1$gaussians[,2])[gauss1_logic];
  value1_sub=(drawInfo1$gaussians[,3])[gauss1_logic];
  axis1_logic=drawInfo1$xAxis>=minMass & drawInfo1$xAxis<=maxMass;
  X1_axis=drawInfo1$xAxis[axis1_logic];
  Y1_axis=drawInfo1$yAxis[axis1_logic];
  if(length(mean1_sub)==0) {print("warning: no peaks in arg1 for this mass range"); return(-1);}
  maxRelVal1=100*max(value1_sub)/max(drawInfo1$gaussians[,3]);
  
  #para arg2
  gauss2_logic=drawInfo2$gaussians[,1]>=minMass & drawInfo2$gaussians[,1]<=maxMass;
  mean2_sub=(drawInfo2$gaussians[,1])[gauss2_logic];
  sigma2_sub=(drawInfo2$gaussians[,2])[gauss2_logic];
  value2_sub=(drawInfo2$gaussians[,3])[gauss2_logic];
  axis2_logic=drawInfo2$xAxis>=minMass & drawInfo2$xAxis<=maxMass;
  X2_axis=drawInfo2$xAxis[axis2_logic];
  Y2_axis=drawInfo2$yAxis[axis2_logic];
  if(length(mean2_sub)==0) {print("warning: no peaks in arg2 for this mass range"); return(-1);}
  maxRelVal2=100*max(value2_sub)/max(drawInfo2$gaussians[,3]);
  
  #para arg3
  gauss3_logic=drawInfo3$gaussians[,1]>=minMass & drawInfo3$gaussians[,1]<=maxMass;
  mean3_sub=(drawInfo3$gaussians[,1])[gauss3_logic];
  sigma3_sub=(drawInfo3$gaussians[,2])[gauss3_logic];
  value3_sub=(drawInfo3$gaussians[,3])[gauss3_logic];
  axis3_logic=drawInfo3$xAxis>=minMass & drawInfo3$xAxis<=maxMass;
  X3_axis=drawInfo3$xAxis[axis3_logic];
  Y3_axis=drawInfo3$yAxis[axis3_logic];
  if(length(mean3_sub)==0) {print("warning: no peaks in arg3 for this mass range"); return(-1);}
  maxRelVal3=100*max(value3_sub)/max(drawInfo3$gaussians[,3]);
  
  #mensaje informativo
  msg1<-sprintf("arg1 max relative value in range:%.1f %%", maxRelVal1)  
  msg2<-sprintf("arg2 max relative value in range:%.1f %%", maxRelVal2) 
  msg3<-sprintf("arg3 max relative value in range:%.1f %%", maxRelVal3) 
  print(msg1);
  print(msg2);
  print(msg3);
  
  #normalización de los datos en rango de 0..100
  value1_sub=value1_sub*100/max(Y1_axis);
  value2_sub=value2_sub*100/max(Y2_axis);
  value3_sub=value3_sub*100/max(Y3_axis);
  Y1Factor=100/max(Y1_axis);
  Y2Factor=100/max(Y2_axis);
  Y3Factor=100/max(Y3_axis);
  Y1_axis=Y1_axis*Y1Factor; #100/max(Y1_axis);
  Y2_axis=Y2_axis*Y2Factor; #100/max(Y2_axis);
  Y3_axis=Y3_axis*Y3Factor; #100/max(Y2_axis);
  
  X=c(minMass, maxMass);  #extremos del eje X
  Y=c(0, 325);      #extremos del eje Y
  
  #marco de imagen: sólo con eje X
  plot(X, Y, type="p", col="white", cex=0.01, axes = FALSE, main="Deconvolution", xlab="mz(Da)", ylab="relative average concentration", las=1, col.axis="black");
  axis(1); ##visualiza el eje X
  legendY1=sprintf("top factor=%.4f", Y1Factor);
  legendY2=sprintf("central factor=%.4f", Y2Factor);
  legendY3=sprintf("bottom factor=%.4f", Y3Factor);
  #  print(legendY1, legendY2);
  #  legend("topleft", legend=c(legendY1, legendY2));
  mtext(legendY1, side=3, adj=1);
  mtext(legendY2, side=3);
  mtext(legendY3, side=3, adj=0);
  
  #resolución en X ajustada a 1/5 de la sigma inferior
  minSigma1<-min(sigma1_sub); #required to determine the increments in X.
  minSigma2<-min(sigma2_sub); 
  minSigma3<-min(sigma3_sub); 
  minSigma=min(c(minSigma1, minSigma2, minSigma3));
  deltaX<-minSigma/5; #5 points fit within the lowest sigma (resolution in X).
  if(deltaX<1e-6) deltaX=1e-6;
  
  #para cada muestra
  for(res in 1:3)
  {
    if(res==1) #low resolution
    {
      mean_t=mean1_sub; #info de gausianas
      sigma_t=sigma1_sub;
      value_t=value1_sub;
      nGauss=length(mean_t); #número de gausianas
      Y=Y1_axis;        #eje Y
      X=X1_axis;        #eje X
    }
    if(res==2) #median resolution
    {
      mean_t=mean2_sub; #info de gausianas
      sigma_t=sigma2_sub;
      value_t=value2_sub;
      nGauss=length(mean_t); #número de gausianas
      Y=Y2_axis;        #eje Y
      X=X2_axis;
    }
    if(res==3) #high resolution
    {
      mean_t=mean3_sub; #info de gausianas
      sigma_t=sigma3_sub;
      value_t=value3_sub;
      nGauss=length(mean_t); #número de gausianas
      Y=Y3_axis;        #eje Y
      X=X3_axis;
    }
    #Blue circles are drawn corresponding to the original data.
    points(X, Y+115*(res-1), pch=21, col="blue",lwd=0.75, cex=0.75);  
    
    X<-seq(minMass, maxMass, deltaX); # X axis común
    totalY=0; #used to represent the sum curve.
    
    #para cada gausiana de cada pico
    for(i in 1 : nGauss) 
    {
      mean =mean_t[i]; #gaussians parameters
      sigma=sigma_t[i];
      value=value_t[i];
      
      Y=value*exp(-((X-mean)^2)/(2*(sigma^2))); #magnitude of the current Gaussian.
      totalY=totalY+Y; #the sum is updated.
      color=((i-1)%%8)+1;
      lines(X, Y+115*(res-1), type="l", col=color, lwd=1); #the current Gaussian is drawn.
    }
    lines(X, totalY+115*(res-1), type="l", col="red", lwd=2); #the sum curve is drawn.
  }
  
  #Lineas discontinuas verticales: mean of gaussians from mean spectrum 3
  for(i in 1:length(mean3_sub))
  {
    yHigh=value3_sub[i]+230
    Y<-1:yHigh;
    X<-rep(mean3_sub[i], times=yHigh)
    lines(X, Y, type="l", col="green", lwd=0.75, lty=2); 
  }
  #centroides from rMSI2
  if(length(rMSI2_peaks1)>1)
  {
    logicMass<-rMSI2_peaks1>=minMass & rMSI2_peaks1<=maxMass;
    rMSI2_mass<-rMSI2_peaks1[logicMass];
    lineLength=5;
    Y<-1:lineLength;
    for(i in 1:length(rMSI2_mass))
    {
      X<-rep(rMSI2_mass[i], times=lineLength)
      lines(X, Y, type="l", col="blue", lwd=3)
    }
  }
  if(length(rMSI2_peaks2)>1)
  {
    logicMass<-rMSI2_peaks2>=minMass & rMSI2_peaks2<=maxMass;
    rMSI2_mass<-rMSI2_peaks2[logicMass];
    for(i in 1:length(rMSI2_mass))
    {
      X<-rep(rMSI2_mass[i], times=lineLength)
      lines(X, Y+115, type="l", col="blue", lwd=3)
    }
  }
  if(length(rMSI2_peaks3)>1)
  {
    logicMass<-rMSI2_peaks3>=minMass & rMSI2_peaks3<=maxMass;
    rMSI2_mass<-rMSI2_peaks3[logicMass];
    for(i in 1:length(rMSI2_mass))
    {
      X<-rep(rMSI2_mass[i], times=lineLength)
      lines(X, Y+230, type="l", col="blue", lwd=3)
    }
  }
  return(0);
}


#' rPlotDeconv2
#' Presents an image with the Gaussians of two peaks,
#' along with the curve resulting from the sum of them (in red).
#' The blue circles indicate the original concentration data to be adjusted.
#' The drawInfo parameters come from rGetGaussians()
#' 
#' @param drawInfo1: gaussians info from low  resolution peak
#' @param drawInfo2: gaussians info from high resolution peak
#' @param minMass: minimun mz
#' @param maxMass: maximun mz
#' @param rMSI2_peaks1: low  resolution peaks from rMSI peak matrix
#' @param rMSI2_peaks1: high resolution peaks from rMSI peak matrix
#' drawInfo$gaussians:
#' column 1 -> mean values of each Gaussian.
#' column 2 -> standard deviation.
#' column 3 -> factor associated with each Gaussian.
#' drawInfo$axis:
#' column 1 -> Y axis (concentration data to adjust).
#' column 2 -> X axis (in Daltons).
#'
#' @return -1 if KO; 0 if OK
#' @export
#' 
rPlotDeconv2<-function(drawInfo1, drawInfo2, minMass=0, maxMass=0, rMSI2_peaks1=-1, rMSI2_peaks2=-1)
{
  #Los tres ejes X se unifican 
  minMass1=min(drawInfo1$xAxis);
  minMass2=min(drawInfo2$xAxis);
  maxMass1=max(drawInfo1$xAxis);
  maxMass2=max(drawInfo2$xAxis);
  
  if(minMass==0) {minMass=min(c(minMass1, minMass2));} #eje completo
  else 
    {
    minMass=minMass; 
    if(minMass<minMass1 & minMass<minMass2)
      {minMass=min(c(minMass1, minMass2));
      print("warning: the low mass was update");}
    }
  if(maxMass==0) {maxMass=max(c(maxMass1, maxMass2));}
  else 
  {
    maxMass=maxMass; 
    if(maxMass>maxMass1 & maxMass>maxMass2)
      {maxMass=max(c(maxMass1, maxMass2));
    print("warning: the high mass was update");}
  }
  
  if(maxMass<=minMass) {print("mass range is wrong"); return(-1);}

  #mz extremas de cada eje
  x1Range=c(minMass1, maxMass1);
  x2Range=c(minMass2, maxMass2);
  
  #se extrae un subconjunto de datos delimitado por minMz y maxMz
  #para arg1
  gauss1_logic=drawInfo1$gaussians[,1]>=minMass & drawInfo1$gaussians[,1]<=maxMass;
   mean1_sub=(drawInfo1$gaussians[,1])[gauss1_logic];
  sigma1_sub=(drawInfo1$gaussians[,2])[gauss1_logic];
  value1_sub=(drawInfo1$gaussians[,3])[gauss1_logic];
  axis1_logic=drawInfo1$xAxis>=minMass & drawInfo1$xAxis<=maxMass;
  X1_axis=drawInfo1$xAxis[axis1_logic];
  Y1_axis=drawInfo1$yAxis[axis1_logic];
  if(length(mean1_sub)==0) {print("warning: no peaks in arg1 for this mass range"); return(-1);}
  maxRelVal1=100*max(value1_sub)/max(drawInfo1$gaussians[,3]);
  
  #para arg2
  gauss2_logic=drawInfo2$gaussians[,1]>=minMass & drawInfo2$gaussians[,1]<=maxMass;
   mean2_sub=(drawInfo2$gaussians[,1])[gauss2_logic];
  sigma2_sub=(drawInfo2$gaussians[,2])[gauss2_logic];
  value2_sub=(drawInfo2$gaussians[,3])[gauss2_logic];
  axis2_logic=drawInfo2$xAxis>=minMass & drawInfo2$xAxis<=maxMass;
  X2_axis=drawInfo2$xAxis[axis2_logic];
  Y2_axis=drawInfo2$yAxis[axis2_logic];
  if(length(mean2_sub)==0) {print("warning: no peaks in arg2 for this mass range"); return(-1);}
  maxRelVal2=100*max(value2_sub)/max(drawInfo2$gaussians[,3]);

  #mensaje informativo
  msg1<-sprintf("arg1 max relative value in range:%.1f %%", maxRelVal1)  
  msg2<-sprintf("arg2 max relative value in range:%.1f %%", maxRelVal2) 
  print(msg1);
  print(msg2);
  
  #normalización de los datos en rango de 0..100
  value1_sub=value1_sub*100/max(Y1_axis);
  value2_sub=value2_sub*100/max(Y2_axis);
  Y1Factor=100/max(Y1_axis);
  Y2Factor=100/max(Y2_axis);
  Y1_axis=Y1_axis*Y1Factor; #100/max(Y1_axis);
  Y2_axis=Y2_axis*Y2Factor; #100/max(Y2_axis);
  
  X=c(minMass, maxMass);  #extremos del eje X
  Y=c(0, 225);      #extremos del eje Y
  
  #marco de imagen: sólo con eje X
  plot(X, Y, type="p", col="white", cex=0.01, axes = FALSE, main="Deconvolution", xlab="mz(Da)", ylab="relative average concentration", las=1, col.axis="black");
  axis(1); ##visualiza el eje X
  legendY1=sprintf("top factor=%.4f", Y1Factor);
  legendY2=sprintf("bottom factor=%.4f", Y2Factor);
#  print(legendY1, legendY2);
#  legend("topleft", legend=c(legendY1, legendY2));
  mtext(legendY1, side=3, adj=1);
  mtext(legendY2, side=3, adj=0);
  
  #resolución en X ajustada a 1/5 de la sigma inferior
  minSigma1<-min(sigma1_sub); #required to determine the increments in X.
  minSigma2<-min(sigma2_sub); 
  minSigma=min(c(minSigma1, minSigma2));
  deltaX<-minSigma/5; #5 points fit within the lowest sigma (resolution in X).
  if(deltaX<1e-6) deltaX=1e-6;
  
  #para cada muestra
  for(res in 1:2)
  {
    if(res==1) #low resolution
    {
      mean_t=mean1_sub; #info de gausianas
      sigma_t=sigma1_sub;
      value_t=value1_sub;
      nGauss=length(mean_t); #número de gausianas
      Y=Y1_axis;        #eje Y
      X=X1_axis;        #eje X
    }
    if(res==2) #median resolution
    {
      mean_t=mean2_sub; #info de gausianas
      sigma_t=sigma2_sub;
      value_t=value2_sub;
      nGauss=length(mean_t); #número de gausianas
      Y=Y2_axis;        #eje Y
      X=X2_axis;
    }
    #Blue circles are drawn corresponding to the original data.
    points(X, Y+115*(res-1), pch=21, col="blue",lwd=0.75, cex=0.75);  
    
    X<-seq(minMass, maxMass, deltaX); # X axis común
    totalY=0; #used to represent the sum curve.
    
    #para cada gausiana de cada pico
    for(i in 1 : nGauss) 
    {
      mean =mean_t[i]; #gaussians parameters
      sigma=sigma_t[i];
      value=value_t[i];
      
      Y=value*exp(-((X-mean)^2)/(2*(sigma^2))); #magnitude of the current Gaussian.
      totalY=totalY+Y; #the sum is updated.
      color=((i-1)%%8)+1;
      lines(X, Y+115*(res-1), type="l", col=color, lwd=1); #the current Gaussian is drawn.
    }
    lines(X, totalY+115*(res-1), type="l", col="red", lwd=2); #the sum curve is drawn.
  }
  
  #Lineas discontinuas verticales: mean of gaussians from mean spectrum 3
  for(i in 1:length(mean2_sub))
  {
    yHigh=value2_sub[i]+115
    Y<-1:yHigh;
    X<-rep(mean2_sub[i], times=yHigh)
    lines(X, Y, type="l", col="green", lwd=0.75, lty=2); 
  }
  #centroides from rMSI2
  if(length(rMSI2_peaks1)>1)
  {
    logicMass<-rMSI2_peaks1>=minMass & rMSI2_peaks1<=maxMass;
    rMSI2_mass<-rMSI2_peaks1[logicMass];
    lineLength=5;
    Y<-1:lineLength;
    for(i in 1:length(rMSI2_mass))
    {
      X<-rep(rMSI2_mass[i], times=lineLength)
      lines(X, Y, type="l", col="blue", lwd=3)
    }
  }
  if(length(rMSI2_peaks2)>1)
  {
    logicMass<-rMSI2_peaks2>=minMass & rMSI2_peaks2<=maxMass;
    rMSI2_mass<-rMSI2_peaks2[logicMass];
    for(i in 1:length(rMSI2_mass))
    {
      X<-rep(rMSI2_mass[i], times=lineLength)
      lines(X, Y+115, type="l", col="blue", lwd=3)
    }
  }
  return(0);
}


#' rPlotPeaks
#' Presents an image with the centroids of two data set into a mass range,
#' The parameters come from rGetGaussians() and rMSI2::LoadPeakMatrix()
#' 
#' @param gaussians -> gaussians info matrix from return of rGetGaussians()
#' @param MSI2peak  -> list from rMSI2::LoadPeakMatrix()
#' @param minMass -> low  mass to analyze (Da)
#' @param minMass -> high mass to analyze (Da)
#' siremPeak$gaussians:
#' column 1 -> mean values of each Gaussian.
#' column 2 -> standard deviation.
#' column 3 -> factor associated with each Gaussian.
#' MSI2peak$mass: mass vector
#' MSI2peak$intensity:
#' columns -> concentration for each mass.
#' rows    -> concentration for each pixel.
#'
#' @return -1 if KO; 0 if OK
#' @export
#' 
rPlotPeaks<-function(gaussians, MSI2peak, minMass, maxMass) 
{
  logicSiremMass<-gaussians[,1]>=minMass & gaussians[,1]<=maxMass;
  logicMSImass  <-MSI2peak$mass>=minMass & MSI2peak$mass<=maxMass;
  siremMass<-(gaussians[,1])[logicSiremMass];
  MSI2mass  <-MSI2peak$mass[logicMSImass];

  #Los ejes X se unifican 
  x1Size<-length(siremMass);
  x2Size<-length(MSI2mass);

  if(x1Size<2 | x2Size<2) {print("no data to draw"); return(-1);}
  
  mzLowArray =c(siremMass[1], MSI2mass[1]);
  mzHighArray=c(siremMass[x1Size], MSI2mass[x2Size]);
  minX=max(mzLowArray); #valor mínimo del eje unificado
  maxX=min(mzHighArray);#valor máximo del eje unificado
  
  logicSiremMass<-gaussians[,1]>=minX & gaussians[,1]<=maxX;
  logicMSImass  <-MSI2peak$mass>=minX & MSI2peak$mass<=maxX;
  siremMass<-(gaussians[,1])[logicSiremMass];
  siremMag <-(gaussians[,3])[logicSiremMass];
  MSI2mass <-MSI2peak$mass[logicMSImass];
  
  meanMag=rep(0, times=length(MSI2peak$mass));
  
  for(i in 1:length(MSI2peak$mass))
      meanMag[i]<-mean(MSI2peak$intensity[,i]);
      
  MSI2mag<-meanMag[logicMSImass];
  
#  MSI2mag  <-(MSI2peak$intensity[1,])[logicMSImass];
  
  maxMSI2mag=max(MSI2mag);
  MSI2mag<-100.0*MSI2mag/maxMSI2mag;
  
  X=c(minX, maxX);  #extremos del eje X
  Y=c(0, 220);      #extremos del eje Y
  
  #marco de imagen: sólo con eje X
  plot(X, Y, type="p", col="white", cex=0.01, axes = FALSE, main="Centroids", xlab="mz(Da)", ylab="relative average concentration", las=1, col.axis="black");
  axis(1); ##visualiza el eje X
  
  #para cada gausiana de cada pico
  for(i in 1 : length(siremMass)) 
    {
      X=c(siremMass[i], siremMass[i]); #gaussians mean
      Y=c(0, siremMag[i])
      lines(X, Y, type="l", col="blue", lwd=1); #the current Gaussian is drawn.
    }
  for(i in 1 : length(MSI2mass))
  {
    X=c(MSI2mass[i], MSI2mass[i]); #gaussians mean
    Y=c(0, MSI2mag[i])
    lines(X, Y+105, type="l", col="red", lwd=1); #the current Gaussian is drawn.
  }
  

  return(0);
}



#' rPlotSirem
#' Presents an image with peak concentration and sirem information.
#' In red the magnitude curve and in green the sirem curve.
#' Blue circles indicate data.
#'
#' @param peakInfo: list returned from rGetSiremPeaks().
#' @param minMz: minimun mz
#' @param maxMz: maximun mz
#'
#' @export
#' 
rPlotSirem<-function(peaksInfo, minMz, maxMz)
  {
  xAxisLogic=peaksInfo$massAxis>=minMz & peaksInfo$massAxis<=maxMz;
  xAxis=peaksInfo$massAxis[xAxisLogic];
  yAxis=peaksInfo$siremPeaks$scanMean[xAxisLogic];
  sirem=peaksInfo$siremPeaks$sirem[xAxisLogic];
  yFactor=100/max(yAxis);
  yAxis=yAxis*yFactor;
  
  if(length(xAxis)<2) {print("no data to draw"); return(-1);}
  
  nPoints<-length(xAxis); #number of original data on each axis.
  minX<-xAxis[1];         #minimum value on the X axis.
  maxX<-xAxis[nPoints];   #maximum value on the X axis.
  if(nPoints<2) {print("no data to draw."); return(-1);}
  
  #Blue circles are drawn corresponding to the original data.
  plot(xAxis, yAxis, type="p",col="blue",lwd=1, main="Deconvolution", xlab="mz(Daltons)", 
       ylab="Concentration/sirem", las=1, col.axis="black");  
  totalY=0; #used to represent the sum curve.

  maxSirem<-max(sirem);
  sFactor=40/maxSirem;# maximum height of 40%
  
  legendY1=sprintf("concentration factor(red)=%.4f", yFactor);
  legendY2=sprintf("sirem factor(green)=%.4f", sFactor);
  #  print(legendY1, legendY2);
  #  legend("topleft", legend=c(legendY1, legendY2));
  mtext(legendY1, side=3, adj=0);
  mtext(legendY2, side=3, adj=1);
  
  lines(xAxis, yAxis, type="l",col="red",lwd=2); #concentration
  
  lines(xAxis, sirem*sFactor, type="p",col="blue",lwd=1); #sirem
  lines(xAxis, sirem*sFactor, type="l",col="green",lwd=2);
}
