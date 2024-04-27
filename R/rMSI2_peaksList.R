
#genera dos matrices con info extraida del fihero pasado
#magMatrix : cada fila mantiene las  magnitudes de los picos en un pixel
#massMatrix: cada fila mantiene las  masas      de los picos en un pixel
#' @export
rGetPeaksList<-function(fileName, params, pxLow, pxHigh, mzLow, mzHigh)
{
  #carga de valores
  X<-pxLow : pxHigh;
  minMass=1e32;
  maxMass=0;
  minInt=1e32;
  maxInt=0;
  
  imzML_fname <- path.expand(fileName)
  dirName<-dirname(imzML_fname) #hasta último '/'
  baseName <- basename(imzML_fname) #file
  peaksFile<-paste0(dirName, "/PeakPicking/", baseName);
  
  nC<-nchar(peaksFile);
  nC<-nC-6; #restamos .imzML

  peaksFile<-sprintf("%s-peaks.imzML",substr(peaksFile, 1, nC));
  
  myDataInt<-rMSI2::LoadMsiData(imzML_fname);
  #gaussians from average spectrum
  siremPeaks<-rGetSiremPeaks(myDataInt, params, mzLow, mzHigh);
  gaussInfoMean<-rGetGaussians(myDataInt, siremPeaks, 0, params$minMeanPxMag); 
  
  if(0)
  for(px in pxLow : pxHigh)
  {
    pxPeaksList<-rMSI2::readimzML_singlePixelPeakList(peaksFile, px);
    logicIndex<-pxPeaksList$mass>=mzLow & pxPeaksList$mass<=mzHigh;
    mass=pxPeaksList$mass[logicIndex];
    nMass<-length(mass);
    intensity=pxPeaksList$intensity[logicIndex];
    if(nMass>0)
      for(i in 1 : nMass)
      {
        if(mass[i]<minMass) minMass=mass[i];
        if(mass[i]>maxMass) maxMass=mass[i];
        if(intensity[i]<minInt) minInt=intensity[i];
        if(intensity[i]>maxInt) maxInt=intensity[i];
      }
  }
  #yRange=c(minMass, maxMass);
  #xRange=c(minInt, maxInt);
  #yRange=c(mzLow, mzHigh);
  nGauss=length(gaussInfoMean$gaussians[,1])
  minMz=gaussInfoMean$gaussians[1,1]-0.02;
  maxMz=gaussInfoMean$gaussians[nGauss,1]+0.02;
  yRange=c(minMz, maxMz);
  xRange=c(pxLow, pxHigh);
  plot(xRange, yRange, type="p", col="white", cex=0.01, main="Deconvolution", xlab="pixels", ylab="mz(Da)", las=1, col.axis="black");
#  plot(X, Y, type="p",col="blue",lwd=1, main="Deconvolution", xlab="mz(Daltons)", 
#       ylab="Concentration", las=1, col.axis="black");  
  
  for(px in pxLow : pxHigh)
  {
    pxPeaksList<-rMSI2::readimzML_singlePixelPeakList(peaksFile, px);
    logicIndex<-pxPeaksList$mass>=mzLow & pxPeaksList$mass<=mzHigh;
    nMass<-length(pxPeaksList$mass[logicIndex]);
      
#      X<-pxPeaksList$intensity[logicIndex];
      Y<-pxPeaksList$mass     [logicIndex];
      X<-rep(px, times=nMass);
      #X<-pxPeaksList$intensity[logicIndex];
      if(nMass>0)
        lines(X, Y, type="p", col="blue", cex=0.2);  
      
      
  }
  
  xyLow<-myDataInt$pos[pxLow,];
  xyHigh<-myDataInt$pos[pxHigh,];
  sprintf("low(X/Y):(%d/%d) high(x/y):(%d/%d)", xyLow[1], xyLow[2], xyHigh[1], xyHigh[2]);
  
  
  for(px in pxLow : pxHigh)
  {
    gaussInfo<-rGetGaussians(myDataInt, siremPeaks, px, params$minMeanPxMag);
    nGauss<-nrow(gaussInfo$gaussians); #number of Gaussians to represent.
    Y<-gaussInfo$gaussians[,1];
    X<-rep(px, times=nGauss);
    if(nGauss>0)
      lines(X, Y, type="p", col="red", cex=0.2);  
    
  }
  
  #Lineas discontinuas: mean of gaussians from mean spectrum
  pxList<-pxLow:pxHigh
  for(i in 1:length(gaussInfoMean$gaussians[,1]))
  {
    Y<-rep(gaussInfoMean$gaussians[i, 1], times=length(pxList))
    lines(pxList, Y, type="l", col="green", lwd=2, lty=2); 
  }
  
  return(0);
}

#genera dos matrices con info extraida del fihero pasado
#magMatrix : cada fila mantiene las  magnitudes de los picos en un pixel
#massMatrix: cada fila mantiene las  masas      de los picos en un pixel
#' @export
rGetPeaksListFromROI<-function(fileName, params, roi, mzLow, mzHigh)
{
  #carga de valores
#  X<-pxLow : pxHigh;
  minMass=1e32;
  maxMass=0;
  minInt=1e32;
  maxInt=0;
  
  imzML_fname <- path.expand(fileName)
  dirName<-dirname(imzML_fname) #hasta último '/'
  baseName <- basename(imzML_fname) #file
  peaksFile<-paste0(dirName, "/PeakPicking/", baseName);
  
  nC<-nchar(peaksFile);
  nC<-nC-6; #restamos .imzML
  
  peaksFile<-sprintf("%s-peaks.imzML",substr(peaksFile, 1, nC));
  
  myDataInt<-rMSI2::LoadMsiData(imzML_fname);
  

  #se genera la lista de pixeles del roi
  xSize=roi[2,1]-roi[1,1];
  ySize=roi[2,2]-roi[1,2];
  pxList<-array(1, xSize*ySize);
  z<-1;
  for(i in 1:length(myDataInt$pos[,1]))
  {
    if(myDataInt$pos[i,1]>=roi[1,1] & myDataInt$pos[i,1]<=roi[2,1]
       & myDataInt$pos[i,2]>=roi[1,2] & myDataInt$pos[i,2]<=roi[2,2]) 
          {pxList[z]<-i; z<-z+1;}
  }
  
  xRange=c(1, xSize*ySize);
  yRange=c(mzLow, mzHigh);
  
  plot(xRange, yRange, type="p", col="white", cex=0.01, main="Deconvolution", xlab="pixels", ylab="mz(Da)", las=1, col.axis="black");
  #  plot(X, Y, type="p",col="blue",lwd=1, main="Deconvolution", xlab="mz(Daltons)", 
  #       ylab="Concentration", las=1, col.axis="black");  
  
  for(px in pxList)
  {
    pxPeaksList<-rMSI2::readimzML_singlePixelPeakList(peaksFile, px);
    logicIndex<-pxPeaksList$mass>=mzLow & pxPeaksList$mass<=mzHigh;
    nMass<-length(pxPeaksList$mass[logicIndex]);
    
    #      X<-pxPeaksList$intensity[logicIndex];
    Y<-pxPeaksList$mass     [logicIndex];
    X<-rep(px, times=nMass);
    #X<-pxPeaksList$intensity[logicIndex];
    if(nMass>0)
      lines(X-pxList[1], Y, type="p", col="blue", cex=0.2);  
    
    
  }
  
  xyLow<-myDataInt$pos[pxList[1],];
  xyHigh<-myDataInt$pos[pxList[xSize*ySize],];
  sprintf("low(X/Y):(%d/%d) high(x/y):(%d/%d)", xyLow[1], xyLow[2], xyHigh[1], xyHigh[2]);
  
  siremPeaks<-rGetSiremPeaks(myDataInt, params, mzLow, mzHigh);
  myList<-list();
  index=1;
  maxValue=0;
  minValue=1e10;
  for(px in pxList)
  {
    gaussInfo<-rGetGaussians(myDataInt, siremPeaks, px, params$minMeanPxMag);
    nGauss<-nrow(gaussInfo$gaussians); #number of Gaussians to represent.
    Y<-gaussInfo$gaussians[,1];
#    X<-rep(px, times=nGauss);
#    if(nGauss>0)
#      lines(X-pxList[1], Y, type="p", col="red", cex=0.2);  
    if(nGauss>0)
      {
      myList[[index]]<-Y;
      }
    else 
      {myList[[index]]<--1;}
    index<-index+1;
  }
  
  index=1;
  for(px in pxList)
  {
      if(myList[index]>=0)
      {
        X<-rep(px, times=length(myList[index]));
        lines(X-pxList[1], myList[index], type="p", col="red", cex=0.2);  
      }
    
  }
  #gaussians from average spectrum
  gaussInfo<-rGetGaussians(myDataInt, siremPeaks, 0, params$minMeanPxMag); #sp promediado
  
  #Lineas discontinuas: mean of gaussians from mean spectrum
  for(i in 1:length(gaussInfo$gaussians[,1]))
  {
    Y<-rep(gaussInfo$gaussians[i, 1], times=length(pxList))
    lines(pxList-pxList[1], Y, type="l", col="green", lwd=2, lty=2); 
  }
  
  return(0);
}
