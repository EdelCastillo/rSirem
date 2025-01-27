#########################################################################
#     rSirem - R package for MSI data deconvolution
#     Copyright (C) marzo 2024, Esteban del Castillo Pérez
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

#funciones para validar la eficacia de Sirem
#------------------------------------------------------------------------------------------------------
#
#' fitQualitySiremDeconv()
#' Segrega los picos compuestos en los que existe deconvolución, indica los picos simples que lo componen
#' e identifica el pico de alta resolución con el que tiene el mejor ajuste
#' 
#' fitQualitySirem()
#' Identifica el pico de alta resolución más próximo al pico de baja resolución e indica la desviación
#' Lo hace sobre todos los picos. Para info de Sirem
#' 
#' fitQualityPere()
#' Identifica el pico de alta resolución más próximo al pico de baja resolución e indica la desviación
#' Lo hace sobre todos los picos. Para info de rMSI2 (Pere)
#' 
#' peaksDeconvolved()
#' Para cada pico de magnitud compuesto que contiene a picos simples de sirem,
#' se anotan los picos de sirem incluidos.
#' 
#' getUpeakSiremMassesFromUpeak()
#' retorna un array con las masas de sirem asociadas a un pico de magnitud
#' 
#' getUpeakGaussFromMass()
#' retorna un array con las masas de sirem contenidas en el mismo pico compuesto que la masa dada
#' 
#' siremPeaksFilter()
#' Retorna un array indicando si el pico compuesto debe ser considerado o no
#' Sólo se consideran aquellos que contienen a un pico de rMSI2
#' 
#' sirem_vs_rMSI2()
#' genera datos estadísticos con las desviaciones de ambos algoritmos respecto a una muestra patrón


#functions code
#------------------------------------------------------------------------------------------------------

#' fitQualitySirem
#' Identifica el pico de alta resolución más próximo al pico de baja resolución e indica la desviación
#' Lo hace sobre todos los picos. Para info de Sirem
#' 
#' @param reference  -> gaussians info from rGetGaussians() (high resolution)
#' @param testSirem  -> sirem     info from rGetGaussians() (low resolution)
#' @param testGauss  -> gaussians info from rGetGaussians() (low resolution)
#' @param refMinMag  -> mínima magnitud de picos para su consideración (high resolution)
#' @param testMinMag -> mínima magnitud de picos para su consideración (low resolution)
#' @param minMass    -> low  mass to analyze (Da)
#' @param minMass    -> high mass to analyze (Da)
#'
#' @return a matrix with columns: "mzTest", "mzRef", "ppm", "maxDev", "repe", "deconv" 
#'  mzTest -> mz from low  resolution sample
#'  mzRef  -> mz from high resolution sample
#'  ppm    -> deviation from low to high mz
#'  maxDev -> '1' if greater than ppm of one scan
#'  repe   -> '1' if mzRef is the same for two consecutives mzTest 
#'  deconv -> '1' if is a deconvolved peak
#' @export
#' 
fitQualitySirem<-function(reference, testSirem, testGauss, refMinMag=1e-6, testMinMag=1e-6, minMass=0, maxMass=0, histo=FALSE)
{
  fail<-matrix(nrow=2, ncol=5);
  nSingletestGaussPks=length(testGauss$gaussians[,1]);
  nSingleRefPks=length(reference$gaussians[,1]);
  
  deviation <-matrix(nrow = nSingletestGaussPks, ncol = 6);
  deviation[,1]=rep(0, times=nSingletestGaussPks);
  deviation[,2]=rep(0, times=nSingletestGaussPks);
  deviation[,3]=rep(0, times=nSingletestGaussPks);
  deviation[,4]=rep(0, times=nSingletestGaussPks);
  deviation[,5]=rep(0, times=nSingletestGaussPks);
  deviation[,6]=rep(0, times=nSingletestGaussPks);
  pksDeconv<-peaksDeconvolved(testSirem, testGauss);
  
  maxRefMag=max(reference$gaussians[,3]);
  if(refMinMag>1e-6) refMinMag=maxRefMag*refMinMag/100;
  maxTestMag=max(testGauss$gaussians[,3]);
  if(testMinMag>1e-6) testMinMag=maxTestMag*testMinMag/100;
  
  minMassRef=reference$gaussians[1,1];
  maxMassRef=reference$gaussians[nSingleRefPks,1];
  minMassSirem=testGauss$gaussians[1,1];
  maxMassSirem=testGauss$gaussians[nSingletestGaussPks,1];

  if(minMass==0 | minMass<minMassSirem) minMass<-minMassSirem;
  if(maxMass==0 | maxMass>maxMassSirem) maxMass<-maxMassSirem;
  
  for(iPk in 1:nSingletestGaussPks) #para cada pico del test
  {
    if(testGauss$gaussians[iPk, 1]==0) next;
    offset=0;
    iMass=rGetIndexFromMass(testGauss$xAxis[iPk], testGauss$xAxis);
    if(iPk==nSingletestGaussPks)
      {testPPM=1e6*(testGauss$xAxis[iMass]-testGauss$xAxis[iMass-1])/testGauss$xAxis[iMass-1];}
    else
      {testPPM=1e6*(testGauss$xAxis[iMass+1]-testGauss$xAxis[iMass])/testGauss$xAxis[iMass];}
    
    testMass=testGauss$gaussians[iPk, 1];
    retMass<-nearestValue(testMass, reference$gaussians[, 1]);
    if(retMass==-1) {offset=-1;}
    else {offset<-abs(testMass-retMass)}
    ppm<-1e6*offset/testMass;
    
    deviation[iPk, 1]=testMass;
    deviation[iPk, 2]=retMass;
    deviation[iPk, 3]=ppm;
    
    if(ppm>1.5*testPPM) #desviación > 1.5*resolución de masa (low resolution)
      {deviation[iPk, 4]=2;}# deviation[iPk, 3]=0;}
    else if(ppm>testPPM) #desviación >1 && <= 1.5*resolución 
      {deviation[iPk, 4]=1;}
    else                #desviación <=1*resolución 
      {deviation[iPk, 4]=0;}
    
    deviation[iPk, 5]=0;
    if(iPk>1)
      if(deviation[iPk, 2]==deviation[iPk-1, 2]) #comparten misma masa de ref.
      {deviation[iPk, 5]=1;}
    
    deviation[iPk, 6]=0;    
    if(length(pksDeconv$gaussMassList)>0)
    {
      retMass<-nearestValue(testMass, pksDeconv$gaussMassList);
      if(retMass==testMass)
        deviation[iPk, 6]=1; #pico deconvolucionado
    }
  }
  if(histo==TRUE)
    {
    hist(deviation[, 3], main="Histogram of rSirem deviations", xlab="ppm", ylab="Frequency");
    #  legend("topright", legend="SNR=1");
    }
  colnames(deviation)<-c("mzTest", "mzRef", "ppm", "maxDev", "repe", "deconv")
  return(deviation);
}

#' fitQualitySiremDeconv
#' Identifica el pico de alta resolución más próximo al pico de baja resolución e indica la desviación
#' Lo hace exclusivamente sobre los picos deconvolucionados. Para info de Sirem
#' 
#' @param reference  -> gaussians info from rGetGaussians() (high resolution)
#' @param testSirem  -> sirem     info from rGetGaussians() (low resolution)
#' @param testGauss  -> gaussians info from rGetGaussians() (low resolution)
#' @param refMinMag  -> mínima magnitud de picos para su consideración (high resolution)
#' @param testMinMag -> mínima magnitud de picos para su consideración (low resolution)
#' @param minMass    -> low  mass to analyze (Da)
#' @param minMass    -> high mass to analyze (Da)
#'
#' @return a matrix with columns: "mzTest", "mzRef", "ppm", "maxDev", "repe" 
#'  mzTest -> mz from low  resolution sample
#'  mzRef  -> mz from high resolution sample
#'  ppm    -> deviation from low to high mz
#'  maxDev -> '1' if greater than ppm of a scan
#'  repe   -> '1' if mzRef is the same for two consecutives mzTest 
#' @export
#' 
fitQualitySiremDeconv<-function(reference, testSirem, testGauss, refMinMag=1e-6, testMinMag=1e-6, minMass=0, maxMass=0, histo=FALSE)
{
  fail<-matrix(nrow=2, ncol=5);
  
  pksDeconv<-peaksDeconvolved(testSirem, testGauss);
  
  nSingleRefPks=length(reference$gaussians[,1]);
  nSingletestGaussPks=length(pksDeconv$gaussMassList);
  if(nSingletestGaussPks==0)
    {print("Warning: there are not deconvolved peaks"); return(fail);}
    
  deviation <-matrix(nrow = nSingletestGaussPks, ncol=5);
  
  minMassRef=reference$gaussians[1,1];
  maxMassRef=reference$gaussians[nSingleRefPks,1];
  minMassSirem=pksDeconv$gaussMassList[1];
  maxMassSirem=pksDeconv$gaussMassList[nSingletestGaussPks];
  
  if(minMass==0 | minMass<minMassSirem) minMass<-minMassSirem;
  if(maxMass==0 | maxMass>maxMassSirem) maxMass<-maxMassSirem;
  
  for(iPk in 1:nSingletestGaussPks) #para cada pico del test
  {
    offset=0;
    iMass=rGetIndexFromMass(pksDeconv$gaussMassList[iPk], testGauss$xAxis);
    if(iPk==nSingletestGaussPks)
    {testPPM=1e6*(testGauss$xAxis[iMass]-testGauss$xAxis[iMass-1])/testGauss$xAxis[iMass-1];}
    else
    {testPPM=1e6*(testGauss$xAxis[iMass+1]-testGauss$xAxis[iMass])/testGauss$xAxis[iMass];}
    
    testMass=pksDeconv$gaussMassList[iPk];
    retMass<-nearestValue(testMass, reference$gaussians[, 1]);
    if(retMass==-1) {offset=-1;}
    else {offset<-abs(testMass-retMass)}
    ppm<-1e6*offset/testMass;
    
    deviation[iPk, 1]=testMass;
    deviation[iPk, 2]=retMass;
    deviation[iPk, 3]=ppm;
    
    if(ppm>1.5*testPPM) 
      {deviation[iPk, 4]=2;}# deviation[iPk, 3]=0;}
    else if(ppm>testPPM) 
      {deviation[iPk, 4]=1;}
    else
      {deviation[iPk, 4]=0;}
    
    deviation[iPk, 5]=0;
    if(iPk>1)
      if(deviation[iPk, 2]==deviation[iPk-1, 2])
      {deviation[iPk, 5]=1;}
    
  }
  if(histo==TRUE)
    {
    hist(deviation[, 3], main="Histogram of deconvolved rSirem deviations", xlab="ppm", ylab="Frequency");
  #  legend("topright", legend="SNR=1");
    }
  colnames(deviation)<-c("mzTest", "mzRef", "ppm", "maxDev", "repe")
  return(deviation);
}


#' fitQualityPere
#' Identifica el pico de alta resolución más próximo al pico de baja resolución e indica la desviación
#' Lo hace sobre todos los picos. Para info de rMSI2 (Pere)
#' 
#' @param reference  -> gaussians info from rGetGaussians() (high resolution)
#' @param testPere   -> list from  rMSI2::LoadPeakMatrix() from rMSI2::processWizard()
#' @param testSirem  -> gaussians info from rGetGaussians() used for mass resolution
#' @param refMinMag  -> mínima magnitud de picos para su consideración (high resolution)
#' @param testMinMag -> mínima magnitud de picos para su consideración (low resolution)
#' @param minMass    -> low  mass to analyze (Da)
#' @param minMass    -> high mass to analyze (Da)
#'
#' @return a matrix with columns: "mzTest", "mzRef", "ppm", "maxDev", "repe" 
#'  mzTest -> mz from low  resolution sample
#'  mzRef  -> mz fron high resolution sample
#'  ppm    -> deviation from low to high mz
#'  maxDev -> '1' if greater than ppm of a scan
#'  repe   -> '1' if mzRef is the same for two consecutives mzTest 
#' @export
#' 
fitQualityPere<-function(reference, testPere, testSirem, refMinMag=1e-6, testMinMag=1e-6, minMass=0, maxMass=0, histo=FALSE)
{
  fail<-matrix(nrow=2, ncol=5);
  nSingleTestPerePks=length(testPere$mass);
  nSingleRefPks=length(reference$gaussians[,1]);
  
  colMeansPere<-colMeans(testPere$intensity);
  
  maxRefMag=max(reference$gaussians[,3]);
  if(refMinMag>1e-6) refMinMag=maxRefMag*refMinMag/100;
  maxTestMag=max(colMeansPere);
  if(testMinMag>1e-6) testMinMag=maxTestMag*testMinMag/100;
  
  minMassRef=reference$gaussians[1,1];
  maxMassRef=reference$gaussians[nSingleRefPks,1];
  
  pereMassLogic=testPere$mass>=minMassRef & testPere$mass<=maxMassRef;
  pereMass=testPere$mass[pereMassLogic];
  minMassPere=pereMass[1];
  iMaxPere=length(pereMass);
  maxMassPere=pereMass[iMaxPere];
  
  if(minMass==0 | minMass<minMassPere) minMass<-minMassPere;
  if(maxMass==0 | maxMass>maxMassPere) maxMass<-maxMassPere;

  refMassLogic=reference$gaussians[,1]>=minMass & reference$gaussians[,1]<=maxMass;
  refMass=reference$gaussians[,1][refMassLogic];
  
  siremMassLogic=testSirem$xAxis>=minMass & testSirem$xAxis<=maxMass;
  siremMass=testSirem$xAxis[siremMassLogic];
  iMaxSirem=length(siremMass);
  
  iMinMass=rGetIndexFromMass(minMass, testPere$mass);
  iMaxMass=rGetIndexFromMass(maxMass, testPere$mass);
  
  deviation <-matrix(nrow = iMaxPere, ncol = 5);
  
  for(iPk in 1:iMaxPere) #para cada pico del test
  {
    offset=0;
    nearMassSirem<-nearestValue(pereMass[iPk], siremMass);
    iNearMass=rGetIndexFromMass(nearMassSirem, siremMass);
    
    if(iNearMass>=iMaxSirem)
    {testPPM=1e6*(siremMass[iMaxSirem]-siremMass[iMaxSirem-1])/siremMass[iMaxSirem-1];}
    else
    {testPPM=1e6*(siremMass[iNearMass+1]-siremMass[iNearMass])/siremMass[iNearMass];}

    testMass=pereMass[iPk];
    retMass<-nearestValue(testMass, refMass);
    if(retMass==-1) {offset=-1;}
    else {offset<-abs(testMass-retMass)}
    ppm<-1e6*offset/testMass;
    
    deviation[iPk, 1]=testMass;
    deviation[iPk, 2]=retMass;
    deviation[iPk, 3]=ppm;
    
    if(ppm>1.5*testPPM) 
      {deviation[iPk, 4]=2;}# deviation[iPk, 3]=0;}
    else if(ppm>testPPM) 
      {deviation[iPk, 4]=1;}
    else
      {deviation[iPk, 4]=0;}
    
    if(iPk>1)
    {
      if(deviation[iPk, 2]==deviation[iPk-1, 2]) {deviation[iPk, 5]=1;}
      else {deviation[iPk, 5]=0;}
    }
    else {deviation[iPk, 5]=0;}
  }
  if(histo==TRUE)
    {
    hist(deviation[, 3], main="Histogram of rMSI2 deviations", xlab="ppm", ylab="Frequency");
#  title("rMSI2 deviations");
#  legend("topright", legend="SNR=1");
    }
  colnames(deviation)<-c("mzTest", "mzRef", "ppm", "maxDev", "repe")
  return(deviation);
}

#' fitQualityPereSirem
#' Aisla aquellos picos de rMSI2 y rSirem que comparten al mismo pico de alta resolución como más próximo
#' Para ellos, se obtienen las desviaciones respecto al pico de alta resolución
#' 
#' @param reference  -> gaussians info from rGetGaussians() (high resolution)
#' @param testPere   -> list from  rMSI2::LoadPeakMatrix() from rMSI2::processWizard()
#' @param testSirem  -> gaussians info from rGetGaussians() (low resolution)
#' @param refMinMag  -> mínima magnitud de picos para su consideración (high resolution)
#' @param testMinMag -> mínima magnitud de picos para su consideración (low resolution)
#' @param minMass    -> low  mass to analyze (Da)
#' @param minMass    -> high mass to analyze (Da)
#'
#' @return a matrix with columns: "mzrMSI2", "mzRef", "ppm", "maxDev", "mzSirem", "mzRef", "ppm", "maxDev" 
#'  mzrMSI2 -> mz from rMSI2 sample
#'  mzRef   -> mz from high resolution sample
#'  ppm     -> deviation from rMSI2 to high mz
#'  maxDev -> '1' if greater than ppm of a scan; 2 if greater than 1.5 ppm of a scan
#'  mzSirem -> mz from Sirem sample
#'  mzRef   -> mz from high resolution sample
#'  ppm     -> deviation from SIrem to high mz
#'  maxDev -> '1' if greater than ppm of a scan; 2 if greater than 1.5 ppm of a scan
#'  
#' @export
#' 
fitQualityPereSirem<-function(reference, testPere, testSirem, refMinMag=1e-6, testMinMag=1e-6, minMass=0, maxMass=0, histo=FALSE)
{
  fail<-matrix(nrow=2, ncol=5);
  nSingleTestPerePks=length(testPere$mass);
  nSingleRefPks=length(reference$gaussians[,1]);
  
  colMeansPere<-colMeans(testPere$intensity);
  
  maxRefMag=max(reference$gaussians[,3]);
  if(refMinMag>1e-6) refMinMag=maxRefMag*refMinMag/100;
  maxTestMag=max(colMeansPere);
  if(testMinMag>1e-6) testMinMag=maxTestMag*testMinMag/100;
  
  minMassPere=min(testPere$mass);
  maxMassPere=max(testPere$mass);
  
  if(minMass==0 | minMass<minMassPere) minMass<-minMassPere;
  if(maxMass==0 | maxMass>maxMassPere) maxMass<-maxMassPere;
  
  refGaussLogic=reference$gaussians[,1]>=minMass & reference$gaussians[,1]<=maxMass;
  refGauss=reference$gaussians[,1][refGaussLogic];
  
  refMassLogic=reference$xAxis>=minMass & reference$xAxis<=maxMass;
  refMass=reference$xAxis[refMassLogic];
  
  siremMassLogic=testSirem$xAxis>=minMass & testSirem$xAxis<=maxMass;
  siremMass=testSirem$xAxis[siremMassLogic];
  iMaxSirem=length(siremMass);
  
  siremGaussLogic=testSirem$gaussians[,1]>=minMass & testSirem$gaussians[,1]<=maxMass;
  siremGauss=testSirem$gaussians[,1][siremGaussLogic];
  
  minMassRef=min(refMass);
  maxMassRef=max(refMass);
  
  pereMass=testPere$mass;  
  iMinMass=rGetIndexFromMass(minMassRef, testPere$mass);
  iMaxMass=rGetIndexFromMass(maxMassRef, testPere$mass);
  
  deviation <-matrix(nrow = iMaxMass-iMinMass+1, ncol = 8);
  iMatrix=1;
#  for(iPk in 1:iMaxPere) #para cada pico del test
  for(iPk in iMinMass:iMaxMass) #para cada pico del test
    {
    offset=0;
    nearMassSirem<-nearestValue(pereMass[iPk], siremMass);
    iNearMass=rGetIndexFromMass(nearMassSirem, siremMass);
    
    if(iNearMass>=iMaxSirem)
    {testPPM=1e6*(siremMass[iMaxSirem]-siremMass[iMaxSirem-1])/siremMass[iMaxSirem-1];}
    else
    {testPPM=1e6*(siremMass[iNearMass+1]-siremMass[iNearMass])/siremMass[iNearMass];}
    
    rMSI2Mass=pereMass[iPk];
    highResMassNearestPere<-nearestValue(rMSI2Mass, refGauss); #gaussiana highRes más próxima a Pere
    
    lowResMassNearestMSI2 <-nearestValue(rMSI2Mass, siremGauss); #gaussiana lowRes más próxima a Pere
    highResMassNearestSirem<-nearestValue(lowResMassNearestMSI2, refGauss); #gaussiana highRes más próxima a Sirem
    if(highResMassNearestSirem != highResMassNearestPere) next; #si se trata de gaussianas distintas
    
    lowResMassNearestHighResMass <-nearestValue(highResMassNearestPere, siremGauss);
    
    if(highResMassNearestPere==-1 | lowResMassNearestHighResMass==-1) {offset=-1;}
    else 
    {
        offsetPere<-abs(rMSI2Mass-highResMassNearestPere);
        offsetSirem<-abs(lowResMassNearestHighResMass-highResMassNearestPere)
        }
    ppmPere<-1e6*offsetPere/rMSI2Mass;
    ppmSirem<-1e6*offsetSirem/lowResMassNearestHighResMass;
    
    deviation[iMatrix, 1]=rMSI2Mass;
    deviation[iMatrix, 2]=highResMassNearestPere;
    deviation[iMatrix, 3]=ppmPere;
    
    if(ppmPere>1.5*testPPM) 
      {deviation[iMatrix, 4]=2;}# deviation[iMatrix, 3]=0;}
    else if(ppmPere>testPPM) 
      {deviation[iMatrix, 4]=1;}
    else
      {deviation[iMatrix, 4]=0;}

    deviation[iMatrix, 5]=lowResMassNearestHighResMass;
    deviation[iMatrix, 6]=highResMassNearestPere;
    deviation[iMatrix, 7]=ppmSirem;
    
    if(ppmSirem>1.5*testPPM) 
      {deviation[iMatrix, 8]=2; deviation[iMatrix, 7]=0;}
    else if(ppmSirem>testPPM) 
      {deviation[iMatrix, 8]=1;}
    else
      {deviation[iMatrix, 8]=0;}
  iMatrix=iMatrix+1;
  }
  if(histo==TRUE)
    {
    hist(deviation[, 3], main="Histogram of rMSI2 deviations", xlab="ppm", ylab="Frequency");
    hist(deviation[, 7], main="Histogram of rSIREM deviations", xlab="ppm", ylab="Frequency");
  #  title("rMSI2 deviations");
  #  legend("topright", legend="SNR=1");
    }
  colnames(deviation)<-c("mzrMSI2", "mzRef", "ppm", "maxDev", "mzSirem", "mzRef", "ppm", "maxDev")
  return(deviation[1:iMatrix-1,]);
}



#' nearestValue
#' Return the nearest value in data
#' Algoritmo de aproximaciones sucesivas
#' 
#' @param value -> reference value
#' @param data  -> array of sort values
#'
#' @return nearest value in data to value; -1 if value es out of range
#' @export
#' 
nearestValue<-function(value, data)
{
  indexLow<-1;
  indexHigh<-length(data);
  
  if(indexHigh==indexLow) return(data[1]);
  if(indexHigh==indexLow+1)
  {
    if(value-data[indexLow] <= data[indexHigh]-value) {return(data[indexLow]);}
    else {return(data[indexHigh]);}
  }
  
  if(value<data[indexLow])       {return(data[indexLow]);}
  else if(value>data[indexHigh]) {return(data[indexHigh]);}
  if(value==data[indexHigh])     return(value);
  if(value==data[indexLow])      return(value);
  
  while(1)
    {
    indexCenter<-round((indexHigh+indexLow)/2);
    if(value==data[indexCenter]) return(value);
    if(value<data[indexCenter]) {indexHigh<-indexCenter; }
    else {indexLow <-indexCenter;}
    if(indexHigh==indexLow+1)
    {
      if(value-data[indexLow] <= data[indexHigh]-value) {return(data[indexLow]);}
      else {return(data[indexHigh]);}
    }
    
  }
}

#' peaksDeconvolved
#' Para cada pico de magnitud compuesto que contiene a más de un pico simples de sirem,
#' se anotan los picos de sirem incluidos y las gausianas asociadas
#' 
#' @param peaksInfo -> sirem peaks, from rGetSiremPeaks()
#' @param gaussInfo -> gaussians,   from rGetGaussians()
#' @param uPeakList -> lista de picos compuestos a considerar, esten o no deconvolucionados
#' 
#' @return a list: siremMassMatrix, gaussMassMatrix, gaussMagMatrix, gaussMassList
#'  siremMassMatrix mantiene las masas asociadas a los picos de sirem
#'  gaussMassMatrix mantiene las masas asociadas a las gausianas luego del ajuste con el valor promediado de magnitudes
#'  gaussMagMatrix  mantiene las magnitudes asociadas a cada gausiana  
#'  gaussList es un array con las gausianas de gaussMassMatrix
#'  siremMassMatrix & gaussMassMatrix and gaussMagMatrix have the same format:
#'          rows=united magPeaks deconvolved;
#'          col[] =central mass of each peak
#'          col[1]=index to magnitude united Peak
#'          col[2]=low  mass to magnitude united peak
#'          col[3]=high mass to magnitude united peak
#'          col[4]=number of single peaks into magnitude peak
#'          col[5...]=list of masses into magnitude united peaks
#'          
#'  #attn: en peakInfo, si en un pico de  magnitud no se encuentra ningún máximo de sirem,
#          viene insertado un pico de sirem coincidente con el pico de magnitud.
#' @export
#' 
peaksDeconvolved<-function(peaksInfo, gaussInfo, uPeakList=c())
{
  fail<-list(siremMassMatrix=0, gaussMassMatrix=0, gaussMagMatrix=0, gaussMassList=c());
  if(length(uPeakList)==0) uPeakList=c(-1);
  
    nScans=length(peaksInfo$siremPeaks$sirem);
    nMagPeaks=length(peaksInfo$siremPeaks$magnitudePeaks[,1]);
    nSiremPeaks=length(peaksInfo$siremPeaks$siremPeaks[,1]);
    nUnitedMagPeaks=length(peaksInfo$siremPeaks$unitedMagnitudePeaks[,1]);
           
    retRow=0;
    sPeaks=rep(0, times=100);
    iSirem=1;
    maxSiremPeaks=0;
    uPeaksCount=0;
    for(united in 1:nUnitedMagPeaks)
    {
      uPeaksLow =peaksInfo$siremPeaks$unitedMagnitudePeaks[united, 1]+1;
      if(uPeaksLow==0) next; #descartado
      uPeaksHigh=peaksInfo$siremPeaks$unitedMagnitudePeaks[united, 2]+1;
      nUpeaks=uPeaksHigh-uPeaksLow+1;
      k=0;
      iLowMag <-peaksInfo$siremPeaks$magnitudePeaks[uPeaksLow,  1]+1;#indice a mz low
      iHighMag<-peaksInfo$siremPeaks$magnitudePeaks[uPeaksHigh, 3]+1;#indice a mz high
      for(j in iSirem:nSiremPeaks)
      {
        iMaxSirem<-peaksInfo$siremPeaks$siremPeaks[j, 2]+1;
        if(iMaxSirem>iHighMag) {iSirem=j; break;}
        if(iMaxSirem>=iLowMag & iMaxSirem<=iHighMag)
        {
          k=k+1;
        }
      }
      if(k>nUpeaks | !is.na(match(united, uPeakList)))
      {
        uPeaksCount=uPeaksCount+1;
        if(k>maxSiremPeaks) maxSiremPeaks=k;
      }
    }
    #retM mantiene los índices a las masas centrales de cada pico de sirem
    #retM[,1]=referencia al pico compuesto
    #retM[,2]=índice a la masa mínima del pico de magnitud compuesto
    #retM[,3]=índice a la masa máxima del pico de magnitud compuesto
    #retM[,4]=número de picos de sirem dentro del pico compuesto
    #retM[,5]=secuencia ordenada de las masas centrales de cada pico de sirem
    retM<-matrix(nrow=uPeaksCount, ncol=maxSiremPeaks+4);
    
    iSirem=1;
    #se buscan picos de magnitud que contengan a picos de sirem
    #attn: si en un pico de  magnitud no se encuentra ningún máximo de sirem,
    #se inserta un pico de sirem coincidente con el pico de magnitud
    for(united in 1:nUnitedMagPeaks)
    {
      uPeaksLow =peaksInfo$siremPeaks$unitedMagnitudePeaks[united, 1]+1;
      if(uPeaksLow==0) next; #descartado
      uPeaksHigh=peaksInfo$siremPeaks$unitedMagnitudePeaks[united, 2]+1;
      nUpeaks=uPeaksHigh-uPeaksLow+1;
      k=0;
      iLowMag <-peaksInfo$siremPeaks$magnitudePeaks[uPeaksLow,  1]+1;#indice a mz low
      iHighMag<-peaksInfo$siremPeaks$magnitudePeaks[uPeaksHigh, 3]+1;#indice a mz high
      for(j in iSirem:nSiremPeaks)
      {
        iMaxSirem<-peaksInfo$siremPeaks$siremPeaks[j, 2]+1;
        if(iMaxSirem>iHighMag) {iSirem=j; break;}
        if(iMaxSirem>=iLowMag & iMaxSirem<=iHighMag)
        {
          k=k+1;
          sPeaks[k]=j;
        }
      }
    
      #deben haber más picos de sirem que de magnitud para que haya deconv
      #interesan los picos deconvolucionados->más de un pico de sirem
      #k>nUpeaks -> uPeaks deconvolved
      #!is.na(match(united, uPeakList) -> true if united is into uPeakList
      if(k>nUpeaks | !is.na(match(united, uPeakList)))
      {
        retRow=retRow+1;
        retM[retRow,]=0;
        retM[retRow, 1]<-united;        #united mag peak index
        retM[retRow, 2]<-iLowMag;       #mag peak low index
        retM[retRow, 3]<-iHighMag;      #mag peak high index
        retM[retRow, 4]<-k;             #number of gaussians in united peaks
        for(n in 1:k)
        {
        retM[retRow, n+4]<-peaksInfo$siremPeaks$siremPeaks[sPeaks[n],2]; #sirem peak index ini
        }
      }
    }
    
    if(retRow==0){print("warning: there are not peaks deconvolved"); return(fail);}
    
    #massM mantiene las masas centrales de cada pico de sirem
    #massM[,1]=referencia al pico compuesto
    #massM[,2]=índice a la masa mínima del pico de magnitud compuesto
    #massM[,3]=índice a la masa máxima del pico de magnitud compuesto
    #massM[,4]=número de picos de sirem dentro del pico compuesto
    #massM[,5]=secuencia ordenada de las masas centrales de cada pico de sirem
    massM<-matrix(nrow=uPeaksCount, ncol=maxSiremPeaks+4);
    for(i in 1:retRow)
    {
      massM[i,]=0;
      massM[i,1]=retM[i,1]; #masa al máximo del pico de mag
      massM[i,2]=peaksInfo$massAxis[retM[i,2]];
      massM[i,3]=peaksInfo$massAxis[retM[i,3]];
      massM[i,4]=retM[i,4]; #numero de picos de sirem
      
      for(j in 1:retM[i,4])
        {
        massM[i, 4+j]=peaksInfo$massAxis[retM[i, 4+j]];
        }
    }
    
    #matriz de masas de las gausianas
    #Se seleccionan los valores medios de cada gausiana incluidos en el pico compuesto
    massGausM<-matrix(nrow=uPeaksCount, ncol=maxSiremPeaks+4);
    gaussMeanPeaks<-gaussInfo$gaussians[,1];
    massList<-rep(0, maxSiremPeaks);
    iMassList=1;
    
    for(uPeak in 1:retRow)
    {
      massGausM[uPeak,]=0;
      massGausM[uPeak,1]=massM[uPeak, 1];
      massGausM[uPeak,2]=massM[uPeak, 2];
      massGausM[uPeak,3]=massM[uPeak, 3];
      massGausM[uPeak,4]=massM[uPeak, 4];
      massGausLogic<-(gaussMeanPeaks>=massM[uPeak, 2] & gaussMeanPeaks<=massM[uPeak, 3]);
      massGaussA<-gaussMeanPeaks[massGausLogic];
      for(j in 1:massM[uPeak, 4])
      {
#        nearestMass=nearestValue(massM[uPeak, 2+j], gaussMeanPeaks);
        massGausM[uPeak, 4+j]=massGaussA[j];
        massList[iMassList]=massGaussA[j];
        iMassList=iMassList+1;
      }
    }
    
    #matriz de magnitudes de las gausianas
    #Se seleccionan los valores de magnitud de cada gausiana incluidos en el pico compuesto
    magGausM<-matrix(nrow=uPeaksCount, ncol=maxSiremPeaks+4);
    gaussMagPeaks<-gaussInfo$gaussians[,3];
    magList<-rep(0, maxSiremPeaks);
    iMagList=1;

    for(uPeak in 1:retRow)
    {
      magGausM[uPeak,]=0;
      magGausM[uPeak,1]=massM[uPeak, 1];
      magGausM[uPeak,2]=massM[uPeak, 2];
      magGausM[uPeak,3]=massM[uPeak, 3];
      magGausM[uPeak,4]=massM[uPeak, 4];
      massGausLogic<-(gaussMeanPeaks>=massM[uPeak, 2] & gaussMeanPeaks<=massM[uPeak, 3]);
      magGaussA<-gaussMagPeaks[massGausLogic];
      for(j in 1:massM[uPeak, 4])
      {
        magGausM[uPeak, 4+j]=magGaussA[j];
        iMagList=iMagList+1;
      }
    }
    mColNames=c("uPks", "mzInit", "mzEnd", "nPks")
    for(i in 1:(ncol(massGausM)-4))
    {
      newStr<-sprintf("mzPk_%d", i);
      mColNames<- c(mColNames, newStr);
    }
    colnames(massM)     <-mColNames;
    colnames(massGausM) <-mColNames;
    colnames(magGausM)  <-mColNames;
    ret<-list(siremMassMatrix=massM, gaussMassMatrix=massGausM, gaussMagMatrix=magGausM, gaussMassList=massList);
    return(ret);
}


#' getUpeakSiremMassesFromUpeak
#' retorna un array con las masas de sirem asociadas a un pico de magnitud
#' 
#' @param siremPeaks -> sirem peaks info, from rGetSiremPeaks()
#' @param uPeak      -> united peaks
#' @return array con las masas asociadas al pico dado
#' 
getUpeakSiremMassesFromUpeak<-function(siremPeaks, uPeak)
{
  if(uPeak<0 | uPeak>length(siremPeaks$siremPeaks$unitedMagnitudePeaks[,1]))
  {print("warning: uPeak out of range"); return(-1);}
  
  sPeaks=siremPeaks$siremPeaks$siremPeaks[,1];
  #se buscanlos  picos de sirem contenidos en el pico de magnitud dado
  #attn: si en un pico de  magnitud no se encuentra ningún máximo de sirem,
  #se inserta un pico de sirem coincidente con el pico de magnitud
    uPeaksLow =siremPeaks$siremPeaks$unitedMagnitudePeaks[uPeak, 1]+1; #pico inferior unido
    uPeaksHigh=siremPeaks$siremPeaks$unitedMagnitudePeaks[uPeak, 2]+1; #pico superior unido
    
    iLowMag <-siremPeaks$siremPeaks$magnitudePeaks[uPeaksLow,  1];#indice a mz low
    iHighMag<-siremPeaks$siremPeaks$magnitudePeaks[uPeaksHigh, 3];#indice a mz high
    siremPksLogic=sPeaks>=iLowMag & sPeaks<=iHighMag;
    sPeaks=sPeaks[siremPksLogic];
    maxSiremUpeaks=rep(0, times=length(sPeaks));
    for(i in 1:length(sPeaks))
    {
      #índice a la masa que identifica al máximo del pico de sirem
      iMaxSiremMz =siremPeaks$siremPeaks$siremPeaks[sPeaks[i], 2];
      
      #masa que identifica al máximos del pico de sirem
      #attn: los índices en peaksInfo son relativos a cero!!
      #R tine al '1' como índice mínimo
      maxSiremUpeaks[i]=siremPeaks$massAxis[iMaxSiremMz+1];
    }
  return(maxSiremUpeaks);
  }

#' getUpeakGaussFromMass
#' retorna un array con las masas de sirem contenidas en el mismo pico compuesto que la masa dada
#' 
#' @param siremPeaks -> sirem peaks, from rGetSiremPeaks()
#' @param mass       -> masa 
#' @return array con las masas asociadas
#' 
getUpeakGaussFromMass<-function(siremPeaks, gaussInfo, mass)
{
  nUpeaks<-length(siremPeaks$siremPeaks$unitedMagnitudePeaks[,1]);
  
  #se buscan los  picos de sirem contenidos en el pico de magnitud dado
  #attn: si en un pico de  magnitud no se encuentra ningún máximo de sirem,
  #se inserta un pico de sirem coincidente con el pico de magnitud
  for(i in 1:nUpeaks)
  {
    uPeaksLow =siremPeaks$siremPeaks$unitedMagnitudePeaks[i, 1]+1; #pico inferior unido
    uPeaksHigh=siremPeaks$siremPeaks$unitedMagnitudePeaks[i, 2]+1; #pico superior unido
    iMzPkLow =siremPeaks$siremPeaks$magnitudePeaks[uPeaksLow, 1]+1;#índice a masa inferior
    iMzPkHigh=siremPeaks$siremPeaks$magnitudePeaks[uPeaksHigh, 3]+1;#índice a masa superior
    
    if(mass>=siremPeaks$massAxis[iMzPkLow] & mass<=siremPeaks$massAxis[iMzPkHigh])
      {
      magUPeakLowMass =siremPeaks$massAxis[iMzPkLow];
      magUPeakHighMass=siremPeaks$massAxis[iMzPkHigh];
      break;
      }
  }
  gaussList<-gaussInfo$gaussians[,1];
  gaussLogic<-gaussList>=magUPeakLowMass & gaussList<=magUPeakHighMass;
  
#  maxSiremUpeaks<-getUpeakSiremMassesFromUpeak(siremPeaks, uPeak);
  return(gaussList[gaussLogic]);
}

#' siremPeaksFilter
#' Retorna un array indicando si el pico compuesto incluye a picos de rMSI2
#' Sólo se consideran aquellos que contienen a un pico de rMSI
#' @param siremPeaks -> sirem peaks, from rGetSiremPeaks()
#' @param rMSIPeaks  -> rMSI2 peaks, from  rMSI2::LoadPeakMatrix() from rMSI2::processWizard()
#' @return array indicando si el pico compuesto incluye a picos de rMSI2
#' @export
#' 
siremPeaksFilter<-function(siremPeaks, rMSIPeaks)
{
  nUpeaks=length(siremPeaks$siremPeaks$unitedMagnitudePeaks[,1]);
  uPks=rep(FALSE, times=nUpeaks);
  
  for(up in 1:nUpeaks) #para cada pico compuesto
  {
    uPkLow =siremPeaks$siremPeaks$unitedMagnitudePeaks[up, 1]+1; #índice a pico simple inferior
    uPkHigh=siremPeaks$siremPeaks$unitedMagnitudePeaks[up, 2]+1; #índice a pico simple superior
    #$magnitudePeaks[uPkLow,  1]+1=índice a masa (scan) inferior del pico simple
    #$magnitudePeaks[uPkLow,  1]+1=índice a masa (scan) superior del pico simple
    uLowPksDa =scans2Daltons(siremPeaks$siremPeaks$magnitudePeaks[uPkLow,  1]+1, siremPeaks$massAxis); #mz inferior
    uHighPksDa=scans2Daltons(siremPeaks$siremPeaks$magnitudePeaks[uPkHigh, 3]+1, siremPeaks$massAxis); #mz superior

    lowNearMass =nearestValue(uLowPksDa,  rMSIPeaks$mass); #mz de rMSI2 más próxima a mz inferior del pico compuesto
    highNearMass=nearestValue(uHighPksDa, rMSIPeaks$mass); #mz de rMSI2 más próxima a mz superior del pico compuesto
    #si mz de rMSI2 está contenida en el pico compuesto, se anota al pico compuesto
    if((lowNearMass  >=uLowPksDa & lowNearMass  <=uHighPksDa) | 
       (highNearMass >=uLowPksDa & highNearMass <=uHighPksDa))
      {uPks[up]=TRUE;}
  }
  return(uPks);
}

#' siremPeaksFilterDeconv
#' Retorna una matriz con los centroides de rMSI2 contenidos en los picos compuestos deconvolucionados de rSIREM
#' 
#' @param siremPeaks -> sirem peaks, from rGetSiremPeaks()
#' @param rMSIPeaks  -> rMSI2 peaks, from rMSI2::LoadPeakMatrix() from rMSI2::processWizard()
#' @return una matriz con los centroides de rMSI2 contenidos en los picos compuestos de rSIREM
#' columna 1   referencia al pico compuesto de rSIREM
#' columna 2.. centroides de rMSI2
#' @export
#' 
#' 
rMSI2PeaksFilterDeconv<-function(siremPeaks, gaussInfo, rMSIPeaks)
{
  peaksDec<-peaksDeconvolved(siremPeaks, gaussInfo);
  nUpeaks=length(peaksDec$gaussMassMatrix[,1]);
  uPks=rep(FALSE, times=nUpeaks);
  nRMSIpeaks=0; maxLength=0;
  for(up in 1:nUpeaks) #para cada pico compuesto
  {
    uLowPksDa =peaksDec$gaussMassMatrix[up, 2];
    uHighPksDa=peaksDec$gaussMassMatrix[up, 3];
    lowNearMass =nearestValue(uLowPksDa,  rMSIPeaks$mass); #mz de rMSI2 más próxima a mz inferior del pico compuesto
    highNearMass=nearestValue(uHighPksDa, rMSIPeaks$mass); #mz de rMSI2 más próxima a mz superior del pico compuesto
    #si mz de rMSI2 está contenida en el pico compuesto, se anota al pico compuesto
    if((lowNearMass  >=uLowPksDa & lowNearMass  <=uHighPksDa) | 
       (highNearMass >=uLowPksDa & highNearMass <=uHighPksDa))
    {
      uPks[up]=TRUE;
      logic=rMSIPeaks$mass>=lowNearMass & rMSIPeaks$mass<=highNearMass;
      rMSImz=rMSIPeaks$mass[logic];
      if(length(rMSImz)>maxLength) maxLength=length(rMSImz);
      nRMSIpeaks=nRMSIpeaks+1;
    }
  }
  matrixRMSI<-matrix(nrow=nRMSIpeaks, ncol=maxLength+2);
  index=1;
  for(up in 1:nUpeaks) #para cada pico compuesto
  {
    if(uPks[up]!=TRUE) next;
    uLowPksDa =peaksDec$gaussMassMatrix[up, 2];
    uHighPksDa=peaksDec$gaussMassMatrix[up, 3];
    
    lowNearMass =nearestValue(uLowPksDa,  rMSIPeaks$mass); #mz de rMSI2 más próxima a mz inferior del pico compuesto
    highNearMass=nearestValue(uHighPksDa, rMSIPeaks$mass); #mz de rMSI2 más próxima a mz superior del pico compuesto
    
    #Se determinan los picos de rMSI2 dentro del pico compuesto de rSIREM
    logic=rMSIPeaks$mass>=lowNearMass & rMSIPeaks$mass<=highNearMass;
    rMSImz=rMSIPeaks$mass[logic];
    matrixRMSI[index, 1]=peaksDec$gaussMassMatrix[up, 1];
    matrixRMSI[index, 2]=length(rMSImz);
    for(i in 1:length(rMSImz)) 
      matrixRMSI[index, i+2]=rMSImz[i];
    index=index+1;
  }
  return(matrixRMSI);
}

#' rMSI2PeaksFilter
#' Retorna una matriz con los centroides de rMSI2 contenidos en los picos compuestos de rSIREM
#' 
#' @param siremPeaks -> sirem peaks, from rGetSiremPeaks()
#' @param rMSIPeaks  -> rMSI2 peaks, from rMSI2::LoadPeakMatrix() from rMSI2::processWizard()
#' @return una matriz con los centroides de rMSI2 contenidos en los picos compuestos de rSIREM
#' columna 1   referencia al pico compuesto de rSIREM
#' columna 2.. centroides de rMSI2
#' @export
#' 
#' 
rMSI2PeaksFilter<-function(siremPeaks, rMSIPeaks)
{
  nUpeaks=length(siremPeaks$siremPeaks$unitedMagnitudePeaks[,1]);
  uPks=rep(FALSE, times=nUpeaks);
  nRMSIpeaks=0; maxLength=0;
  for(up in 1:nUpeaks) #para cada pico compuesto
  {
    uPkLow =siremPeaks$siremPeaks$unitedMagnitudePeaks[up, 1]+1; #índice a pico simple inferior
    uPkHigh=siremPeaks$siremPeaks$unitedMagnitudePeaks[up, 2]+1; #índice a pico simple superior
    #$magnitudePeaks[uPkLow,  1]+1=índice a masa (scan) inferior del pico simple
    #$magnitudePeaks[uPkLow,  1]+1=índice a masa (scan) superior del pico simple
    uLowPksDa =scans2Daltons(siremPeaks$siremPeaks$magnitudePeaks[uPkLow,  1]+1, siremPeaks$massAxis); #mz inferior
    uHighPksDa=scans2Daltons(siremPeaks$siremPeaks$magnitudePeaks[uPkHigh, 3]+1, siremPeaks$massAxis); #mz superior
    
    lowNearMass =nearestValue(uLowPksDa,  rMSIPeaks$mass); #mz de rMSI2 más próxima a mz inferior del pico compuesto
    highNearMass=nearestValue(uHighPksDa, rMSIPeaks$mass); #mz de rMSI2 más próxima a mz superior del pico compuesto
    #si mz de rMSI2 está contenida en el pico compuesto, se anota al pico compuesto
    if((lowNearMass  >=uLowPksDa & lowNearMass  <=uHighPksDa) | 
       (highNearMass >=uLowPksDa & highNearMass <=uHighPksDa))
      {
      uPks[up]=TRUE;
      logic=rMSIPeaks$mass>=lowNearMass & rMSIPeaks$mass<=highNearMass;
      rMSImz=rMSIPeaks$mass[logic];
      if(length(rMSImz)>maxLength) maxLength=length(rMSImz);
      nRMSIpeaks=nRMSIpeaks+1;
      }
  }
  matrixRMSI<-matrix(nrow=nRMSIpeaks, ncol=maxLength+2);
  index=1;
  for(up in 1:nUpeaks) #para cada pico compuesto
  {
    if(uPks[up]!=TRUE) next; #este pico compuesto de rSIREM no contiene picos de rMSI2
    
    uPkLow =siremPeaks$siremPeaks$unitedMagnitudePeaks[up, 1]+1; #índice a pico simple inferior
    uPkHigh=siremPeaks$siremPeaks$unitedMagnitudePeaks[up, 2]+1; #índice a pico simple superior
    #$magnitudePeaks[uPkLow,  1]+1=índice a masa (scan) inferior del pico simple
    #$magnitudePeaks[uPkLow,  1]+1=índice a masa (scan) superior del pico simple
    uLowPksDa =scans2Daltons(siremPeaks$siremPeaks$magnitudePeaks[uPkLow,  1]+1, siremPeaks$massAxis); #mz inferior
    uHighPksDa=scans2Daltons(siremPeaks$siremPeaks$magnitudePeaks[uPkHigh, 3]+1, siremPeaks$massAxis); #mz superior
    
    lowNearMass =nearestValue(uLowPksDa,  rMSIPeaks$mass); #mz de rMSI2 más próxima a mz inferior del pico compuesto
    highNearMass=nearestValue(uHighPksDa, rMSIPeaks$mass); #mz de rMSI2 más próxima a mz superior del pico compuesto
    
    #Se determinan los picos de rMSI2 dentro del pico compuesto de rSIREM
    logic=rMSIPeaks$mass>=lowNearMass & rMSIPeaks$mass<=highNearMass;
    rMSImz=rMSIPeaks$mass[logic];
    matrixRMSI[index, 1]=up;
    matrixRMSI[index, 2]=length(rMSImz);
    for(i in 1:length(rMSImz)) 
      matrixRMSI[index, i+2]=rMSImz[i];
    index=index+1;
  }
  return(matrixRMSI);
}
  

#' scans2Daltons
#' convierte un array de scans a Daltons
#' @param scans -> array de valores correspondientes a m/z en massAxis; scans[1]=0...
#'    ejemplo: scans=5.5 -> massAxis[5] + (massAxis[6]-massAxis[5])/2
#' @param massAxis -> eje de masas
#' @return array en Daltons correspondientes al array en scans
#' 
scans2Daltons<-function(scans, massAxis)
{
  ret=rep(0, times=length(scans));
  for(scan in 1:length(scans))
  {
  tmp=floor(scans[scan]);
  #interval in Daltons between two consecutive neighboring scans of mean.
  if(tmp+1 < length(massAxis)) #if it is within range.
    {delta=massAxis[tmp+1]-massAxis[tmp];} #posterior delta.
  else
    {delta=massAxis[tmp]-massAxis[tmp-1];} #previous delta.
  
  offset=delta*(scans[scan]-tmp); #displacement with respect to the origin of the compound peak.
  ret[scan]=massAxis[tmp]+offset; #converted value.
  }
  return(ret); 
}

#' sirem_vs_rMSI2
#' genera datos estadísticos con las desviaciones de ambos algoritmos respecto a una muestra patrón
#' 
#' @param sample -> sample. Valids: "C30k" y "C60k"
#' @param SNR    -> signal to noise ratio. Valids: 1, 2, 3, 5, 7
#' @param histo  -> TRUE si se quiere visualizar los histogramas
#'
#' @return nothing
#'
#' @export
#' 
sirem_vs_rMSI2<-function(sample, SNR, histo=FALSE)
{
  library(rSirem)

  if(SNR!=1 & SNR!=2 & SNR!=3 & SNR!=5 & SNR!=7) {print("Warning: unknown SNR; expected:1,2,3,5,7");  return();}
  load("/home/esteban/MALDI/rSirem_local/C120_all.RData")
  if(sample=="C60k")
    {
    #carga ficheros con los rangos de masas para la muestra de 60k
    load("/home/esteban/MALDI/rSirem_local/C60_all.RData")
    #carga info sobre la muestra
    myData<-rMSI2::LoadMsiData("/mnt2/MALDI/Cerebellum_30_60_120k/NoAlineado/C_60k/231211_Au_P_MBr_cblm_60k.imzML");
    #carga la matriz de picos generada por rMSI2
    peakMatrixPath<-sprintf("/mnt2/MALDI/Cerebellum_30_60_120k/Alineado/C_60k_SNR%d", SNR)
    rMSI2_snr <- rMSI2::LoadPeakMatrix(file.path(peakMatrixPath, "merged-peakmatrix.pkmat"))
    }
  else if(sample=="C30k")
    {
    load("/home/esteban/MALDI/rSirem_local/C30_all.RData")
    myData<-rMSI2::LoadMsiData("/mnt2/MALDI/Cerebellum_30_60_120k/NoAlineado/C_30k/231211_Au_P_MBr_cblm_30k.imzML");
    peakMatrixPath<-sprintf("/mnt2/MALDI/Cerebellum_30_60_120k/Alineado/C_30k_SNR%d", SNR)
    rMSI2_snr <- rMSI2::LoadPeakMatrix(file.path(peakMatrixPath, "merged-peakmatrix.pkmat"))
    }
  else {print("Warning: unknown sample");  return();}
  sprintf("", sample)
  
  #para el rango de masa de 700 a 900 Da
  if(sample=="C30k")
    {siremPeaks<-siremPeaks30_700_900n10s0;} #fichero de picos de sirem
  else if(sample=="C60k")
    {siremPeaks<-siremPeaks60_700_900n10s0;}
  
  #centroides de rMSI2 contenidos en picos compuestos deconvolucionados de rSIREM
  gaussInfo<-rGetGaussians(myData, siremPeaks, 0, 0.1); #gausianas (luego se alteran)
  rMSI2_rSIREM_700_900_snr <-rMSI2PeaksFilterDeconv(siremPeaks, gaussInfo, rMSI2_snr)
  
  #Se extraen los picos unidos de sirem que incluyen a algún centroide de rMSI2
  #si un pico compuesto no contiene a ningún pico de rMSI2 se marca e interpreta como ruidoso
  goodUpeaks=siremPeaksFilter(siremPeaks, rMSI2_snr); #picos compuestos que contienen picos de rMSI2
  #Se marcan los picos unidos que no incluyen a ningún centroide de rMSI2, para su descarte
  siremPeaks$siremPeaks$unitedMagnitudePeaks[,1][!goodUpeaks]=-1; #marcado
  #Se obtienen las gaussianas de los picos no descartados
  gaussInfo<-rGetGaussians(myData, siremPeaks, 0, 0.1); #no se consideran los picos marcados
  gaussInfoA=gaussInfo;
  
  #desviaciones de masa de cada pico sobre el patrón
  fq_700_900_rMSI_snr  <-fitQualityPere       (gaussInfo120_700_900n10ns, rMSI2_snr,  gaussInfo, histo)
  fq_700_900_rMSI_SIREM_snr  <-fitQualityPereSirem (gaussInfo120_700_900n10ns, rMSI2_snr,  gaussInfo, histo)
  fq_700_900_sirem_snr <-fitQualitySirem      (gaussInfo120_700_900n10ns, siremPeaks, gaussInfo, histo)
  fq_700_900_sirem_snrD<-fitQualitySiremDeconv(gaussInfo120_700_900n10ns, siremPeaks, gaussInfo, histo)
  #se eliminan de la matriz los picos deconvolucionados
  logic<-fq_700_900_sirem_snr[, 6]==0; 
  rIndex=1:length(logic);
  rIndex=rIndex[logic];
  #array de desviaciones
  fq_700_900_sirem_snrND=fq_700_900_sirem_snr[rIndex,];
  
  
  #para el rango de masa de 500 a 700 Da
  if(sample=="C30k")
    {siremPeaks<-siremPeaks30_500_700n10s0;}
  else if(sample=="C60k")
    {siremPeaks<-siremPeaks60_500_700n10s0;}
  
  #centroides de rMSI2 contenidos en picos compuestos deconvolucionados de rSIREM
  gaussInfo<-rGetGaussians(myData, siremPeaks, 0, 0.1);
  rMSI2_rSIREM_500_700_snr <-rMSI2PeaksFilterDeconv(siremPeaks, gaussInfo, rMSI2_snr)
  
  goodUpeaks=siremPeaksFilter(siremPeaks, rMSI2_snr); 
  siremPeaks$siremPeaks$unitedMagnitudePeaks[,1][!goodUpeaks]=-1; 
  gaussInfo<-rGetGaussians(myData, siremPeaks, 0, 0.1); 
  gaussInfoB=gaussInfo;
  
  fq_500_700_rMSI_snr  <-fitQualityPere       (gaussInfo120_500_700n10ns, rMSI2_snr,  gaussInfo, histo)
  fq_500_700_rMSI_SIREM_snr  <-fitQualityPereSirem (gaussInfo120_500_700n10ns, rMSI2_snr,  gaussInfo, histo)
  fq_500_700_sirem_snr <-fitQualitySirem      (gaussInfo120_500_700n10ns, siremPeaks, gaussInfo, histo)
  fq_500_700_sirem_snrD<-fitQualitySiremDeconv(gaussInfo120_500_700n10ns, siremPeaks, gaussInfo, histo)
  #se eliminan de la matriz los picos deconvolucionados
  logic<-fq_500_700_sirem_snr[, 6]==0; 
  rIndex=1:length(logic);
  rIndex=rIndex[logic];
  fq_500_700_sirem_snrND=fq_500_700_sirem_snr[rIndex,]; #no deconvolucionados
  
  #para el rango de masa de 300 a 500 Da
  if(sample=="C30k")
    {siremPeaks<-siremPeaks30_300_500n10s0;}
  else if(sample=="C60k")
    {siremPeaks<-siremPeaks60_300_500n10s0;}
  
  #centroides de rMSI2 contenidos en picos compuestos deconvolucionados de rSIREM
  gaussInfo<-rGetGaussians(myData, siremPeaks, 0, 0.1);
  rMSI2_rSIREM_300_500_snr <-rMSI2PeaksFilterDeconv(siremPeaks, gaussInfo, rMSI2_snr)
  
  goodUpeaks=siremPeaksFilter(siremPeaks, rMSI2_snr);
  siremPeaks$siremPeaks$unitedMagnitudePeaks[,1][!goodUpeaks]=-1;
  gaussInfo<-rGetGaussians(myData, siremPeaks, 0, 0.1);
  gaussInfoC=gaussInfo;
  
  fq_300_500_rMSI_snr  <-fitQualityPere       (gaussInfo120_300_500n10ns, rMSI2_snr,  gaussInfo, histo)
  fq_300_500_rMSI_SIREM_snr  <-fitQualityPereSirem (gaussInfo120_300_500n10ns, rMSI2_snr,  gaussInfo, histo)
  fq_300_500_sirem_snr <-fitQualitySirem      (gaussInfo120_300_500n10ns, siremPeaks, gaussInfo, histo)
  fq_300_500_sirem_snrD<-fitQualitySiremDeconv(gaussInfo120_300_500n10ns, siremPeaks, gaussInfo, histo)
  #se eliminan de la matriz los picos deconvolucionados
  logic<-fq_300_500_sirem_snr[, 6]==0; 
  rIndex=1:length(logic);
  rIndex=rIndex[logic];
  fq_300_500_sirem_snrND=fq_300_500_sirem_snr[rIndex,];
 
  #----------------------------
  # trabajos temporales
  #extrae info de los iones deconvolucionados
  #logicA<-fq_300_500_sirem_snr[, 6]!=0; 
  #D_300_500=fq_300_500_sirem_snr[logicA,1];
  #logicB<-fq_500_700_sirem_snr[, 6]!=0; 
  #D_500_700=fq_500_700_sirem_snr[logicB,1];
  #logicC<-fq_700_900_sirem_snr[, 6]!=0; 
  #D_700_900=fq_700_900_sirem_snr[logicC,1];
  #logicDeconv=c(logicA, logicB, logicC);
  #save (logicDeconv, file="/home/esteban/MALDI/rSirem_local/logicDeconv_30_300_900.RData")
  #----------------------------
  
   
  #Resultados para rMSI2
  print("About rMSI2 peaks:")
  print("  range 300-500 Da:")
  m=mean(fq_300_500_rMSI_snr[,3]); sigma=sd(fq_300_500_rMSI_snr[,3]); md=median(fq_300_500_rMSI_snr[,3]);
  txt=sprintf("    size=%d", length(fq_300_500_rMSI_snr[,1])); print(txt);
  txt=sprintf("    mean=%.5f sigma=%.5f median=%.5f", m, sigma, md); print(txt);
  print("  range 500-700 Da:")
  txt=sprintf("    size=%d", length(fq_500_700_rMSI_snr[,3])); print(txt);
  m=mean(fq_500_700_rMSI_snr[,3]); sigma=sd(fq_500_700_rMSI_snr[,3]); md=median(fq_500_700_rMSI_snr[,3]);
  txt=sprintf("    mean=%.5f sigma=%.5f median=%.5f", m, sigma, md); print(txt);
  print("  range 700-900 Da:")
  txt=sprintf("    size=%d", length(fq_700_900_rMSI_snr[,3])); print(txt);
  m=mean(fq_700_900_rMSI_snr[,3]); sigma=sd(fq_700_900_rMSI_snr[,3]); md=median(fq_700_900_rMSI_snr[,3]);
  txt=sprintf("    mean=%.5f sigma=%.5f median=%.5f", m, sigma, md); print(txt);
  
  #desviaciones para el rango de masas unificado (300:900 Da)  
  totalDiff_rMSI2<-c(fq_300_500_rMSI_snr[,3], fq_500_700_rMSI_snr[,3], fq_700_900_rMSI_snr[,3])
  #Histograma de las desviaciones
  if(histo==TRUE)
    {
    hist(totalDiff_rMSI2, main="Histogram of rMSI2 deviations", xlab="ppm", ylab="Frequency");
    legendTxt<-sprintf("mz=300:900 Da\nsample=%s; SNR=%d", sample, SNR)
    legend("topright", legend=legendTxt);
    }
  print("  range 300-900 Da:")
  txt=sprintf("    size=%d", length(totalDiff_rMSI2)); print(txt);
  m=mean(totalDiff_rMSI2); sigma=sd(totalDiff_rMSI2); md=median(totalDiff_rMSI2);
  txt=sprintf("    mean=%.5f sigma=%.5f median=%.5f", m, sigma, md); print(txt);
  print("");
  
  #Resultados para rMSI2-rSIREM (desviaciones de rMSI2 en ppm)
  #Se analizan picos de rMSI2 y de rSIREM que comparten picos de alta reesolución
  #desviaciones para el rango de masas unificado (300:900 Da) 
  #desviaciones de rMSI2
  totalDiff_rMSI2_SIREM_A<-c(fq_300_500_rMSI_SIREM_snr[,3], fq_500_700_rMSI_SIREM_snr[,3], fq_700_900_rMSI_SIREM_snr[,3])
  #Histograma de las desviaciones
  if(histo==TRUE)
    {
    hist(totalDiff_rMSI2_SIREM_A, main="Histogram of rMSI2 vs rMSI2 A deviations", xlab="ppm", ylab="Frequency");
    legendTxt<-sprintf("mz=300:900 Da\nsample=%s; SNR=%d", sample, SNR)
    legend("topright", legend=legendTxt);
    }
  print("About rMSI2_rSIREM-rMSI2 peaks:")
  print("  range 300-500 Da:")
  m=mean(fq_300_500_rMSI_SIREM_snr[,3]); sigma=sd(fq_300_500_rMSI_SIREM_snr[,3]); md=median(fq_300_500_rMSI_SIREM_snr[,3]);
  txt=sprintf("    size=%d", length(fq_300_500_rMSI_SIREM_snr[,1])); print(txt);
  txt=sprintf("    mean=%.5f sigma=%.5f median=%.5f", m, sigma, md); print(txt);
  print("  range 500-700 Da:")
  txt=sprintf("    size=%d", length(fq_500_700_rMSI_SIREM_snr[,3])); print(txt);
  m=mean(fq_500_700_rMSI_SIREM_snr[,3]); sigma=sd(fq_500_700_rMSI_SIREM_snr[,3]); md=median(fq_500_700_rMSI_SIREM_snr[,3]);
  txt=sprintf("    mean=%.5f sigma=%.5f median=%.5f", m, sigma, md); print(txt);
  print("  range 700-900 Da:")
  txt=sprintf("    size=%d", length(fq_700_900_rMSI_snr[,3])); print(txt);
  m=mean(fq_700_900_rMSI_SIREM_snr[,3]); sigma=sd(fq_700_900_rMSI_SIREM_snr[,3]); md=median(fq_700_900_rMSI_SIREM_snr[,3]);
  txt=sprintf("    mean=%.5f sigma=%.5f median=%.5f", m, sigma, md); print(txt);
  print("  range 300-900 Da:")
  txt=sprintf("    size=%d", length(totalDiff_rMSI2_SIREM_A)); print(txt);
  m=mean(totalDiff_rMSI2_SIREM_A); sigma=sd(totalDiff_rMSI2_SIREM_A); md=median(totalDiff_rMSI2_SIREM_A);
  txt=sprintf("    mean=%.5f sigma=%.5f median=%.5f", m, sigma, md); print(txt);
  print("");
  

  #Resultados para rMSI2-rSIREM (desviaciones de rSIREM en ppm)
  #Se analizan picos de rMSI2 y de rSIREM que comparten picos de alta resolución
  #desviaciones para el rango de masas unificado (300:900 Da)  
  #desviaciones de rSIREM
  totalDiff_rMSI2_SIREM_B<-c(fq_300_500_rMSI_SIREM_snr[,7], fq_500_700_rMSI_SIREM_snr[,7], fq_700_900_rMSI_SIREM_snr[,7])
  #Histograma de las desviaciones
  if(histo==TRUE)
    {
    hist(totalDiff_rMSI2_SIREM_B, main="Histogram of rMSI2_SIREM B deviations", xlab="ppm", ylab="Frequency");
    legendTxt<-sprintf("mz=300:900 Da\nsample=%s; SNR=%d", sample, SNR)
    legend("topright", legend=legendTxt);
    }
  
  print("About rMSI2_rSIREM-rSIREM peaks:")
  print("  range 300-500 Da:")
  m=mean(fq_300_500_rMSI_SIREM_snr[,7]); sigma=sd(fq_300_500_rMSI_SIREM_snr[,7]); md=median(fq_300_500_rMSI_SIREM_snr[,7]);
  txt=sprintf("    size=%d", length(fq_300_500_rMSI_SIREM_snr[,1])); print(txt);
  txt=sprintf("    mean=%.5f sigma=%.5f median=%.5f", m, sigma, md); print(txt);
  print("  range 500-700 Da:")
  txt=sprintf("    size=%d", length(fq_500_700_rMSI_SIREM_snr[,7])); print(txt);
  m=mean(fq_500_700_rMSI_SIREM_snr[,7]); sigma=sd(fq_500_700_rMSI_SIREM_snr[,7]); md=median(fq_500_700_rMSI_SIREM_snr[,7]);
  txt=sprintf("    mean=%.5f sigma=%.5f median=%.5f", m, sigma, md); print(txt);
  print("  range 700-900 Da:")
  txt=sprintf("    size=%d", length(fq_700_900_rMSI_SIREM_snr[,7])); print(txt);
  m=mean(fq_700_900_rMSI_SIREM_snr[,7]); sigma=sd(fq_700_900_rMSI_SIREM_snr[,7]); md=median(fq_700_900_rMSI_SIREM_snr[,7]);
  txt=sprintf("    mean=%.5f sigma=%.5f median=%.5f", m, sigma, md); print(txt);
  print("  range 300-900 Da:")
  txt=sprintf("    size=%d", length(totalDiff_rMSI2_SIREM_B)); print(txt);
  m=mean(totalDiff_rMSI2_SIREM_B); sigma=sd(totalDiff_rMSI2_SIREM_B); md=median(totalDiff_rMSI2_SIREM_B);
  txt=sprintf("    mean=%.5f sigma=%.5f median=%.5f", m, sigma, md); print(txt);
  print("");
  
  #Resultados para los picos de rSirem descartando los picos deconvolucionados
  totalDiff_siremNDA<-c(fq_300_500_sirem_snrND[,1], fq_500_700_sirem_snrND[,1], fq_700_900_sirem_snrND[,1])
  totalDiff_siremND<-c(fq_300_500_sirem_snrND[,3], fq_500_700_sirem_snrND[,3], fq_700_900_sirem_snrND[,3])
  logic=totalDiff_siremNDA>0 & totalDiff_siremND<50; #(<50)excluye deficiencias en la comparación
  totalDiff_siremND=totalDiff_siremND[logic];
  if(histo==TRUE)
    {
    hist(totalDiff_siremND, main="Histogram of rSirem deviations (not deconv)", xlab="ppm", ylab="Frequency");
    legend("topright", legend=legendTxt);
    }
  siremMean=mean(totalDiff_siremND)
  siremSd=sd(totalDiff_siremND)
  msg<-sprintf("total w/o deconv rSirem peaks=%5d; mean=%.4f; sigma=%.4f", length(totalDiff_siremND), siremMean, siremSd);
  print(msg)
  
  #Resultados para todos los picos de rSirem (incluye a los picos deconvolucionados) 
  totalDiff_siremA<-c(fq_300_500_sirem_snr[,1], fq_500_700_sirem_snr[,1], fq_700_900_sirem_snr[,1])
  totalDiff_sirem<-c(fq_300_500_sirem_snr[,3], fq_500_700_sirem_snr[,3], fq_700_900_sirem_snr[,3])
  logic=totalDiff_siremA>0 & totalDiff_sirem<50; #(<50)excluye deficiencias en la comparación
  totalDiff_sirem=totalDiff_sirem[logic];
  
#  totalDiff_siremA=totalDiff_siremA[logic];
#  save(totalDiff_siremA, ascii=TRUE, file=("~/MALDI/rSirem/rSIREM_30k_snr3.txt"))
  
  if(histo==TRUE)
    {
    hist(totalDiff_sirem, main="Histogram of rSirem deviations", xlab="ppm", ylab="Frequency");
    legend("topright", legend=legendTxt);
    }
  siremMean=mean(totalDiff_sirem)
  siremSd=sd(totalDiff_sirem)
  msg<-sprintf("           total rSirem peaks=%5d; mean=%.4f; sigma=%.4f", length(totalDiff_sirem), siremMean, siremSd);
  print(msg)
  
  #Resultados para los picos de rSirem considerando solo los picos deconvolucionados
  totalDiff_siremD<-c(fq_300_500_sirem_snrD[,3], fq_500_700_sirem_snrD[,3], fq_700_900_sirem_snrD[,3])
  if(histo==TRUE)
    {
    hist(totalDiff_siremD, main="Histogram of rSirem deconvolution deviations", xlab="ppm", ylab="Frequency");
    legend("topright", legend=legendTxt);
    }
  siremMean=mean(totalDiff_siremD)
  siremSd=sd(totalDiff_siremD)
  msg<-sprintf("    deconvoluted rSirem peaks=%5d; mean=%.4f; sigma=%.4f", length(totalDiff_siremD), siremMean, siremSd);
  print(msg)
  print("");
  
  print("About rSIREM peaks:")
  print("  range 300-500 Da:")
  logic      =fq_300_500_rMSI_snr[,1]>0; #mz
  mzDiffMSI  =fq_300_500_rMSI_snr[,1][logic];
  logic      =fq_300_500_sirem_snr[,1]>0;
  mzDiffSIREM=fq_300_500_sirem_snr[,1][logic];
  devi       =fq_300_500_sirem_snr[,3][logic]; #deviations
  neIndex= notEqual(mzDiffMSI, mzDiffSIREM, 15); #index to deconvolved peaks
  m=mean(devi); sigma=sd(devi); md=median(devi);
  txt=sprintf("         total: mean=%.5f sigma=%.5f median=%.5f", m, sigma, md);
  print(txt);
  neDevi=devi[neIndex]; #deviations of deconvolved peaks
  m=mean(neDevi); sigma=sd(neDevi); md=median(neDevi);
  txt=sprintf("        deconv: mean=%.5f sigma=%.5f median=%.5f", m, sigma, md);
  print(txt);
  eDevi=devi[!neIndex]; #deviations of not deconvolved peaks
  m=mean(eDevi); sigma=sd(eDevi); md=median(eDevi);
  txt=sprintf("    not deconv: mean=%.5f sigma=%.5f median=%.5f", m, sigma, md);
  print(txt);

  print("  range 500-700 Da:")
  logic      =fq_500_700_rMSI_snr[,1]>0;
  mzDiffMSI  =fq_500_700_rMSI_snr[,1][logic];
  logic      =fq_500_700_sirem_snr[,1]>0;
  mzDiffSIREM=fq_500_700_sirem_snr[,1][logic];
  devi       =fq_500_700_sirem_snr[,3][logic]; #deviations
  neIndex= notEqual(mzDiffMSI, mzDiffSIREM, 15); #index to deconvolved peaks
  m=mean(devi); sigma=sd(devi); md=median(devi);
  txt=sprintf("         total: mean=%.5f sigma=%.5f median=%.5f", m, sigma, md);
  print(txt);
  neDevi=devi[neIndex]; #deviations of deconvolved peaks
  m=mean(neDevi); sigma=sd(neDevi); md=median(neDevi);
  txt=sprintf("        deconv: mean=%.5f sigma=%.5f median=%.5f", m, sigma, md);
  print(txt);
  eDevi=devi[!neIndex]; #deviations of not deconvolved peaks
  m=mean(eDevi); sigma=sd(eDevi); md=median(eDevi);
  txt=sprintf("    not deconv: mean=%.5f sigma=%.5f median=%.5f", m, sigma, md);
  print(txt);
  
  print("  range 700-900 Da:")
  logic      =fq_700_900_rMSI_snr[,1]>0;
  mzDiffMSI  =fq_700_900_rMSI_snr[,1][logic];
  logic      =fq_700_900_sirem_snr[,1]>0;
  mzDiffSIREM=fq_700_900_sirem_snr[,1][logic];
  devi       =fq_700_900_sirem_snr[,3][logic]; #deviations
  neIndex= notEqual(mzDiffMSI, mzDiffSIREM, 15); #index to deconvolved peaks
  m=mean(devi); sigma=sd(devi); md=median(devi);
  txt=sprintf("         total: mean=%.5f sigma=%.5f median=%.5f", m, sigma, md);
  print(txt);
  neDevi=devi[neIndex]; #deviations of deconvolved peaks
  m=mean(neDevi); sigma=sd(neDevi); md=median(neDevi);
  txt=sprintf("        deconv: mean=%.5f sigma=%.5f median=%.5f", m, sigma, md);
  print(txt);
  eDevi=devi[!neIndex]; #deviations of not deconvolved peaks
  m=mean(eDevi); sigma=sd(eDevi); md=median(eDevi);
  txt=sprintf("    not deconv: mean=%.5f sigma=%.5f median=%.5f", m, sigma, md);
  print(txt);
  
  print("  range 300-900 Da")
  mzDiffMSI<-c(fq_300_500_rMSI_snr[,1], fq_500_700_rMSI_snr[,1], fq_700_900_rMSI_snr[,1]);
  logic=mzDiffMSI>0;
  mzDiffMSI=mzDiffMSI[logic];
  mzDiffSIREM<-c(fq_300_500_sirem_snr[,1], fq_500_700_sirem_snr[,1], fq_700_900_sirem_snr[,1])
  logic=mzDiffSIREM>0;
  mzDiffSIREM=mzDiffSIREM[logic];
  devi=c(fq_300_500_sirem_snr[,3], fq_500_700_sirem_snr[,3], fq_700_900_sirem_snr[,3]) #deviations
  devi=devi[logic];
  neIndex= notEqual(mzDiffMSI, mzDiffSIREM, 15); #index to deconvolved peaks
  m=mean(devi); sigma=sd(devi); md=median(devi);
  txt=sprintf("         total: mean=%.5f sigma=%.5f median=%.5f", m, sigma, md);
  print(txt);
  neDevi=devi[neIndex]; #deviations of deconvolved peaks
  m=mean(neDevi); sigma=sd(neDevi); md=median(neDevi);
  txt=sprintf("        deconv: mean=%.5f sigma=%.5f median=%.5f", m, sigma, md);
  print(txt);
  eDevi=devi[!neIndex]; #deviations of not deconvolved peaks
  m=mean(eDevi); sigma=sd(eDevi); md=median(eDevi);
  txt=sprintf("    not deconv: mean=%.5f sigma=%.5f median=%.5f", m, sigma, md);
  print(txt);
  #  print("deconvolved peaks:");
#  print(mzDiffSIREM[neIndex])
#  return(0);  
  
  print("list of rSirem peaks with deviation > x ppm:")
#  totalWrong_mz4<-c(fq_300_500_sirem_snr[,4], fq_500_700_sirem_snr[,4], fq_700_900_sirem_snr[,4])
  totalWrong_mz3<-c(fq_300_500_sirem_snr[,3], fq_500_700_sirem_snr[,3], fq_700_900_sirem_snr[,3])
  totalWrong_mz5<-c(fq_300_500_sirem_snr[,5], fq_500_700_sirem_snr[,5], fq_700_900_sirem_snr[,5])
  logic=totalWrong_mz3 > 5 | totalWrong_mz5 > 0; #0: <= 1 scans; 1:>1 && <=1.5 scans; 2; >1.5 scans
  totalWrong_mz<-c(fq_300_500_sirem_snr[,1], fq_500_700_sirem_snr[,1], fq_700_900_sirem_snr[,1])
  totalWrong_mz<-totalWrong_mz[logic];
  print(totalWrong_mz);

  print("list of rMSI2 peaks with deviation > x ppm:")
#  totalWrong_mz4<-c(fq_300_500_rMSI_snr[,4], fq_500_700_rMSI_snr[,4], fq_700_900_rMSI_snr[,4])
  totalWrong_mz3<-c(fq_300_500_rMSI_snr[,3], fq_500_700_rMSI_snr[,3], fq_700_900_rMSI_snr[,3])
  totalWrong_mz5<-c(fq_300_500_rMSI_snr[,5], fq_500_700_rMSI_snr[,5], fq_700_900_rMSI_snr[,5])
  logic=totalWrong_mz3 > 5 | totalWrong_mz5 > 0; #0: <= 1 scans; 1:>1 && <=1.5 scans; 2; >1.5 scans
  totalWrong_mz<-c(fq_300_500_rMSI_snr[,1], fq_500_700_rMSI_snr[,1], fq_700_900_rMSI_snr[,1])
  totalWrong_mz<-totalWrong_mz[logic];
  print(totalWrong_mz);
  #  return (totalDiff_siremD);
}

#compara dos arrays de longitud posiblemente distinta y retorna un array lógico donde son ciertos los elementos del segundo argumento
#que no están incluidos en el primero, aceptando como iguales aquellos con una diferencia <= al tercer argumento (en ppm)
notEqual<-function(totalDiff_rMSI2, totalDiff_sirem, ppm)
{
  sizeMSI  =length(totalDiff_rMSI2);
  sizeSirem=length(totalDiff_sirem);
  equalSirem=rep(FALSE, sizeSirem);
  isolated=rep(0, sizeMSI);
  equal=rep(0, sizeMSI);
  n=0; k=0;
  
  #búsqueda no optimizada
  for(i in 1:sizeMSI)
  {
    hit=FALSE;
    for(j in 1:sizeSirem)
    {
      diff=abs(totalDiff_rMSI2[i]-totalDiff_sirem[j]);
      diff_ppm=diff*1e6/totalDiff_rMSI2[i];
      if(diff_ppm<=ppm)
        {equal[k]=totalDiff_sirem[j]; hit=TRUE; k=k+1; equalSirem[j]=TRUE; break;}
    }
    if(hit==FALSE)
      {isolated[n]=totalDiff_rMSI2[i]; n=n+1;}
  }
  msg=sprintf("    sizes of rMSI2=%4d; SIREM=%4d, matches=%4d rMSI2_isolated=%4d rSIREM_isolated=%4d", sizeMSI, sizeSirem, k, n, sizeSirem-k)
  print(msg);
#  diffSirem=totalDiff_sirem[!equalSirem];
  return (!equalSirem);
}

#' histo2
#' Histograma adaptado para visualizar los hits en base de datos de aduptos de Jorde Capellades
#' bin size=1; 0>= x <=15; x.column=4
#' 
#' @param x       -> datos (read.csv("Jordi_Capellades_file"))
#' @param title   -> título principal
#'
#' @return nothing
#'
#' @export
#' 
histoJC<-function(x, title)
  {
  logic=x[,4]<=15; 
  data=x[,4][logic]; 
  bins=seq(0,15,1); 
  hist(data, breaks=bins, freq=FALSE, main=title);
  
}

#' violinSirem()
#' imagen violin para las desviaciones de rSIREM y rMSI2
#' @param type: SIREM, SIREM_d, SIREM_nd, MSI2, MSI2_MSI2, MSI2_SIREM
#' SIREM    -> para todos los centroides de rSIREM
#' SIREM_d  -> para los controides exclusivos de rSIREM (no existentes en rMSI2)
#' SIREM_nd -> para los centroides compartidos entre rSIREM y rMSI2
#' MSI2     -> para todos los centroides de rMSI2
#' rMSI2_MSI2  -> rMSI2 sobre los centroides comunes entre rMSI2 y rSIREM
#' rMSI2_SIREM -> SIREM sobre los centroides comunes entre rMSI2 y rSIREM
#' @return a list with two statistical data matrices.
#' @export
#' 
violinSirem<-function(type="SIREM")
{
  if(type=="SIREM")
    base="SIREM"
  else if(type=="SIREM_d")
    base="SIREMD"
  else if(type=="SIREM_nd")
    base="SIREMND"
  else if(type=="MSI2")
    base="MSI2"
  else if(type=="MSI2_MSI2")
    base="rMSI2_SIREM_MSI2"
  else if(type=="MSI2_SIREM")
    base="rMSI2_SIREM_SIREM"
  else
    {
    txt=sprintf("unknow type %s. Valids types are: SIREM, SIREM_d, SIREM_nd, MSI2, MSI2_MSI2, MSI2_SIREM", type);
    print(txt);
    return();
    }
  txt=sprintf("Violin test for %s data", base); print(txt);
  
#  library(rSirem)
  library(vioplot)
  
  #carga de datos desde fichero
  txt=sprintf("/home/esteban/MALDI/rSirem_local/%s_30_60.RData", base)
  load(txt);
  
  result30 <-matrix(nrow = 4, ncol = 4); #para alojar resultados de 30k
  #para 30k y SNR=1
  txt=sprintf("%s_30_1", base); var=eval(parse(text = txt));
  deviation=var;
  result30[1,1]=1;
  result30[1,2]=mean(var);
  result30[1,3]=sd(var);
  result30[1,4]=median(var);
  size30_1=length(var)
  
  #para 30k y SNR=3
  txt=sprintf("%s_30_3", base); var=eval(parse(text = txt));
  var_30_3=var;
  deviation=c(deviation, var);
  result30[2,1]=3;
  result30[2,2]=mean(var);
  result30[2,3]=sd(var);
  result30[2,4]=median(var);
  size30_3=length(var)
  
  #para 30k y SNR=5
  txt=sprintf("%s_30_5", base); var=eval(parse(text = txt));
  deviation=c(deviation, var);
  result30[3,1]=5;
  result30[3,2]=mean(var);
  result30[3,3]=sd(var);
  result30[3,4]=median(var);
  size30_5=length(var)
  
  #para 30k y SNR=7
  txt=sprintf("%s_30_7", base); var=eval(parse(text = txt));
  deviation=c(deviation, var);
  result30[4,1]=7;
  result30[4,2]=mean(var);
  result30[4,3]=sd(var);
  result30[4,4]=median(var);
  size30_7=length(var)
  colnames(result30)<-c("SNR", "mean", "sigma", "median");
  
  txt=sprintf("#peaks for 30k: SNR1=%4d SNR3=%4d SNR5=%4d SNR7=%4d", size30_1, size30_3, size30_5, size30_7)
  print(txt);
  SNR30_1=rep("30k_1", size30_1)
  SNR30_3=rep("30k_3", size30_3)
  SNR30_5=rep("30k_5", size30_5)
  SNR30_7=rep("30k_7", size30_7)

  result60 <-matrix(nrow = 4, ncol = 4);
  #para 60k y SNR=1
  txt=sprintf("%s_60_1", base); var=eval(parse(text = txt));
  result60[1,1]=1;
  result60[1,2]=mean(var);
  result60[1,3]=sd(var);
  result60[1,4]=median(var);
  size60_1=length(var)
  deviation=c(deviation, var);
  
  #para 60k y SNR=3
  txt=sprintf("%s_60_3", base); var=eval(parse(text = txt));
  result60[2,1]=3;
  result60[2,2]=mean(var);
  result60[2,3]=sd(var);
  result60[2,4]=median(var);
  size60_3=length(var)
  deviation=c(deviation, var);
  var_60_3=var;
  
  #para 60k y SNR=5
  txt=sprintf("%s_60_5", base); var=eval(parse(text = txt));
  result60[3,1]=5;
  result60[3,2]=mean(var);
  result60[3,3]=sd(var);
  result60[3,4]=median(var);
  size60_5=length(var)
  deviation=c(deviation, var);
  
  #para 60k y SNR=7
  txt=sprintf("%s_60_7", base); var=eval(parse(text = txt));
  result60[4,1]=7;
  result60[4,2]=mean(var);
  result60[4,3]=sd(var);
  result60[4,4]=median(var);
  size60_7=length(var)
  deviation=c(deviation, var);
  colnames(result60)<-c("SNR", "mean", "sigma", "median");
  
  txt=sprintf("#peaks for 60k: SNR1=%4d SNR3=%4d SNR5=%4d SNR7=%4d", size60_1, size60_3, size60_5, size60_7)
  print(txt);
  #para poder visualizar varios violines juntos
  SNR60_1=rep("60k_1", size60_1)
  SNR60_3=rep("60k_3", size60_3)
  SNR60_5=rep("60k_5", size60_5)
  SNR60_7=rep("60k_7", size60_7)
  
  #Se visualiza un violin para cada SNR en 30k y 60k
  resolution_SNR=c(SNR30_1, SNR30_3, SNR30_5, SNR30_7, SNR60_1, SNR60_3, SNR60_5, SNR60_7);
  vioplot(deviation ~ resolution_SNR, col = 2:9, xlab="massResolution_SNR", ylab="deviation (ppm)")
  
  #Se visualiza un violin para SNR_3 en 30k y 60k
  deviation2=c(var_30_3, var_60_3);
  SNR30b_3=rep("30k", length(var_30_3))
  SNR60b_3=rep("60k", length(var_60_3))
  resolution2_SNR=c(SNR30b_3, SNR60b_3);
  vioplot(deviation2 ~ resolution2_SNR, col = c(3,7), xlab="massResolution (SNR=3)", ylab="deviation (ppm)")

  print("30k"); print(result30); 
  print("60k"); print(result60);
#  ret=list(var_30=var_30_3, var_60=var_60_3)
  ret=list(matrix_30k=result30, matrix_60k=result60);
  return(ret);
  
#  NOTA: para representar la figura SI_Fig_8
# activar: ret=list(var_30=var_30_3, var_60=var_60_3) y resactivar el ret actual
# desde consola
#  var_MSI2=violinSirem(type="MSI2_MSI2")
#  var_MSI2_S=violinSirem(type="MSI2_SIREM")
  # deviation2=c(var_MSI2$var_30, var_MSI2_S$var_30, var_MSI2$var_60, var_MSI2_S$var_60);
  # SNR30_3=rep("30k_rMSI2", length(var_MSI2$var_30)); SNR30b_3=rep("30k_rSIREM", length(var_MSI2_S$var_30));
  # SNR60_3=rep("60k_rMSI2", length(var_MSI2$var_60)); SNR60b_3=rep("60k_rSIREM", length(var_MSI2_S$var_60));
  # resolution2_SNR=c(SNR30_3, SNR30b_3, SNR60_3, SNR60b_3);
  # vioplot(deviation2 ~ resolution2_SNR, col = c(3,7), xlab="massResolution (SNR=3)", ylab="deviation (ppm)")
  
}

#' annotations()
#' Resultados sobre la cantidad de anotaciones en bases de datos.
#' Compara los resultados separando las masas extras de las comunes en rSIREM.
#' son masas comunes las que coinciden o están muy próximas entre rMSI2 y rSIREM
#' son masas extras las no comunes de rSIREM.
#' Presenta 4 violines: dos para masas comunes y dos para masas extras en 30k y 60k

#' @argument ppm -> se consideran mz comunes las que están incluidas en esta desviación
#' @export
#' 
annotations<-function(ppm=10)
{
  library(rSirem)
  library(vioplot)
  #ficheros con resultados
  rMSI_30_3   <- read.csv("/home/esteban/MALDI/Paper/anotaciones/rMSI2_30k_snr3_peaks_matches.csv",  header=TRUE, stringsAsFactors=FALSE);
  rMSI_60_3   <- read.csv("/home/esteban/MALDI/Paper/anotaciones/rMSI2_60k_snr3_peaks_matches.csv",  header=TRUE, stringsAsFactors=FALSE);
  rSIREM_30_3 <- read.csv("/home/esteban/MALDI/Paper/anotaciones/rSIREM_30k_snr3_peaks_matches.csv", header=TRUE, stringsAsFactors=FALSE);
  rSIREM_60_3 <- read.csv("/home/esteban/MALDI/Paper/anotaciones/rSIREM_60k_snr3_peaks_matches.csv", header=TRUE, stringsAsFactors=FALSE);
  
  neIndex_30= notEqual(rMSI_30_3[,3], rSIREM_30_3[,3], ppm); #index to extra peaks 30k SNR=3
  neIndex_60= notEqual(rMSI_60_3[,3], rSIREM_60_3[,3], ppm); #index to extra peaks 60k SNR=3
  
  #separación en 4 arrays
  extraSIREM_30=rSIREM_30_3[,4][ neIndex_30]; #extra 30K
  extraSIREM_60=rSIREM_60_3[,4][ neIndex_60]; #extra 60k
  baseSIREM_30 =rSIREM_30_3[,4][!neIndex_30]; #common 30k
  baseSIREM_60 =rSIREM_60_3[,4][!neIndex_60]; #common 60k
  
  #presentación de 4 violines
  B_30=rep("Common_30k", length(baseSIREM_30));
  B_60=rep("Common_60k", length(baseSIREM_60));
  E_30=rep("Extra_30k", length(extraSIREM_30));
  E_60=rep("Extra_60k", length(extraSIREM_60));
  
  vio=c(B_30, B_60, E_30, E_60); #los violines respetan el orden alfabético
  annoted=c(baseSIREM_30, baseSIREM_60, extraSIREM_30, extraSIREM_60);
  
  vioplot(annoted ~ vio, col = c(2:5), xlab="", ylab="#annotations")
#  legend("topright", legend="C=commons mz; E=extra mz");
  
  #datos estadísticos
  mB30=mean(baseSIREM_30); mB60=mean(baseSIREM_60); mE30=mean(extraSIREM_30); mE60=mean(extraSIREM_60);
  sdB30=sd(baseSIREM_30); sdB60=sd(baseSIREM_60); sdE30=sd(extraSIREM_30); sdE60=sd(extraSIREM_60);
  mdB30=median(baseSIREM_30); mdB60=median(baseSIREM_60); mdE30=median(extraSIREM_30); mdE60=median(extraSIREM_60);
  
  print("violins about numbers of annotations")
  txt=sprintf(" Common_30: mean=%.5f sigma=%.5f median=%.5f", mB30, sdB30, mdB30);
  print(txt);
  txt=sprintf(" Common_60: mean=%.5f sigma=%.5f median=%.5f", mB60, sdB60, mdB60);
  print(txt);
  txt=sprintf("  Extra_30: mean=%.5f sigma=%.5f median=%.5f", mE30, sdE30, mdE30);
  print(txt);
  txt=sprintf("  Extra_60: mean=%.5f sigma=%.5f median=%.5f", mE60, sdE60, mdE60);
  print(txt);
  
  }

#' rPeaksMatrixModif()
#' Modifies a matrix of peaks by making peaks of higher magnitude
#' or less than a cut-off value, become null
#' @param peaksMatrix -> matrix[pixels, ions]
#' @param cutLevel -> cutting intensity.
#' @param up -> TRUE  if peaks exceeding cutLevel are nulled.
#'           -> FALSE if peaks that do not exceed cutLevel are cancelled.
#'
#' @return -> the new matrix [pixels, ions]
#' @export
#' 
rPeaksMatrixModif<-function(peaksMatrix, cutLevel, up=TRUE)
{
  nPixels=nrow(peaksMatrix);
  nIons=ncol(peaksMatrix);
  newMatrix=matrix(nrow=nPixels, ncol=nIons);
  if(up==TRUE)
  {
    for(col in 1:nIons)
    {
      logic=peaksMatrix[,col]>cutLevel;
      newMatrix[,col]=peaksMatrix[,col];
      newMatrix[logic, col]=0;
    }
  }
  else
  {
    for(col in 1:nIons)
    {
      logic=peaksMatrix[,col]<cutLevel;
      newMatrix[,col]=peaksMatrix[,col];
      newMatrix[logic, col]=0;
    }
  }
  
  return (newMatrix);
}

#' rPeaksMatrixRandom()
#' Altera ciertas partes de una matriz
#' @param peaksMatrix matrix original
#' @param logic vector lógico. TRUE si la columna se ha de alterar
#' @param mode="RAMDOM" si la afectación es aleatoria para la columna TRUE.
#'        mode="NULL" si se hacen nulas las columnas TRUE
#' @param modify    porcentaje de columnas afectadas entre las que son TRUE 
#' @param magnitude porcentaje máximo de cambio sobre el valor promediado, descartando los valore no nulos
#' @return la matrix modificada
#' @export
#'  
rPeaksMatrixRandom<-function(peaksMatrix, logic, mode="RANDOM", modify=10, magnitude=10)
{
  nPixels=nrow(peaksMatrix); #filas de la matriz
  nCols=ncol(peaksMatrix); #columnas de la matriz
  countD=0; 
  
  if(nCols!= length(logic))
  {print("Warning: los tamaños de logic y columnas de peaksMatrix no coinciden"); return();}
  for(i in 1:nCols)
    if(logic[i]==TRUE) countD=countD+1; #potenciales columnas afectadas
  countLimit=countD*modify/100.0; #columnas afectadas
  count=0;
  
  newMatrix=peaksMatrix[, logic]; #copia de la matriz original
  nCols=ncol(newMatrix); #columnas de la matriz
  
  if(countLimit<countD)
  {
    set.seed(1); #semilla para random
    rand=runif(countLimit, 1, countD);
  }
  else
  {rand=1:countLimit;}
  
  set.seed(1); #semilla para random
  for(index in 1:length(rand))
  {
    ion=as.integer(rand[index]);
    logic=peaksMatrix[,ion]!=0;
    meanMag=mean(peaksMatrix[logic,ion]);
    meanMag=meanMag*magnitude/100;
    
    #si tipo RANDOM y la columna debe alterarse
    if(mode=="RANDOM" && count<countLimit)
    {
      count=count+1; #contador de columnas alteradas
      #newIntensity=runif(nPixels, 0, max(peaksMatrix[,ion]));#valores aleatorios
      for(px in 1:nPixels) 
      {
        #if(peaksMatrix[px,index]!=0) #solo se alteran los píxeles no nulos
        {
          #delta=newMatrix[px,ion]*magnitude/100;
          #rIntensity=runif(1, -delta, delta)
          rIntensity=runif(1, -meanMag, meanMag)
          newMatrix[px,ion]=newMatrix[px,ion]+rIntensity;
          if(newMatrix[px,ion]<0) newMatrix[px,ion]=-newMatrix[px,ion];
        }
        
      }
    }
    #si tipo NULL y la columna debe alterarse
    else if(mode=="NULL" && count<countLimit)
    {
      count=count+1;
      newMatrix[, ion]=0;
    }
  }
  return(newMatrix);
}

