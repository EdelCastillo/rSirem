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
fitQualitySirem<-function(reference, testSirem, testGauss, refMinMag=1e-6, testMinMag=1e-6, minMass=0, maxMass=0)
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
    
    if(ppm>1.5*testPPM) 
      {deviation[iPk, 4]=2; deviation[iPk, 3]=0;}
    else if(ppm>testPPM) 
      {deviation[iPk, 4]=1;}
    else
      {deviation[iPk, 4]=0;}
    
    deviation[iPk, 5]=0;
    if(iPk>1)
      if(deviation[iPk, 2]==deviation[iPk-1, 2])
      {deviation[iPk, 5]=1;}
    
    deviation[iPk, 6]=0;    
    if(length(pksDeconv$gaussMassList)>0)
    {
      retMass<-nearestValue(testMass, pksDeconv$gaussMassList);
      if(retMass==testMass)
        deviation[iPk, 6]=1;
    }
  }
  hist(deviation[, 3], main="Histogram of rSirem deviations", xlab="ppm", ylab="Frequency");
  #  legend("topright", legend="SNR=1");
  
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
fitQualitySiremDeconv<-function(reference, testSirem, testGauss, refMinMag=1e-6, testMinMag=1e-6, minMass=0, maxMass=0)
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
      {deviation[iPk, 4]=2; deviation[iPk, 3]=0;}
    else if(ppm>testPPM) 
      {deviation[iPk, 4]=1;}
    else
      {deviation[iPk, 4]=0;}
    
    deviation[iPk, 5]=0;
    if(iPk>1)
      if(deviation[iPk, 2]==deviation[iPk-1, 2])
      {deviation[iPk, 5]=1;}
    
  }
  hist(deviation[, 3], main="Histogram of deconvolved rSirem deviations", xlab="ppm", ylab="Frequency");
  #  legend("topright", legend="SNR=1");
  
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
fitQualityPere<-function(reference, testPere, testSirem, refMinMag=1e-6, testMinMag=1e-6, minMass=0, maxMass=0)
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
      {deviation[iPk, 4]=2; deviation[iPk, 3]=0;}
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
  hist(deviation[, 3], main="Histogram of rMSI2 deviations", xlab="ppm", ylab="Frequency");
#  title("rMSI2 deviations");
#  legend("topright", legend="SNR=1");
  
  colnames(deviation)<-c("mzTest", "mzRef", "ppm", "maxDev", "repe")
  return(deviation);
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
#' Para cada pico de magnitud compuesto que contiene a picos simples de sirem,
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
    #retM mantiene los índices a  las masas centrales de cada pico de sirem
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
#' Retorna un array indicando si el pico compuesto debe ser considerado o no
#' Sólo se consideran aquellos que contienen a un pico de rMSI
#' @param siremPeaks -> sirem peaks, from rGetSiremPeaks()
#' @param rMSIPeaks  -> rMSI2 peaks, from  rMSI2::LoadPeakMatrix() from rMSI2::processWizard()
#' @export
#' 
siremPeaksFilter<-function(siremPeaks, rMSIPeaks)
{
  nUpeaks=length(siremPeaks$siremPeaks$unitedMagnitudePeaks[,1]);
  uPks=rep(FALSE, times=nUpeaks);
  
  for(up in 1:nUpeaks)
  {
    uPkLow =siremPeaks$siremPeaks$unitedMagnitudePeaks[up, 1]+1;
    uPkHigh=siremPeaks$siremPeaks$unitedMagnitudePeaks[up, 2]+1;
    uLowPksDa =scans2Daltons(siremPeaks$siremPeaks$magnitudePeaks[uPkLow,  1]+1, siremPeaks$massAxis);
    uHighPksDa=scans2Daltons(siremPeaks$siremPeaks$magnitudePeaks[uPkHigh, 3]+1, siremPeaks$massAxis);

    lowNearMass =nearestValue(uLowPksDa,  rMSIPeaks$mass);
    highNearMass=nearestValue(uHighPksDa, rMSIPeaks$mass);
    if((lowNearMass  >=uLowPksDa & lowNearMass  <=uHighPksDa) | 
       (highNearMass >=uLowPksDa & highNearMass <=uHighPksDa))
        {uPks[up]=TRUE;}
  }
  return(uPks);
  
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
#'
#' @return nothing
#'
#' @export
#' 
sirem_vs_rMSI2<-function(sample, SNR)
{
  library(rSirem)
  if(SNR!=1 & SNR!=2 & SNR!=3 & SNR!=5 & SNR!=7) {print("Warning: unknown SNR; expected:1,2,3,5,7");  return();}
  load("/home/esteban/MALDI/rSirem/C120_all.RData")
  if(sample=="C60k")
    {
    #carga ficheros con los rangos de masas para la muestra de 60k
    load("/home/esteban/MALDI/rSirem/C60_all.RData")
    #carga info sobre la muestra
    myData<-rMSI2::LoadMsiData("/media/esteban/disk_1TB/esteban/MALDI/Samples/Cerebellum_30_60_120k/NoAlineado/C_60k/231211_Au_P_MBr_cblm_60k.imzML");
    #carga la matriz de picos generada por rMSI2
    peakMatrixPath<-sprintf("/media/esteban/disk_1TB/esteban/MALDI/Samples/Cerebellum_30_60_120k/Alineado/C_60k_SNR%d", SNR)
    rMSI2_snr <- rMSI2::LoadPeakMatrix(file.path(peakMatrixPath, "merged-peakmatrix.pkmat"))
    }
  else if(sample=="C30k")
    {
    load("/home/esteban/MALDI/rSirem/C30_all.RData")
    myData<-rMSI2::LoadMsiData("/media/esteban/disk_1TB/esteban/MALDI/Samples/Cerebellum_30_60_120k/NoAlineado/C_30k/231211_Au_P_MBr_cblm_30k.imzML");
    peakMatrixPath<-sprintf("/media/esteban/disk_1TB/esteban/MALDI/Samples/Cerebellum_30_60_120k/Alineado/C_30k_SNR%d", SNR)
    rMSI2_snr <- rMSI2::LoadPeakMatrix(file.path(peakMatrixPath, "merged-peakmatrix.pkmat"))
    }
  else {print("Warning: unknown sample");  return();}
  sprintf("", sample)
  
  #para el rango de masa de 700 a 900 Da
  if(sample=="C30k")
    {siremPeaks<-siremPeaks30_700_900n10s0;} #fichero de picos de sirem
  else if(sample=="C60k")
    {siremPeaks<-siremPeaks60_700_900n10s0;}
  #Se extraen los picos de sirem 'coincidentes' con los de rMSI2
  goodUpeaks=siremPeaksFilter(siremPeaks, rMSI2_snr);
  #Se marcan los que no 'coinciden' para su descarte
  siremPeaks$siremPeaks$unitedMagnitudePeaks[,1][!goodUpeaks]=-1;
  #Se obtienen las gausianas
  gaussInfo<-rGetGaussians(myData, siremPeaks, 0, 0.1);
  #desviaciones de masa de cada pico sobre el patrón
  fq_700_900_rMSI_snr  <-fitQualityPere       (gaussInfo120_700_900n10ns, rMSI2_snr,  gaussInfo)
  fq_700_900_sirem_snr <-fitQualitySirem      (gaussInfo120_700_900n10ns, siremPeaks, gaussInfo)
  fq_700_900_sirem_snrD<-fitQualitySiremDeconv(gaussInfo120_700_900n10ns, siremPeaks, gaussInfo)
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
  goodUpeaks=siremPeaksFilter(siremPeaks, rMSI2_snr);
  siremPeaks$siremPeaks$unitedMagnitudePeaks[,1][!goodUpeaks]=-1;
  gaussInfo<-rGetGaussians(myData, siremPeaks, 0, 0.1);
  fq_500_700_rMSI_snr  <-fitQualityPere       (gaussInfo120_500_700n10ns, rMSI2_snr,  gaussInfo)
  fq_500_700_sirem_snr <-fitQualitySirem      (gaussInfo120_500_700n10ns, siremPeaks, gaussInfo)
  fq_500_700_sirem_snrD<-fitQualitySiremDeconv(gaussInfo120_500_700n10ns, siremPeaks, gaussInfo)
  #se eliminan de la matriz los picos deconvolucionados
  logic<-fq_500_700_sirem_snr[, 6]==0; 
  rIndex=1:length(logic);
  rIndex=rIndex[logic];
  fq_500_700_sirem_snrND=fq_500_700_sirem_snr[rIndex,];
  
  #para el rango de masa de 300 a 500 Da
  if(sample=="C30k")
  {siremPeaks<-siremPeaks30_300_500n10s0;}
  else if(sample=="C60k")
  {siremPeaks<-siremPeaks60_300_500n10s0;}
  goodUpeaks=siremPeaksFilter(siremPeaks, rMSI2_snr);
  siremPeaks$siremPeaks$unitedMagnitudePeaks[,1][!goodUpeaks]=-1;
  gaussInfo<-rGetGaussians(myData, siremPeaks, 0, 0.1);
  fq_300_500_rMSI_snr  <-fitQualityPere       (gaussInfo120_300_500n10ns, rMSI2_snr,  gaussInfo)
  fq_300_500_sirem_snr <-fitQualitySirem      (gaussInfo120_300_500n10ns, siremPeaks, gaussInfo)
  fq_300_500_sirem_snrD<-fitQualitySiremDeconv(gaussInfo120_300_500n10ns, siremPeaks, gaussInfo)
  #se eliminan de la matriz los picos deconvolucionados
  logic<-fq_300_500_sirem_snr[, 6]==0; 
  rIndex=1:length(logic);
  rIndex=rIndex[logic];
  fq_300_500_sirem_snrND=fq_300_500_sirem_snr[rIndex,];
  
  #Resultados para rMSI2
  #desviaciones para el rango de masas unificado (300:900 Da)  
  totalDiff_rMSI2<-c(fq_300_500_rMSI_snr[,3], fq_500_700_rMSI_snr[,3], fq_700_900_rMSI_snr[,3])
  #Histograma de las desviaciones
  hist(totalDiff_rMSI2, main="Histogram of rMSI2 deviations", xlab="ppm", ylab="Frequency");
  legendTxt<-sprintf("mz=300:900 Da\nsample=%s; SNR=%d", sample, SNR)
  legend("topright", legend=legendTxt);
  #valores medios y desviación estándar
  rMSI2Mean=mean(totalDiff_rMSI2)
  rMSI2Sd=sd(totalDiff_rMSI2)
  #presenta resultados
  msg<-sprintf("           total rMSI2  peaks=%5d; mean=%.4f; sigma=%.4f", length(totalDiff_rMSI2), rMSI2Mean, rMSI2Sd);
  print(msg)
  
  #Resultados para los picos de rSirem descartando los picos deconvolucionados
  totalDiff_siremND<-c(fq_300_500_sirem_snrND[,3], fq_500_700_sirem_snrND[,3], fq_700_900_sirem_snrND[,3])
  logic=totalDiff_siremND>0 & totalDiff_siremND<50; #(<50)excluye deficiencias en la comparación
  totalDiff_siremND=totalDiff_siremND[logic];
  hist(totalDiff_siremND, main="Histogram of rSirem deviations (not deconv)", xlab="ppm", ylab="Frequency");
  legend("topright", legend=legendTxt);
  siremMean=mean(totalDiff_siremND)
  siremSd=sd(totalDiff_siremND)
  msg<-sprintf("total w/o deconv rSirem peaks=%5d; mean=%.4f; sigma=%.4f", length(totalDiff_siremND), siremMean, siremSd);
  print(msg)
  
  #Resultados para todos los picos de rSirem (incluye a los picos deconvolucionados) 
  totalDiff_sirem<-c(fq_300_500_sirem_snr[,3], fq_500_700_sirem_snr[,3], fq_700_900_sirem_snr[,3])
  logic=totalDiff_sirem>0 & totalDiff_sirem<50; #(<50)excluye deficiencias en la comparación
  totalDiff_sirem=totalDiff_sirem[logic];
  hist(totalDiff_sirem, main="Histogram of rSirem deviations", xlab="ppm", ylab="Frequency");
  legend("topright", legend=legendTxt);
  siremMean=mean(totalDiff_sirem)
  siremSd=sd(totalDiff_sirem)
  msg<-sprintf("           total rSirem peaks=%5d; mean=%.4f; sigma=%.4f", length(totalDiff_sirem), siremMean, siremSd);
  print(msg)
  
  #Resultados para los picos de rSirem considerando solo los picos deconvolucionados
  totalDiff_siremD<-c(fq_300_500_sirem_snrD[,3], fq_500_700_sirem_snrD[,3], fq_700_900_sirem_snrD[,3])
  hist(totalDiff_siremD, main="Histogram of rSirem deconvolution deviations", xlab="ppm", ylab="Frequency");
  legend("topright", legend=legendTxt);
  siremMean=mean(totalDiff_siremD)
  siremSd=sd(totalDiff_siremD)
  msg<-sprintf("    deconvoluted rSirem peaks=%5d; mean=%.4f; sigma=%.4f", length(totalDiff_siremD), siremMean, siremSd);
  print(msg)
  print("list of deconvolved peaks with deviation >= 1.5 scans:")
  totalWrong_mz<-c(fq_300_500_sirem_snrD[,4], fq_500_700_sirem_snrD[,4], fq_700_900_sirem_snrD[,4])
  logic=totalWrong_mz==1;
  totalWrong_mz<-c(fq_300_500_sirem_snrD[,1], fq_500_700_sirem_snrD[,1], fq_700_900_sirem_snrD[,1])
  totalWrong_mz<-totalWrong_mz[logic];
  totalWrong_mz
}
