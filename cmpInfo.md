---
title: "Complementary information"
author: "E del Castillo"
date: '2024-04-27'
output: html_document
---


## **Graphics tools**:


### rPlotDeconv3<-function(gaussInfo1, gaussInfo2, gaussInfo3, minMass=0, maxMass=0, rMSI2_peaks1=-1, rMSI2_peaks2=-1, rMSI2_peaks3=-1)
> #### Presents an image with the gaussians passed as arguments for tree data.
> **Description of the parameters:**
```
    gaussInfo1: gaussians info from low    resolution peak: rGetGaussians();   
    gaussInfo2: gaussians info from median resolution peak: rGetGaussians();   
    gaussInfo3: gaussians info from high   resolution peak: rGetGaussians();   
       minMass: minimun mz;   
       maxMass: maximun mz;   
  rMSI2_peaks1: low    resolution peaks from rMSI2 peak matrix;   
  rMSI2_peaks2: median resolution peaks from rMSI2 peak matrix;   
  rMSI3_peaks3: high   resolution peaks from rMSI2 peak matrix;  
  NOTE: rMSI2_peaks come from rMSI2::LoadPeakMatrix(), from rMSI2::processWizard()
```

> For **example**:  
> rPlotDeconv3(gaussInfo_A, gaussInfo_B, gaussInfo_C, 769.4, 769.75);   
> rPlotDeconv3(gaussInfo_A, gaussInfo_B, gaussInfo_C, 769.4, 769.75, rMSI2Pks_A\$mass, rMSI2Pks_B\$mass);  


### rPlotDeconv2<-function(gaussInfo1, gaussInfo2, minMass=0, maxMass=0, rMSI2_peaks1=-1, rMSI2_peaks2=-1)
> #### Presents an image with the gaussians passed as arguments for two data

> For **example**:
> rPlotDeconv2(gaussInfo_A, gaussInfo_B, 769.4, 769.75);   
> rPlotDeconv2(gaussInfo_A, gaussInfo_B, 769.4, 769.75, rMSI2Pks_A\$mass);


### rPlotDeconv\<-function(gaussInfo)
>#### Presents an image with the Gaussians passed as arguments,
>
> **Description of the parameters**  
> gaussInfo: information returned by rGetGaussians().

> For **example**:  
> rPlotDeconv(gaussInfo)


### rPlotPeaks<-function(gaussians, rMSI2peaks, minMass, maxMass) 
>#### Presents an image with the gaussians info passed as arguments for one data.

> For **example**:  
> rPlotPeaks(gaussInfo\$gaussians, rMSI2peaks, 769.4, 769.75)


### rPlotSirem<-function(siremInfo, minMz, maxMz)
>#### Presents an image with the concentration and SIREM info passed as arguments.

> **Description of the parameters:**
```
  siremInfo: sirem data and peaks obtained from rGetSiremPeaks();   
    minMass: minimun mz;   
    maxMass: maximun mz;  
```

## **Tests tools:**


### fitQualitySirem<-function(reference, testSirem, testGauss, refMinMag=1e-6, testMinMag=1e-6, minMass=0, maxMass=0)
>#### Identifies the high resolution peak closest to the low resolution peak and indicates the deviation. It does it on all the peaks. For information about Sirem.

> **Description of the parameters**
```
   reference: gaussians info from rGetGaussians() (high resolution);   
   testSirem: sirem     info from rGetGaussians() (low resolution);   
   testGauss: gaussians info from rGetGaussians() (low resolution);   
   refMinMag: minimum magnitude of peaks for consideration. (high resolution);   
  testMinMag: minimum magnitude of peaks for consideration. (low resolution);   
     minMass: low  mass to analyze (Da);   
     minMass: high mass to analyze (Da);   
```

> **return** a matrix with columns: "mzTest", "mzRef", "ppm", "maxDev", "repe", "deconv" 

> For **example**:   
> fqSirem<-fitQualitySirem(gaussInfo_A, siremPeaks_B, gaussInfo_B);



### fitQualitySiremDeconv<-function(reference, testSirem, testGauss, refMinMag=1e-6, testMinMag=1e-6, minMass=0, maxMass=0)
>#### Identifies the high resolution peak closest to the low resolution peak and indicates the deviation. It does so exclusively on the deconvolved peaks. For information about Sirem.

> For **example**:   
> fqSiremD<-fitQualitySiremDeconv(gaussInfoA, siremPeaks_B, gaussInfo_B)



### fitQualityPere<-function(reference, testPere, testSirem, refMinMag=1e-6, testMinMag=1e-6, minMass=0, maxMass=0);
>#### Identifies the high resolution peak closest to the low resolution peak and indicates the deviation. It does it on all the peaks. For info on rMSI2 (Pere).

> **Description of new the parameters:**   
testPere   -> list from  rMSI2::LoadPeakMatrix() from rMSI2::processWizard()

>For **example**:   
> fqrMSI2<-fitQualityPere(gaussInfo_A, perePks, gaussInfo_B)


### peaksDeconvolved<-function(peaksInfo, gaussInfo, uPeakList=c())
>#### For each composite magnitude peak containing simple sirem peaks, the included sirem peaks and associated Gaussians are noted.

> **Description of new the parameters:**   
> uPeakList: lista de picos compuestos a considerar, esten o no deconvolucionados. 

> **return** a list: siremMassMatrix, gaussMassMatrix, gaussMagMatrix, gaussMassList.  
> siremMassMatrix: maintains the masses associated with sirem peaks;  
> gaussMassMatrix: maintains the masses associated with the Gaussians after adjustment with the average value of magnitudes;  
> gaussMagMatrix: maintains the magnitudes associated with each Gaussian;
> gaussList: is an array with the Gaussians of gaussMassMatrix; 

> For **example**:   
> pkD<-peaksDeconvolved(siremPeaks, gaussInfo, c())



### sirem_vs_rMSI2<-function(sample, SNR)
>#### Generates statistical data with the deviations of rSIREM and rMSI2 with respect to a standard sample.
> Valid for a particular case; not generalizable. Requires prior information.

> **Description of the parameters:**
```
  sample: sample. Valids: "C30k" y "C60k";   
  SNR   : signal to noise ratio. Valids: 1, 2, 3, 5, 7;
```

> For **example**:   
> sirem_vs_rMSI2<-function("30k", 1);

# Examples
### Deconvolution of a peak over the m/z range 769.4 to 769.65

> myData<-rMSI2::LoadMsiData("absolute_path_to_file_30k.imzML");

> params<-list("algorithm"=0, "minMeanPxMag"=1, "minSectionDensity"=5, "noiseLevel"=100, "tileSide"=4, "referenceType"=1, "cutLevels"=c(50), 
             "siremSensitivity"=c(0.0001, 0.05), "magSensitivity"=c(0.01, 1), "normalization"=0, "relativeMag"=0);
             
> siremPeaks_30k<-rGetSiremPeaks(myData, params, 769.4, 769.65);

**On the overaged spectrum of all pixels.** 

> gaussInfo_30k<-rGetGaussians(myData, siremPeaks_30k, 0, 1);   
> rPlotDeconv(gaussInfo_30k)

**About the 100 pixel spectrum.**  

> gaussInfo_30k<-rGetGaussians(myData, siremPeaks_30k, 100, 1);   
> rPlotDeconv(gaussInfo_30k)

### Deviations from the pattern

> myData<-rMSI2::LoadMsiData("absolute_path_to_file_120k.imzML");

> params<-list("algorithm"=0, "minMeanPxMag"=1, "minSectionDensity"=5, "noiseLevel"=100, "tileSide"=4, "referenceType"=1, "cutLevels"=c(50), 
             "siremSensitivity"=c(0.0001, 0.05), "magSensitivity"=c(0.01, 1), "normalization"=0, "relativeMag"=0);
             
> siremPeaks_120k<-rGetSiremPeaks(myData, params, 769.4, 769.65);

> gaussInfo_120k<-rGetGaussians(myData, siremPeaks_120k, 0, 1);   

> rPlotDeconv2(gaussInfo_30k, gaussInfo_120k)

 **Include rMSI2 centroids**
 
> After the peak matrix has been generated with rMSI2::processWizard():

> rMSI2_peaks <- rMSI2::LoadPeakMatrix(file.path("absolute_path_to_folder", "merged-peakmatrix.pkmat"))

> rPlotDeconv2(gaussInfo_30k, gaussInfo_120k, 769.4, 769.75, rMSI2_peaks$mass)

### Quality tests

> NOTE: For this test it is appropriate to deconvolve a wide range of masses.

> fq<-fitQualitySiremDeconv(gaussInfo_120k, siremPeaks_30k, gaussInfo_30k)

### Deconvolution information

> pkD<-peaksDeconvolved(siremPeaks_30k, gaussInfo_30k)


