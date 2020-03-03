# Frequently asked questions  

## Table of Contents  

[Where does the bias correction algorithm come from?](#where-does-the-bias-correction-algorithm-come-from)  
[What kind of data can be corrected by rBiasCorrection/BiasCorrector?](#what-kind-of-data-can-be-corrected-by-rBiasCorrectionbiascorrector)   
[Do my input files need to be in a special format?](#do-my-input-files-need-to-be-in-a-special-format)  
[Are there any requirements for naming the files?](#are-there-any-requirements-for-naming-the-files)  
[What is exactly done during rBiasCorrection's/BiasCorrector's data preprocessing?](#what-is-exactly-done-during-rBiasCorrectionsbiascorrectors-data-preprocessing)  
[What are the regression statistics?](#what-are-the-regression-statistics)  
[What are 'substitutions' in my final results?](#what-are-substitutions-in-my-final-results)  


## Where does the bias correction algorithm come from?  

rBiasCorrection/BiasCorrector is the user friendly implementation of the algorithms, described by Moskalev et. al in their article *'Correction of PCR-bias in quantitative DNA methylation studies by means of cubic polynomial regression'*, published 2011 in *Nucleic acids research, Oxford University Press* (DOI: [https://doi.org/10.1093/nar/gkr213](https://doi.org/10.1093/nar/gkr213)).  

### Citation:  
```
@article{10.1093/nar/gkr213,
    author = {Moskalev, Evgeny A. and Zavgorodnij, Mikhail G. and Majorova, Svetlana P. and Vorobjev, Ivan A. and Jandaghi, Pouria and Bure, Irina V. and Hoheisel, Jörg D.},
    title = "{Correction of PCR-bias in quantitative DNA methylation studies by means of cubic polynomial regression}",
    journal = {Nucleic Acids Research},
    volume = {39},
    number = {11},
    pages = {e77-e77},
    year = {2011},
    month = {04},
    abstract = "{DNA methylation profiling has become an important aspect of biomedical molecular analysis. Polymerase chain reaction (PCR) amplification of bisulphite-treated DNA is a processing step that is common to many currently used methods of quantitative methylation analysis. Preferential amplification of unmethylated alleles—known as PCR-bias—may significantly affect the accuracy of quantification. To date, no universal experimental approach has been reported to overcome the problem. This study presents an effective method of correcting biased methylation data. The procedure includes a calibration performed in parallel to the analysis of the samples under investigation. DNA samples with defined degrees of methylation are analysed. The observed deviation of the experimental results from the expected values is used for calculating a regression curve. The equation of the best-fitting curve is then used for correction of the data obtained from the samples of interest. The process can be applied irrespective of the locus interrogated and the number of sites analysed, avoiding an optimization of the amplification conditions for each individual locus.}",
    issn = {0305-1048},
    doi = {10.1093/nar/gkr213},
    url = {https://dx.doi.org/10.1093/nar/gkr213},
    eprint = {http://oup.prod.sis.lan/nar/article-pdf/39/11/e77/16775711/gkr213.pdf},
}
```

## What kind of data can be corrected by rBiasCorrection/BiasCorrector?  

<!-- BiasCorrector can handle two types of input data:  
  
- Type 1: one locus in many samples (e.g. pyrosequencing data)  
- Type 2: many loci in one sample (e.g. next-generation sequencing data or microarray data)
-->
Currently, both R packages, `rBiasCorrection` and `BiasCorrector`, can correct measurement biases in DNA methylation data of the type "one locus in many biological samples". The programme has been tested on data derived by bisulphite pyrosequencing, next-generation sequencing, and oligonucleotide microarrays. A future implementation is planned for correcting data of the type "many loci in one biological sample". However with some effort, the latter can be transformed to data of the first type in order be corrected with `BiasCorrector`. 

## Do my input files need to be in a special format?  

Yes, rBiasCorrection/BiasCorrector places very strict requirements on the file format. Below is a description of the exact requirements. <!-- for the two types of input data, which differ in several aspects.-->  
All uploaded files must  

* be in CSV format [file endings: \*.csv and \*.CSV]  
* contain the column headers in the first row   
* the number of CpG sites per locus ID has to be equal in every provided file  


<!--### Type 1: one locus in many samples (e.g. pyrosequencing data)  -->

* Experimental data:  
  + the first column contains the sample IDs (alpha-numeric characters only)  
  + the other columns contain the results of your methylation analysis with one column for each CpG site  
  + sample IDs may occur more than once, indicating repeated measurements of the same sample (in this case, the mean values of the replicates will be used for bias correction)   
  
* Calibration data:  
  + the first column contains the percentages of the actual methylation of the calibration samples (calibration steps, numeric)  
  + these calibration steps (CS) must be in the range 0 <= CS <= 100  
  + a minimum of four distinct calibration steps are required  
  + the other columns contain the calibration sample's results of the methylation analysis with one column for each CpG site  
  + calibration steps may occur more than once, indicating repeated measurements of the same calibration sample (in this case, the mean values of the repeated measurements will be used for calculation of the calibration curve)   

### Example files  

Example files are available for download, to demonstrate how to preprare files appropriately: 
* calibration data: [Example_calibration.csv](tests/testthat/testdata/cal_type_1.csv)  
* experimental data: [Example_experimental.csv](tests/testthat/testdata/exp_type_1.csv)  

### Template files  

Template files are available, if you want to copy-paste your data. Please note that you might have to adjust the column headers and sample IDs or calibration steps: 
* calibration data: [Template_calibration.csv](inst/template_calibration.csv)  
* experimental data: [Template_experimental.csv](inst/template_experimental.csv)  
  

<!--### Type 2: many loci in one sample (e.g. next-generation sequencing data or microarray data)  

- Experimental data: 

  -- the first column contains the locus IDs (alphanumeric)  
  -- locus IDs may occur more than once (indicating repeated measurements of the same locus; in this case, the mean-values of the repeated measurements will be used for bias correction)
  -- all other columns contain the results of your methylation analysis for each CpG site of the respective locus  
  
- Calibration data:  

  -- the first column contains the locus IDs (alphanumeric)  
  -- locus IDs may occur more than once (indicating repeated measurements of the same locus; in this case, the mean-values of the repeated measurements will be used for bias correction)
  -- all other columns contain the results of your methylation analysis for each CpG site of the respective locus  
  -- for bias correction of type 2 data, you need to provide one separate calibration file for each degree of methylation  
  -- a minimum of four calibration data files (four distinct calibration steps) are provided  
  -- all provided calibration files have equal dimensions (number of rows * number of columns), equal column names and equal locus IDs  
  -- the calibration steps (CS) must be in the range 0 <= CS <= 100  -->


## Are there any requirements for naming the files? 

The filename must not contain additional dots (".") beyond the one in the file ending.

<!--Files of the input data type 1 (one locus in many samples) do not place specific requirements for naming the files.  
On the contrary, files of the input data type 2 (many loci in one sample), and in particular the files containing the calibration data, place very strict requirements on the file naming:  

Every filname must be of the following pattern:  'anyfilename'_CS###.csv  

The suffix '_CS###.csv'  
- must begin with '_CS', otherwise the file is going to be rejected during the data preprocessing step (CS is the short form of 'calibration step')  
- the placeholder '###' must be replaced with the respective degree of methylation of the calibration data contained in the specific file  
-- it can be or an integer number between 0 and 100 (integers < 0 or > 100 will be rejected during the data preprocessing step)  
-- or a numeric number between 0 an 100, indicated by an underscore ('_') as decimal separator (e.g. '12_5' meaning '12.5')   

Example: to upload a file for bias correction of type 2, that contains the calibration data for the calibration step '12.5' (true degree of calibration = 12.5) it need to be named the following:  
  *'my-calibrationfile_CS12_5.csv'* --> 


## What is exactly done during rBiasCorrection's/BiasCorrector's data preprocessing?  

During the preprocessing, all requirements of the input files as stated in [Do my input files need to be in a special format?](#do-my-input-files-need-to-be-in-a-special-format) are checked. Furthermore, the mean methylation percentages of all CpG columns are calculated for every provided file. 

If any of the abovementioned file requirements is not met, an error will occur. For example, an error message will pop up if any calibration step is not within the range of 0 <= CS <= 100 or if you provided less than four calibration steps in your input data. 


## What are the regression statistics?

The regression statistics table shows the regression parameters of the hyperbolic and the cubic polynomial regression. 

- Column 1 presents the CpG site's ID. 
- Column 2 contains the mean of the relative absolute errors for every interrogated CpG site. 
- Columns 3-9 comprise the sum of squared errors of the hyperbolic regression ('SSE [h]') and the coefficients of the hyperbolic equation that describes the hyperbolic regression curves for the respective CpG sites. 
- Columns 10-15 summarise the sum of squared errors of the cubic polynomial regression ('SSE [c]') and the coefficients of the cubic polynomial equations. 
- The rows highlighted with a green background colour indicate the regression method (hyperbolic or cubic polynomial) that is suggested by BiasCorrector for correcting data. This automatic choice of the regression method relies on either minimising the value of SSE (the default setting) or minimising the average relative error as selected by the user in the Settings tab.


## What are 'substitutions' in my final results?

Substitutions occur if no result is found in the range of plausible values between 0 and 100 during the BiasCorrection. A 'border zone' is implemented in the ranges 0 – 10% and 100 + 10%. If a result is in the range -10% < x < 0% or 100% < x < 110% , the value is substituted in the final results with 0% or 100%, respectively. 
Values beyond these border zones will be substituted with a blank value in the final output, as they seem implausible and could indicate substantial errors in the underlying data. 
For a detailed feedback, the substitutions table shows the results of the algorithm 'BiasCorrected value' and the corresponding substitution 'Substituted value' for the respective CpG site. 
