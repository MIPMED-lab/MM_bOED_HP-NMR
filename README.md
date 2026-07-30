# MM_bOED_HP-NMR
Code developed for the paper "A Bayesian Study to Optimally Unveil LDH Kinetic Mechanisms In Vivo" by David Gomez-Cabeza. 

## Directories: ##
- ***BayesianInference***: Directory containing all the scripts used for the Bayesian inference of parameters for the simplistic models (with and without pyruvate transport term) and the Michaelis Menten model. Only contains inferences using the experiment using 3.2 mM pyruvate. 
- ***ComputationalBayesianOED***: Directory containing all the scripts used for the in silico validation of the novel Bayesian Optimal Experimental Design strategy presented in this work compared to random and intuition-derived experiments. 
- ***ExperimentalBayesianOED_NoiseAssesment***: Directory containing notebooks to perform inferences in a pseudo-experiment with different levels of heteroscedastic noise and evaluate the gain in information of the posteriors. 
- ***ExperimentalBayesianOED_Validation***: Directory containing all the scripts used to perform the in vivo validation of the novel Bayesian Optimal Experimental Design strategy presented in this work. Structure is parallel to the one found in ComputationalBayesianOED.
- ***ExperimentalBayesianOED_Y0Measured***: Directory containing all the scripts used to employ our bOED approach to fine tune a mathematical model (competitive respresion by substrate) using informative priors for y0s thanks to experimental measurements. 
- ***InitalParameterFitPriorDef***: Directory containing scripts for parameter estimation (Maximum Likelihood Estimation) tests for all investigated models, used to adjust resonable uninformative priors for the model parameters. 
-***LDHActivityAssey***: Directory containing all the scritps used to process the experimental data for the LDH activity assay using increasing pyruvate concentrations and scripts used to predict/understant the ex vivo LDH activity and to infer model parameters using the extracted data. 
- ***LDHQuantification_Processing***: Directory containing all the scripts used to extract, process and plot all the flow cytometry data generated to quantify the LDH concentration in HepG2 cells. This includes antibody concentration selection curves, particle calibration to convert fluorescence to number of molecules and computation of number of LDH molecules for cells. 
- ***NADHQuantification_Processing***: Directory containing all the scripts used to extract and process the data from the NADH/NAD quantitation colorimetric assays (https://www.sigmaaldrich.com/ES/es/product/sigma/mak037?srsltid=AfmBOoriA5bwgg8jkIj1ytxyyDjLnLMD-kyeRT5F7tT0mke56WEb-sLl).
- ***NoCO2Asseys***: Directory containing all the scripts used to process and plot all the experimental data acquired to estudy the cause on the lower pyruvate production when removing CO2 in culture conditions overnight. This includes flow cytometry data (LDH and MCT1 estimations), colorimetric data (NADH) and pH-meter data (pH).
- ***ProcessedData***: Directory containing CSV files for all experiments with the estimated mean and standard deviation from experimental traces for HP pyruvate and lactate. 
- ***PyrBuffer_ViabilityAssayCellCount***: Directory containing the scripts used to extract, process and plot the data for the cell viability assays performed in this work (AlamarBlue HS and trypan blue staining). 
- ***StanCopy***: Directory containing the version copy of cmd-stan used in this work (https://mc-stan.org/). 
- ***StructuralIdentifiability***: Directory containing all the MATLAB scripts necessary to perform structural identifiability analysis of the models used in this work using GenSSI (https://github.com/genssi-developer/GenSSI) or Strike-Goldd (https://github.com/afvillaverde/strike-goldd). 


## Scripts: ##
-  ***GenerateAverageData***: Notebook to generate average data from experiment traces for the experiments considering only 3.2 mM pyruvate with changing conditions (number of cells, membrane integrity and CO2 availability). The code also extrapolates the HP-Pyruvate concentration at t=0 (mixing of pyruvate with cells) from its T1 decay value.
-  ***GenerateAverageData_DiffPyr***: Notebook to generate average data from experiment traces for the experiments considering 3 million cells with changing pyruvate concentrations. The code also extrapolates the HP-Pyruvate concentration at t=0 (mixing of pyruvate with cells) from its T1 decay value.
- ***VisualiseControlsHP***: Notebook to extract and plot the mean and standard deviation for the control experiments (0 mM pyruvate and 0 cells at 3.2 mM pyruvate).
- ***EntropyEstimationPosterior***: Julia script with the functions to estimate Entropy from multivariate priors/posterior distributions. See "Information content analysis reveals desirable aspects of in vivo experiments of a synthetic circuit". 
- ***GenerateExamplePlots_Figure1***: Notebook to generate example distribution and experiment plots. Used to generate panels for Figure 1 of the paper. 

## HELP ##
For any inquiry, please contact Dr. David Gomez-Cabeza at dgomez@ibecbarcelona.eu or dgomezcabeza@outlook.com. 
