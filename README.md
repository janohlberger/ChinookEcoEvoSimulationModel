# ChinookEcoEvoSimulationModel
The code in this repository allows reconstruction of the results published in Ohlberger et al. 2019 PNAS.

Ohlberger J, Schindler DE, Ward EJ, Walsworth TE, Essington TE. Resurgence of an apex marine predator and the decline in prey body size. Proc National Acad Sci. 2019; 116(52):26682–9. 

The article DOI is https://doi.org/10.1073/pnas.1910930116, and the  ZENODO DOI is https://doi.org/10.5281/zenodo.6324335.

The repository consists of three files: 
1. ChinookEcoEvoCode.R contains the main model function ChinookEcoEvoModel()
2. ChinookEcoEvoParameters.R contains the input parameters to the model
3. ChinookEcoEvoRun.R sources the other two files to run the analyses

Simply open the file ChinookEcoEvoRun.R and run the few lines of code. These simulations do not require data input, though parameters were based on literature values, as described in Ohlberger et al. 2019. The analyses were run in R (v.3.5.1). Running the simulations using original parameter values may take several hours to days, depending on the machine/cluster.
