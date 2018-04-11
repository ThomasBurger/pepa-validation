README

This directory contains 3 folders:

- "preproc_data", which contains the datasets
	Exp1_R2_pep (also termed UPSPepx2Check-preproc.Rdata) 
	Exp1_R25_pep (also termed UPSpep25-preproc.Rdata) 

- "code", which contains all the necessary codes to run the experiments:
	addSharedPep.R		>> used to artificialy introduce shared peptides in the real datasets
	dataGenerator-v7.R	>> generates simulated data
	groupttest.R		>> PeptideModel
	Sim2MSnSet-corrected.R	>> Convert simulated data into the MSnSet format, for DAPAR compliancy
	test-one-simu.R		>> 
	test-one-UPSx2.R	>>
	test-one-UPSx25.R	>> 3 similar scripts, one for each type of data. Called by test-several-XXX.R
	test-several-simu.R	>>
	test-several-UPSx2.R	>>
	test-several-UPSx25.R	>> Main script - 3 similar scripts, one for each type of data. 

- "comparaisons" (the french word for "comparisons") which is empty, but where the PDFs of the plots will be automatically saved along the experiments.

1/ set the working directory to the folder named "codes"
setwd("./20170301/codes");

2/ Depending on the dataset to test, run one of the following commands
	source("test-several-simu.R");
	source("test-several-UPSx25.R");
	source("test-several-UPSx2.R");
Possibly,the 3 scripts can be run in a row.

3/ Wait for a while...(depending on your machine as well as on the number of iterations)

4/ Browse the results in the folder named "comparaisons"

5/ To modify the parameters (shared peptides, number of iterations, etc.) open the corresponding test-several-XXX.R file and go back to step 2.