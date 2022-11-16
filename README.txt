Data and code for the article "A role for population heterogeneity in Plasmodium vivid recurrence".
There are two folders, one for each of the two used data sets:
 
 1. Papua New Guinea data
 2. Thailand-Myanmar data

Papua New Guinea data folder:
 - Data folder:
	All "Albinama" files are from Robinson et al. (2015) and the original source https://datadryad.org/stash/dataset/doi:10.5061/dryad.m1n03
	Combined_data.RData: combination of data from the different Albinama files
 - Matlab functions:
	modelfit_1rec_1w.m:
	modelfit_1rec_vill.m:
	modelfit.m:
 - Model_fit_1_recurrence.m:
	Fitting to all data (Fig. 2, Fig. S2)
	Confidence intervals for fitting to all data (Tables S1-4)
	Fitting to data by village (Fig. S3, Fig. S4)
	Confidence intervals for fitting to data by village (Tables S5-8)
 - PNG_data_analysis.R:
	Survival curves (Fig. 1B)
	Weekly incidence rate (Fig. 1D)
 	Transmission and relapse rate (from model fits) (Fig. 5F)
 - Results folder:
	b-CI(11opt)m1-PNG-all.mat: data to compute the 95% CIs for the parameters from the fit to the all data and model 1 (also files for models 2-4)
	b-CI(11opt)m1-PNG-vill.mat: data to compute the 95% CIs for the parameters from the fit to the data by village and model 1 (also files for models 2-4)


Thailand-Myanmar data folder:
 - Analysis of simulated data.R:
	Correlation between first and second recurrence times (Table 1, Table S12)
	Survival curves (Fig. 4C and D, Fig. S12)
	Average number of recurrences per individual and year by time to first recurrence (Fig. S13)
 - Data folder:
	Combined_Time_Data.mat: the file "Combined_Time_Event.RData" converted to at MATLAB data file.
	Combined_Time_Event.RData: data from Taylor et al. (2019) from the original source https://github.com/jwatowatson/RecurrentVivax/blob/master/RData/TimingModel/Combined_Time_Event.RData
 - Data analysis and plots.R:
	Contribution to recurrences by the number of recurrences (Table S10)
	First and second recurrence time: correlations (Table S11, Table 1, Fig. S5, Table S12), Cox regression (Table S13) and likelihood-ratio test for comparing models 3 & 4 fit to first 2 recurrences
	Number of recurrences by time to first recurrence (Fig. S13)
 - Matlab functions:
	modelfit_1rec.m: fits a model to data for one recurrence (grouped by drug and study)
  	modelfit_2rec.m: fits a model to data for the first and second recurrence (grouped by drug and study)
 	modelfit.m: computes the fraction of uninfected individuals at different time points
 - Model_fit_1st_recurrence.m:
	This script uses the functions modelfit.m and modelfit_1rec.m.
	Fit to a part of the 1st recurrence data, PMQ+ vs blood-stage (Fig. S9)
	Fitting to data grouped by drug and study (Fig. 2, Fig. S10)
	Confidence intervals for fitting to data grouped by drug and study (Tables S14-17)
 - Model_fit_2_recurrences.m:
	Fitting to the 1st and 2nd recurrence data (Maximum likelihood estimates in Tables S18-21)
	Confidence intervals via bootstrapping (CIs in Tables S18-21)
	Model fit results: parameter estimates for different follow-up schemes, 1 or 2 infection rates, and different numbers of relapse risk groups for model 3
 	Different plots (Figures S11, 17-21 and Tables S23-25)
 - Results:
	bootstrap_CI_data(11opt)m1(1st-rec).mat: data to compute the 95% CIs for the parameters from the fit to the first recurrence and model 1 (also files for models 2-4)
	bootstrap_CI_data(11opt)m1.mat: data to compute the 95% CIs for the parameters from the fit to the first and second recurrence and model 1 (also files for models 2-4)
	surv_curve_data_by_drug_and_study.mat: survival curve data to plot the data next to the model estimates (Fig. S17-21)
 - Survival curves.R:
	Time to first recurrence: blood-stage vs primaquine + blood-stage (Fig. 1A, Fig 1C)
	Survival curves for time to 1st recurrence by drug and study (Fig. S8)
	Survival curves for time to 2nd recurrence by drug and study (Fig. 3 confidence region)
	Survival curves for time to 2nd recurrence by time to first recurrence quartiles (Fig. S7)
 - Simulations_1year.m:
	Model simulations of 1,000 individuals for 1 year using models 1-3


