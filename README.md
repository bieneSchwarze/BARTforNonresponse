# BARTforNonresponse
Analyzing Nonresponse in Longitudinal Surveys Using Bayesian Additive Regression Trees: A Nonparametric Event History Analysis
Date: 18.12.2019 (Version 2.0.0, revision 1)

This is the online material for the publication (currently under review): Analyzing Nonresponse in Longitudinal Surveys Using Bayesian Additive Regression Trees:
A Nonparametric Event History Analysis

By Sabine Zinn (1,2) & Timo Gnambs (2,3)

(1) German Institute for Economic Research, 
(2) Leibniz Institute for Educational Trajectories, 
(3) Johannes Kepler University Linz

Abstract: 	
Increasing nonresponse rates are a pressing issue for many longitudinal panel studies. Respondents frequently either refuse participation in single survey waves (temporary dropout) or discontinue participation altogether (permanent dropout). Contemporary statistical methods that are used to elucidate predictors of survey nonresponse are typically limited to small variable sets and ignore complex interaction patterns. The innovative approach of Bayesian additive regression trees (BART) is an elegant way to overcome these limitations because it does not specify a parametric form for the relationship between the outcome and its predictors. We present a BART event history analysis that allows identifying predictors for different types of nonresponse to anticipate response rates for upcoming survey waves. We apply our novel method to data from the German National Educational Panel study including N = 4,559 students in grade 5 that observed nonresponse rates of up to 36% across five waves. A cross-validation and comparison with logistic regression models with LASSO (least absolute shrinkage and selection operator) penalization underline the advantages of the approach. Our results highlight the potential of Bayesian discrete time event modeling for the long-term projection of panel stability across multiple survey waves. Finally, potential applications of this approach for operational use in survey management are outlined.

We use data from the National Educational Panel Study (NEPS): Starting Cohort Grade 5, doi:10.5157/NEPS:SC3:8.0.0. 
From 2008 to 2013, NEPS data was collected as part of the Framework Program for the Promotion of Empirical Educational Research funded by the German Federal Ministry of Education and Research (BMBF). 
As of 2014, NEPS is carried out by the Leibniz Institute for Educational Trajectories (LIfBi) at the University of Bamberg in cooperation with a nationwide network. 
Data access requires the conclusion of a Data Use Agreement with the Leibniz Institute for Educational Trajectories (LIfBi), see https://www.neps-data.de/en-us/datacenter/dataaccess/datauseagreements.aspx.

This GitHub project contains the source code for data preparation and analysis (with BART and logistic regression with LASSO penalization), as well as for validation and plotting the results.
1. loadData.R (editing and preparing data using the NEPS SUF files)

2. modelPermanentDropout_BART.R (BART model for passing over to permanent drop out; in our case this means leaving the NEPS school context)

3. modelTempDropout_BART.R (BART model for temporary dropout in the school context)

4. combineResults_getJointProbabilities.R (predicting participation probabilities in the distinct survey waves for new samples)

A1.	modelPermanentDropout_BART_crossValidation.R (cross validation for BART model for permanent dropout)

A2.	modelPermanentDropout_BART_validationLeaveOneWaveOut.R (validation 'leave last wave out' for BART model for permanent dropout)

A3.	modelPermanentDropout_logitLasso.R (logistic regression with LASSO for permanent dropout, for model comparision)
A4.	modelPermanentDropout_logitLasso_crossValidation.R (cross validation for LASSO logistic regression model for permanent dropout)

A5.	modelPermanentDropout_logitLasso_validationLeaveOneWaveOut.R (validation 'leave last wave out' for LASSO logistic regression model for permanent dropout)

B1.	modelTempDropout _BART_crossValidation.R (cross validation for BART model for temporary dropout)

B2.	modelTempDropout _BART_validationLeaveOneWaveOut.R (validation 'leave last wave out' for BART model for temporary dropout)

B3.	modelTempDropout_logitLasso.R (logistic regression with LASSO for temporary dropout, for model comparision)

B4.	modelTempDropout _logitLasso_crossValidation.R (cross validation for LASSO logistic regression model for temporary dropout)

B5.	modelTempDropout _logitLasso_validationLeaveOneWaveOut.R (validation 'leave last wave out' for LASSO logistic regression model for temporary dropout)

