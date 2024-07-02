# MultiOmicPDA

This is the repository of the "A time-dependent explainable radiomic analysis from the multi-omic cohort of CPTAC-Pancreatic Ductal Adenocarcinoma" paper.

The code is written in R and Python and it is composed by the following  parts:
1. _data_ folder, that should contains all data in a csv format and the pyradiomics configuration *.yaml* file in
2. _python_scripts_ folder, containing all pre-processing operations for merging different data, extraction and validation of the radiomic signature
3. _R_scripts_ folder, containing the feature selection in UV, feature dicotomizing and survival models training and test 

Preprocessing part: Data Reading and formatting, dropping of non relevant features, categorization, statistical analysis and feature selection.
Analysis part: Survival Curves, Models Training and Validation, XAI(t).
The python libraries requested are:

-_pandas, torch, pycox, reticulate, numpy, seaborn, scipy, statsmodels, scikit-learn_

The R libraries requested are:

-_survival, surviminer, survex, randomForestSRC, gbm, survivalsvm, ggsurvfit, ggplot2, pec, caret, SurvMetrics, mlr3proba, mlr3extralearners, mlr3pipelines, survivalmodels, reticulate_
