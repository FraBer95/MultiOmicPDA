# MultiOmicPDA

This is the repository of the "A time-dependent explainable radiomic analysis from the multi-omic cohort of CPTAC-Pancreatic Ductal Adenocarcinoma" paper.

The code is written in R and Python and it is composed by two main parts:

Preprocessing part: Data Reading and formatting, dropping of non relevant features, categorization, statistical analysis and feature selection.
Analysis part: Survival Curves, Models Training and Validation, XAI(t).
The python libraries requested are:

-_pandas, torch, pycox, reticulate, numpy, seaborn, scipy, statsmodels, scikit-learn_

The R libraries requested are:

-_survival, surviminer, survex, randomForestSRC, gbm, survivalsvm, ggsurvfit, ggplot2, pec, caret, SurvMetrics, mlr3proba, mlr3extralearners, mlr3pipelines, survivalmodels, reticulate_
