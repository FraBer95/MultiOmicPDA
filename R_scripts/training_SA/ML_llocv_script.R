library(survival)
library(survminer)
library(survex)
library(ggsurvfit)
library(ggplot2)
library(pec)
library(caret)
library(SurvMetrics)
library(mlr3extralearners)
library(mlr3pipelines)
library(paradox)
library(mlr3tuning)
library(survivalmodels)
library(reticulate)
library(progressr)
library(survivalsvm)
library(SurvMetrics)
library(gbm)
library(mlr3proba)

source("./Radiomics/utils.R")

use_condaenv("E:\\Users\\Berloco\\anaconda3\\envs\\survEnv", required = TRUE)

py_set_seed(1969)
set.seed(1969)


OS <- FALSE
REC <- TRUE

if (OS) {
  data <- read.csv("./Radiomics/data_OS_HPC.csv")
} else {data <- read.csv("./Radiomics/data_REC_HPC.csv")}


RADIOMICS <- c(TRUE, FALSE)
CLINICAL <- c(TRUE, FALSE)
GENOMIC <- c(TRUE, FALSE)

flag_combs <- expand.grid(RADIOMICS, CLINICAL, GENOMIC)

gbm <- lrn("surv.gbm", n.trees=100, interaction.depth = 10, bag.fraction=0.9, shrinkage=0.001)
svm <- lrn("surv.svm", type = "vanbelle2", diff.meth = "makediff3", gamma.mu = 0.2)

learners <- list()
learners[[1]] <- svm
learners[[2]] <- gbm

#dataframe with overall perfomance
df_perf <- data.frame(
    config = character(),
    model = character(),
    C_Index = numeric(),
    Brier_score = numeric(),
    stringsAsFactors = FALSE
  )

for (learner in learners) {

  for (i in 1:nrow(flag_combs)) {

    current_flags <- flag_combs[i, ]

    if (any(current_flags)){

      RADIOMICS <- current_flags[[1]]
      CLINICAL <- current_flags[[2]]
      GENOMIC <- current_flags[[3]]

    #print(paste("Radiomics", RADIOMICS, "Clinical", CLINICAL, "Genomic", GENOMIC, sep=" "))
    if (GENOMIC & !CLINICAL & !RADIOMICS & learner$id == "surv.gbm")
    {
      print("Set minimun objs in Node = 5, for genomic train in GBM for OS task")
      learner <- lrn("surv.gbm", n.trees=100, interaction.depth = 10, bag.fraction=0.9, shrinkage=0.1,  n.minobsinnode=5)
    }


    #returns in a list (data, configuration)
    if (OS) {
      dta_conf <- feature_sel_OS(data, RADIOMICS, CLINICAL, GENOMIC)
    } else{ dta_conf <- feature_sel_REC(data, RADIOMICS, CLINICAL, GENOMIC)}

    dta <- dta_conf[[1]]
    conf <- dta_conf[[2]]

    perfomance_list <- list()
    test_probs <- list()
    surv_func_list <- list()



    paste("Training Model:", learner$id, sep = " ")

    avg_CI <- 0
    avg_AUC <- 0
    avg_Brier <- 0



    for (i in 1:dim(dta)[1]) {

        start_time <- Sys.time()
        print(paste("Iteration n.",i,"of", dim(dta)[1], sep = " "))


        dta_LOOCV <- dta[-i,]
        time <- dta_LOOCV$time
        status <- dta_LOOCV$status

        train_task <- as_task_surv(x = dta_LOOCV,
                                 time = "time",
                                 event = "status")


        composite_learner <- as_learner(ppl(
          "distrcompositor",
          learner = learner,
          estimator = "kaplan",
          form = "ph"))

        composite_learner$train(train_task)

        class(composite_learner) <- c(class(composite_learner), "LearnerSurv")

        if (learner$id == "surv.svm"){

                model_name <- "survSVM"

                svm_predict <- function(model, newdata, times) {
                    if (nrow(newdata) == 1){
                        newdata <- rbind(newdata, newdata)
                        t(model$predict_newdata(newdata)$distr$survival(times))[1, , drop=FALSE]
                    }
                    else{
                        t(model$predict_newdata(newdata)$distr$survival(times))
                    }
                }

                model_explainer <- survex::explain(composite_learner,
                                          data = dta_LOOCV[, -c(1,2)],
                                          y = Surv(dta_LOOCV$time, dta_LOOCV$status),
                                          predict_survival_function = svm_predict,
                                          label = paste0(model_name, "Explainer", sep = " "),
                                          verbose = FALSE)

        } else {

          model_name <- "survGBM"
          model_explainer <- survex::explain(composite_learner,
                                          data = dta_LOOCV[, -c(1,2)],
                                          y = Surv(dta_LOOCV$time, dta_LOOCV$status),
                                          label = paste0(model_name, "Explainer", sep = " "),
                                          verbose = FALSE)

          }

        #Performance on training
        perf <- model_performance(model_explainer)

        print(paste("CI", perf$result[[1]], "AUC", perf$result[[2]], "Brier", perf$result[[4]] ))
        avg_CI <- avg_CI+perf$result[[1]]
        avg_AUC <- avg_AUC+perf$result[[2]]
        avg_Brier <- avg_Brier+perf$result[[4]]



        print(paste0("Computing survshap..."))
        survshap <- predict_parts(model_explainer,
                                  new_observation = dta[i, -c(1,2)],
                                  aggregation_method = "integral",
                                  calculation_method = "kernelshap",
                                  type="survshap")

        if (i==1) {

            #creating only for having data structures, the values will be updated on each iteration
            print("Creating global survshap for all test explanations...")

            eval_times <- survshap$eval_times
            surv_f_mtrx <- model_explainer$predict_survival_function(composite_learner, dta[, -c(1,2)], eval_times)

            glob_surv <- predict_parts(model_explainer,
                                  new_observation = dta[1:dim(dta)[1], -c(1,2)],
                                  aggregation_method = "integral",
                                  calculation_method = "kernelshap",
                                  type="survshap")
        }

        if (model_name == 'survSVM'){
          new_data <- rbind(dta[i, -c(1,2)], dta[i, -c(1,2)])
          test_probs <- append(test_probs, composite_learner$predict_newdata(new_data)$crank[1])

        } else {
          test_probs <- append(test_probs, composite_learner$predict_newdata(dta[i, -c(1,2)])$crank)
        }

        surv_f_mtrx[i,] <- model_explainer$predict_survival_function(composite_learner, dta[i, -c(1,2)], eval_times)


        print(paste("Aggregating explanation of sample number", i, "...", sep = " "))
        glob_surv[["aggregate"]][[i]] <- survshap[["aggregate"]]

        end_time <- Sys.time()
        time_elapsed <- end_time - start_time
        print(paste("Elapsed time for iteration", i, ":", time_elapsed, units(time_elapsed), sep = " "))
      }


      test_CI <- c_index(survival::Surv(dta$time, dta$status), risk = unlist(test_probs))
      test_IBS <- integrated_brier_score(survival::Surv(dta$time, dta$status), surv = surv_f_mtrx, times = eval_times)


      avg_CI <- avg_CI/dim(dta)[1]
      avg_AUC <- avg_AUC/dim(dta)[1]
      avg_Brier <- avg_Brier/dim(dta)[1]

      avg_metrics <- list(avg_CI, avg_AUC, avg_Brier)


      print(paste("Evaluation Metrics for", model_name, "...", sep= " "))
      print(paste("Average C-Index", avg_CI, "Average C/D AUC", avg_AUC, "Average Brier", avg_Brier, sep= " "))



      plot(glob_surv, subtitle = paste("Created for", model_name, "after LOOCV (m models X n patients, 60 X 60)", sep = " "))
      plot(glob_surv, geom = "beeswarm", subtitle = paste("Created for", model_name, "after LOOCV (m models X n patients, 60 X 60)", sep = " "))

      save_ws(model_name, RADIOMICS, CLINICAL, GENOMIC, OS)
      df_current <- data.frame(
                  config = conf,
                  model = model_name,
                  C_Index = test_CI,
                  Brier_score = test_IBS,
                  stringsAsFactors = FALSE
      )
      df_perf <- rbind(df_perf, df_current)
    }
  }
}

if (OS) {
  write.csv(df_perf, './Radiomics_result/OS/Performances_OS.csv', row.names = FALSE)
} else { write.csv(df_perf, './Radiomics_result/REC/Performances_REC.csv', row.names = FALSE)}



