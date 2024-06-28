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
library(randomForestSRC)

source("./Radiomics_result/utils.R")

use_condaenv("E:\\Users\\Berloco\\anaconda3\\envs\\survEnv", required = TRUE)

#py_set_seed(1969)
#set.seed(1969)
set_seed(1969, 1969, 1969)


OS <- FALSE
REC <- TRUE

if (OS) {
  data <- read.csv("./Radiomics/data_OS_HPC.csv")
} else {data <- read.csv("./Radiomics/data_REC_HPC.csv")}


RADIOMICS <- c(TRUE, FALSE)
CLINICAL <- c(TRUE, FALSE)
GENOMIC <- c(TRUE, FALSE)

flag_combs <- expand.grid(RADIOMICS, CLINICAL, GENOMIC)



gbm <- lrn("surv.gbm", n.trees=100, interaction.depth = 10, bag.fraction=0.9, shrinkage=0.001, n.minobsinnode=2)
#svm <- lrn("surv.svm", gamma.mu = 0.2)
svm <- lrn("surv.svm", type = "vanbelle2", diff.meth = "makediff3", gamma.mu = 0.2)

cox <- coxph(Surv(data$time,data$status)~ ., data)
cox$id <- "coxph"

formula <- Surv(time, status) ~ .
rfs <- rfsrc(formula, data = data)
rfs$id <- "rfs"


learners <- list()
learners[[1]] <- cox
learners[[2]] <- rfs
learners[[3]] <- svm
learners[[4]] <- gbm


#dataframe with overall perfomance
df_perf <- data.frame(
    config = character(),
    model = character(),
    C_Index = numeric(),
    Brier_score = numeric(),
    stringsAsFactors = FALSE
  )

boot_iters <- 100
max_resamp <- 900
boot_indices500 <- createResample(data$status, times = max_resamp)
boot_indices <- list()

for (k in 1:max_resamp){

    train_indices <- boot_indices500[[k]]

    train_data <- data[train_indices, ]
    if (OS) {result <- check_binary_columns(train_data[, -c(1)])}
    else {result <- check_binary_columns(train_data[, -c(11)])}

    if (all(result)) {
      boot_indices <- append(boot_indices, list(train_indices))
    }
  if (length(boot_indices)== 100){
    break
  }
}


for (learner in learners) {

  for (i in 1:nrow(flag_combs)) {

    current_flags <- flag_combs[i, ]

    if (any(current_flags)){

      RADIOMICS <- current_flags[[1]]
      CLINICAL <- current_flags[[2]]
      GENOMIC <- current_flags[[3]]

    #print(paste("Radiomics", RADIOMICS, "Clinical", CLINICAL, "Genomic", GENOMIC, sep=" "))
    #if (GENOMIC & !CLINICAL & !RADIOMICS & learner$id == "surv.gbm")
   # {
    #  print("Set minimun objs in Node = 3, for genomic train in GBM for OS task")
    #  learner <- lrn("surv.gbm", n.trees=100, interaction.depth = 10, bag.fraction=0.9, shrinkage=0.1,  n.minobsinnode=3)
    #}


    #returns in a list (data, configuration)
    if (OS) {
      dta_conf <- feature_sel_OS(data, RADIOMICS, CLINICAL, GENOMIC)
    } else{ dta_conf <- feature_sel_REC(data, RADIOMICS, CLINICAL, GENOMIC)}

    dta <- dta_conf[[1]]
    conf <- dta_conf[[2]]

    perfomance_list <- list()
    test_probs <- list()
    surv_func_list <- list()

    CI_list <- list()
    IBS_list <- list()

    paste("Training Model:", learner$id, sep = " ")

    avg_CI <- 0
    avg_AUC <- 0
    avg_Brier <- 0

    for (i in 1: boot_iters)
    {

      train_indices <- boot_indices[[i]]
      test_indices <- setdiff(seq_len(nrow(data)), train_indices)


      train_data <- dta[train_indices, ]
      test_data <- dta[-train_indices, ]



      start_time <- Sys.time()
      print(paste("Iteration n.", i, "of", boot_iters, sep = " "))


      if (learner$id == 'coxph')
        {
          model_name <- "coxPH"
          result <- withCallingHandlers(
            {
              learner <- coxph(Surv(train_data$time, train_data$status)~ ., train_data, x=TRUE, model=TRUE)
              TRUE
            },
            warning = function(w) {
              message("Warning caught: ", conditionMessage(w), "skipping model...")
              invokeRestart("muffleWarning")
              FALSE
            }
            )

          if(!result){
            next
          }


          learner$id <- "coxph"

          #risks <- predict(learner, type = "lp", new_data = test_data[, -c(1,2)])
          #test_CI <- c_index(survival::Surv(test_data$time, test_data$status), risk =risks)
          cox_explainer_test <- survex::explain(learner, data = test_data[,-c(1,2)],
                                        y = survival::Surv(test_data$time, test_data$status),
                                        verbose=FALSE)
          perf_success <- tryCatch(
                  {
                    perfList <- model_performance(cox_explainer_test)
                    TRUE
                  },
                    error = function(cond) {
                      message("Error Occurred during inference!")
                      message(conditionMessage(cond))
                      return(FALSE)
                    }
                  )
          if (!perf_success){
            next
          }

          CI_list <- append(CI_list, perfList$result[[1]])
          IBS_list <- append(IBS_list, perfList$result[[4]])
        }

        else if (learner$id == 'rfs'){

          model_name <- "survRF"
          formula <- Surv(time, status) ~ .
          learner <- rfsrc(formula, data = train_data,  ntree = 500)

          learner$id <- "rfs"
          rf_explainer_test <- survex::explain(learner,data = test_data[,-c(1,2)],
                                        y = survival::Surv(test_data$time, test_data$status),
                                        verbose=FALSE)

          perfList <- model_performance(rf_explainer_test)
          CI_list <- append(CI_list, perfList$result[[1]])
          IBS_list <- append(IBS_list, perfList$result[[4]])

        }

      else {
          train_task <- as_task_surv(x = train_data,
                                   time = "time",
                                   event = "status")


          composite_learner <- as_learner(ppl(
            "distrcompositor",
            learner = learner,
            estimator = "kaplan",
            form = "ph"))

          train_success <- tryCatch(
          {
            composite_learner$train(train_task)
            TRUE
          },
            error = function(cond) {
              message("Error Occurred during training!")
              message(conditionMessage(cond))
              return(FALSE)
            }
          )

          if (!train_success) {
          next
          }



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
                                            data = test_data[,-c(1,2)],
                                            y = survival::Surv(test_data$time, test_data$status),
                                            predict_survival_function = svm_predict,
                                            label = paste0(model_name, "Explainer", sep = " "),
                                            verbose=FALSE)

                  perfList <- model_performance(model_explainer)
                  CI_list <- append(CI_list, perfList$result[[1]])
                  IBS_list <- append(IBS_list, perfList$result[[4]])

          } else {

              model_name <- "survGBM"
              model_explainer <- survex::explain(composite_learner,
                                              data = test_data[,-c(1,2)],
                                              y = survival::Surv(test_data$time, test_data$status),
                                              label = paste0(model_name, "Explainer", sep = " "),
                                              verbose=FALSE)

              perfList <- model_performance(model_explainer)
              CI_list <- append(CI_list, perfList$result[[1]])
              IBS_list <- append(IBS_list, perfList$result[[4]])

            }

          #Performance on training
          #perf <- model_performance(model_explainer)

          #print(paste("CI", perf$result[[1]], "AUC", perf$result[[2]], "Brier", perf$result[[4]] ))
          #avg_CI <- avg_CI+perf$result[[1]]
          #avg_AUC <- avg_AUC+perf$result[[2]]
          #avg_Brier <- avg_Brier+perf$result[[4]]



          # print(paste0("Computing survshap..."))
          # survshap <- predict_parts(model_explainer,
          #                           new_observation = dta[i, -c(1,2)],
          #                           aggregation_method = "integral",
          #                           calculation_method = "kernelshap",
          #                           type="survshap")
          #
          # if (i==1) {
          #
          #     #creating only for having data structures, the values will be updated on each iteration
          #     print("Creating global survshap for all test explanations...")
          #
          #     eval_times <- survshap$eval_times
          #     surv_f_mtrx <- model_explainer$predict_survival_function(composite_learner, dta[, -c(1,2)], eval_times)
          #
          #     glob_surv <- predict_parts(model_explainer,
          #                           new_observation = dta[1:dim(dta)[1], -c(1,2)],
          #                           aggregation_method = "integral",
          #                           calculation_method = "kernelshap",
          #                           type="survshap")
          # }

          # if (model_name == 'survSVM'){
          #   new_data <- rbind(dta[i, -c(1,2)], dta[i, -c(1,2)])
          #   test_probs <- append(test_probs, composite_learner$predict_newdata(new_data)$crank[1])
          #
          # } else {
          #   test_probs <- append(test_probs, composite_learner$predict_newdata(dta[i, -c(1,2)])$crank)
          # }
          #
          # surv_f_mtrx[i,] <- model_explainer$predict_survival_function(composite_learner, dta[i, -c(1,2)], eval_times)


          #print(paste("Aggregating explanation of sample number", i, "...", sep = " "))
          #glob_surv[["aggregate"]][[i]] <- survshap[["aggregate"]]

          end_time <- Sys.time()
          time_elapsed <- end_time - start_time
          print(paste("Elapsed time for iteration", i, ":", time_elapsed, units(time_elapsed), sep = " "))
        }


        #test_CI <- c_index(survival::Surv(dta$time, dta$status), risk = unlist(test_probs))
        #test_IBS <- integrated_brier_score(survival::Surv(dta$time, dta$status), surv = surv_f_mtrx, times = eval_times)

    }
      avg_CI <- mean(unlist(CI_list))
      avg_IBS <- mean(unlist(IBS_list))

      avg_metrics <- list(avg_CI, avg_Brier)

      print(paste("Evaluation Metrics for", model_name, "...", sep= " "))
      print(paste("Average C-Index", avg_CI, "Average Brier", avg_Brier, sep= " "))



        #plot(glob_surv, subtitle = paste("Created for", model_name, "after LOOCV (m models X n patients, 60 X 60)", sep = " "))
        #plot(glob_surv, geom = "beeswarm", subtitle = paste("Created for", model_name, "after LOOCV (m models X n patients, 60 X 60)", sep = " "))

      save_ws(model_name, RADIOMICS, CLINICAL, GENOMIC, OS)
      df_current <- data.frame(
                  config = conf,
                  model = model_name,
                  C_Index = avg_CI,
                  Brier_score = avg_IBS,
                  stringsAsFactors = FALSE
      )
      df_perf <- rbind(df_perf, df_current)
      }
    }
  }

  if (OS) {
    write.csv(df_perf, './Radiomics_result_bootstrap/OS/Bootstrap_Performances_OS.csv', row.names = FALSE)
  } else { write.csv(df_perf, './Radiomics_result_bootstrap/REC/Bootstrap_Performances_REC.csv', row.names = FALSE)}





