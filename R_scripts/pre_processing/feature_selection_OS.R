# # # ##################################################################################################################################################
# # # # PACKAGES LOADING
# # # ##################################################################################################################################################

suppressPackageStartupMessages({
  # Load the installed Package
  library(survival)
  library(survex)
  library(randomForestSRC)
  library(survivalmodels)
  library(reticulate)
  library(survminer)
  library(stringr)
  library(tableone)
  library(forestmodel)
  library(patchwork)
  library(sjPlot)
})

#
#set path
setwd("/Users/gianmariazaccaria/Documents/R_studio/CPTAC-PDA/Integrated_analysis/test_8_RADIOsig+CLIN+MUT_250624_R")

# # # ##################################################################################################################################################
# # # # DATA PREPARATION
# # # ##################################################################################################################################################
# 
# #Loading of multi-modal dataset
# discovery <- read.csv(gsub(" ", "", paste(getwd(),"/multimodal_dataset_IVStages.csv")))
# 
# #Loading of the clinical dataset for adding OS time and other clinical features
# load("integrated_dataset_070224.Rda")
# clinical <- df_train[,c("Barcode","df_clinical.OS","df_follow_up.bmi_categorical","path_detail.lymph_nodes_positive","path_detail.tumor_largest_dimension_diameter_categorical","df_follow_up.ecog_performance_status","df_follow_up.karnofsky_performance_status","df_follow_up.disease_response")];
# remove(df_train)
# 
# #Include novel features in discovery-set
# discovery$time <- 0
# discovery$BMI <- "NULL"
# discovery$LN_pos <- 0
# discovery$tumor_largest_dia <- "NULL"
# discovery$ecogps <- "NULL"
# discovery$karn <- "NULL"
# discovery$clin_resp <- "NULL"
# 
# for(index_i in seq(1,dim(discovery)[1],1)){
#   for (index_j in seq(1,dim(clinical)[1],1)) {
#     if(str_equal(discovery$patient_id[index_i],clinical$Barcode[index_j])==TRUE){
#       discovery$time[index_i] <- clinical$df_clinical.OS[index_j]
#       discovery$BMI[index_i] <- clinical$df_follow_up.bmi_categorical[index_j]
#       discovery$LN_pos[index_i] <- clinical$path_detail.lymph_nodes_positive[index_j]
#       discovery$tumor_largest_dia[index_i] <- clinical$path_detail.tumor_largest_dimension_diameter_categorical[index_j]
#       discovery$ecogps[index_i] <- clinical$df_follow_up.ecog_performance_status[index_j]
#       discovery$karn[index_i] <- clinical$df_follow_up.karnofsky_performance_status[index_j]
#       discovery$clin_resp[index_i] <- clinical$df_follow_up.disease_response[index_j]
#     }
#   }
# }
# 
# #Remove clinical df
# remove(clinical)
# 
# #Preprocessing, transforming all categorical features in categorical numeric
# table(discovery$BMI)
# discovery$BMI_reclass = 0
# discovery$BMI_reclass[discovery$BMI == "OW"] = 1
# discovery$BMI_reclass[discovery$BMI == "NA"] = NA
# table(discovery$BMI_reclass)
# 
# table(discovery$LN_pos)
# discovery$LN_pos_reclass = 0
# discovery$LN_pos_reclass[discovery$LN_pos >= 1] = 1 #Verificare questo
# discovery$LN_pos_reclass[is.na(discovery$LN_pos)] = NA
# table(discovery$LN_pos_reclass)
# 
# table(discovery$tumor_largest_dia)
# discovery$tumor_largest_dia_reclass = 0
# discovery$tumor_largest_dia_reclass[discovery$tumor_largest_dia == ">=3.5cm"] = 1
# discovery$tumor_largest_dia_reclass[discovery$tumor_largest_dia == "NA"] = NA
# table(discovery$tumor_largest_dia_reclass)
# 
# table(discovery$ecogps)
# discovery$ecogps_reclass = 0
# discovery$ecogps_reclass[discovery$ecogps == "3" | discovery$ecogps == "4" | discovery$ecogps == "5"] = 1
# discovery$ecogps_reclass[discovery$ecogps == "Not Reported" | discovery$ecogps == "Unknown"] = NA
# table(discovery$ecogps_reclass)
# 
# table(discovery$karn)
# discovery$karn_reclass = 0
# discovery$karn_reclass[discovery$karn == "30" | discovery$karn == "40" | discovery$karn == "50" | discovery$karn == "60" | discovery$karn == "70" | discovery$karn == "80" | discovery$karn == "90" | discovery$karn == "100"] = 1
# discovery$karn_reclass[discovery$karn == "Not Reported" | discovery$karn == "Unknown"] = NA
# table(discovery$karn_reclass)
# 
# table(discovery$clin_resp)
# discovery$clin_resp_reclass = 0
# discovery$clin_resp_reclass[discovery$clin_resp == "CR-Complete Response" | discovery$clin_resp == "PR-Partial Response"] = 1
# discovery$clin_resp_reclass[discovery$clin_resp == "Not Reported" | discovery$clin_resp == "Unknown"] = NA
# table(discovery$clin_resp_reclass)
# 
# # ##################################################################################################################################################
# # # SURVIVAL ANALYSIS  - OS
# # ##################################################################################################################################################
# 
# #Standardizing outcome
# time_str = "time"
# event_str = "status"
# time = discovery$time
# event = discovery$status
# 
# ##################################################################################################################################################
# # RADIOMICS FEATURES - SURVIVAL ANALYSIS
# ##################################################################################################################################################
# 
# # Create dataframe of results
# df_res_radio = data.frame(matrix(nrow = 16, ncol = 7))
# colnames(df_res_radio) <- c("variables","pval","method","pval.txt","cut_off","C-index","std")
# df_res_radio[,1] <- colnames(discovery[,2:17]);
# 
# for (i in colnames(discovery[,2:17])) {
# 
#   #Cut-off and survival analysis
#   res.cut <- surv_cutpoint(discovery, time = time_str, event = event_str, variables = i); cutoff <- res.cut$cutpoint[1,1]
#   summary(res.cut)
# 
#   # pdf(file = gsub(" ", "", paste("/Users/gianmariazaccaria/Documents/R_studio/CPTAC-PDA/Integrated_analysis/test_4_entire_dataset_progression_020224/results/UV_analysis/radiomics/",i,"_MSRS.pdf")))
#   # plot(res.cut,i, palette = "npg")
#   # dev.off()
# 
#   res.cat <- surv_categorize(res.cut)
# 
#   #Aggiungo valore categorizzato della variabile al dataframe
#   string <- gsub(" ", "", paste(i,"_cate"))    #Creo nome colonna
#   discovery[string] <- res.cat[,3]
# 
#   fit <- survfit(Surv(time, event) ~res.cat[,3], data = res.cat); fit
#   pp = ggsurvplot(fit, data=res.cat, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                   palette= c("red","blue"), legend.labs = c("high","low"),
#                   legend="top", risk.table.title="N. at risk", title=i, font.x=12,font.y=12,
#                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2);
# 
#   pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/radiomics/OS/",i,"_KM.pdf")))
#   print(pp, newpage = FALSE)
#   dev.off()
# 
#   cox <- coxph(formula=Surv(time,event) ~ res.cat[,3], data=res.cat); summary(cox)
# 
#   df_res_radio[df_res_radio$variables == i,2:4] <- surv_pvalue(fit)[,2:4]
#   df_res_radio[df_res_radio$variables == i,5] <- cutoff
#   df_res_radio[df_res_radio$variables == i,6] <- round(cox$concordance["concordance"], 2);
#   df_res_radio[df_res_radio$variables == i,7] <- round(cox$concordance["std"], 2);
# }
# 
# #Salvo nomi variabili categorizzate
# radiomics_sign <- df_res_radio$variables[df_res_radio$pval < 0.15]; radiomics_sign
# 
##################################################################################################################################################
# CLINICAL FEATURES - SURVIVAL ANALYSIS
##################################################################################################################################################
# 
# # # #Definition of lists for each panel construction
# # # splots_clinical <- list(); splots_exposure <- list();
# #
# # df_clinical.gender
# table(discovery$df_clinical.gender)
# name <- ("Gender")
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(df_clinical.gender),data=discovery); sopravv
# pp = ggsurvplot(sopravv, data=discovery, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                 palette= c("red","blue"), legend.labs = c("female","male"),
#                                 legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time, event) ~ as.factor(df_clinical.gender), data=discovery); summary(cox)
# 
# pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/OS/",name,"_KM.pdf")))
# print(pp, newpage = FALSE)
# dev.off()
# 
# # # discovery$df_clinical.age_categorical
# table(discovery$df_clinical.age_categorical)
# name <- ("Age at diag")
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(discovery$df_clinical.age_categorical), data=discovery); sopravv
# pp=ggsurvplot(sopravv, data=discovery, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                 palette= c("red","blue"), legend.labs = c("Low Age","High Age"),
#                                 legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time, event) ~ as.factor(discovery$df_clinical.age_categorical), data=discovery); summary(cox)
# 
# pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/OS/",name,"_KM.pdf")))
# print(pp, newpage = FALSE)
# dev.off()
# 
# # # BMI
# table(discovery$BMI_reclass)
# df_2 <- subset(discovery, !is.na(BMI_reclass))
# table(df_2$BMI_reclass)
# time2 = df_2$time
# event2 = df_2$status
# 
# name <- ("BMI")
# sopravv <- survfit(Surv(time2, as.numeric(event2))~as.factor(df_2$BMI_reclass), data=df_2); sopravv
# pp=ggsurvplot(sopravv, data=df_2, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                 palette= c("red","blue"), legend.labs = c("Under Weight","Over Weight"),
#                                 legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time2, event2) ~ as.factor(df_2$BMI_reclass), data=df_2); summary(cox)
# 
# pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/OS/",name,"_KM.pdf")))
# print(pp, newpage = FALSE)
# dev.off()
# 
# discovery$df_clinical.ajcc_pathologic_stage
# table(discovery$df_clinical.ajcc_pathologic_stage)
# name <- ("Stage")
# 
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(df_clinical.ajcc_pathologic_stage_reclass_2),data=discovery); sopravv
# pp=ggsurvplot(sopravv, data=discovery, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#               risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#               palette= c("red","blue"), legend.labs = c("Stage IB/IIA/IIB","Stage III/IIIB/IV"),
#               legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#               font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# 
# discovery$df_clinical.ajcc_pathologic_stage
# table(discovery$df_clinical.ajcc_pathologic_stage_reclass_2)
# name <- ("Stage")
# 
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(df_clinical.ajcc_pathologic_stage_reclass_2),data=discovery); sopravv
# pp=ggsurvplot(sopravv, data=discovery, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                 palette= c("red","blue"), legend.labs = c("Stage IB/IIA/IIB","Stage III/IIIB/IV"),
#                                 legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time, event) ~ as.factor(df_clinical.ajcc_pathologic_stage_reclass_2), data=discovery); summary(cox)
# 
# pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/OS/",name,"_KM.pdf")))
# print(pp, newpage = FALSE)
# dev.off()
# 
# df_clinical.tumor_grade
# table(discovery$df_clinical.tumor_grade)
# name <- ("Tumor Grade")
# 
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(discovery$df_clinical.tumor_grade_reclass),data=discovery); sopravv
# pp=ggsurvplot(sopravv, data=discovery, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#               risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#               palette= c("red","blue"), legend.labs = c("G1/G2","G3/G4"),
#               legend="top", risk.table.title="N. at risk", title=name, font.x=16,font.y=16,
#               font.tickslab=14, font.legend=16, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time, event) ~ discovery$df_clinical.tumor_grade_reclass, data=discovery); summary(cox)
# 
# pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/OS/",name,"_KM.pdf")))
# print(pp, newpage = FALSE)
# dev.off()
# 
# # discovery$LN_pos_reclass
# table(discovery$LN_pos_reclass);
# name <- ("Lymphnode Involvement")
# 
# df_2 <- subset(discovery, !is.na(discovery$LN_pos_reclass))
# 
# time2 = df_2$time
# event2 = df_2$status
# 
# sopravv <- survfit(Surv(time2, as.numeric(event2))~as.factor(df_2$LN_pos_reclass),data=df_2); sopravv
# pp=ggsurvplot(sopravv, data=df_2, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                 palette= c("red","blue"), legend.labs = c("0", "at least 1"),
#                                 legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time2, event2) ~ as.factor(df_2$LN_pos_reclass), data=df_2); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/OS/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # discovery$path_detail.tumor_largest_dimension_diameter_categorical
# table(discovery$tumor_largest_dia_reclass)
# name <- ("Tumor diameter")
# 
# df_2 <- subset(discovery, !is.na(discovery$tumor_largest_dia_reclass))
# time2 = df_2$time
# event2 = df_2$status
# 
# sopravv <- survfit(Surv(time2, as.numeric(event2))~as.factor(df_2$tumor_largest_dia_reclass),data=df_2); sopravv
# pp=ggsurvplot(sopravv, data=df_2, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                 palette= c("red","blue"), legend.labs = c("lower 3.5", "higher 3.5"),
#                                 legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time2, event2) ~ as.factor(df_2$tumor_largest_dia_reclass), data=df_2); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/OS/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # df_follow_up.ecog_performance_status
# table(discovery$ecogps_reclass); #freq(discovery$df_follow_up.ecog_performance_status)
# name <- ("ECOGps")
# 
# df_2 <- subset(discovery, !is.na(discovery$ecogps_reclass))
# time2 = df_2$time
# event2 = df_2$status
# 
# sopravv <- survfit(Surv(time2, as.numeric(event2))~as.factor(df_2$ecogps_reclass),data=df_2); sopravv
# pp=ggsurvplot(sopravv, data=df_2, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                 palette= c("red","blue"), legend.labs = c("0/1/2/3", "4/5"),
#                                 legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time2, event2) ~ as.factor(df_2$ecogps_reclass), data=df_2); summary(cox)
# 
# pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/OS/",name,"_KM.pdf")))
# print(pp, newpage = FALSE)
# dev.off()
# 
# df_follow_up.karnofsky_performance_status
# table(discovery$karn_reclass);
# name <- ("Karnofsky")
# 
# df_2 <- subset(discovery, !is.na(discovery$karn_reclass))
# time2 = df_2$time
# event2 = df_2$status
# 
# sopravv <- survfit(Surv(time2, as.numeric(event2))~as.factor(df_2$karn_reclass),data=df_2); sopravv
# pp=ggsurvplot(sopravv, data=df_2, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                 palette= c("red","blue"), legend.labs = c( "0", "higher then 0%"),
#                                 legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time2, event2) ~ as.factor(df_2$karn_reclass), data=df_2); summary(cox)
# 
# pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/OS/",name,"_KM.pdf")))
# print(pp, newpage = FALSE)
# dev.off()
# 
# discovery$df_follow_up.disease_response
# table(discovery$clin_resp_reclass)
# name <- ("Clinical Response")
# 
# df_2 <- subset(discovery, !is.na(discovery$clin_resp_reclass))
# time2 = df_2$time
# event2 = df_2$status
# 
# sopravv <- survfit(Surv(time2, as.numeric(event2))~as.factor(df_2$clin_resp_reclass),data=df_2); sopravv
# pp <- ggsurvplot(sopravv, data=df_2, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                     risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                     palette= c("red","blue"), legend.labs = c("Not Responders", "Responders"),
#                                     legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                     font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time2, event2) ~ as.factor(df_2$clin_resp_reclass), data=df_2); summary(cox)
# 
# pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/OS/",name,"_KM.pdf")))
# print(pp, newpage = FALSE)
# dev.off()
# 
# df_clinical.residual_disease
# table(discovery$df_clinical.residual_disease_reclass)
# 
# name <- "Residual Disease Reclass"
# sopravv <- survfit(Surv(time, status)~as.factor(discovery$df_clinical.residual_disease_reclass),data=discovery); sopravv
# pp=ggsurvplot(sopravv, data=discovery, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                  risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                  palette= c("red","blue","darkgreen"), legend.labs = c("R1/R2", "R0","RX"),
#                                  legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                  font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time, status) ~ as.factor(discovery$df_clinical.residual_disease_reclass), data=discovery); summary(cox)
# 
# pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/OS/",name,"_KM.pdf")))
# print(pp, newpage = FALSE)
# dev.off()
# 
# # ##################################################################################################################################################
# # # MUTATIONAL FEATURES - SURVIVAL ANALYSIS
# # ##################################################################################################################################################
# 
# #Definition of lists for each panel construction
# # splots_mutational <- list();
# #
# # # KRAS
# table(discovery$KRAS)
# name <- "KRAS"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(KRAS),data=discovery); sopravv
# pp=ggsurvplot(sopravv, data=discovery, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                   palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                   legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(KRAS), data=discovery); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/OS/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # # TP53
# table(discovery$TP53)
# name <- "TP53"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(TP53),data=discovery); sopravv
# pp=ggsurvplot(sopravv, data=discovery, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                   palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                   legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(TP53), data=discovery); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/OS/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # CDKN2A
# table(discovery$CDKN2A)
# name <- "CDKN2A"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(CDKN2A),data=discovery); sopravv
# pp=ggsurvplot(sopravv, data=discovery, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                   palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                   legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(CDKN2A), data=discovery); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/OS/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # SMAD4
# table(discovery$SMAD4)
# name <- "SMAD4"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(SMAD4),data=discovery); sopravv
# pp=ggsurvplot(sopravv, data=discovery, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                   palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                   legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(SMAD4), data=discovery); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/OS/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # TTN
# table(discovery$TTN)
# name <- "TTN"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(TTN),data=discovery); sopravv
# pp=ggsurvplot(sopravv, data=discovery, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                   palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                   legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(TTN), data=discovery); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/OS/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# #MUC16
# table(discovery$MUC16)
# name <- "MUC16"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(MUC16),data=discovery); sopravv
# pp=ggsurvplot(sopravv, data=discovery, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                   palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                   legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(MUC16), data=discovery); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/OS/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# #KMT2D
# table(discovery$KMT2D)
# name <- "KMT2D"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(KMT2D),data=discovery); sopravv
# pp=ggsurvplot(sopravv, data=discovery, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                   palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                   legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(KMT2D), data=discovery); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/OS/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # CSMD1
# table(discovery$CSMD1)
# name <- "CSMD1"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(CSMD1),data=discovery); sopravv
# pp=ggsurvplot(sopravv, data=discovery, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                   palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                   legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(CSMD1), data=discovery); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/OS/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # RYR2
# table(discovery$RYR2)
# name <- "RYR2"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(RYR2),data=discovery); sopravv
# pp=ggsurvplot(sopravv, data=discovery, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                   palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                   legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(RYR2), data=discovery); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/OS/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # RYR1
# table(discovery$RYR1)
# name <- "RYR1"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(RYR1),data=discovery); sopravv
# pp=ggsurvplot(sopravv, data=discovery, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability OS",
#                                    risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                    palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                    legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                    font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(RYR1), data=discovery); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/OS/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # ##################################################################################################################################################
# # # SAVING DATAFRAMES
# # ##################################################################################################################################################
# 
# save(discovery,file="discovery_set_OS.Rda")
# save(df_res_radio,file="df_res_radio_discovery_set_OS.Rda")
# 
# ##################################################################################################################################################
# # LOADINING DATAFRAMES
# ##################################################################################################################################################

load(gsub(" ", "", paste(getwd(),"/results/discovery_set_OS.Rda")))
load(gsub(" ", "", paste(getwd(),"/results/df_res_radio_discovery_set_OS.Rda")))

# # # ##################################################################################################################################################
# # # # PATIENTS' CHARACTERISTICS
# # # ##################################################################################################################################################
# 
# #Including features
# discovery_patchar <- discovery[,c("df_clinical.age_categorical", "df_clinical.gender", "df_clinical.ajcc_pathologic_stage_reclass_2","df_clinical.tumor_grade_reclass",
#                    "KRAS","TP53", "CDKN2A", "SMAD4", "TTN", "MUC16", "KMT2D", "CSMD1", "RYR2","RYR1", "status", "BMI_reclass", "LN_pos_reclass",
#                    "tumor_largest_dia_reclass", "ecogps_reclass", "karn_reclass", "clin_resp_reclass", "original_firstorder_Kurtosis_cate",
#                    "original_glcm_Correlation_cate", "log.sigma.1.0.mm.3D_gldm_DependenceVariance_cate", "original_shape_Elongation_cate",
#                    "original_glcm_JointEnergy_cate", "log.sigma.1.0.mm.3D_firstorder_Skewness_cate", "log.sigma.2.0.mm.3D_glszm_SmallAreaLowGrayLevelEmphasis_cate",
#                    "wavelet.LLH_glcm_ClusterTendency_cate", "wavelet.HLH_glcm_Imc1_cate", "original_shape_LeastAxisLength_cate", "wavelet.HHL_firstorder_Mean_cate",
#                    "wavelet.HHL_glcm_Imc1_cate", "wavelet.HHL_glszm_LargeAreaEmphasis_cate", "log.sigma.2.0.mm.3D_firstorder_Median_cate", "wavelet.HLL_glcm_InverseVariance_cate",
#                    "original_firstorder_90Percentile_cate", "df_clinical.residual_disease_reclass")]
# 
# #Rename features
# colnames(discovery_patchar)[1] <- c("Age at diagnosis")
# colnames(discovery_patchar)[2] <- c("Gender")
# colnames(discovery_patchar)[3] <- c("Stage")
# colnames(discovery_patchar)[4] <- c("Grade")
# colnames(discovery_patchar)[15] <- c("OS status")
# colnames(discovery_patchar)[16] <- c("BMI")
# colnames(discovery_patchar)[17] <- c("Lymphnode involvement")
# colnames(discovery_patchar)[18] <- c("Tumor largest diameter")
# colnames(discovery_patchar)[19] <- c("ECOGps")
# colnames(discovery_patchar)[20] <- c("Karnofsky")
# colnames(discovery_patchar)[21] <- c("Clinical response")
# colnames(discovery_patchar)[22] <- c("Original firstorder Kurtosis")
# colnames(discovery_patchar)[37] <- c("Original firstorder 90 Percentile")
# colnames(discovery_patchar)[23] <- c("Original glcm Correlation")
# colnames(discovery_patchar)[26] <- c("Original glcm Joint Energy")
# colnames(discovery_patchar)[25] <- c("Original shape Elongation")
# colnames(discovery_patchar)[31] <- c("Original shape Least Axis Length")
# colnames(discovery_patchar)[24] <- c("loG.sigma.1.0.mm.3D gldm Dependence Variance")
# colnames(discovery_patchar)[27] <- c("loG.sigma.1.0.mm.3D firstorder Skewness")
# colnames(discovery_patchar)[28] <- c("loG.sigma.2.0.mm.3D glszm Small Area Low Gray Level Emphasis")
# colnames(discovery_patchar)[35] <- c("loG.sigma.2.0.mm.3D firstorder Median")
# colnames(discovery_patchar)[29] <- c("Wavelet.LLH glcm Cluster Tendency")
# colnames(discovery_patchar)[30] <- c("Wavelet.HLH glcm Imc1")
# colnames(discovery_patchar)[36] <- c("Wavelet.HLL glcm Inverse Variance")
# colnames(discovery_patchar)[32] <- c("Wavelet.HHL firstorder Mean")
# colnames(discovery_patchar)[33] <- c("Wavelet.HHL glcm Imc1")
# colnames(discovery_patchar)[34] <- c("Wavelet.HHL glszm Large Area Emphasis")
# colnames(discovery_patchar)[38] <- c("Residual Disease")
# 
# summary(discovery_patchar)
# 
# #Rename classes and rearrange levels
# table(discovery_patchar$`Age at diagnosis`)
# discovery_patchar$`Age at diagnosis`[discovery_patchar$`Age at diagnosis` == 0] <- "<65 ys"
# discovery_patchar$`Age at diagnosis`[discovery_patchar$`Age at diagnosis` == 1] <- ">=65 ys"
# discovery_patchar$`Age at diagnosis` <- as.factor(discovery_patchar$`Age at diagnosis`)
# summary(discovery_patchar$`Age at diagnosis`)
# 
# table(discovery_patchar$Gender)
# discovery_patchar$Gender[discovery_patchar$Gender == 0] <- "Female"
# discovery_patchar$Gender[discovery_patchar$Gender == 1] <- "Male"
# discovery_patchar$Gender <- as.factor(discovery_patchar$Gender)
# discovery_patchar$Gender <- factor(discovery_patchar$Gender, levels=c('Male', 'Female'))
# summary(discovery_patchar$Gender)
# 
# table(discovery_patchar$Stage)
# discovery_patchar$Stage[discovery_patchar$Stage == 0] <- "IB/IIA/IIB"
# discovery_patchar$Stage[discovery_patchar$Stage == 1] <- "III/IIIB/IV"
# discovery_patchar$Stage <- as.factor(discovery_patchar$Stage)
# summary(discovery_patchar$Stage)
# 
# table(discovery_patchar$Grade)
# discovery_patchar$Grade[discovery_patchar$Grade == 0] <- "G1/G2"
# discovery_patchar$Grade[discovery_patchar$Grade == 1] <- "G3/G4"
# discovery_patchar$Grade <- as.factor(discovery_patchar$Grade)
# summary(discovery_patchar$Grade)
# 
# table(discovery_patchar$`OS status`)
# discovery_patchar$`OS status`[discovery_patchar$`OS status` == 0] <- "Alive"
# discovery_patchar$`OS status`[discovery_patchar$`OS status` == 1] <- "Dead"
# discovery_patchar$`OS status` <- as.factor(discovery_patchar$`OS status`)
# summary(discovery_patchar$`OS status`)
# 
# table(discovery_patchar$BMI)
# discovery_patchar$BMI[discovery_patchar$BMI == 0] <- "<25"
# discovery_patchar$BMI[discovery_patchar$BMI == 1] <- ">=25"
# discovery_patchar$BMI[is.na(discovery_patchar$BMI)] = "NA"
# discovery_patchar$BMI <- as.factor(discovery_patchar$BMI)
# summary(discovery_patchar$BMI)
# 
# table(discovery_patchar$`Lymphnode involvement`)
# discovery_patchar$`Lymphnode involvement`[discovery_patchar$`Lymphnode involvement` == 0] <- "0 nodes"
# discovery_patchar$`Lymphnode involvement`[discovery_patchar$`Lymphnode involvement` == 1] <- ">0 nodes"
# discovery_patchar$`Lymphnode involvement`[is.na(discovery_patchar$`Lymphnode involvement`)] = "NA"
# discovery_patchar$`Lymphnode involvement` <- as.factor(discovery_patchar$`Lymphnode involvement`)
# summary(discovery_patchar$`Lymphnode involvement`)
# 
# table(discovery_patchar$`Tumor largest diameter`)
# discovery_patchar$`Tumor largest diameter`[discovery_patchar$`Tumor largest diameter` == 0] <- "<3.5 cm"
# discovery_patchar$`Tumor largest diameter`[discovery_patchar$`Tumor largest diameter` == 1] <- ">=3.5 cm"
# discovery_patchar$`Tumor largest diameter`[is.na(discovery_patchar$`Tumor largest diameter`)] = "NA"
# discovery_patchar$`Tumor largest diameter` <- as.factor(discovery_patchar$`Tumor largest diameter`)
# summary(discovery_patchar$`Tumor largest diameter`)
# 
# table(discovery_patchar$ECOGps)
# discovery_patchar$ECOGps[discovery_patchar$ECOGps == 0] <- "0-3"
# discovery_patchar$ECOGps[discovery_patchar$ECOGps == 1] <- "4-5"
# discovery_patchar$ECOGps[is.na(discovery_patchar$ECOGps)] = "NA"
# discovery_patchar$ECOGps <- as.factor(discovery_patchar$ECOGps)
# summary(discovery_patchar$ECOGps)
# 
# table(discovery_patchar$Karnofsky)
# discovery_patchar$Karnofsky[discovery_patchar$Karnofsky == 0] <- "0 %"
# discovery_patchar$Karnofsky[discovery_patchar$Karnofsky == 1] <- ">0 %"
# discovery_patchar$Karnofsky[is.na(discovery_patchar$Karnofsky)] = "NA"
# discovery_patchar$Karnofsky <- as.factor(discovery_patchar$Karnofsky)
# summary(discovery_patchar$Karnofsky)
# 
# table(discovery_patchar$`Clinical response`)
# discovery_patchar$`Clinical response`[discovery_patchar$`Clinical response` == 0] <- "Not responders"
# discovery_patchar$`Clinical response`[discovery_patchar$`Clinical response` == 1] <- "Responders"
# discovery_patchar$`Clinical response`[is.na(discovery_patchar$`Clinical response`)] = "NA"
# discovery_patchar$`Clinical response` <- factor(discovery_patchar$`Clinical response`, levels=c('Responders', 'Not responders','NA'))
# summary(discovery_patchar$`Clinical response`)
# 
# table(discovery_patchar$`Residual Disease`)
# discovery_patchar$`Residual Disease`[discovery_patchar$`Residual Disease` == 0] <- "R0"
# discovery_patchar$`Residual Disease`[discovery_patchar$`Residual Disease` == 1] <- "R1/R2"
# discovery_patchar$`Residual Disease`[discovery_patchar$`Residual Disease` == 3] <- "RX"
# discovery_patchar$`Residual Disease`[is.na(discovery_patchar$`Residual Disease`)] = "NA"
# discovery_patchar$`Residual Disease` <- as.factor(discovery_patchar$`Residual Disease`)
# table(discovery_patchar$`Residual Disease`); summary(discovery_patchar$`Residual Disease`)
# 
# discovery_patchar[discovery_patchar == "high"] = "High"
# discovery_patchar[discovery_patchar == "low"] = "Low"
# discovery_patchar[is.na(discovery_patchar)] = "NA"
# 
# discovery_patchar[discovery_patchar == 0] = "Not mutated"
# discovery_patchar[discovery_patchar == 1] = "Mutated"
# 
# table(discovery_patchar$KRAS)
# discovery_patchar$KRAS <- factor(discovery_patchar$KRAS, levels=c('Not mutated', 'Mutated'))
# table(discovery_patchar$KRAS)
# 
# table(discovery_patchar$TP53)
# discovery_patchar$TP53 <- factor(discovery_patchar$TP53, levels=c('Not mutated', 'Mutated'))
# table(discovery_patchar$TP53)
# 
# table(discovery_patchar$CDKN2A)
# discovery_patchar$CDKN2A <- factor(discovery_patchar$CDKN2A, levels=c('Not mutated', 'Mutated'))
# table(discovery_patchar$CDKN2A)
# 
# table(discovery_patchar$SMAD4)
# discovery_patchar$SMAD4 <- factor(discovery_patchar$SMAD4, levels=c('Not mutated', 'Mutated'))
# table(discovery_patchar$SMAD4)
# 
# table(discovery_patchar$TTN)
# discovery_patchar$TTN <- factor(discovery_patchar$TTN, levels=c('Not mutated', 'Mutated'))
# table(discovery_patchar$TTN)
# 
# table(discovery_patchar$MUC16)
# discovery_patchar$MUC16 <- factor(discovery_patchar$MUC16, levels=c('Not mutated', 'Mutated'))
# table(discovery_patchar$MUC16)
# 
# table(discovery_patchar$KMT2D)
# discovery_patchar$KMT2D <- factor(discovery_patchar$KMT2D, levels=c('Not mutated', 'Mutated'))
# table(discovery_patchar$KMT2D)
# 
# table(discovery_patchar$CSMD1)
# discovery_patchar$CSMD1 <- factor(discovery_patchar$CSMD1, levels=c('Not mutated', 'Mutated'))
# table(discovery_patchar$CSMD1)
# 
# table(discovery_patchar$RYR2)
# discovery_patchar$RYR2 <- factor(discovery_patchar$RYR2, levels=c('Not mutated', 'Mutated'))
# table(discovery_patchar$RYR2)
# 
# table(discovery_patchar$RYR1)
# discovery_patchar$RYR1 <- factor(discovery_patchar$RYR1, levels=c('Not mutated', 'Mutated'))
# table(discovery_patchar$RYR1)
# 
# table(discovery_patchar$`Original firstorder Kurtosis`)
# discovery_patchar$`Original firstorder Kurtosis` <- factor(discovery_patchar$`Original firstorder Kurtosis`, levels=c('Low', 'High'))
# table(discovery_patchar$`Original firstorder Kurtosis`)
# 
# table(discovery_patchar$`Original firstorder 90 Percentile`)
# discovery_patchar$`Original firstorder 90 Percentile` <- factor(discovery_patchar$`Original firstorder 90 Percentile`, levels=c('Low', 'High'))
# table(discovery_patchar$`Original firstorder 90 Percentile`)
# 
# table(discovery_patchar$`Original glcm Correlation`)
# discovery_patchar$`Original glcm Correlation` <- factor(discovery_patchar$`Original glcm Correlation`, levels=c('Low', 'High'))
# table(discovery_patchar$`Original glcm Correlation`)
# 
# table(discovery_patchar$`Original shape Elongation`)
# discovery_patchar$`Original shape Elongation` <- factor(discovery_patchar$`Original shape Elongation`, levels=c('Low', 'High'))
# table(discovery_patchar$`Original shape Elongation`)
# 
# table(discovery_patchar$`Original glcm Joint Energy`)
# discovery_patchar$`Original glcm Joint Energy` <- factor(discovery_patchar$`Original glcm Joint Energy`, levels=c('Low', 'High'))
# table(discovery_patchar$`Original glcm Joint Energy`)
# 
# table(discovery_patchar$`loG.sigma.1.0.mm.3D gldm Dependence Variance`)
# discovery_patchar$`loG.sigma.1.0.mm.3D gldm Dependence Variance` <- factor(discovery_patchar$`loG.sigma.1.0.mm.3D gldm Dependence Variance`, levels=c('Low', 'High'))
# table(discovery_patchar$`loG.sigma.1.0.mm.3D gldm Dependence Variance`)
# 
# table(discovery_patchar$`loG.sigma.1.0.mm.3D firstorder Skewness`)
# discovery_patchar$`loG.sigma.1.0.mm.3D firstorder Skewness` <- factor(discovery_patchar$`loG.sigma.1.0.mm.3D firstorder Skewness`, levels=c('Low', 'High'))
# table(discovery_patchar$`loG.sigma.1.0.mm.3D firstorder Skewness`)
# 
# table(discovery_patchar$`loG.sigma.2.0.mm.3D glszm Small Area Low Gray Level Emphasis`)
# discovery_patchar$`loG.sigma.2.0.mm.3D glszm Small Area Low Gray Level Emphasis` <- factor(discovery_patchar$`loG.sigma.2.0.mm.3D glszm Small Area Low Gray Level Emphasis`, levels=c('Low', 'High'))
# table(discovery_patchar$`loG.sigma.2.0.mm.3D glszm Small Area Low Gray Level Emphasis`)
# 
# table(discovery_patchar$`Original shape Least Axis Length`)
# discovery_patchar$`Original shape Least Axis Length` <- factor(discovery_patchar$`Original shape Least Axis Length`, levels=c('Low', 'High'))
# table(discovery_patchar$`Original shape Least Axis Length`)
# 
# table(discovery_patchar$`loG.sigma.2.0.mm.3D firstorder Median`)
# discovery_patchar$`loG.sigma.2.0.mm.3D firstorder Median` <- factor(discovery_patchar$`loG.sigma.2.0.mm.3D firstorder Median`, levels=c('Low', 'High'))
# table(discovery_patchar$`loG.sigma.2.0.mm.3D firstorder Median`)
# 
# table(discovery_patchar$`Wavelet.LLH glcm Cluster Tendency`)
# discovery_patchar$`Wavelet.LLH glcm Cluster Tendency` <- factor(discovery_patchar$`Wavelet.LLH glcm Cluster Tendency`, levels=c('Low', 'High'))
# table(discovery_patchar$`Wavelet.LLH glcm Cluster Tendency`)
# 
# table(discovery_patchar$`Wavelet.HLL glcm Inverse Variance`)
# discovery_patchar$`Wavelet.HLL glcm Inverse Variance` <- factor(discovery_patchar$`Wavelet.HLL glcm Inverse Variance`, levels=c('Low', 'High'))
# table(discovery_patchar$`Wavelet.HLL glcm Inverse Variance`)
# 
# table(discovery_patchar$`Wavelet.HLH glcm Imc1`)
# discovery_patchar$`Wavelet.HLH glcm Imc1` <- factor(discovery_patchar$`Wavelet.HLH glcm Imc1`, levels=c('Low', 'High'))
# table(discovery_patchar$`Wavelet.HLH glcm Imc1`)
# 
# table(discovery_patchar$`Wavelet.HHL firstorder Mean`)
# discovery_patchar$`Wavelet.HHL firstorder Mean` <- factor(discovery_patchar$`Wavelet.HHL firstorder Mean`, levels=c('Low', 'High'))
# table(discovery_patchar$`Wavelet.HHL firstorder Mean`)
# 
# table(discovery_patchar$`Wavelet.HHL glcm Imc1`)
# discovery_patchar$`Wavelet.HHL glcm Imc1` <- factor(discovery_patchar$`Wavelet.HHL glcm Imc1`, levels=c('Low', 'High'))
# table(discovery_patchar$`Wavelet.HHL glcm Imc1`)
# 
# table(discovery_patchar$`Wavelet.HHL glszm Large Area Emphasis`)
# discovery_patchar$`Wavelet.HHL glszm Large Area Emphasis` <- factor(discovery_patchar$`Wavelet.HHL glszm Large Area Emphasis`, levels=c('Low', 'High'))
# table(discovery_patchar$`Wavelet.HHL glszm Large Area Emphasis`)
# 
# #Reordering columns
# discovery_patchar <- discovery_patchar[, c(22, 37, 23, 25, 26, 31, 24, 27, 28, 35, 29, 30, 36, 32, 33, 34, 1, 2, 3, 4, 16, 17, 18, 19, 20, 21, 38, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)]
# 
# ## Vector of variables to summarize
# myVars <- dput(names(discovery_patchar))
# 
# ## Vector of categorical variables that need transformation
# catVars <- dput(names(discovery_patchar))
# 
# ## Create a TableOne object
# tab2 <- CreateTableOne(vars = myVars, data = discovery_patchar)
# tab2
# summary(tab2)
# tab2Mat <- print(tab2, showAllLevels = TRUE, formatOptions = list(big.mark = ","))
# 
# ## Save to a CSV file
# write.csv(tab2Mat, file = gsub(" ","",paste(getwd(),"/results/pat_char/pat_char_overall_discovery.csv")))
# 
# #Stratified groups
# tab3 <- CreateTableOne(vars = myVars, strata = "OS status" , data = discovery_patchar)
# print(tab3, formatOptions = list(big.mark = ","))
# 
# #Exporting tableone
# tab3Mat <- print(tab3, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

# ## Save to a CSV file
# write.csv(tab3Mat, file = gsub(" ","",paste(getwd(),"/results/pat_char/pat_char_OS_discovery.csv")))
#
##################################################################################################################################################
# UNIVARIATE KAPLAN-MEYER IMAGES EXPORT
##################################################################################################################################################

centr = theme_classic(base_family = 'Helvetica') + theme(plot.title = element_text(hjust = 0.5, size = 14))

splots_OS <- list();
# # original_firstorder_Kurtosis_cate
# sopravv <- survfit(Surv(time,status)~as.factor(original_firstorder_Kurtosis_cate),data=discovery); sopravv
# splots_radiomics_OS[[1]] <- ggsurvplot(sopravv, data=discovery, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability OS",
#                              risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
#                              palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE,
#                              title="Original firstorder Kurtosis",
#                              legend="none", font.x=12,font.y=12, pval.coord = c(0,0.1), font.title = 12, font.tickslab=12,pval.size=4,risk.table.height = FALSE, ggtheme = centr); splots_radiomics_OS[[1]]

# original_firstorder_90Percentile_cate
sopravv <- survfit(Surv(time,status)~as.factor(original_firstorder_90Percentile_cate),data=discovery); sopravv
splots_OS[[1]] <- ggsurvplot(sopravv, data=discovery, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE, ylab="Probability OS",
                                       risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
                                       palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE,
                                       title="Original Fo 90P", font.title=18,
                                       legend="top", legend.labs=c("high","low"),
                                       font.x=14,font.y=14, pval.coord = c(0,0.1), font.tickslab=14,font.legend=12,pval.size=5, risk.table.height = FALSE, ggtheme = centr); splots_OS[[1]]

# original_shape_Elongation_cate
sopravv <- survfit(Surv(time,status)~as.factor(original_shape_Elongation_cate),data=discovery); sopravv
splots_OS[[2]] <- ggsurvplot(sopravv, data=discovery, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability OS",
                                       risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
                                       palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE,
                                       title="Original Shape Elong", font.title=18,
                                       legend="top", legend.labs=c("high","low"),
                                       font.x=14,font.y=14, pval.coord = c(0,0.1), font.tickslab=14,font.legend=12,pval.size=5, risk.table.height = FALSE, ggtheme = centr); splots_OS[[2]]

# original_glcm_JointEnergy_cate
sopravv <- survfit(Surv(time,status)~as.factor(original_glcm_JointEnergy_cate),data=discovery); sopravv
splots_OS[[3]] <- ggsurvplot(sopravv, data=discovery, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability OS",
                                       risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
                                       palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE,
                                       title="Original GLCM J En", font.title=18,
                                       legend="top", legend.labs=c("high","low"),
                                       font.x=14,font.y=14, pval.coord = c(0,0.1), font.tickslab=14,font.legend=12,pval.size=5, risk.table.height = FALSE, ggtheme = centr); splots_OS[[3]]

# # log.sigma.1.0.mm.3D_gldm_DependenceVariance_cate
# sopravv <- survfit(Surv(time,status)~as.factor(log.sigma.1.0.mm.3D_gldm_DependenceVariance_cate),data=discovery); sopravv
# splots_radiomics_OS[[5]] <- ggsurvplot(sopravv, data=discovery, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability OS",
#                                        risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
#                                        palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE,
#                                        title="loG.sigma.1.0.mm.3D gldm D V",
#                                        legend="none", font.x=12,font.y=12, pval.coord = c(0,0.1), font.title = 12, font.tickslab=12,pval.size=4,risk.table.height = FALSE, ggtheme = centr); splots_radiomics_OS[[5]]

# log.sigma.1.0.mm.3D_firstorder_Skewness_cate
sopravv <- survfit(Surv(time,status)~as.factor(log.sigma.1.0.mm.3D_firstorder_Skewness_cate),data=discovery); sopravv
splots_OS[[4]] <- ggsurvplot(sopravv, data=discovery, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability OS",
                                       risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
                                       palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE,
                                       title="LoG.sig.1.mm Fo Skew", font.title=18,
                                       legend="top", legend.labs=c("high","low"),
                                       font.x=14,font.y=14, pval.coord = c(0,0.1), font.tickslab=14,font.legend=12,pval.size=5, risk.table.height = FALSE, ggtheme = centr); splots_OS[[4]]

# log.sigma.2.0.mm.3D_firstorder_Median_cate
sopravv <- survfit(Surv(time,status)~as.factor(log.sigma.2.0.mm.3D_firstorder_Median_cate),data=discovery); sopravv
splots_OS[[5]] <- ggsurvplot(sopravv, data=discovery, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability OS",
                                       risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
                                       palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE,
                                       title="LoG.sig.2.mm Fo Med", font.title=18,
                                       legend="top", legend.labs=c("high","low"),
                                       font.x=14,font.y=14, pval.coord = c(0,0.1), font.tickslab=14,font.legend=12,pval.size=5, risk.table.height = FALSE, ggtheme = centr); splots_OS[[5]]

# # wavelet.LLH_glcm_ClusterTendency_cate
# sopravv <- survfit(Surv(time,status)~as.factor(wavelet.LLH_glcm_ClusterTendency_cate),data=discovery); sopravv
# splots_radiomics_OS[[8]] <- ggsurvplot(sopravv, data=discovery, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability OS",
#                                       risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
#                                       palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE,
#                                       title="wavelet.LLH glcm Cl Te",
#                                       legend="none", font.x=12,font.y=12, pval.coord = c(0,0.1), font.title = 12, font.tickslab=12,pval.size=4,risk.table.height = FALSE, ggtheme = centr); splots_radiomics_OS[[8]]
# 
# wavelet.HLH_glcm_Imc1_cate
sopravv <- survfit(Surv(time,status)~as.factor(wavelet.HLH_glcm_Imc1_cate),data=discovery); sopravv
splots_OS[[6]] <- ggsurvplot(sopravv, data=discovery, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability OS",
                                       risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
                                       palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE,
                                       title="Wave.HLH GLCM IMC1", font.title=18,
                                       legend="top", legend.labs=c("high","low"),
                                       font.x=14,font.y=14, pval.coord = c(0,0.1), font.tickslab=14,font.legend=12,pval.size=5, risk.table.height = FALSE, ggtheme = centr); splots_OS[[6]]

# # wavelet.HLL_glcm_InverseVariance_cate
# sopravv <- survfit(Surv(time,status)~as.factor(wavelet.HLL_glcm_InverseVariance_cate),data=discovery); sopravv
# splots_radiomics_OS[[10]] <- ggsurvplot(sopravv, data=discovery, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability OS",
#                                        risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
#                                        palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE,
#                                        title="wavelet.HLL glcm Inv Var",
#                                        legend="none", font.x=12,font.y=12, pval.coord = c(0,0.1), font.title = 12, font.tickslab=12,pval.size=4,risk.table.height = FALSE, ggtheme = centr); splots_radiomics_OS[[10]]

# # wavelet.HHL_firstorder_Mean_cate
# sopravv <- survfit(Surv(time,status)~as.factor(wavelet.HHL_firstorder_Mean_cate),data=discovery); sopravv
# splots_radiomics_OS[[11]] <- ggsurvplot(sopravv, data=discovery, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability OS",
#                                        risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
#                                        palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE,
#                                        title="wavelet.HHL firstorder Mean",
#                                        legend="none", legend.labs = c("High","Low"), font.legend=16,
#                                        font.x=12,font.y=12, pval.coord = c(0,0.1), font.title = 12, font.tickslab=12,pval.size=4,risk.table.height = FALSE, ggtheme = centr); splots_radiomics_OS[[11]]

# # wavelet.HHL_glcm_Imc1_cate
# sopravv <- survfit(Surv(time,status)~as.factor(wavelet.HHL_glcm_Imc1_cate),data=discovery); sopravv
# splots_radiomics_OS[[12]] <- ggsurvplot(sopravv, data=discovery, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability OS",
#                                        risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
#                                        palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE,
#                                        title="wavelet.HHL glcm Imc1",
#                                        legend="none", font.x=12,font.y=12, pval.coord = c(0,0.1), font.title = 12, font.tickslab=12,pval.size=4,risk.table.height = FALSE, ggtheme = centr); splots_radiomics_OS[[12]]

# # wavelet.HHL_glszm_LargeAreaEmphasis_cate
# sopravv <- survfit(Surv(time,status)~as.factor(wavelet.HHL_glszm_LargeAreaEmphasis_cate),data=discovery); sopravv
# splots_radiomics_OS[[13]] <- ggsurvplot(sopravv, data=discovery, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability OS",
#                                        risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
#                                        palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE,
#                                        title="wavelet.HHL glszm L AE",
#                                        legend="none", font.x=12,font.y=12, pval.coord = c(0,0.1), font.title = 12, font.tickslab=12,pval.size=4,risk.table.height = FALSE, ggtheme = centr); splots_radiomics_OS[[13]]
# Gender
sopravv <- survfit(Surv(time,status)~as.factor(df_clinical.gender),data=discovery); sopravv
splots_OS[[7]] <- ggsurvplot(sopravv, data=discovery, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability OS",
                             risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
                             palette= c("#B3C99C","#06948E"), legend.labs = c("Females","Males"),
                             legend="top", risk.table.title=FALSE,
                             title="Gender", font.title=18,
                             font.x=14,font.y=14, pval.coord = c(0,0.1),
                             font.tickslab=14,font.legend=12,pval.size=5,risk.table.height = FALSE, ggtheme = centr); splots_OS[[7]]

# Grade
sopravv <- survfit(Surv(time,status)~as.factor(df_clinical.tumor_grade_reclass),data=discovery); sopravv
splots_OS[[8]] <- ggsurvplot(sopravv, data=discovery, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability OS",
                             risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
                             palette= c("#06948E","#B3C99C"), legend.labs = c("G1/G2","G3"),
                             legend="top", risk.table.title=FALSE,
                             title="Grade", font.title=18,
                             font.x=14,font.y=14, pval.coord = c(0,0.1),
                             font.tickslab=14,font.legend=12,pval.size=5,risk.table.height = FALSE, ggtheme = centr); splots_OS[[8]]

# Residual disease
sopravv <- survfit(Surv(time,status)~as.factor(df_clinical.residual_disease_reclass),data=discovery); sopravv
splots_OS[[9]] <- ggsurvplot(sopravv, data=discovery, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability OS",
                             risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
                             palette= c("#06948E","#B3C99C","darkgray"), legend.labs = c("R1/R2","R0","RX"),
                             legend="top", risk.table.title=FALSE,
                             title="Residual Disease", font.title=18,
                             font.x=14,font.y=14, pval.coord = c(0,0.1),
                             font.tickslab=14,font.legend=12,pval.size=5,risk.table.height = FALSE, ggtheme = centr); splots_OS[[9]]

# KRAS
sopravv <- survfit(Surv(time,status)~as.factor(KRAS),data=discovery); sopravv
splots_OS[[10]] <- ggsurvplot(sopravv, data=discovery, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability OS",
                             risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
                             palette= c("blue","#9BABB8"), legend.labs = c("Not mut","Mut"),
                             legend="top", risk.table.title=FALSE,
                             title="KRAS", font.title=18,
                             font.x=14,font.y=14, pval.coord = c(0,0.1),
                             font.tickslab=14,font.legend=12,pval.size=5,risk.table.height = FALSE, ggtheme = centr); splots_OS[[10]]

# TTN
sopravv <- survfit(Surv(time,status)~as.factor(TTN),data=discovery); sopravv
splots_OS[[11]] <- ggsurvplot(sopravv, data=discovery, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability OS",
                             risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
                             palette= c("blue","#9BABB8"), legend.labs = c("Not mut","Mut"),
                             legend="top", risk.table.title=FALSE,
                             title="TTN", font.title=18,
                             font.x=14,font.y=14, pval.coord = c(0,0.1),
                             font.tickslab=14,font.legend=12,pval.size=5,risk.table.height = FALSE, ggtheme = centr); splots_OS[[10]]

res <- arrange_ggsurvplots(splots_OS, print = TRUE, nrow = 3, ncol = 4); res;
ggsave(gsub(" ","",paste(getwd(),"/results/UV_analysis/Discovery/","UV_KM_panel_OS_R.png")), res, width = 4000, height = 3000, units = "px", dpi = 300)

# ##################################################################################################################################################
# # RADIOMICS - MULTIVARIATE SURVIVAL ANALYSIS
# ##################################################################################################################################################

#Modello di cox multiva riato
cox_radio <- coxph(formula=Surv(time, status) ~
               original_firstorder_Kurtosis_cate +
               original_firstorder_90Percentile_cate +
               # original_glcm_Correlation_cate +
               original_shape_Elongation_cate +
               original_glcm_JointEnergy_cate +
               # original_shape_LeastAxisLength_cate +
               log.sigma.1.0.mm.3D_gldm_DependenceVariance_cate +
               log.sigma.1.0.mm.3D_firstorder_Skewness_cate +
               # log.sigma.2.0.mm.3D_glszm_SmallAreaLowGrayLevelEmphasis_cate +
               log.sigma.2.0.mm.3D_firstorder_Median_cate +
               wavelet.LLH_glcm_ClusterTendency_cate +
               wavelet.HLL_glcm_InverseVariance_cate +
               wavelet.HLH_glcm_Imc1_cate +
               wavelet.HHL_firstorder_Mean_cate +
               wavelet.HHL_glcm_Imc1_cate +
               wavelet.HHL_glszm_LargeAreaEmphasis_cate,
             data=discovery, x =TRUE); summary(cox_radio)
# ggforest(cox, data=discovery, main = "Probability OS - Train Set")

# # ##################################################################################################################################################
# # # CLINICAL - MULTIVARIATE SURVIVAL ANALYSIS
# # ##################################################################################################################################################

#Modello di cox multivariato
cox_clin <- coxph(formula=Surv(time, status) ~
             df_clinical.gender +
             df_clinical.tumor_grade_reclass +
             df_clinical.residual_disease_reclass,
             data=discovery, x = TRUE); summary(cox_clin)

# ##################################################################################################################################################
# # MUTATIONS - MULTIVARIATE SURVIVAL ANALYSIS
# ##################################################################################################################################################

# Modello di cox multivariato
cox_mut <- coxph(formula=Surv(time, status) ~
               KRAS +
               TTN,
             data=discovery, x = TRUE); summary(cox_mut)

# # ##################################################################################################################################################
# # #INTEGRATED FOREST PLOTS
# # ##################################################################################################################################################
# 
# all.models <- list()
# all.models[[1]] <- cox_radio; summary(all.models[[1]])
# all.models[[2]] <- cox_clin; summary(all.models[[2]])
# all.models[[3]] <-cox_mut; summary(all.models[[3]])
# 
# # forest_model(model_list = list("Radiomics" = all.models[[1]], "Clinical" = all.models[[2]],"Mutational" = all.models[[3]]))
# 
# p <- plot_models(all.models,
#                  grid = F,
#                  axis.title = "HR (Probability OS)",
#                  legend.title = NULL,
#                  m.labels = c("Radiomics","Clinical","Mutational"),
#                  axis.labels = c("TTN, mutated vs. not",
#                                  "KRAS, mutated vs. not",
#                                  "Residual Disease, R0 vs. R1/R2",
#                                  "Grade, G3 vs. G1/G2",
#                                  "Gender, M vs. F",
#                                  "Wavelet.HHL GLSZM LArea Emphasis, low vs. high",
#                                  "Wavelet.HHL GLCM Imc1, low vs. high",
#                                  "Wavelet.HHL Fo Mean, low vs. high",
#                                  "Wavelet.HLH GLCM Imc1, low vs. high",
#                                  "Wavelet.HLL GLCM Inv Variance, low vs. high",
#                                  "Wavelet.LLH GLCM Cl Tendency, low vs. high",
#                                  "LoG.sigma.2.0.mm.3D Fo Median, low vs. high",
#                                  "LoG.sigma.1.0.mm.3D Fo Skewness, low vs. high",
#                                  "LoG.sigma.1.0.mm.3D GLDM DV, low vs. high",
#                                  "Original GLCM Jo Energy, low vs. high",
#                                  "Original Shape Elongation, low vs. high",
#                                  "Original Fo 90P, low vs. high", 
#                                  "Original Fo Kurtosis, low vs. high"),
#                  p.threshold = c(0.10, 0.05, 0.01),
#                  wrap.labels = 200,
#                  vline.color = "black",
#                  line.size = 3,
#                  dot.size = 4,
#                  axis.lim = c(0.01,20),
#                  show.values = TRUE,
#                  value.size = 4,
#                  spacing = 0.7,
#                  colors = c("#2D5F8A","#06948E","#e66101"))
# p + set_theme(base = theme_transparent(base_family = 'Helvetica'),
#               legend.pos = "left",
#               legend.size = 1.5,
#               legend.item.size = 1.2,
#               legend.item.backcol = "white",
#               axis.title.size = 1.5,
#               axis.textsize.y = 1.2,
#               axis.textsize.x = 1.2,
#               axis.title.color = "black",
#               axis.textcolor = "black",
# )
# 
# ggsave(gsub(" ","",paste(getwd(),"/results/UV_analysis/Discovery/","Forest_plot_OS.png")), p, width = 3100, height = 4000, units = "px", dpi = 300)
# 
# 
