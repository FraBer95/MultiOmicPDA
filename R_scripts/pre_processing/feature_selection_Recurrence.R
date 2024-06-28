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
# # # ##################################################################################################################################################
# # # # PRE-PROCESSING
# # # ##################################################################################################################################################

#set path
setwd("path_to_generic_folder")

# #Loading of multi-modal dataset
# discovery <- read.csv(gsub(" ", "", paste(getwd(),"/multimodal_dataset.csv")))
# 
# #Loading of the clinical dataset for adding Recurrence time and other clinical features
# load("integrated_dataset_070224.Rda")
# clinical <- df_train[,c("Barcode","df_clinical.months_to_prog_rec","df_clinical.prog_rec_status","df_follow_up.bmi_categorical","path_detail.lymph_nodes_positive","path_detail.tumor_largest_dimension_diameter_categorical","df_follow_up.ecog_performance_status","df_follow_up.karnofsky_performance_status","df_follow_up.disease_response")];
# remove(df_train)
# 
# #Include novel features in discovery-set
# discovery$time <- 0
# discovery$status <- 0
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
#       discovery$time[index_i] <- clinical$df_clinical.months_to_prog_rec[index_j]
#       discovery$status[index_i] <- clinical$df_clinical.prog_rec_status[index_j]
#       discovery$BMI[index_i] <- clinical$df_follow_up.bmi_categorical[index_j]
#       discovery$LN_pos[index_i] <- clinical$path_detail.lymph_nodes_positive[index_j]
#       discovery$tumor_largest_dia[index_i] <- clinical$path_detail.tumor_largest_dimension_diameter_categorical[index_j]
#       discovery$ecogps[index_i] <- clinical$df_follow_up.ecog_performance_status[index_j]
#       discovery$karn[index_i] <- clinical$df_follow_up.karnofsky_performance_status[index_j]
#       discovery$clin_resp[index_i] <- clinical$df_follow_up.disease_response[index_j]
#     }
#   }
# }

#Remove clinical df
remove(clinical)
# 
# # Subsetting patients without Recurrence
# discovery_sub <- subset(discovery,status != "NA")
# remove(discovery)
# discovery_sub$status <- as.numeric(discovery_sub$status)
# 
# #Preprocessing, transforming all categorical features in categorical numeric
# table(discovery_sub$BMI)
# discovery_sub$BMI_reclass = 0
# discovery_sub$BMI_reclass[discovery_sub$BMI == "OW"] = 1
# discovery_sub$BMI_reclass[discovery_sub$BMI == "NA"] = NA
# table(discovery_sub$BMI_reclass)
# 
# table(discovery_sub$LN_pos)
# discovery_sub$LN_pos_reclass = 0
# discovery_sub$LN_pos_reclass[discovery_sub$LN_pos >= 1] = 1
# discovery_sub$LN_pos_reclass[is.na(discovery_sub$LN_pos)] = NA
# table(discovery_sub$LN_pos_reclass)
# 
# table(discovery_sub$tumor_largest_dia)
# discovery_sub$tumor_largest_dia_reclass = 0
# discovery_sub$tumor_largest_dia_reclass[discovery_sub$tumor_largest_dia == ">=3.5cm"] = 1
# discovery_sub$tumor_largest_dia_reclass[discovery_sub$tumor_largest_dia == "NA"] = NA
# table(discovery_sub$tumor_largest_dia_reclass)
# 
# table(discovery_sub$ecogps)
# discovery_sub$ecogps_reclass = 0
# discovery_sub$ecogps_reclass[discovery_sub$ecogps == "3" | discovery_sub$ecogps == "4" | discovery_sub$ecogps == "5"] = 1
# discovery_sub$ecogps_reclass[discovery_sub$ecogps == "Not Reported" | discovery_sub$ecogps == "Unknown"] = NA
# table(discovery_sub$ecogps_reclass)
# 
# table(discovery_sub$karn)
# discovery_sub$karn_reclass = 0
# discovery_sub$karn_reclass[discovery_sub$karn == "30" | discovery_sub$karn == "40" | discovery_sub$karn == "50" | discovery_sub$karn == "60" | discovery_sub$karn == "70" | discovery_sub$karn == "80" | discovery_sub$karn == "90" | discovery_sub$karn == "100"] = 1
# discovery_sub$karn_reclass[discovery_sub$karn == "Not Reported" | discovery_sub$karn == "Unknown"] = NA
# table(discovery_sub$karn_reclass)
# 
# table(discovery_sub$clin_resp)
# discovery_sub$clin_resp_reclass = 0
# discovery_sub$clin_resp_reclass[discovery_sub$clin_resp == "CR-Complete Response" | discovery_sub$clin_resp == "PR-Partial Response"] = 1
# discovery_sub$clin_resp_reclass[discovery_sub$clin_resp == "Not Reported" | discovery_sub$clin_resp == "Unknown"] = NA
# table(discovery_sub$clin_resp_reclass)
# 
# # # ##################################################################################################################################################
# # # # SURVIVAL ANALYSIS  - Recurrence
# # # ##################################################################################################################################################
# 
# #Standardizing outcome
# time_str = "time"
# event_str = "status"
# time = discovery_sub$time
# event = discovery_sub$status
# 
# ##################################################################################################################################################
# # RADIOMICS FEATURES - SURVIVAL ANALYSIS
# ##################################################################################################################################################
# 
# # Create dataframe of results
# df_res_radio = data.frame(matrix(nrow = 14, ncol = 7))
# colnames(df_res_radio) <- c("variables","pval","method","pval.txt","cut_off","C-index","std")
# df_res_radio[,1] <- colnames(discovery_sub[,2:15]);
# 
# for (i in colnames(discovery_sub[,2:17])) {
# 
#   #Cut-off and survival analysis
#   res.cut <- surv_cutpoint(discovery_sub, time = time_str, event = event_str, variables = i); cutoff <- res.cut$cutpoint[1,1]
#   summary(res.cut)
# 
#   # pdf(file = gsub(" ", "", paste("/Users/gianmariazaccaria/Documents/R_studio/CPTAC-PDA/Integrated_analysis/test_4_entire_dataset_progression_020224/results/UV_analysis/radiomics/",i,"_MSRS.pdf")))
#   # plot(res.cut,i, palette = "npg")
#   # dev.off()
#   #
#   res.cat <- surv_categorize(res.cut)
# 
#   #Aggiungo valore categorizzato della variabile al dataframe
#   string <- gsub(" ", "", paste(i,"_cate"))    #Creo nome colonna
#   discovery_sub[string] <- res.cat[,3]
# 
#   fit <- survfit(Surv(time, event) ~res.cat[,3], data = res.cat); fit
#   pp = ggsurvplot(fit, data=res.cat, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                   palette= c("red","blue"), legend.labs = c("high","low"),
#                   legend="top", risk.table.title="N. at risk", title=i, font.x=12,font.y=12,
#                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2);
# 
#   # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/radiomics/progression/",i,"_KM.pdf")))
#   # print(pp, newpage = FALSE)
#   # dev.off()
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
# ##################################################################################################################################################
# # CLINICAL FEATURES - SURVIVAL ANALYSIS
# ##################################################################################################################################################
# 
# # # # #Definition of lists for each panel construction
# # # # splots_clinical <- list(); splots_exposure <- list();
# # #
# # # df_clinical.gender
# table(discovery_sub$df_clinical.gender)
# name <- ("Gender")
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(df_clinical.gender),data=discovery_sub); sopravv
# pp = ggsurvplot(sopravv, data=discovery_sub, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                 palette= c("red","blue"), legend.labs = c("female","male"),
#                                 legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time, event) ~ as.factor(df_clinical.gender), data=discovery_sub); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # # discovery_sub$df_clinical.age_categorical
# table(discovery_sub$df_clinical.age_categorical)
# name <- ("Age at diag")
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(discovery_sub$df_clinical.age_categorical), data=discovery_sub); sopravv
# pp=ggsurvplot(sopravv, data=discovery_sub, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                 palette= c("red","blue"), legend.labs = c("Low Age","High Age"),
#                                 legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time, event) ~ as.factor(discovery_sub$df_clinical.age_categorical), data=discovery_sub); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # # # BMI
# table(discovery_sub$BMI_reclass)
# df_2 <- subset(discovery_sub, !is.na(BMI_reclass))
# table(df_2$BMI_reclass)
# time2 = df_2$time
# event2 = df_2$status
# 
# name <- ("BMI")
# sopravv <- survfit(Surv(time2, as.numeric(event2))~as.factor(df_2$BMI_reclass), data=df_2); sopravv
# pp=ggsurvplot(sopravv, data=df_2, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                 palette= c("red","blue"), legend.labs = c("Under Weight","Over Weight"),
#                                 legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time2, event2) ~ as.factor(df_2$BMI_reclass), data=df_2); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # discovery_sub$df_clinical.ajcc_pathologic_stage
# table(discovery_sub$df_clinical.ajcc_pathologic_stage_reclass_2)
# name <- ("Stage")
# 
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(df_clinical.ajcc_pathologic_stage_reclass_2),data=discovery_sub); sopravv
# pp=ggsurvplot(sopravv, data=discovery_sub, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                 palette= c("red","blue"), legend.labs = c("Stage IB/IIA/IIB","Stage III/IIIB/IV"),
#                                 legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time, event) ~ as.factor(df_clinical.ajcc_pathologic_stage_reclass_2), data=discovery_sub); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # df_clinical.tumor_grade
# table(discovery_sub$df_clinical.tumor_grade)
# name <- ("Tumor Grade")
# 
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(discovery_sub$df_clinical.tumor_grade_reclass),data=discovery_sub); sopravv
# pp=ggsurvplot(sopravv, data=discovery_sub, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#               risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#               palette= c("red","blue"), legend.labs = c("G1/G2","G3/G4"),
#               legend="top", risk.table.title="N. at risk", title=name, font.x=16,font.y=16,
#               font.tickslab=14, font.legend=16, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time, event) ~ discovery_sub$df_clinical.tumor_grade_reclass, data=discovery_sub); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # discovery_sub$LN_pos_reclass
# table(discovery_sub$LN_pos_reclass);
# name <- ("Lymphnode Involvement")
# 
# df_2 <- subset(discovery_sub, !is.na(discovery_sub$LN_pos_reclass))
# 
# time2 = df_2$time
# event2 = df_2$status
# 
# sopravv <- survfit(Surv(time2, as.numeric(event2))~as.factor(df_2$LN_pos_reclass),data=df_2); sopravv
# pp=ggsurvplot(sopravv, data=df_2, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                 palette= c("red","blue"), legend.labs = c("0", "at least 1"),
#                                 legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time2, event2) ~ as.factor(df_2$LN_pos_reclass), data=df_2); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # discovery_sub$path_detail.tumor_largest_dimension_diameter_categorical
# table(discovery_sub$tumor_largest_dia_reclass)
# name <- ("Tumor diameter")
# 
# df_2 <- subset(discovery_sub, !is.na(discovery_sub$tumor_largest_dia_reclass))
# time2 = df_2$time
# event2 = df_2$status
# 
# sopravv <- survfit(Surv(time2, as.numeric(event2))~as.factor(df_2$tumor_largest_dia_reclass),data=df_2); sopravv
# pp=ggsurvplot(sopravv, data=df_2, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                 palette= c("red","blue"), legend.labs = c("lower 3.5", "higher 3.5"),
#                                 legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time2, event2) ~ as.factor(df_2$tumor_largest_dia_reclass), data=df_2); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # df_follow_up.ecog_performance_status
# table(discovery_sub$ecogps_reclass); #freq(discovery_sub$df_follow_up.ecog_performance_status)
# name <- ("ECOGps")
# 
# df_2 <- subset(discovery_sub, !is.na(discovery_sub$ecogps_reclass))
# time2 = df_2$time
# event2 = df_2$status
# 
# sopravv <- survfit(Surv(time2, as.numeric(event2))~as.factor(df_2$ecogps_reclass),data=df_2); sopravv
# pp=ggsurvplot(sopravv, data=df_2, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                 palette= c("red","blue"), legend.labs = c("0/1/2/3", "4/5"),
#                                 legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time2, event2) ~ as.factor(df_2$ecogps_reclass), data=df_2); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # df_follow_up.karnofsky_performance_status
# table(discovery_sub$karn_reclass);
# name <- ("Karnofsky")
# 
# df_2 <- subset(discovery_sub, !is.na(discovery_sub$karn_reclass))
# time2 = df_2$time
# event2 = df_2$status
# 
# sopravv <- survfit(Surv(time2, as.numeric(event2))~as.factor(df_2$karn_reclass),data=df_2); sopravv
# pp=ggsurvplot(sopravv, data=df_2, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                 risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                 palette= c("red","blue"), legend.labs = c( "0", "higher then 0%"),
#                                 legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                 font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time2, event2) ~ as.factor(df_2$karn_reclass), data=df_2); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # discovery_sub$df_follow_up.disease_response
# table(discovery_sub$clin_resp_reclass)
# name <- ("Clinical Response")
# 
# df_2 <- subset(discovery_sub, !is.na(discovery_sub$clin_resp_reclass))
# time2 = df_2$time
# event2 = df_2$status
# 
# sopravv <- survfit(Surv(time2, as.numeric(event2))~as.factor(df_2$clin_resp_reclass),data=df_2); sopravv
# pp <- ggsurvplot(sopravv, data=df_2, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                     risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                     palette= c("red","blue"), legend.labs = c("Not Responders", "Responders"),
#                                     legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                     font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time2, event2) ~ as.factor(df_2$clin_resp_reclass), data=df_2); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # df_clinical.residual_disease
# table(discovery_sub$df_clinical.residual_disease_reclass)
# 
# name <- "Residual Disease Reclass"
# sopravv <- survfit(Surv(time, status)~as.factor(df_clinical.residual_disease_reclass),data=discovery_sub); sopravv
# pp=ggsurvplot(sopravv, data=discovery_sub, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                  risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                  palette= c("red","blue","darkgreen"), legend.labs = c("R0", "R1/R2","RX"),
#                                  legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                  font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# 
# cox <- coxph(formula=Surv(time, status) ~ as.factor(discovery_sub$df_clinical.residual_disease), data=discovery_sub); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/clinical/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # ##################################################################################################################################################
# # # MUTATIONAL FEATURES - SURVIVAL ANALYSIS
# # ##################################################################################################################################################
# 
# # # KRAS
# table(discovery_sub$KRAS)
# name <- "KRAS"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(KRAS),data=discovery_sub); sopravv
# pp=ggsurvplot(sopravv, data=discovery_sub, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                   palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                   legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(KRAS), data=discovery_sub); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # # TP53
# table(discovery_sub$TP53)
# name <- "TP53"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(TP53),data=discovery_sub); sopravv
# pp=ggsurvplot(sopravv, data=discovery_sub, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                   palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                   legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(TP53), data=discovery_sub); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # CDKN2A
# table(discovery_sub$CDKN2A)
# name <- "CDKN2A"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(CDKN2A),data=discovery_sub); sopravv
# pp=ggsurvplot(sopravv, data=discovery_sub, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                   palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                   legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(CDKN2A), data=discovery_sub); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # SMAD4
# table(discovery_sub$SMAD4)
# name <- "SMAD4"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(SMAD4),data=discovery_sub); sopravv
# pp=ggsurvplot(sopravv, data=discovery_sub, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                   palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                   legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(SMAD4), data=discovery_sub); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # TTN
# table(discovery_sub$TTN)
# name <- "TTN"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(TTN),data=discovery_sub); sopravv
# pp=ggsurvplot(sopravv, data=discovery_sub, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                   palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                   legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(TTN), data=discovery_sub); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# #MUC16
# table(discovery_sub$MUC16)
# name <- "MUC16"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(MUC16),data=discovery_sub); sopravv
# pp=ggsurvplot(sopravv, data=discovery_sub, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                   palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                   legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(MUC16), data=discovery_sub); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# #KMT2D
# table(discovery_sub$KMT2D)
# name <- "KMT2D"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(KMT2D),data=discovery_sub); sopravv
# pp=ggsurvplot(sopravv, data=discovery_sub, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                   palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                   legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(KMT2D), data=discovery_sub); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # CSMD1
# table(discovery_sub$CSMD1)
# name <- "CSMD1"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(CSMD1),data=discovery_sub); sopravv
# pp=ggsurvplot(sopravv, data=discovery_sub, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                   palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                   legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(CSMD1), data=discovery_sub); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # RYR2
# table(discovery_sub$RYR2)
# name <- "RYR2"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(RYR2),data=discovery_sub); sopravv
# pp=ggsurvplot(sopravv, data=discovery_sub, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                   risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                   palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                   legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                   font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(RYR2), data=discovery_sub); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# 
# # RYR1
# table(discovery_sub$RYR1)
# name <- "RYR1"
# sopravv <- survfit(Surv(time, as.numeric(event))~as.factor(RYR1),data=discovery_sub); sopravv
# pp=ggsurvplot(sopravv, data=discovery_sub, risk.table=TRUE, conf.int=TRUE, xlim=c(0,70), pval=TRUE,ylab="Probability Recurrence",
#                                    risk.table.y.col=TRUE, risk.table.y.text=TRUE, xlab="Time (months)", break.time.by=10,
#                                    palette= c("red","blue"), legend.labs = c("No Mut","Mut"),
#                                    legend="top", risk.table.title="N. at risk", title=name, font.x=12,font.y=12,
#                                    font.tickslab=12, font.legend=14, pval.size=6.5, risk.table.height = 0.2); pp
# cox <- coxph(formula=Surv(time, event) ~ as.factor(RYR1), data=discovery_sub); summary(cox)
# 
# # pdf(file = gsub(" ", "", paste(getwd(),"/results/UV_analysis/Discovery/mutations/progression/",name,"_KM.pdf")))
# # print(pp, newpage = FALSE)
# # dev.off()
# # 
# # ##################################################################################################################################################
# # # SAVING DATAFRAMES 
# # ##################################################################################################################################################
# 
# save(discovery_sub,file="discovery_set_REC.Rda")
# save(df_res_radio,file="df_res_radio_discovery_set_REC.Rda")

# ##################################################################################################################################################
# # LOADINING DATAFRAMES
# ##################################################################################################################################################

load(gsub(" ", "", paste(getwd(),"/results/discovery_set_REC.Rda")))
load(gsub(" ", "", paste(getwd(),"/results/df_res_radio_discovery_set_REC.Rda")))

# # # ##################################################################################################################################################
# # # # PATIENTS' CHARACTERISTICS
# # # ##################################################################################################################################################
# 
# #Including features
# discovery_patchar <- discovery_sub[,c("df_clinical.age_categorical", "df_clinical.gender", "df_clinical.ajcc_pathologic_stage_reclass_2","df_clinical.tumor_grade_reclass",
#                                   "KRAS","TP53", "CDKN2A", "SMAD4", "TTN", "MUC16", "KMT2D", "CSMD1", "RYR2","RYR1", "status", "BMI_reclass", "LN_pos_reclass",
#                                   "tumor_largest_dia_reclass", "ecogps_reclass", "karn_reclass", "clin_resp_reclass",
#                                    "original_shape_Elongation_cate","original_firstorder_90Percentile_cate","log.sigma.2.0.mm.3D_glszm_SmallAreaLowGrayLevelEmphasis_cate", 
#                                   "wavelet.LHL_glcm_InverseVariance_cate", "wavelet.LLH_firstorder_Variance_cate", "wavelet.HHL_firstorder_Mean_cate",
#                                   "wavelet.HHL_glszm_LargeAreaEmphasis_cate", "log.sigma.4.0.mm.3D_glszm_SmallAreaLowGrayLevelEmphasis_cate",
#                                   "log.sigma.1.0.mm.3D_glcm_ClusterShade_cate", "original_firstorder_Kurtosis_cate", "log.sigma.2.0.mm.3D_firstorder_Median_cate",
#                                   "original_shape_LeastAxisLength_cate", "original_glcm_JointEnergy_cate", "original_glcm_Correlation_cate" ,"df_clinical.residual_disease_reclass")]
# 
# #Rename features
# colnames(discovery_patchar)[1] <- c("Age at diagnosis")
# colnames(discovery_patchar)[2] <- c("Gender")
# colnames(discovery_patchar)[3] <- c("Stage")
# colnames(discovery_patchar)[4] <- c("Grade")
# colnames(discovery_patchar)[15] <- c("Recurrence status")
# colnames(discovery_patchar)[16] <- c("BMI")
# colnames(discovery_patchar)[17] <- c("Lymphnode involvement")
# colnames(discovery_patchar)[18] <- c("Tumor largest diameter")
# colnames(discovery_patchar)[19] <- c("ECOGps")
# colnames(discovery_patchar)[20] <- c("Karnofsky")
# colnames(discovery_patchar)[21] <- c("Clinical response")
# colnames(discovery_patchar)[22] <- c("Original shape Elongation")
# colnames(discovery_patchar)[23] <- c("Original firstorder 90Percentile")
# colnames(discovery_patchar)[31] <- c("Original firstorder Kurtosis")
# colnames(discovery_patchar)[33] <- c("Original shape Least Axis Length")
# colnames(discovery_patchar)[34] <- c("Original glcm Joint Energy")
# colnames(discovery_patchar)[35] <- c("Original glcm Correlation")
# colnames(discovery_patchar)[24] <- c("loG.sigma.2.0.mm.3D glszm Small Area Low Gray Level Emphasis")
# colnames(discovery_patchar)[29] <- c("loG.sigma.4.0.mm.3D glszm Small Area Low Gray Level Emphasis")
# colnames(discovery_patchar)[30] <- c("loG.sigma.1.0.mm.3D glcm Cluster Shade")
# colnames(discovery_patchar)[32] <- c("loG.sigma.2.0.mm.3D firstorder Median")
# colnames(discovery_patchar)[25] <- c("Wavelet.LHL glcm Inverse Variance")
# colnames(discovery_patchar)[26] <- c("Wavelet.LLH firstorder Variance")
# colnames(discovery_patchar)[27] <- c("Wavelet.HHL firstorder Mean")
# colnames(discovery_patchar)[28] <- c("Wavelet.HHL glszm Large Area Emphasis")
# colnames(discovery_patchar)[36] <- c("Residual Disease")
# # summary(discovery_patchar)
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
# table(discovery_patchar$`Recurrence status`)
# discovery_patchar$`Recurrence status`[discovery_patchar$`Recurrence status` == 0] <- "Not Recurrence"
# discovery_patchar$`Recurrence status`[discovery_patchar$`Recurrence status` == 1] <- "Recurrence"
# discovery_patchar$`Recurrence status` <- as.factor(discovery_patchar$`Recurrence status`)
# summary(discovery_patchar$`Recurrence status`) 
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
# table(discovery_patchar$`Original shape Elongation`)
# discovery_patchar$`Original shape Elongation` <- factor(discovery_patchar$`Original shape Elongation`, levels=c('Low', 'High'))
# table(discovery_patchar$`Original shape Elongation`)
# 
# table(discovery_patchar$`Original firstorder 90Percentile`)
# discovery_patchar$`Original firstorder 90Percentile` <- factor(discovery_patchar$`Original firstorder 90Percentile`, levels=c('Low', 'High'))
# table(discovery_patchar$`Original firstorder 90Percentile`)
# 
# table(discovery_patchar$`loG.sigma.2.0.mm.3D glszm Small Area Low Gray Level Emphasis`)
# discovery_patchar$`loG.sigma.2.0.mm.3D glszm Small Area Low Gray Level Emphasis` <- factor(discovery_patchar$`loG.sigma.2.0.mm.3D glszm Small Area Low Gray Level Emphasis`, levels=c('Low', 'High'))
# table(discovery_patchar$`loG.sigma.2.0.mm.3D glszm Small Area Low Gray Level Emphasis`)
# 
# table(discovery_patchar$`Wavelet.LHL glcm Inverse Variance`)
# discovery_patchar$`Wavelet.LHL glcm Inverse Variance` <- factor(discovery_patchar$`Wavelet.LHL glcm Inverse Variance`, levels=c('Low', 'High'))
# table(discovery_patchar$`Wavelet.LHL glcm Inverse Variance`)
# 
# table(discovery_patchar$`Wavelet.LLH firstorder Variance`)
# discovery_patchar$`Wavelet.LLH firstorder Variance` <- factor(discovery_patchar$`Wavelet.LLH firstorder Variance`, levels=c('Low', 'High'))
# table(discovery_patchar$`Wavelet.LLH firstorder Variance`)
# 
# table(discovery_patchar$`Wavelet.HHL firstorder Mean`)
# discovery_patchar$`Wavelet.HHL firstorder Mean` <- factor(discovery_patchar$`Wavelet.HHL firstorder Mean`, levels=c('Low', 'High'))
# table(discovery_patchar$`Wavelet.HHL firstorder Mean`)
# 
# table(discovery_patchar$`Wavelet.HHL glszm Large Area Emphasis`)
# discovery_patchar$`Wavelet.HHL glszm Large Area Emphasis` <- factor(discovery_patchar$`Wavelet.HHL glszm Large Area Emphasis`, levels=c('Low', 'High'))
# table(discovery_patchar$`Wavelet.HHL glszm Large Area Emphasis`)
# 
# table(discovery_patchar$`loG.sigma.4.0.mm.3D glszm Small Area Low Gray Level Emphasis`)
# discovery_patchar$`loG.sigma.4.0.mm.3D glszm Small Area Low Gray Level Emphasis` <- factor(discovery_patchar$`loG.sigma.4.0.mm.3D glszm Small Area Low Gray Level Emphasis`, levels=c('Low', 'High'))
# table(discovery_patchar$`loG.sigma.4.0.mm.3D glszm Small Area Low Gray Level Emphasis`)
# 
# table(discovery_patchar$`loG.sigma.1.0.mm.3D glcm Cluster Shade`)
# discovery_patchar$`loG.sigma.1.0.mm.3D glcm Cluster Shade` <- factor(discovery_patchar$`loG.sigma.1.0.mm.3D glcm Cluster Shade`, levels=c('Low', 'High'))
# table(discovery_patchar$`loG.sigma.1.0.mm.3D glcm Cluster Shade`)
# 
# table(discovery_patchar$`Original firstorder Kurtosis`)
# discovery_patchar$`Original firstorder Kurtosis` <- factor(discovery_patchar$`Original firstorder Kurtosis`, levels=c('Low', 'High'))
# table(discovery_patchar$`Original firstorder Kurtosis`)
# 
# table(discovery_patchar$`loG.sigma.2.0.mm.3D firstorder Median`)
# discovery_patchar$`loG.sigma.2.0.mm.3D firstorder Median` <- factor(discovery_patchar$`loG.sigma.2.0.mm.3D firstorder Median`, levels=c('Low', 'High'))
# table(discovery_patchar$`loG.sigma.2.0.mm.3D firstorder Median`)
# 
# table(discovery_patchar$`Original shape Least Axis Length`)
# discovery_patchar$`Original shape Least Axis Length` <- factor(discovery_patchar$`Original shape Least Axis Length`, levels=c('Low', 'High'))
# table(discovery_patchar$`Original shape Least Axis Length`)
# 
# table(discovery_patchar$`Original glcm Joint Energy`)
# discovery_patchar$`Original glcm Joint Energy` <- factor(discovery_patchar$`Original glcm Joint Energy`, levels=c('Low', 'High'))
# table(discovery_patchar$`Original glcm Joint Energy`)
# 
# table(discovery_patchar$`Original glcm Correlation`)
# discovery_patchar$`Original glcm Correlation` <- factor(discovery_patchar$`Original glcm Correlation`, levels=c('Low', 'High'))
# table(discovery_patchar$`Original glcm Correlation`)
# 
# #Reordering columns
# discovery_patchar <- discovery_patchar[, c(31, 23, 35, 22, 34, 33, 30, 24, 32, 29, 26, 25, 27, 28, 1, 2, 3, 4, 16, 17, 18, 19, 20, 21, 36, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)]
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
# print(tab2, showAllLevels = TRUE, formatOptions = list(big.mark = ","))
# # summary(tab2)
# 
# #Stratified groups
# tab3 <- CreateTableOne(vars = myVars, strata = "Recurrence status" , data = discovery_patchar)
# print(tab3, formatOptions = list(big.mark = ","))
# 
# #Exporting tableone
# tab3Mat <- print(tab3, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
# ## Save to a CSV file
# # write.csv(tab3Mat, file = gsub(" ","",paste(getwd(),"/results/pat_char/pat_char_REC_discovery.csv")))

# ##################################################################################################################################################
# # UNIVARIATE KAPLAN-MEYER IMAGES EXPORT
# ##################################################################################################################################################
# RECURRENCE
#Style
centr = theme_classic(base_family = 'Helvetica') + theme(plot.title = element_text(hjust = 0.5, size = 14))
# Initialize list
splots_REC <- list();
# 
# # Gender
# sopravv <- survfit(Surv(time,status)~as.factor(df_clinical.gender),data=discovery_sub); sopravv
# splots_REC[[1]] <- ggsurvplot(sopravv, data=discovery_sub, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability Recurrence",
#                               risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
#                               palette= c("#B3C99C","#008837"), legend.labs = c("Females","Males"),
#                               legend="top", risk.table.title=FALSE, title="Gender", font.x=14,font.y=14, pval.coord = c(0,0.1),
#                               font.tickslab=14,font.legend=12,pval.size=5,risk.table.height = FALSE, ggtheme = centr); splots_REC[[1]]
# # Grade
# sopravv <- survfit(Surv(time,status)~as.factor(df_clinical.tumor_grade_reclass),data=discovery_sub); sopravv
# splots_REC[[2]] <- ggsurvplot(sopravv, data=discovery_sub, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability Recurrence",
#                               risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
#                               palette= c("#008837","#B3C99C"), legend.labs = c("G1/G2","G3/G4"),
#                               legend="top", risk.table.title=FALSE, title="Grade", font.x=14,font.y=14, pval.coord = c(0,0.1),
#                               font.tickslab=14,font.legend=12,pval.size=5,risk.table.height = FALSE, ggtheme = centr); splots_REC[[2]]
# 
# ggforest(splots_REC[[1]], data=discovery_sub, main = "Probability Recurrence - Whole Set")

# # original_firstorder_Kurtosis_cate
# sopravv <- survfit(Surv(time,status)~as.factor(original_firstorder_Kurtosis_cate),data=discovery_sub); sopravv
# splots_radiomics_REC[[1]] <- ggsurvplot(sopravv, data=discovery_sub, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability Rec",
#                                        risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
#                                        palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE, 
#                                        title="Original firstorder Kurtosis", 
#                                        legend="none", font.x=12,font.y=12, pval.coord = c(0,0.1), font.title = 12, font.tickslab=12,pval.size=4,risk.table.height = FALSE, ggtheme = centr); splots_radiomics_REC[[1]]
# 
# # original_glcm_Correlation_cate
# sopravv <- survfit(Surv(time,status)~as.factor(original_glcm_Correlation_cate),data=discovery_sub); sopravv
# splots_radiomics_REC[[2]] <- ggsurvplot(sopravv, data=discovery_sub, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability Rec",
#                                        risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
#                                        palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE, 
#                                        title="Original glcm Corr", 
#                                        legend="none", font.x=12,font.y=12, pval.coord = c(0,0.1), font.title = 12, font.tickslab=12,pval.size=4,risk.table.height = FALSE, ggtheme = centr); splots_radiomics_REC[[2]]
# 
# original_shape_Elongation_cate
sopravv <- survfit(Surv(time,status)~as.factor(original_shape_Elongation_cate),data=discovery_sub); sopravv
splots_REC[[1]] <- ggsurvplot(sopravv, data=discovery_sub, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability REC",
                                       risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
                                       palette= c("#FFC8C8","#e66101"), legend.labs = c("high","low"),
                                       legend="top", risk.table.title=FALSE,
                                       title="Original Shape Elong",font.title=18,
                                       font.x=14,font.y=14, pval.coord = c(0,0.1),
                                       font.tickslab=14,font.legend=12,pval.size=5,risk.table.height = FALSE, ggtheme = centr); splots_REC[[1]]
# 
# # original_glcm_JointEnergy_cate
# sopravv <- survfit(Surv(time,status)~as.factor(original_glcm_JointEnergy_cate),data=discovery_sub); sopravv
# splots_radiomics_REC[[4]] <- ggsurvplot(sopravv, data=discovery_sub, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability Rec",
#                                        risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
#                                        palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE, 
#                                        title="Original glcm J E",
#                                        legend="none", font.x=12,font.y=12, pval.coord = c(0,0.1), font.title = 12, font.tickslab=12,pval.size=4,risk.table.height = FALSE, ggtheme = centr); splots_radiomics_REC[[4]]
# 
# # original_shape_LeastAxisLength_cate
# sopravv <- survfit(Surv(time,status)~as.factor(original_shape_LeastAxisLength_cate),data=discovery_sub); sopravv
# splots_radiomics_REC[[5]] <- ggsurvplot(sopravv, data=discovery_sub, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability Rec",
#                                         risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
#                                         palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE, 
#                                         title="Original shape Least AL",
#                                         legend="none", font.x=12,font.y=12, pval.coord = c(0,0.1), font.title = 12, font.tickslab=12,pval.size=4,risk.table.height = FALSE, ggtheme = centr); splots_radiomics_REC[[5]]
# 
# # log.sigma.2.0.mm.3D_glszm_SmallAreaLowGrayLevelEmphasis_cate
# sopravv <- survfit(Surv(time,status)~as.factor(log.sigma.2.0.mm.3D_glszm_SmallAreaLowGrayLevelEmphasis_cate),data=discovery_sub); sopravv
# splots_radiomics_REC[[6]] <- ggsurvplot(sopravv, data=discovery_sub, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability Rec",
#                                        risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
#                                        palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE, 
#                                        title="loG.sigma.2.0.mm.3D glszm SALGLE", 
#                                        legend="none", font.x=12,font.y=12, pval.coord = c(0,0.1), font.title = 12, font.tickslab=12,pval.size=4,risk.table.height = FALSE, ggtheme = centr); splots_radiomics_REC[[6]]
# 
# log.sigma.2.0.mm.3D_firstorder_Median_cate
sopravv <- survfit(Surv(time,status)~as.factor(log.sigma.2.0.mm.3D_firstorder_Median_cate),data=discovery_sub); sopravv
splots_REC[[2]] <- ggsurvplot(sopravv, data=discovery_sub, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability REC",
                                       risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
                                       palette= c("#FFC8C8","#e66101"), legend.labs = c("high","low"),
                                       legend="top", risk.table.title=FALSE,
                                       title="LoG.sig.2.0.mm Fo Med", font.title=18,
                                       font.x=14,font.y=14, pval.coord = c(0,0.1),
                                       font.tickslab=14,font.legend=12,pval.size=5,risk.table.height = FALSE, ggtheme = centr); splots_REC[[2]]
 
# Residual disease
sopravv <- survfit(Surv(time,status)~as.factor(df_clinical.residual_disease_reclass),data=discovery_sub); sopravv
splots_REC[[3]] <- ggsurvplot(sopravv, data=discovery_sub, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability REC",
                              risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
                              palette= c("#008837","#B3C99C","darkgray"), legend.labs = c("R1/R2","R0","RX"),
                              legend="top", risk.table.title=FALSE,
                              title="Residual Disease", font.title=18,
                              font.x=14,font.y=14, pval.coord = c(0,0.1),
                              font.tickslab=14,font.legend=12,pval.size=5,risk.table.height = FALSE, ggtheme = centr); splots_REC[[3]]


# # wavelet.LLH_firstorder_Variance_cate
# sopravv <- survfit(Surv(time,status)~as.factor(wavelet.LLH_firstorder_Variance_cate),data=discovery_sub); sopravv
# splots_radiomics_REC[[8]] <- ggsurvplot(sopravv, data=discovery_sub, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability Rec",
#                                        risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
#                                        palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE,
#                                        title="wavelet.LLH firstorder Var", 
#                                        legend="none", font.x=12,font.y=12, pval.coord = c(0,0.1), font.title = 12, font.tickslab=12,pval.size=4,risk.table.height = FALSE, ggtheme = centr); splots_radiomics_REC[[8]]
# 
# # wavelet.HHL_firstorder_Mean_cate
# sopravv <- survfit(Surv(time,status)~as.factor(wavelet.HHL_firstorder_Mean_cate),data=discovery_sub); sopravv
# splots_radiomics_REC[[9]] <- ggsurvplot(sopravv, data=discovery_sub, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability Rec",
#                                        risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
#                                        palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE, 
#                                        title="wavelet.HHL f Mean", 
#                                        legend="none", font.x=12,font.y=12, pval.coord = c(0,0.1), font.title = 12, font.tickslab=12,pval.size=4,risk.table.height = FALSE, ggtheme = centr); splots_radiomics_REC[[9]]
# 
# # wavelet.HHL_glszm_LargeAreaEmphasis_cate
# sopravv <- survfit(Surv(time,status)~as.factor(wavelet.HHL_glszm_LargeAreaEmphasis_cate),data=discovery_sub); sopravv
# splots_radiomics_REC[[10]] <- ggsurvplot(sopravv, data=discovery_sub, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability Rec",
#                                         risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
#                                         palette= c("#FFC8C8","#e66101"), risk.table.title=FALSE, 
#                                         title="wavelet.HHL glszm LAE", 
#                                         legend="none", legend.labs = c("High","Low"), font.legend=16,
#                                         font.x=12,font.y=12, pval.coord = c(0,0.1), font.title = 12, font.tickslab=12,pval.size=4,risk.table.height = FALSE, ggtheme = centr); splots_radiomics_REC[[10]]

# KRAS
sopravv <- survfit(Surv(time,status)~as.factor(KRAS),data=discovery_sub); sopravv
splots_REC[[4]] <- ggsurvplot(sopravv, data=discovery_sub, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability REC",
                              risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
                              palette= c("blue","#9BABB8"), legend.labs = c("Not mut","Mut"),
                              legend="top", risk.table.title=FALSE,
                              title="KRAS", font.title=18,
                              font.x=14,font.y=14, pval.coord = c(0,0.1),
                              font.tickslab=14,font.legend=12,pval.size=5,risk.table.height = FALSE, ggtheme = centr); splots_REC[[4]]

# # SMAD4
# sopravv <- survfit(Surv(time,status)~as.factor(SMAD4),data=discovery_sub); sopravv
# splots_REC[[5]] <- ggsurvplot(sopravv, data=discovery_sub, risk.table=FALSE, conf.int=TRUE, xlim=c(0,60), pval=TRUE,ylab="Probability REC",
#                               risk.table.y.col=FALSE, risk.table.y.text=FALSE, xlab="Time (months)", break.time.by=10,
#                               palette= c("blue","#9BABB8"), legend.labs = c("Not mut","Mut"),
#                               legend="top", risk.table.title=FALSE,
#                               title="SMAD4", font.title=18,
#                               font.x=14,font.y=14, pval.coord = c(0,0.1),
#                               font.tickslab=14,font.legend=12,pval.size=5,risk.table.height = FALSE, ggtheme = centr); splots_REC[[5]]

res <- arrange_ggsurvplots(splots_REC, print = F, nrow = 1, ncol = 4); res
ggsave(gsub(" ","",paste(getwd(),"/results/UV_analysis/Discovery/","UV_KM_radio_panel_REC.png")), res, width = 4000, height = 1000, units = "px", dpi = 300)

# ##################################################################################################################################################
# # RADIOMICS - MULTIVARIATE SURVIVAL ANALYSIS
# ##################################################################################################################################################

# Modello di cox multivariato
cox_radio <- coxph(formula=Surv(time, status) ~
               original_firstorder_Kurtosis_cate +
               # original_firstorder_90Percentile_cate +
               original_glcm_Correlation_cate +
               original_shape_Elongation_cate +
               original_glcm_JointEnergy_cate +
               original_shape_LeastAxisLength_cate +
               # log.sigma.1.0.mm.3D_glcm_ClusterShade_cate +
               log.sigma.2.0.mm.3D_glszm_SmallAreaLowGrayLevelEmphasis_cate +
               log.sigma.2.0.mm.3D_firstorder_Median_cate +
               # log.sigma.4.0.mm.3D_glszm_SmallAreaLowGrayLevelEmphasis_cate +
               wavelet.LLH_firstorder_Variance_cate +
               # wavelet.LHL_glcm_InverseVariance_cate +
               wavelet.HHL_firstorder_Mean_cate +
               wavelet.HHL_glszm_LargeAreaEmphasis_cate,
             data=discovery_sub, x =TRUE); summary(cox_radio)
# ggforest(cox, data=discovery_sub, main = "Probability Recurrence - Train Set")

# ##################################################################################################################################################
# # CLINICAL - MULTIVARIATE SURVIVAL ANALYSIS
# ##################################################################################################################################################

# Modello di cox multivariato
cox_clin <- coxph(formula=Surv(time, status) ~
             df_clinical.gender +
             df_clinical.tumor_grade_reclass +
             df_clinical.residual_disease_reclass,
             data=discovery_sub, x = TRUE); summary(cox_clin)
# ggforest(cox_clin, data=discovery_sub, main = "Probability Recurrence - Whole Set")

# ##################################################################################################################################################
# # MUTATIONS - MULTIVARIATE SURVIVAL ANALYSIS
# ##################################################################################################################################################

# Modello di cox multivariato
cox_mut <- coxph(formula=Surv(time, status) ~
               KRAS +
               SMAD4,
             data=discovery_sub, x = TRUE); summary(cox_mut)

# ##################################################################################################################################################
# #INTEGRATED FOREST PLOTS
# ##################################################################################################################################################

all.models <- list()
all.models[[1]] <- cox_radio; summary(all.models[[1]])
all.models[[2]] <- cox_clin; summary(all.models[[2]])
all.models[[3]] <-cox_mut; summary(all.models[[3]])

p <- plot_models(all.models,
                 grid = F,
                 axis.title = "HR (Probability REC)",
                 legend.title = NULL,
                 m.labels = c("Radiomics","Clinical","Mutational"),
                 axis.labels = c("SMAD4, mutated vs. not",
                                 "KRAS, mutated vs. not",
                                 "Residual Disease, R0 vs. R1/R2",
                                 "Grade, G3 vs. G1/G2",
                                 "Gender, M vs. F",
                                 "Wavelet.HHL GLSZM LArea Emphasis, low vs. high",
                                 "Wavelet.HHL Fo Mean, low vs. high",
                                 "Wavelet.LLH Fo Var, low vs. high",
                                 "LoG.sigma.2.0.mm.3D Fo Median, low vs. high",
                                 "LoG.sigma.2.0.mm.3D GLSZM SALGLE, low vs. high",
                                 "Original Shape LAL, low vs. high",
                                 "Original GLCM Jo Energy, low vs. high",
                                 "Original Shape Elongation, low vs. high",
                                 "Original GLCM Correlation, low vs. high", 
                                 "Original Fo Kurtosis, low vs. high"),
                 p.threshold = c(0.10, 0.05, 0.01),
                 wrap.labels = 200,
                 vline.color = "black",
                 line.size = 3,
                 dot.size = 4,
                 axis.lim = c(0.01,20),
                 show.values = TRUE,
                 value.size = 4,
                 spacing = 0.7,
                 colors = c("#2D5F8A","#06948E","#e66101"))
p + set_theme(base = theme_transparent(base_family = 'Helvetica'),
              legend.pos = "left",
              legend.size = 1.5,
              legend.item.size = 1.2,
              legend.item.backcol = "white",
              axis.title.size = 1.5,
              axis.textsize.y = 1.2,
              axis.textsize.x = 1.2,
              axis.title.color = "black",
              axis.textcolor = "black",
)

ggsave(gsub(" ","",paste(getwd(),"/results/UV_analysis/Discovery/","Forest_plot_REC.png")), p, width = 3100, height = 4000, units = "px", dpi = 300)


