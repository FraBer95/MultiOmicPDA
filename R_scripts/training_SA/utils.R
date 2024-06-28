
feature_sel_OS <- function(df, R, C, G){
   if (R) {
      if (C) {
        if (G) {
             dta <- data[, c("time", "status", "original_firstorder_90Percentile_cate","original_shape_Elongation_cate","original_glcm_JointEnergy_cate","log.sigma.1.0.mm.3D_firstorder_Skewness_cate",
          "log.sigma.2.0.mm.3D_firstorder_Median_cate","wavelet.HLH_glcm_Imc1_cate",
          "df_clinical.gender","df_clinical.tumor_grade_reclass","df_clinical.residual_disease_reclass", "TTN")]
             return(list(dta, "Radiomics_Clinical_Genomic"))
        } else {
           dta <- df[, c("time", "status","original_firstorder_90Percentile_cate","original_shape_Elongation_cate","original_glcm_JointEnergy_cate",
                         "log.sigma.1.0.mm.3D_firstorder_Skewness_cate", "log.sigma.2.0.mm.3D_firstorder_Median_cate","wavelet.HLH_glcm_Imc1_cate",
                          "df_clinical.gender","df_clinical.tumor_grade_reclass",
                         "df_clinical.residual_disease_reclass")]
          return(list(dta, "Radiomics_Clinical"))
        }
      } else if (G) {
             dta <- df[, c("time", "status", "original_firstorder_90Percentile_cate","original_shape_Elongation_cate",
                           "original_glcm_JointEnergy_cate","log.sigma.1.0.mm.3D_firstorder_Skewness_cate",
                           "log.sigma.2.0.mm.3D_firstorder_Median_cate","wavelet.HLH_glcm_Imc1_cate","TTN")]
          return(list(dta, "Radiomics_Genomic"))
      } else {
          dta <- df[, c("time", "status", "original_firstorder_Kurtosis_cate","log.sigma.1.0.mm.3D_gldm_DependenceVariance_cate","original_shape_Elongation_cate",
          "original_glcm_JointEnergy_cate","log.sigma.1.0.mm.3D_firstorder_Skewness_cate","wavelet.LLH_glcm_ClusterTendency_cate",
            "wavelet.HLH_glcm_Imc1_cate","wavelet.HHL_firstorder_Mean_cate","wavelet.HHL_glcm_Imc1_cate","wavelet.HHL_glszm_LargeAreaEmphasis_cate",
                  "log.sigma.2.0.mm.3D_firstorder_Median_cate","wavelet.HLL_glcm_InverseVariance_cate","original_firstorder_90Percentile_cate" )]
        return(list(dta, "Radiomics"))
      }
    } else if (C) {
      if (G) {
             dta <- df[, c("time", "status","df_clinical.gender","df_clinical.tumor_grade_reclass",
                           "df_clinical.residual_disease_reclass","TTN")]
              return(list(dta, "Clinical_Genomic"))
      } else {
            dta <- df[, c("time", "status", "df_clinical.gender","df_clinical.tumor_grade_reclass","df_clinical.residual_disease_reclass")]
            return(list(dta, "Clinical"))
      }
    } else if (G) {
       dta <- df[, c("time", "status", "KRAS", "TTN")]
      return(list(dta, "Genomic"))

    }
}

check_binary_columns <- function(df) {
  results <- sapply(df, function(col) {
    num_zeros <- sum(col == 0)
    num_ones <- sum(col == 1)
    return(num_zeros >= 2 && num_ones >= 2)
  })

  return(results)
}



feature_sel_REC <- function(df, R, C, G){
   if (R) {
      if (C) {
        if (G) {
             dta <- data[, c("time", "status", "original_shape_Elongation_cate","log.sigma.2.0.mm.3D_firstorder_Median_cate",
                             "df_clinical.residual_disease_reclass", "KRAS","SMAD4")]
             return(list(dta, "Radiomics_Clinical_Genomic"))
        } else {
           dta <- df[, c("time", "status","original_shape_Elongation_cate",
                         "log.sigma.2.0.mm.3D_firstorder_Median_cate", "df_clinical.residual_disease_reclass")]
          return(list(dta, "Radiomics_Clinical"))
        }
      } else if (G) {
             dta <- df[, c("time", "status", "original_shape_Elongation_cate",
                           "log.sigma.2.0.mm.3D_firstorder_Median_cate", "KRAS","SMAD4")]
          return(list(dta, "Radiomics_Genomic"))
      } else {
          dta <- df[, c("time", "status","original_firstorder_Kurtosis_cate","original_glcm_Correlation_cate","original_shape_Elongation_cate",
                        "original_glcm_JointEnergy_cate","original_shape_LeastAxisLength_cate", "log.sigma.2.0.mm.3D_glszm_SmallAreaLowGrayLevelEmphasis_cate",
                        "log.sigma.2.0.mm.3D_firstorder_Median_cate", "wavelet.LLH_firstorder_Variance_cate","wavelet.HHL_firstorder_Mean_cate","wavelet.HHL_glszm_LargeAreaEmphasis_cate")]
        return(list(dta, "Radiomics"))
      }
    } else if (C) {
      if (G) {
             dta <- df[, c("time", "status", "df_clinical.residual_disease_reclass","KRAS","SMAD4")]
              return(list(dta, "Clinical_Genomic"))
      } else {
            dta <- df[, c("time", "status", "df_clinical.gender","df_clinical.tumor_grade_reclass","df_clinical.residual_disease_reclass")]
            return(list(dta, "Clinical"))
      }
    } else if (G) {
       dta <- df[, c("time", "status","KRAS","SMAD4")]
      return(list(dta, "Genomic"))

    }
}






save_ws <- function(model, RADIOMICS, CLINICAL, GENOMIC, OS) {

   if (OS) { base_path <- file.path('./Radiomics_result', 'OS')}
    else {base_path <- file.path('./Radiomics_result', 'REC')}

  if (model == 'survGBM'){
      posizione <- file.path(base_path, "GBM")
    } else if (model == 'coxPH'){
      posizione <- file.path(base_path, "coxPH")
    } else if (model == 'survRF'){
      posizione <- file.path(base_path, "survRF")
    } else {posizione <- file.path(base_path,"SVM") }

    if (RADIOMICS) {
      if (CLINICAL) {
        if (GENOMIC) {
          workspace_nome <- "Radiomics_Clinical_Genomic_workspace.RData"
        } else {
          workspace_nome <- "Radiomics_Clinical_workspace.RData"
        }
      } else if (GENOMIC) {
        workspace_nome <- "Radiomics_Genomic_workspace.RData"
      } else {
        workspace_nome <- "Radiomics_workspace.RData"
      }
    } else if (CLINICAL) {
      if (GENOMIC) {
        workspace_nome <- "Clinical_Genomic_workspace.RData"
      } else {
        workspace_nome <- "Clinical_workspace.RData"
      }
    } else if (GENOMIC) {
      workspace_nome <- "Genomic_workspace.RData"
    } else {

      workspace_nome <- "Default_workspace.RData"
    }

    save.image(file.path(posizione, workspace_nome))

    print(paste("Workspace salvato come:", workspace_nome))
}