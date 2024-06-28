import pandas as pd

from pre_processing.utils import *
import os
from sklearn.model_selection import train_test_split
from pre_processing.statistical_analysis import statistical_test

from pre_processing.features_correlation import pearson_corr

from pre_processing.plot_roc_curves import plot_roc_test


if __name__=='__main__':

    #experiment settings
    IV_stages = False
    multimodal = True
    loocv = True
    split_wrt_variable = "KRAS_TP53"
    label_list_full=["status", "tumor_grade_reclass"]
    target_list=["status"] #for multiple training at once

    radiomics_path = '../data/radiomics_data.csv' #put your radiomics data in data folder
    clinical_path = '../data/clinical_data.csv' #put your clinical data in data folder

    df_origin, df_with_dup = join_data(radiomics_path, clinical_path, IV_stages) #match between radiomics and clinical dataset

    print(f"Possible labels: {label_list_full} \n"
          f"Selected labels as target: {target_list}")

    for label in target_list:
        print(f"Launching experiment w.r.t. {label}")
        #create directory for each label
        path_with_label = os.path.join('../logs', label)
        if not IV_stages: path_with_label = os.path.join(path_with_label, "no_4stages")
        else: path_with_label = os.path.join(path_with_label, "4stages")
        os.makedirs(path_with_label, exist_ok=True)

        # for el in label_list_full:
        #     print("Values distribution for Training {}: {}".format(el, df_train[el].value_counts()))
        #     print("Values distribution for Test {}: {}".format(el, df_test[el].value_counts()))

        other_target_labels = [item for item in label_list_full if item != label]


        df_train_full = df_origin.drop(columns=other_target_labels)


        #statistical_test(df_train_full, label, path_with_label)

        training_features = pearson_corr(df_train_full, path_with_label, IV_stages)
        training_features.append("patient_id")
        training_features.append("status")
