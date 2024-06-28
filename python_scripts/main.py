import pandas as pd

from pre_processing.utils import *
import os
from sklearn.model_selection import train_test_split
from pre_processing.statistical_analysis import statistical_test

from pre_processing.features_correlation import pearson_corr
from classification.classification import training_models
from pre_processing.plot_roc_curves import plot_roc_test


if __name__=='__main__':

    #experiment settings
    IV_stages = False
    multimodal = True
    loocv = True
    split_wrt_variable = "KRAS_TP53"
    label_list_full=["status", "tumor_grade_reclass"]
    target_list=["status"]

    radiomics_path = '../data/radiomics_data.csv' #put your radiomics data in data folder
    clinical_path = '../data/clinical_data.csv' #put your clinical data in data folder

    df_origin, df_with_dup = join_data(radiomics_path, clinical_path, IV_stages) #match between radiomics and clinical dataset

    print(f"Possible labels: {label_list_full} \n"
          f"Selected labels as target: {target_list}")

    for label in target_list:
        print(f"Launching experiment w.r.t. {label}")
        # create directory for each label
        path_with_label = os.path.join('../logs', label)
        if not IV_stages: path_with_label = os.path.join(path_with_label, "no_4stages")
        else: path_with_label = os.path.join(path_with_label, "4stages")
        os.makedirs(path_with_label, exist_ok=True)

        #stratified train test split

        features = df_origin.drop(columns=[split_wrt_variable], axis=1)

        X_train, X_test, y_train, y_test = train_test_split(features, df_origin[split_wrt_variable], test_size=0.20, stratify=df_origin[label], random_state=42)



        df_train = pd.concat([X_train, y_train], axis=1)
        df_test = pd.concat([X_test, y_test], axis=1)

        # for el in label_list_full:
        #     print("Values distribution for Training {}: {}".format(el, df_train[el].value_counts()))
        #     print("Values distribution for Test {}: {}".format(el, df_test[el].value_counts()))

        ids_train = df_train['patient_id']
        ids_test = df_test['patient_id']
        # df_train.to_csv(os.path.join(path_with_label, 'train.csv'), index=False)
        # df_test.to_csv(os.path.join(path_with_label, 'test.csv'), index=False)
        #values = df_test[label].value_counts()
        df_train_full = add_duplicates(df_with_dup, ids_train)
        y_train_full = df_train_full[label]

        if loocv:
            df_test = add_duplicates(df_with_dup, ids_test)
            y_test_full = df_test[label]
            df_train_full = pd.concat([df_train_full, df_test], axis=0)
            y_train_full = pd.concat([y_train_full, y_test_full], axis=0)

        other_target_labels = [item for item in label_list_full if item != label]

        y_train_full = df_origin[label]
        df_train_full = df_origin.drop(columns=other_target_labels)


        #statistical_test(df_train_full, label, path_with_label)



        feature_AUC = None
        training_features = pearson_corr(df_train_full, feature_AUC, path_with_label, IV_stages)
        training_features.append("patient_id")
        training_features.append("status")


        if multimodal:
            print("Loading mutations and clinical datasets...")
            clin_data = pd.read_csv('/Volumes/ExtremeSSD/pycharmProg/pancreas2/RadiomicsPanc/data/clinical/clin_data.csv')
            mut_data = pd.read_csv('/Volumes/ExtremeSSD/pycharmProg/pancreas2/RadiomicsPanc/data/RNA-Seq/mutations.csv')
            df_train = merge_datasets(df_train_full[training_features], clin_data, mut_data, IV_stages)
            # desired_column = df_train["patient_id"]
            #
            #
            # df_train.drop(columns=["patient_id"], inplace=True)
            #
            #
            # df_train.insert(0, "patient_id", desired_column)
            # desired_column = df_train["status"]
            # df_train.drop(columns=["status"], inplace=True)

            # df_train["status"] = desired_column

            #df_train.to_csv('../data/multimodal_dataset_IVStages.csv', index=False)
            y_train_full = df_train["status"]
            df_train.drop(columns=["status"], inplace=True)
            path_with_label = os.path.join(path_with_label, "clin_gen")
            os.makedirs(path_with_label, exist_ok=True)
            classifiers_list = training_models(df_train, y_train_full, loocv, path_with_label,
                                               radiomic_training=False)

        else: classifiers_list = training_models(df_train_full[training_features], y_train_full, loocv, path_with_label,
                                               radiomic_training=True)

        if not loocv:
            plot_roc_test(classifiers_list, df_test, training_features, label)