import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import json
from scipy.cluster import hierarchy

def correlogramma_withClustering(corr_mtx):
    linkage_matrix = hierarchy.ward(corr_mtx)
    mtx_size = corr_mtx.shape[0]
    fig_size = min(24, 16 * (mtx_size / 1200))

    # Plotta il dendrogramma del clustering gerarchico
    plt.figure(figsize=(fig_size, fig_size))
    dendrogram = hierarchy.dendrogram(linkage_matrix, labels=corr_mtx.columns, leaf_rotation=90)
    plt.title('Dendrogramma di clustering gerarchico')
    plt.xlabel('Variabili')
    plt.ylabel('Distanza')
    plt.show()
    plt.close()

    font_dict = {'family': 'Helvetica', 'size': 14}
    # Plotta il correlogramma con l'ordinamento basato sul clustering gerarchico
    # sns.clustermap(corr_mtx, method='ward', cmap='coolwarm', figsize=(fig_size, fig_size), row_cluster=False, col_cluster=False, linewidths=2)
    # plt.title('Correlogramma con clustering gerarchico')
    #
    # plt.show()
    # plt.close()
    #sns.clustermap(corr_mtx, method='ward', cmap='coolwarm', figsize=(fig_size, fig_size))
    #plt.title('Correlogramma con clustering gerarchico')
    #plt.savefig('/Volumes/ExtremeSSD/pycharmProg/pancreas2/RadiomicsPanc/correlogramma/CORRELATION.png', dpi=400)

def pearson_corr(df, features_auc, path, IV_stages):
    #df = df_t.copy()
    df = df.drop(columns=["patient_id", "Series.description_class", "tumor_grade"])
    corr_mat = df.corr()
    #corr_mat.to_csv(os.path.join(path, "corr_mat_df.csv"), index=False)
    print(corr_mat)
    correlogramma_withClustering(corr_mat)
    # Save Correlation Matrix
    mtx_size = corr_mat.shape[0]
    fig_size = min(24, 16 * (mtx_size / 1200))
    #plt.subplots(figsize=(fig_size, fig_size))
    plt.figure(figsize=(fig_size, fig_size))
    # sns.heatmap(corr_mat, annot=True, cmap='coolwarm', fmt=".2f", xticklabels=True, yticklabels=True)
    # plt.title("Correlation Matrix")
    heatmap = plt.imshow(corr_mat, cmap='coolwarm', interpolation='nearest')
    # plt.xticks(range(len(corr_mat.columns)), corr_mat.columns, rotation=90)
    # plt.yticks(range(len(corr_mat.columns)), corr_mat.columns)

    # Aggiungi la barra dei colori
    plt.colorbar(heatmap)
    print("Saving Correlation matrix...")
    #plt.show()
    #if features_auc is not None:
    #    plt.savefig(os.path.join(path, "correlation_matrix.png"))
    #else:  plt.savefig(os.path.join(path, "correlation_matrix_full.png"))

    # Simple Correlation
    if features_auc is not None:
        correlation_threshold = 0.5
    else: correlation_threshold = 0.3

    high_correlation_pairs = {
        "Feature1": [],
        "Feature2": [],
        "Correlation": [],

    }

    feature_dict = {"radiomic_signature" : []}
    for i in range(len(corr_mat.columns)):
        for j in range(i + 1, len(corr_mat.columns)):
            if abs(corr_mat.iloc[i, j]) > correlation_threshold:
                v1 = corr_mat.columns[i]
                v2 = corr_mat.columns[j]
                correlation_value = corr_mat.iloc[i, j]
                high_correlation_pairs["Feature1"].append(v1)
                high_correlation_pairs["Feature2"].append(v2)
                high_correlation_pairs["Correlation"].append(correlation_value)


    f_deleted = []
    feature_selected = []
    for f1, f2, corr in zip(high_correlation_pairs["Feature1"], high_correlation_pairs["Feature2"], high_correlation_pairs["Correlation"]):

        print("Coppia {} - {}: Correlation value = {}".format(f1, f2,corr))

        #se la coppia risulta correlata, prendo quella con la AUC più alta

        if f1 in f_deleted or f2 in f_deleted: #check if feature already deleted
            continue
        if features_auc is not None:
            auc_0 = round(features_auc.loc[features_auc["Features"] == f1, 'AUC'].values[0], 3)
            auc_1 = round(features_auc.loc[features_auc["Features"] == f2, 'AUC'].values[0], 3)

            if auc_0 >= auc_1:
                index = features_auc.index[features_auc['Features'] == f2].to_list()
                features_auc = features_auc.drop(index)
                f_deleted.append(f2)

            else:
                index = features_auc.index[features_auc['Features'] == f1].tolist()
                features_auc = features_auc.drop(index)
                f_deleted.append(f1)

        else:
            feature_selected.append(f1)
            f_deleted.append(f2)

    if features_auc is not None: return  features_auc['Features'].to_list()

    else:
        feature_selected = list(set(feature_selected))
        feature_dict["radiomic_signature"].append(feature_selected)
        feature_df = pd.DataFrame(feature_dict)
        if IV_stages:
            #feature_df.to_csv(os.path.join(path, 'radiomic_signature_with4stages.csv'), index=False)
            with open(os.path.join(path, 'radiomic_signature_with4stages.json'), "w") as json_file:
                json.dump(feature_dict, json_file)
        else:
            #feature_df.to_csv(os.path.join(path, 'radiomic_signature_no4stages.csv'), index=False)
            with open(os.path.join(path, 'radiomic_signature_no4stages.json'), "w") as json_file:
                json.dump(feature_dict, json_file)
        return feature_selected

