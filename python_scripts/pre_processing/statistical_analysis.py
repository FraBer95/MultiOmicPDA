import os
import pandas as pd
from scipy.stats import chi2_contingency, mannwhitneyu, kruskal
from statsmodels.sandbox.stats.multicomp import multipletests
import seaborn as sns
import matplotlib.pyplot as plt



def statistical_test(df, label, path):

    labels_values_dict = {}
    labels_values_dict[label] = df[label].unique()
    n_label = len(df[label].unique())
    #df.drop(columns=[label], inplace=True)
    groups_dict = {}
    for index, (el, val) in enumerate(labels_values_dict.items()):
        groups_dict[str(index)] = []
        for i in range(len(val)):
            mask = df[el] == val[i]
            df_val = df[mask].copy()
            df_val.drop(columns=[el], inplace=True)
            groups_dict[str(index)].append(df_val)

    features = groups_dict['0'][0].columns

    p_val_dict = {"Features": [],
                  "P-Value": []
                  }

    for feature in features:
        if n_label == 2:
            stat, p = mannwhitneyu(groups_dict['0'][0][feature], groups_dict['0'][1][feature])
            p_val_dict["Features"].append(feature)
            p_val_dict["P-Value"].append(p)
            test = 'mann_whit_' + label
        elif n_label == 3:
            test = 'kruskal_' + label
            stat, p = kruskal(groups_dict['0'][0][feature], groups_dict['0'][1][feature],
                              groups_dict['0'][2][feature])
            p_val_dict["Features"].append(feature)
            p_val_dict["P-Value"].append(p)
        else:
            print("Number of values equal {} in {} not supported!".format(n_label, label))
            raise Exception

        if p < 0.05:
            print(f"Test for {feature} with respect to Target:")
            print(f"Statistic Value: {stat}")
            print(f"P-value: {p}")
            print("----------------")

        p_val_df = pd.DataFrame(p_val_dict)
        significant_f = p_val_df[p_val_df['P-Value'] < 0.05]
        print("Number of significant features: ", significant_f.shape[0])

        corrected_pvalue = multipletests(p_val_df['P-Value'], method='holm-sidak')[1]
        p_val_df['Corrected P-Value'] = corrected_pvalue[:len(p_val_df)]
        significant_con = p_val_df[p_val_df['Corrected P-Value'] < 0.05]
        not_sign_con = p_val_df[p_val_df['Corrected P-Value'] >= 0.05]
        if not significant_con.empty:
            print("Features Continue Significative:")
            print(significant_con)
            print("\n" + "-" * 50 + "\n")
        if not not_sign_con.empty:
            print("Features Continue Non Significative:")
            print(not_sign_con)
            print("\n" + "-" * 50 + "\n")

        significant_f = p_val_df[p_val_df['Corrected P-Value'] < 0.05]
        if significant_f.empty:
            print("There aren't any significant features for label {}!!".format(label))
        else:
            significant_f.to_csv(os.path.join(path, 'p_values_{}.csv'.format(test)), index=False)



