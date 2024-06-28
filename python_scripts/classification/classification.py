import pandas as pd
from sklearn.metrics import RocCurveDisplay, auc, roc_curve
from sklearn.model_selection import StratifiedKFold, LeaveOneOut
import numpy as np
import matplotlib.pyplot as plt
import xgboost as xgb
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
import os

def training_models(X_train, y_train, loocv, path, radiomic_training=True):

    random_state = 42
    rf_classifier = RandomForestClassifier(n_estimators=500, random_state=random_state)
    xgb_classifier = xgb.XGBClassifier(n_estimators=500)
    classifiers = [rf_classifier, xgb_classifier]
    trained_classifiers = []

    if not loocv: #if not leave-one-out
        n_splits = 10
        cv = StratifiedKFold(n_splits=n_splits)

        for classifier in classifiers:
            classifier_name = classifier.__class__.__name__
            print("Training classifier: {} with {}-fold CV...".format(classifier_name, n_splits))
            tprs = []
            aucs = []
            mean_fpr = np.linspace(0, 1, 100)

            fig, ax = plt.subplots(figsize=(6, 6))
            for fold, (train, val) in enumerate(cv.split(X_train, y_train)):

                if isinstance(X_train, pd.DataFrame):
                    X_train_f, y_train_f = X_train.iloc[train].values, y_train.iloc[train].values
                    X_val, y_val = X_train.iloc[val].values, y_train.iloc[val].values

                classifier.fit(X_train_f, y_train_f.ravel())
                viz = RocCurveDisplay.from_estimator(
                    classifier,
                    X_val,
                    y_val,
                    name=f"ROC fold {fold}",
                    alpha=0.3,
                    lw=1,
                    ax=ax,
                    plot_chance_level=(fold == n_splits - 1),
                )

                interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
                interp_tpr[0] = 0.0
                tprs.append(interp_tpr)
                aucs.append(viz.roc_auc)

            mean_tpr = np.mean(tprs, axis=0)
            mean_tpr[-1] = 1.0
            mean_auc = auc(mean_fpr, mean_tpr)
            std_auc = np.std(aucs)
            ax.plot(
                mean_fpr,
                mean_tpr,
                color="b",
                label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
                lw=2,
                alpha=0.8,
            )

            std_tpr = np.std(tprs, axis=0)
            tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
            tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
            ax.fill_between(
                mean_fpr,
                tprs_lower,
                tprs_upper,
                color="grey",
                alpha=0.2,
                label=r"$\pm$ 1 std. dev.",
            )

            ax.set(
                xlim=[-0.05, 1.05],
                ylim=[-0.05, 1.05],
                xlabel="False Positive Rate",
                ylabel="True Positive Rate",
                title=f"Mean ROC curve with variability\n (Positive label 'status') \n for {classifier_name}"
    )
            ax.axis("square")
            ax.legend(loc="lower right")
            plt.show()
            #plt.savefig(os.path.join(path, classifier_name+"_crossVal.png"))

            trained_classifiers.append(classifier)

        return trained_classifiers

    else:

        loo = LeaveOneOut()
        trained_classifiers = []

        for classifier in classifiers:
            classifier_name = classifier.__class__.__name__
            print("Training classifier: {}, with LOOCV...".format(classifier_name))
            prob_list = []
            gt_list = []


            for fold, (train_index, val_index) in enumerate(loo.split(X_train, y_train)):

                # X_training_init = X_train.iloc[train_index]
                y_tr = y_train.iloc[train_index]
                full_df = pd.concat([X_train.iloc[train_index], y_tr], axis=1)

                X_val = X_train.iloc[val_index]
                y_val = y_train.iloc[val_index]
                # Find duplicates in validation set within training set
                val_dup = full_df[full_df["patient_id"].isin(X_val["patient_id"])]

                if not val_dup.empty: #se ho dei duplicati
                    X_val = pd.concat([X_val, val_dup], axis=0) #inserisco i duplicati del validation
                    X_training = full_df[~full_df["patient_id"].isin(val_dup["patient_id"])] #elimino i duplicati dal training

                    # Estrai le righe corrispondenti dagli indici comuni
                    y_tr = X_training[full_df.columns[-1]] #elimino i duplicati dalle labels del training set
                    n_reps = val_dup.shape[0]
                    y_ext = pd.concat([y_val] * n_reps, axis=0) #estendo la label del validation
                    y_val = pd.concat([y_val, y_ext], axis=0)
                    X_training = X_training.iloc[:, :-2]
                    X_val = X_val.iloc[:, :-2]

                if val_dup.empty:
                    X_val.drop(columns=["patient_id"], inplace=True)
                    if radiomic_training:
                        X_training = full_df.drop(columns=["patient_id", "status"])
                    else:
                        X_training = full_df.drop(columns=["patient_id", "status"])


                x = X_training.values

                classifier.fit(x, y_tr)
                prob = classifier.predict_proba(X_val.values)[:, 1]
                prob_list.append(prob)
                gt_list.append(y_val.values)

            y_true = np.concatenate(gt_list)
            predictions = np.concatenate(prob_list)
            fpr, tpr, thresholds = roc_curve(y_true, predictions)

            # Calcola l'area sotto la curva ROC (AUC)
            roc_auc = auc(fpr, tpr)

            # Plotta la curva ROC
            plt.figure()
            plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
            plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title(f'Receiver Operating Characteristic (ROC) Curve for {classifier_name}')
            plt.legend(loc="lower right")
            plt.show()
            #plt.savefig(os.path.join(path, classifier_name+"_loocv.png"))

            trained_classifiers.append(classifier)

        return trained_classifiers
