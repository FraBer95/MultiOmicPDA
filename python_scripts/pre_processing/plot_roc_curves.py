import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve
from sklearn.metrics import roc_auc_score
import pandas as pd
def plot_roc_test(classifiers, test, features, label):

    for classifier in classifiers:
        classifier_name= classifier.__class__.__name__
        X_test = test[features].values
        y_test = test[[label]].values

        fpr = dict()
        tpr = dict()
        roc_auc = dict()

        n_classes = 2

        y_score = classifier.predict_proba(X_test)

        for i in range(n_classes):
            if i > 0:
                fpr[i], tpr[i], _ = roc_curve(y_test[:, 0], y_score[:, i])
                roc_auc[i] = auc(fpr[i], tpr[i])

                plt.plot(
                    fpr[i],
                    tpr[i],
                    lw=2,
                    label="ROC curve for class {} (area = {:.2f})".format(i, roc_auc[i]),
                )

        plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title("Receiver Operating Characteristic for Each Class with {} model".format(classifier_name))
        plt.legend(loc="lower right")
        plt.show()
        plt.close()

        for i in range(n_classes):
            precision, recall, _ = precision_recall_curve(y_test[:, 0], y_score[:, i])
            pr_auc = auc(recall, precision)
            plt.plot(
                recall,
                precision,
                lw=2,
                label="PR curve for class {} (PR AUC = {:.2f})".format(i, pr_auc),
            )

        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.title("Precision-Recall Curve")
        plt.legend(loc="lower right")

        plt.tight_layout()
        plt.show()
