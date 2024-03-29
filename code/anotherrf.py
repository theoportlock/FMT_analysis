#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Random forest analysis of metadata
Theo Portlock
'''
import functions
from scipy import stats
from scipy.spatial import distance
from scipy.stats import spearmanr
from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (accuracy_score, confusion_matrix, classification_report)
from sklearn.metrics import plot_confusion_matrix
from sklearn.metrics import plot_precision_recall_curve
from sklearn.metrics import plot_roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import seaborn as sns
import shap
import skbio

taxo = pd.read_csv('../data/gutTaxo.csv', index_col=0)
msp = pd.read_csv('../data/vect_atlas.csv', index_col=0).T
allmeta = pd.read_csv('../metaAtlasFinish/final/basicMetadata.Theo.SL.2022.06.12.txt', sep='\t',  index_col=0)
hgmi = pd.read_csv('../metaAtlasFinish/finaltheo/hgmiedges.csv', index_col=1)
meta = hgmi.join(allmeta).reset_index().set_index('hgmi_id')

var = 'Disease'

shaps = pd.DataFrame()
#scores = pd.Series()
hgmaid = meta.index.unique()[0]
for hgmaid in meta.index.unique():
    df = msp.join(meta.loc[hgmaid].set_index('sample.ID')[var],how='inner').set_index(var)
    df = pd.DataFrame(StandardScaler().fit_transform(df), index=df.index, columns=df.columns).reset_index()
    model = RandomForestClassifier(
        bootstrap=True,
        ccp_alpha=0.0,
        class_weight='balanced_subsample',
        criterion='gini',
        max_depth=7,
        max_features='sqrt',
        max_leaf_nodes=None,
        max_samples=None,
        min_impurity_decrease=0.0005,
        min_samples_leaf=3,
        min_samples_split=10,
        min_weight_fraction_leaf=0.0,
        n_estimators=300,
        n_jobs=-1,
        oob_score=False,
        random_state=2,
        verbose=0,
        warm_start=False)
    X = df.drop(var, axis=1)
    y = pd.get_dummies(df.xs(var, axis=1))['Healthy']
    X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y,random_state = 1)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    explainer = shap.TreeExplainer(model)
    shaps_values = explainer(X)
    meanabsshap = pd.Series(np.abs(shaps_values.values).mean(axis=0)[:, 0], index=X.columns)
    corrs = [spearmanr(shaps_values.values[:, :, 0][:, x], X.iloc[:,x])[0] for x in range(len(X.columns))]
    final = meanabsshap * np.sign(corrs)
    final.fillna(0, inplace=True)
    shaps[hgmaid] = final
    #scores[hgmaid] = classification_report(y_true=y_test, y_pred=y_pred) + '\n\nAUCROC = ' + str(roc_auc_score(y_test, model.predict_proba(X_test)[:,1]))
    with open(f"../results/{hgmaid}.txt", "w") as f: f.write(classification_report(y_true=y_test, y_pred=y_pred) + '\n\nAUCROC = ' + str(roc_auc_score(y_test, model.predict_proba(X_test)[:,1])))

final.to_frame().join(taxo).sort_values(0).set_index("species").tail(20)[0].plot.barh()
plt.tight_layout(); plt.show()

plotdf = shaps[(np.abs(shaps) > 0.01).any(axis=1)]

plotdf = plotdf.join(taxo['species']).set_index('species').T.join(meta.loc[meta.Disease != "Healthy"].groupby(level=0).first()['Disease']).set_index('Disease')
functions.clustermap(plotdf)






from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import auc
from sklearn.metrics import roc_curve
from sklearn.model_selection import KFold
from sklearn.model_selection import LeaveOneOut 
from sklearn.model_selection import cross_validate
import matplotlib.pyplot as plt
import numpy as np

kf = KFold(10)
#kf = LeaveOneOut()

tprs = []
aucs = []
base_fpr = np.linspace(0, 1, 101)
plt.figure(figsize=(5, 5))
plt.axes().set_aspect('equal', 'datalim')
for i, (train, test) in enumerate(kf.split(X,y)):
    results = model.fit(X.iloc[train], y.iloc[train])
    y_score = model.predict_proba(X.iloc[test])
    fpr, tpr, _ = roc_curve(y.iloc[test], y_score[:, 1])
    plt.plot(fpr, tpr, 'b', alpha=0.15)
    tpr = np.interp(base_fpr, fpr, tpr)
    tpr[0] = 0.0
    aucs.append(auc(fpr, tpr))
    tprs.append(tpr)

tprs = np.array(tprs)
mean_tprs = tprs.mean(axis=0)
std = tprs.std(axis=0)

mean_auc = auc(base_fpr, mean_tpr)
std_auc = np.std(aucs)

tprs_upper = np.minimum(mean_tprs + std, 1)
tprs_lower = mean_tprs - std


plt.plot(base_fpr, mean_tprs, 'b', label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc))
plt.fill_between(base_fpr, tprs_lower, tprs_upper, color='grey', alpha=0.3, label=r"$\pm$ 1 std. dev.",)

plt.plot([0, 1], [0, 1],'r--')
plt.xlim([-0.01, 1.01])
plt.ylim([-0.01, 1.01])
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.show()










import matplotlib.pyplot as plt
from sklearn import svm
from sklearn.metrics import auc
from sklearn.metrics import RocCurveDisplay
from sklearn.model_selection import StratifiedKFold

# Run classifier with cross-validation and plot ROC curves
cv = StratifiedKFold(n_splits=6)
classifier = model
tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)
fig, ax = plt.subplots()
for i, (train, test) in enumerate(cv.split(X, y)):
    classifier.fit(X.iloc[train], y.iloc[train])
    viz = RocCurveDisplay.from_estimator(
        classifier,
        X.iloc[test],
        y.iloc[test],
        name="ROC fold {}".format(i),
        alpha=0.3,
        lw=1,
        ax=ax,
    )
    interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
    interp_tpr[0] = 0.0
    tprs.append(interp_tpr)
    aucs.append(viz.roc_auc)
ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)
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
    title="Receiver operating characteristic example",
)
ax.legend(loc="lower right")
plt.show()

