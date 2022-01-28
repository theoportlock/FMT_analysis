#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Random forest analysis of metadata
Theo Portlock
'''
%autoindent
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.spatial import distance
import skbio

def evaluate(model, test_features, test_labels):
    predictions = model.predict(test_features)
    errors = abs(predictions - test_labels)
    mape = 100 * np.mean(errors / test_labels)
    accuracy = 100 - mape
    print('Model Performance')
    print('Average Error: {:0.4f} degrees.'.format(np.mean(errors)))
    print('Accuracy = {:0.2f}%.'.format(accuracy))
    return accuracy

variable = 'MELD'
samples_metadata = pd.read_csv('metadata.csv', index_col=0).dropna(axis=0)
#samples_gutmsp = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0).T
msp_samples = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)
#samples_gutmeta = samples_gutmsp.join(samples_metadata[variable]).dropna()
msp_taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
taxaType='species'
taxonomy_samples = msp_samples.join(msp_taxonomy[taxaType], how='inner').groupby(taxaType).sum()
samples_taxonomy = taxonomy_samples.T
#box_pairs = list(itertools.combinations([(a,b) for a in samples_taxonomy.columns for b in samples_taxonomy.columns],2))

oralmsp_samples = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0)
oralmsp_taxonomy = pd.read_csv("../oral_downstream_data/oraltaxo.csv", index_col=0)
samples_oraltaxonomy = oralmsp_samples.join(oralmsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
omspdf = samples_oraltaxonomy.join(samples_metadata[variable], how='inner').astype('float')

mspdf = samples_taxonomy.join(samples_metadata[variable], how='inner').astype('float')
#df = pd.get_dummies(samples_metadata)

samples_sm = pd.read_csv("antismashNorm.csv", index_col=0).T
smdf = samples_sm.join(samples_metadata[variable]).dropna()

samples_patric = pd.read_csv("patricNorm.csv", index_col=0).T
padf = samples_patric.join(samples_metadata[variable]).dropna()

samples_pfam = pd.read_csv("pfamNorm.csv", index_col=0).T
pfdf = samples_pfam.join(samples_metadata[variable]).dropna()

samples_card = pd.read_csv("cardNorm.csv", index_col=0).T
cadf = samples_card.join(samples_metadata[variable]).dropna()

samples_cazyme = pd.read_csv("cazymeNorm.csv", index_col=0).T
czdf = samples_cazyme.join(samples_metadata[variable]).dropna()

'''
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
pd.DataFrame(columns=df.columns, data=LabelEncoder().fit_transform(df.values.flatten()).reshape(df.shape))
samples_gutmeta['Type'] = LabelEncoder().fit_transform(samples_gutmeta.Type)
ohe = OneHotEncoder()
meta_ohe = ohe.fit_transform(df)
#clf = RandomForestClassifier()
'''

def runmodel(df):
    #reg = RandomForestRegressor()
    reg = ExtraTreesRegressor()
    X = df.drop(variable, axis=1)
    y = df.xs(variable, axis=1)
    X_train, X_test, y_train, y_test = train_test_split(X, y)
    reg.fit(X_train, y_train)
    return evaluate(reg, X_test, y_test)

final = pd.Series()
final['Gut Species'] = runmodel(mspdf)
final['Oral Species'] = runmodel(omspdf)
final['Secondary Metabolites'] = runmodel(smdf)
final['Virulence Factors'] = runmodel(padf)
final['Protein Family'] = runmodel(pfdf)
final['Antimircobial Resistance'] = runmodel(cadf)
final['Cazymes'] = runmodel(czdf)
final.sort_values(ascending=False, inplace=True)
sns.barplot(x=final, y=final.index)
plt.xlabel('Model Accuracy (%)')
plt.xlim(0,100)
plt.ylabel('Dataset used for MELD prediction')
plt.show()
plt.tight_layout()
plt.savefig('results/OvG_randomfores_MELD.pdf')

'''
y_pred = clf.predict(X_test)

#This is jaccard
print("Accuracy RFC:",metrics.accuracy_score(y_test, y_pred))
print("Accuracy MSE:",metrics.mean_squared_error(y_test, y_pred))
print(metrics.confusion_matrix(y_test, y_pred))

evaluate(clf, X_test, y_test)

from sklearn.model_selection import GridSearchCV
# Create the parameter grid based on the results of random search
param_grid = {
    'bootstrap': [True],
    'max_depth': [80, 90, 100, 110],
    'max_features': [2, 3],
    'min_samples_leaf': [3, 4, 5],
    'min_samples_split': [8, 10, 12],
    'n_estimators': [100, 200, 300, 1000]
}
# Create a based model
rf = RandomForestRegressor()
# Instantiate the grid search model
grid_search = GridSearchCV(estimator = rf, param_grid = param_grid,
                          cv = 3, n_jobs = -1, verbose = 2)

grid_search.fit(X_train, y_train)
best_grid = grid_search.best_estimator_
grid_accuracy = evaluate(best_grid, X_test, y_test)
feature_imp = pd.Series(best_grid.feature_importances_,index=X.columns).sort_values(ascending=False)
sns.barplot(x=feature_imp, y=feature_imp.index)
plt.xlabel('Feature Importance Score')
#plt.ylabel('Metadata')
plt.savefig("meta_MELD_prediction.pdf")
'''

reg = RandomForestRegressor()
X = czdf.drop(variable, axis=1)
y = czdf.xs(variable, axis=1)
X_train, X_test, y_train, y_test = train_test_split(X, y)
reg.fit(X_train, y_train)
feature_imp = pd.Series(reg.feature_importances_,index=X.columns).sort_values(ascending=False)
sns.barplot(x=feature_imp, y=feature_imp.index)

'''
Ar_dist = distance.squareform(distance.pdist(czdf.drop('MELD', axis=1), metric="braycurtis"))
DM_dist = skbio.stats.distance.DistanceMatrix(Ar_dist)
PCoA = skbio.stats.ordination.pcoa(DM_dist)
plot = czdf.join(PCoA.samples[['PC1','PC2']].set_index(cadf.index))
sns.scatterplot(data=plot, x='PC1', y='PC2', hue='MELD', palette='coolwarm')

scaledData = StandardScaler().fit_transform(czdf.drop('MELD', axis=1))
pca = PCA(n_components=2)
result = pca.fit_transform(scaledData)
print(result)
czdf[['PC1', 'PC2']] = result
sns.scatterplot(data=czdf, x='PC1', y='PC2', hue='MELD', palette='coolwarm')

correlation = mspdf.corr(method='spearman')
columns = correlation.nlargest(10, 'MELD').index
correlation_map = np.corrcoef(mspdf[columns].values.T)
sns.set(font_scale=1.0)
heatmap = sns.heatmap(correlation_map, cbar=True, annot=True, square=True, fmt='.2f', yticklabels=columns.values, xticklabels=columns.values)
plt.show()

from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Lasso
from sklearn.linear_model import ElasticNet
from sklearn.tree import DecisionTreeRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.feature_selection import RFE
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
pipelines = []
pipelines.append(('ScaledLR', Pipeline([('Scaler', StandardScaler()),('LR',LinearRegression())])))
pipelines.append(('ScaledLASSO', Pipeline([('Scaler', StandardScaler()),('LASSO', Lasso())])))
pipelines.append(('ScaledEN', Pipeline([('Scaler', StandardScaler()),('EN', ElasticNet())])))
pipelines.append(('ScaledKNN', Pipeline([('Scaler', StandardScaler()),('KNN', KNeighborsRegressor())])))
pipelines.append(('ScaledCART', Pipeline([('Scaler', StandardScaler()),('CART', DecisionTreeRegressor())])))
pipelines.append(('ScaledGBM', Pipeline([('Scaler', StandardScaler()),('GBM', GradientBoostingRegressor())])))
pipelines.append(('ScaledRFR', Pipeline([('Scaler', StandardScaler()),('RFR', RandomForestRegressor())])))
results = []
names = []
for name, model in pipelines:
    kfold = KFold(n_splits=10)
    cv_results = cross_val_score(model, X_train, y_train, cv=kfold, scoring='neg_mean_squared_error')
    results.append(cv_results)
    names.append(name)
    msg = "%s: %f (%f)" % (name, cv_results.mean(), cv_results.std())
    print(msg)
scaler = StandardScaler().fit(X_train)
rescaledX = scaler.transform(X_train)
param_grid = dict(n_neighbors=np.array([2,3,4,5,6]))
model = KNeighborsRegressor()
kfold = KFold(n_splits=10)
grid = GridSearchCV(estimator=model, param_grid=param_grid, scoring='neg_mean_squared_error', cv=kfold)
grid_result = grid.fit(rescaledX, y_train)
means = grid_result.cv_results_['mean_test_score']
stds = grid_result.cv_results_['std_test_score']
params = grid_result.cv_results_['params']
for mean, stdev, param in zip(means, stds, params):
    print("%f (%f) with: %r" % (mean, stdev, param))
print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_))
from sklearn.metrics import mean_squared_error
scaler = StandardScaler().fit(X_train)
rescaled_X_train = scaler.transform(X_train)
model = KNeighborsRegressor(n_neighbors=5)
model.fit(rescaled_X_train, y_train)
# transform the validation dataset
rescaled_X_test = scaler.transform(X_test)
predictions = model.predict(rescaled_X_test)
print (mean_squared_error(y_test, predictions))
compare = pd.DataFrame({'Prediction': predictions, 'Test Data' : y_test})
mcompare = compare.reset_index().melt(id_vars=['index'])

'''
import pyforest
from lazypredict.Supervised import LazyRegressor
from pandas.plotting import scatter_matrix
# Scikit-learn packages
from sklearn.linear_model import LinearRegression
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import ExtraTreesRegressor
from sklearn import metrics
from sklearn.metrics import mean_squared_error
df = mspdf.copy()
X = df.drop(variable, axis=1)
y = df.xs(variable, axis=1)
X_train, X_test, y_train, y_test = train_test_split(X, y)
reg = LazyRegressor(ignore_warnings=False, custom_metric=None)
models, predictions = reg.fit(X_train, X_test, y_train, y_test)
print(models)

def rmse(model, y_test, y_pred, X_train, y_train):
    r_squared = model.score(X_test, y_test)
    mse = mean_squared_error(y_test, y_pred)
    rmse = np.sqrt(mse)
    print('R-squared: ' + str(r_squared))
    print('Mean Squared Error: '+ str(rmse))

def scatter_plot(y_test, y_pred, model_name):
    plt.figure(figsize=(10,6))
    sns.residplot(y_test, y_pred, lowess=True, color='#4682b4',
              line_kws={'lw': 2, 'color': 'r'})
    plt.title(str('MELD vs Residuals for '+ model_name))
    plt.xlabel('MELD',fontsize=16)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.show()

rmse(hist, y_test, y_pred, X_train, y_train)
scatter_plot(y_test, y_pred, 'Histogram-based Extra Trees Regressor')

'''
# first neural network with keras tutorial
from numpy import loadtxt
from keras.models import Sequential
from keras.layers import Dense
# load the dataset
dataset = loadtxt('pima-indians-diabetes.csv', delimiter=',')
# split into input (X) and output (y) variables
X = dataset[:,0:8]
y = dataset[:,8]
# define the keras model
model = Sequential()
model.add(Dense(12, input_dim=8, activation='relu'))
model.add(Dense(8, activation='relu'))
model.add(Dense(1, activation='sigmoid'))
# compile the keras model
model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
# fit the keras model on the dataset
model.fit(X, y, epochs=150, batch_size=10)
# evaluate the keras model
_, accuracy = model.evaluate(X, y)
print('Accuracy: %.2f' % (accuracy*100))
'''
