import matplotlib
matplotlib.use('Agg')
import numpy
from numpy import arange
from matplotlib import pyplot
from pandas import read_csv
from pandas import set_option
from pandas.plotting import scatter_matrix
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Lasso
from sklearn.linear_model import Ridge
from sklearn.linear_model import ElasticNet
from sklearn.tree import DecisionTreeRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.svm import SVR
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.metrics import mean_squared_error
from xgboost import XGBRegressor
from lightgbm import LGBMRegressor
from catboost import CatBoostRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import cross_validate
import re
import sys
import os,glob
import numpy as np
import os
import sys
import argparse
import random
import re
from sklearn.metrics import r2_score
import pandas as pd
import numpy as np
from sklearn import datasets, linear_model
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from scipy import stats
from scipy.stats import spearmanr
from scipy import optimize as op
from sklearn.inspection import permutation_importance

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('inn1', help='input file')
    parser.add_argument('outt1', help='output file')
    opts = parser.parse_args()
    outp = open(opts.outt1, 'w')

    filename = opts.inn1
    names2 = ["RNAPII","DNase","H2A.Z","H3K4me3","H3K4me2","H3K27ac","H3K9ac","H3K36me3","H3K27me3","H3K9me1","H3K9me3","H4K20me1","H3K4me1","H3K79me2","Rloop"]
    dataset2 = read_csv(filename,sep='\t', header=0)
    dataset=dataset2.iloc[:,0:15]
#Correlation matrix plot for the bin
    fig = pyplot.figure(figsize=(12.0, 9.0))
    ax = fig.add_subplot(111)
    cax = ax.matshow(dataset.corr(), vmin=-1, vmax=1, interpolation='none')
    fig.colorbar(cax)
    ticks = numpy.arange(0, 15, 1)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_yticklabels(names2,fontsize=20)
    pyplot.savefig(opts.inn1+"multibin-1.png")

    array = dataset.values
    X = array[:,0:14]
    Y = array[:,14]
    validation_size = 0.20
    seed = 7
    X_train, X_validation, Y_train, Y_validation = train_test_split(X, Y,test_size=validation_size, random_state=seed)

    scoring = "r2", "neg_mean_squared_error","neg_root_mean_squared_error"

#standardize the data and build the model
    pipelines = []
    pipelines.append(('ScaledLR', Pipeline([('Scaler', StandardScaler()),('LR',
    LinearRegression())])))
    pipelines.append(('ScaledLASSO', Pipeline([('Scaler', StandardScaler()),('LASSO',
    Lasso())])))
    pipelines.append(('ScaledEN', Pipeline([('Scaler', StandardScaler()),('EN',
    ElasticNet())])))
    pipelines.append(('ScaledKNN', Pipeline([('Scaler', StandardScaler()),('KNN',
    KNeighborsRegressor())])))
    pipelines.append(('ScaledCART', Pipeline([('Scaler', StandardScaler()),('CART',
    DecisionTreeRegressor())])))
    pipelines.append(('ScaledXGB', Pipeline([('Scaler', StandardScaler()), ('XGB', XGBRegressor())])))
    pipelines.append(('ScaledCATBoost', Pipeline([('Scaler', StandardScaler()), ('CATBoost', CatBoostRegressor())])))
    pipelines.append(('ScaledMLP', Pipeline([('Scaler', StandardScaler()), ('MLP', MLPRegressor())])))
    pipelines.append(('ScaledAB', Pipeline([('Scaler', StandardScaler()),('AB',AdaBoostRegressor())])))
    pipelines.append(('ScaledGBM', Pipeline([('Scaler', StandardScaler()),('GBM',GradientBoostingRegressor())])))
    pipelines.append(('ScaledRF', Pipeline([('Scaler', StandardScaler()),('RF',RandomForestRegressor())])))
    pipelines.append(('ScaledET', Pipeline([('Scaler', StandardScaler()),('ET',ExtraTreesRegressor())])))
    results = []
    results2 = []
    results3 = []
    names = []
    for name, model in pipelines:
        cv_results = cross_validate(model, X_train, Y_train, cv=5, scoring=("r2", "neg_mean_squared_error","neg_mean_absolute_error"))
        results.append(cv_results['test_r2'])
        results2.append(cv_results['test_neg_mean_squared_error'])
        results3.append(cv_results['test_neg_mean_absolute_error'])
        names.append(name)
        msg = "%s: %s: %f (%f)" % (name, "r2", cv_results['test_r2'].mean(), cv_results['test_r2'].std())
        msg2 = "%s:  %s: %f (%f)" % (name, "mse", cv_results['test_neg_mean_squared_error'].mean(), cv_results['test_neg_mean_squared_error'].std())
        msg3 = "%s:  %s: %f (%f)" % (name, "mrse", cv_results['test_neg_mean_absolute_error'].mean(), cv_results['test_neg_mean_absolute_error'].std())
#print r-squared, negative mean squared error, and negative root mean squared error
        outp.write(opts.inn1+': '+msg+'\n')
        outp.write(opts.inn1+': '+msg2 + '\n')
        outp.write(opts.inn1+': '+msg3 + '\n')
#scaled algorithm comparison based on mean squared error
    fig = pyplot.figure(figsize=(12.0, 9.0))
    fig.suptitle('Scaled Algorithm Comparison based on mean_squared_error',fontsize=20)
    ax = fig.add_subplot(111)
    pyplot.boxplot(results2)
    ax.set_xticklabels(names, rotation=40,fontsize=12)
    pyplot.savefig(opts.inn1+"multibin-2.png")
#scaled algorithm comparison based on root mean squared error
    fig = pyplot.figure(figsize=(12.0, 9.0))
    fig.suptitle('Scaled Algorithm Comparison based on root_mean_squared_error',fontsize=20)
    ax = fig.add_subplot(111)
    pyplot.boxplot(results3)
    ax.set_xticklabels(names, rotation=40,fontsize=12)
    pyplot.savefig(opts.inn1+"multibin-3.png")


#GridSearchCV implement for extratreesregressor models
    scaler = StandardScaler().fit(X_train)
    rescaledX = scaler.transform(X_train)
    param_grid = dict(n_estimators=numpy.array([50,100,200,400,800, 1000,1500,2000,3000,4000]))
    model = ExtraTreesRegressor(random_state=seed)
    grid = GridSearchCV(estimator=model, param_grid=param_grid, scoring='neg_mean_squared_error', cv=5, iid=True)
    grid_result = grid.fit(rescaledX, Y_train)
    print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_))
    means = grid_result.cv_results_['mean_test_score']
    stds = grid_result.cv_results_['std_test_score']
    params = grid_result.cv_results_['params']
    for mean, stdev, param in zip(means, stds, params):
        print("%f (%f) with: %r" % (mean, stdev, param))
        if mean==grid_result.best_score_:
            outp.write(opts.inn1+': '+'mean_test_score: '+str(mean) + '\n')
            outp.write(opts.inn1+': '+'std_test_score: ' + str(stdev) + '\n')
            outp.write(opts.inn1+': '+'n_estimators: ' + str(grid_result.best_params_['n_estimators']) + '\n')
#make predictions
    scaler = StandardScaler().fit(X_train)
    rescaledX = scaler.transform(X_train)
    model = ExtraTreesRegressor(random_state=seed, n_estimators=grid_result.best_params_["n_estimators"])
    model.fit(rescaledX, Y_train)
    rescaledValidationX = scaler.transform(X_validation)
    predictions = model.predict(rescaledValidationX)
    importance=model.feature_importances_


#output the relative importance of epigenetics marks and R squared
    count=0
    names3=[]
    importance3=[]
    for i in list(np.argsort(importance)):
        count=count+1
        outp.write(opts.inn1+'\t'+names2[i]+'\t'+str(importance[i])+'\t'+str(count)+'\n')
        names3.append(names2[i])
        importance3.append(importance[i])
    fig = pyplot.figure(figsize=(12.0, 9.0))
    fig.suptitle('The relative importance of epigenetics marks',fontsize=30)
    ax = fig.add_subplot(111)
    pyplot.bar(names3, importance3)
    ax.set_xticklabels(names3, rotation=40,fontsize=14)
    pyplot.savefig(opts.inn1+"multibin-4.png")

    est = sm.OLS(Y_validation, predictions)
    est2 = est.fit()
    outp.write(opts.inn1+': '+'R squared: ' + str(est2.rsquared) + '\n')
    fig=pyplot.figure(figsize=(12.0, 9.0))
    ax=fig.add_subplot(111)
    ax.scatter(Y_validation, predictions)
    ax.set_xlabel('Observed R-loop')
    ax.set_ylabel('Predicted R-loop')
    ax.annotate("R-squared: %f" % (est2.rsquared),(0,predictions.max()))
    pyplot.savefig(opts.inn1+"multibin-5.png")


