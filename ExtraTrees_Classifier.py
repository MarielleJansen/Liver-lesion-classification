# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 15:53:16 2017

@author: Marielle
"""

import xlrd
import numpy as np
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.model_selection import LeaveOneOut
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.svm import SVC
from sklearn.feature_selection import SelectKBest, f_classif,  SelectPercentile

book = xlrd.open_workbook('C:\Users\user\Documents\Characterization\Features_5lesions_splitfeatures.xlsx')


# get the right worksheet
sheet = book.sheet_by_name('T2W+DCE')

Features = []

num_rows = sheet.nrows - 1
num_cols = sheet.ncols
curr_row = 0
while curr_row < num_rows:
    curr_row += 1
    data = []
    data = [sheet.cell(curr_row, row_index).value for row_index in xrange(2, num_cols)]
    if curr_row == 1:
        Features = [data]
    else:
        Features = Features + [data]


lesiontype = [sheet.cell(row_index, 1).value for row_index in xrange(1, num_rows+1)]
sub = [sheet.cell(row_index, 0).value for row_index in xrange(1, num_rows+1)]
sub = np.asarray(sub)

Features = np.asarray(Features)

lesiontype = np.asarray(lesiontype)


Labels = list(lesiontype)


clf = ExtraTreesClassifier(n_estimators=700, max_depth=None,
                           n_jobs=-1, max_features="auto",
                           min_samples_split=4, random_state=None,
                           class_weight='balanced_subsample',
                           min_samples_leaf=5)



X = np.array(xrange(213))
loo = LeaveOneOut()
looPredict = np.zeros((213, 5), dtype=np.float)
result1 = np.zeros((213, 1))
predictions1 = np.zeros((213, 1))



for train, test in loo.split(X):
    train = np.delete(train, np.where(sub == sub[test]))
    xFeatures = Features[train, :]
    y = lesiontype[train]
    
    selection = SelectKBest(score_func=f_classif, k=19)
    BestFeatures = selection.fit(xFeatures, y).transform(xFeatures)
    BestFeaturesIndices = selection.get_support(indices=True)
    FeaturesTest = selection.transform(Features[test])

    
    x = BestFeatures
    clf.fit(x, y)
    testF = FeaturesTest
    testL = lesiontype[test]
    looPredict[test, :] = np.array(clf.predict_proba(testF))
    result1[test] = clf.score(testF, testL)
    predictions1[test] = clf.predict(testF)



print("score of classification: %0.4f " % (sum(result1)/float(213)))



