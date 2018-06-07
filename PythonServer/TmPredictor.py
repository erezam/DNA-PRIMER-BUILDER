#!/usr/bin/python
# -*- coding: utf-8 -*-

import json
import matplotlib.pyplot as plt
import numpy as np
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score


class TmPredictor(object):
    def __init__(self):
        # Create linear regression object
        global regr
        regr = linear_model.LinearRegression()
        # read data from json
        data = json.load(open("Data.json"))
        pairs_vectors = []
        tm = []
        for primer in data["Primers"]:
            pairs_vectors.append(get_pairs_vector(primer["seq"]))
            tm.append(float(primer["tm"]))

        pairs_vectors_sums = []
        for vector in pairs_vectors:
            vector = list(vector.values())
            vector.append(1)
            pairs_vectors_sums.append(vector)

        # Train the model using the training sets
        regr.fit(pairs_vectors_sums, tm)


def get_pairs_vector(seq):
    pairs = {'AA': 0, 'AC': 0, 'AG': 0, 'AT': 0, 'CC': 0, 'CG': 0, 'CT': 0, 'GG': 0, 'GT': 0, 'TT': 0}
    for index, l in enumerate(seq):
        if index + 2 <= len(seq):
            list_pair = sorted(seq[index:index + 2])
            str_pair = ''.join(list_pair)
            pairs[str_pair] = pairs.get(str_pair) + 1
    return pairs



# Make predictions using the testing set
#y_pred = regr.predict([[1, 3, 2, 2, 1, 3, 7, 0, 1, 0, 1]])
#print (y_pred[0])
# The coefficients
#print('Coefficients: \n', regr.coef_)
# The mean squared error
#print("Mean squared error: %.2f" % mean_squared_error(y_test, y_pred))
# Explained variance score: 1 is perfect prediction
#print('Variance score: %.2f' % r2_score(y_test, y_pred))



