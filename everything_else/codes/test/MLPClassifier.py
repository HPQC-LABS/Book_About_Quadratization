# -*- coding: utf-8 -*-
"""
Created on Sat Aug 24 15:02:34 2019

@author: A_Soteriou
"""

import scipy.io as sio
import numpy as np
from random import shuffle
from sklearn.neural_network import MLPClassifier

train_data_percentage = 0.9

# data preprocessing
fgood = open("good.txt","r")
fbad  = open( "bad.txt","r")

lines_good = fgood.readlines()
lines_bad  =  fbad.readlines()

fgood.close()
fbad.close()

data_good = [[int(x) for x in line.split()] for line in lines_good]
data_bad  = [[int(x) for x in line.split()] for line in lines_bad ]
data = np.concatenate((data_good,data_bad))
data_good_len = len(data_good)
data_bad_len  = len(data_bad)

shuffle(data_good)
shuffle(data_bad)

num_train_data = int(np.floor(train_data_percentage*len(data_good)))
num_test_data  = len(data_good) - num_train_data

train_data_good = data_good[:num_train_data]
train_data_bad  =  data_bad[:num_train_data]
test_data_good  = data_good[num_train_data:num_train_data + num_test_data]
test_data_bad   =  data_bad[num_train_data:num_train_data + num_test_data]

train_data = np.concatenate((train_data_good, train_data_bad))
test_data  = np.concatenate((test_data_good, test_data_bad))
labels = [1]*num_train_data +[0]*num_train_data;

# set model parameters
hidden_layers = (30, 20, 20)
max_iter = 200
learning_rate = 0.0004
model = MLPClassifier(hidden_layer_sizes = hidden_layers,max_iter = max_iter,learning_rate_init = learning_rate,verbose = 'true')

# train model
model.fit(train_data,labels)

# save the weights for MATLAB access
i = 1
parameters = {}
for W in model.coefs_:
    parameters['W'+str(i)] = W
    i = i + 1
i = 1
for b in model.intercepts_:
    parameters['b'+str(i)] = b
    i = i + 1
sio.savemat('parameters.mat', parameters)

# calculate accuracy
true_labels = [1]*num_test_data + [0]*num_test_data
score = model.score(test_data,true_labels)
print('accuracy: %.2f%%' % (score*100) )

'''
prob_threshold = 0.5
pred = model.predict_proba(temp)
pred = (pred[:,1] >= prob_threshold)

true_labels = [1]*data_good_len + [0]*data_bad_len
print('overall accuracy: %.2f%%' % (100*(1 - np.mean(abs(true_labels - pred )))) )

print('true positives accuracy: %.2f%%' % (100*(1 - np.mean(abs(true_labels[:data_good_len] - pred[:data_good_len] )))) )

print('true negatives accuracy: %.2f%%' % (100*(1 - np.mean(abs(true_labels[data_good_len:] - pred[data_good_len:] )))) )
'''


