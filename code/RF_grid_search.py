################################################################
#                                                              # 
#   Set up Random Forest -- grid search for protein folding    #
#   system                                                     #
################################################################

import os
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report,confusion_matrix,accuracy_score
from random import shuffle
from timeit import default_timer as timer

direct2 = 'path of final csv files'
filenames = [i for i in os.listdir(direct2) if '.csv' in i]
shuffle(filenames)
print (filenames)

trainfiles, testfiles = train_test_split(filenames, train_size = 0.8, random_state = 42) 
print (testfiles)
x_train = pd.DataFrame()
y_train = pd.DataFrame()
x_test = pd.DataFrame()
y_test = pd.DataFrame()
train_data = pd.DataFrame()
test_data = pd.DataFrame()
for tfile in trainfiles:
    train_data = train_data.append(pd.read_csv(direct2+tfile)[1:])
x_train = train_data.drop(['16029'], axis=1).iloc[:,1:]
y_train = train_data['16029']
y_train=y_train.astype('float')
print(x_train)
print(y_train)

for tfile in testfiles:
    test_data = test_data.append(pd.read_csv(direct2+tfile)[1:])
    temp = pd.read_csv(direct2+tfile)[1:]
    for row in temp.index:
        if float(temp.loc[row, '16029']) != 0 and float(temp.loc[row, '16029']) != 1:
            print (tfile, row, temp.loc[row, '16029'])
x_test = test_data.drop(['16029'], axis=1).iloc[:,1:]
y_test = test_data['16029']
print(x_test)
print(y_test)
y_test=y_test.astype('float')



clf =  RandomForestClassifier(criterion='gini', class_weight = {0:0.5, 1:0.5})
parameters = {
       'n_estimators':(1000,5000),
       'max_depth':(5,100),
       'min_samples_split':(2,5),
       'min_samples_leaf':(1,2) }
grid_search = GridSearchCV(clf, parameters, cv=5,verbose=1,scoring='accuracy')
grid_search.fit(x_train, y_train)
print ('Best Training score: %0.3f' % grid_search.best_score_)
print ('Best parameters set:')
best_parameters = grid_search.best_estimator_.get_params()
for param_name in sorted(parameters.keys()):
    print ('\t%s: %r' % (param_name, best_parameters[param_name]))
predictions = grid_search.predict(x_test)
print ("Testing accuracy:",round(accuracy_score(y_test,predictions),4))
print ("\nComplete report of Testingdata\n",classification_report(y_test, predictions))
print ("\n\nRandom Forest Grid Search- Test Confusion Matrix\n\n",pd.crosstab( y_test, predictions,rownames = ["Actuall"],colnames = ["Predicted"]))
