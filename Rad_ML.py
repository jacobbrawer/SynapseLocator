
import numpy as np
import pandas as pd
import os
from math import sin,cos, pi
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import sklearn

from neuron_morphology.validation.result import InvalidMorphology
from neuron_morphology.validation.result import NodeValidationError

path = 'C:\Users\jacob.brawer\Downloads\Test_Files'

def load_in_swc(path, filename):
    output = []
    r = []
    with open(os.path.join(path,filename), "r") as f:
        r = []
        for line in f:
            try:
                if line.lstrip().startswith("#"):
                    continue
                toks = line.split()
                vals = [
                            int(toks[0]),            # node value
                            int(toks[1]),            # type
                            float(toks[2]),          # x coordinate
                            float(toks[3]),          # y coordinate
                            float(toks[4]),          # z coordinate
                            float(toks[5]),          # radius
                            int(toks[6].rstrip())   # parent node
                        ]
                r.append(vals)
                output.append(vals)
            except IndexError:
                message = "File is not recognized as a valid swc file. One of the columns is missing a value"
                raise InvalidMorphology([NodeValidationError(message, line[0], "Error")])
    return output


def re_org_swc(SWC):
    """Reorganizes an SWC input, which divides by cell, and instead creates seperate lists (columns), for each variable"""
    output = []
    nodenum = []
    proj_type = []
    x = []
    y = []
    z = []
    radius = []
    parent = []
    for col in SWC:
        nodenum.append(col[0])
        proj_type.append(col[1])
        x.append(col[2])
        y.append(col[3])
        z.append(col[4])
        radius.append(col[5])
        parent.append(col[6])
    output.extend([nodenum,proj_type,x,y,z,radius,parent])
    return output

def remove_un(line):
    """Remove unnecessary characters in the string file name"""
    line = line.replace("[", "")
    line = line.replace("]", "")
    line = line.replace("'", "")
    return line

def basal_branch_order(path, filename):
    """Establishes branch order for all basal dendrites within an SWC"""
    SWC = load_in_swc(path, filename)
    basal_info = []
    soma = SWC[0]
    SWC = SWC[1:]
    counter = 0
    branch_index, branch_order, current_order = 1, 1, 1
    branch_change_dict = {}
    branch_dict = {}
    main_branches = []
    for node in SWC:        #finds the node locations for the start of each main branch
        if node[1] == 3:
            branch_index = node[-1]
            node_index = node[0]
            if branch_index == 1:
                main_branches.append(node_index)
                
    branch_order = 1
    for node in SWC:
        if node[1] == 3:
               #searches all projections that branch from the main branch and establishes their node location
            branch_index = node[-1]
            node_index = node[0]
            if node_index in main_branches:
                basal_info.append([node_index, 1, branch_index])
                branch_change_dict[node_index] = branch_order   
            elif branch_index in main_branches:
                branch_order = 2 
                basal_info.append([node_index, branch_order, branch_index])
                branch_dict[node_index] = branch_order
            elif branch_index in branch_dict:
                branch_order = branch_dict[branch_index]
                basal_info.append([node_index, branch_order, branch_index])
                branch_dict[node_index] = branch_order

            elif node_index != branch_index - 1:
                branch_order = branch_dict[branch_index] + 1
                branch_change_dict[branch_index] = branch_order 
                basal_info.append([node_index, branch_order, branch_index])

    return basal_info 

def apical_branch_order(path, filename):
    """Establishes branch order for all apical dendrites within an SWC"""
    SWC = load_in_swc(path, filename)
    apical_info = []
    soma = SWC[0]
    SWC = SWC[1:]
    counter = 0
    branch_index, branch_order, current_order = 1, 1, 1
    branch_change_dict = {}
    branch_dict = {}
    main_branches = []
    for node in SWC:        #finds the node locations for the start of each main branch
        if node[1] == 4:
            branch_index = node[-1]
            node_index = node[0]
            if branch_index == 1:
                main_branches.append(node_index)  
    branch_order = 1
    for node in SWC:
        if node[1] == 4:
               #searches all projections that branch from the main branch and establishes their node location
            branch_index = node[-1]
            node_index = node[0]
            if node_index in main_branches:
                apical_info.append([node_index, 1, branch_index])
                branch_change_dict[node_index] = branch_order   
            elif branch_index in main_branches:
                branch_order = 2 
                apical_info.append([node_index, branch_order, branch_index])
                branch_dict[node_index] = branch_order
            elif branch_index in branch_dict:
                branch_order = branch_dict[branch_index]
                apical_info.append([node_index, branch_order, branch_index])
                branch_dict[node_index] = branch_order

            elif node_index != branch_index - 1:
                branch_order = branch_dict[branch_index] + 1
                branch_change_dict[branch_index] = branch_order 
                apical_info.append([node_index, branch_order, branch_index])

    return apical_info 

# def axon_branch_order(path, filename):
#     """Establishes branch order for all axons within an SWC"""
#     SWC = load_in_swc(path, filename)
#     axon_info = []
#     soma = SWC[0]
#     SWC = SWC[1:]
#     branch_index, branch_order, current_order = 1, 1, 1
#     branch_dict = {}
#     main_branches = []
#     for node in SWC:            #finds the node locations for the start of each main branch
#         if node[1] == 2:
#             branch_index = node[-1]
#             if branch_index == 1:
#                 node_index = node[0]
#                 main_branches.append(node_index)
#                 axon_info.append([node_index, 1, node[-2], branch_index])

#     for main in main_branches:      #searches all projections that branch from the main branch and establishes their node location
#         branch_order = 1
#         for node in SWC:
#             if node[1] == 2:
#                 branch_index = node[-1]
#                 node_index = node[0]
#                 if node_index == main:
#                     branch_order += 1
#                     axon_info.append([node_index, branch_order, node[-2]])
#                     branch_dict[branch_index] = branch_order
#                 elif node_index in branch_dict:
#                     branch_order = branch_dict[node_index] + 1
#                     branch_dict[branch_index] = branch_order
#                     axon_info.append([node_index, branch_order, node[-2], branch_index])
#                 else:
#                     axon_info.append([node_index, branch_order, node[-2], branch_index])
#     return axon_info 

# def basal_branch_order(path, filename):
#     """Establishes branch order for all basal dendrites within an SWC"""
#     SWC = load_in_swc(path, filename)
#     basal_info = []
#     soma = SWC[0]
#     SWC = SWC[1:]
#     branch_index, branch_order, current_order = 1, 1, 1
#     branch_dict = {}
#     main_branches = []
#     for node in SWC:        #finds the node locations for the start of each main branch
#         if node[1] == 3:
#             branch_index = node[-1]
#             if branch_index == 1:
#                 node_index = node[0]
#                 main_branches.append(node_index)
#                 basal_info.append([node_index, 1, node[-2]])

#     for main in main_branches:      #searches all projections that branch from the main branch and establishes their node location
#         branch_order = 1
#         for node in SWC:
#             if node[1] == 3:
#                 branch_index = node[-1]
#                 node_index = node[0]
#                 if node_index == main:
#                     branch_order += 1
#                     basal_info.append([node_index, branch_order, node[-2]])
#                     branch_dict[branch_index] = branch_order
#                 elif node_index in branch_dict:
#                     branch_order = branch_dict[node_index] + 1
#                     branch_dict[branch_index] = branch_order
#                     basal_info.append([node_index, branch_order, node[-2]])
#                 else:
#                     basal_info.append([node_index, branch_order, node[-2]])
#     return basal_info 


def path_length(path, filename, proj_type):
    SWC = load_in_swc(path, filename)
    soma = np.array([SWC[0][2], SWC[0][3], SWC[0][4]])
    SWC = re_org_swc(SWC[1:])
    current_node = SWC[0]
    node_type = SWC[1]
    parent_node = SWC[-1]
    x = SWC[2]
    y = SWC[3]
    z = SWC[4]
    path_dict = {}
    for i in range(len(node_type)):
        if node_type[i] == proj_type and parent_node[i] == 1:
            node = np.array([x[i],y[i],z[i]])
            path_dict[current_node[i]] = [node, np.linalg.norm(node - soma)]
        elif parent_node[i] in path_dict:
            node = np.array([x[i],y[i],z[i]])
            prev_node = path_dict[parent_node[i]][0]
            prev_dist = path_dict[parent_node[i]][1]
            path_dict[current_node[i]] = [node, prev_dist + np.linalg.norm(node - prev_node)]
    return path_dict
    
def coords_calc(location,x,y,z):
    coords = np.array([x[location], y[location], z[location]]) #returns the coordinates based on a given index within the SWC dataset
    return coords

def bdendrite_projections(path, filename):
    """Returns a list of all basal dendrite nodes and their coordinates based on cell number and which dendrite projection it belongs to"""
    SWC = load_in_swc(path, filename)
    SWC = re_org_swc(SWC)
    branch_list = basal_branch_order(path, filename)
    path_dict = path_length(path, filename, 3)
    current_node = SWC[0]
    node_type = SWC[1]
    rad = SWC[5]
    parent_node = SWC[-1]
    x = SWC[2]
    y = SWC[3]
    z = SWC[4]
    basal_coords = []
    for i in range(len(node_type)):
        if node_type[i] == 3:
            for basal in branch_list:
                if basal[0] == current_node[i] and basal[-1] == parent_node[i]:
                    branch_order = basal[1]
                    proj_len = path_dict[current_node[i]][1]
                    coords = coords_calc(i, x, y, z)
                    proj = [branch_order, rad[i], coords, proj_len]
                    basal_coords.append(proj)
    return basal_coords  # returns in the form (cell #, projection type, projection #, coordinates)

def adendrite_projections(path, filename):
    """Returns a list of all apical dendrite nodes and their coordinates based on cell number and which dendrite projection it belongs to"""
    SWC = load_in_swc(path, filename)
    SWC = re_org_swc(SWC)
    branch_list = apical_branch_order(path, filename)
    path_dict = path_length(path, filename, 4)
    current_node = SWC[0]
    node_type = SWC[1]
    rad = SWC[5]
    parent_node = SWC[-1]
    x = SWC[2]
    y = SWC[3]
    z = SWC[4]
    apical_coords = []
    for i in range(len(node_type)):
        if node_type[i] == 4:
            for apical in branch_list:
                if apical[0] == current_node[i] and apical[-1] == parent_node[i]:
                    branch_order = apical[1]
                    proj_len = path_dict[current_node[i]][1]
                    coords = coords_calc(i, x, y, z)
                    proj = [branch_order, rad[i], coords, proj_len]
                    apical_coords.append(proj)
    return apical_coords  # returns in the form (cell #, projection type, projection #, coordinates)


def bdendrite_compiler(path):
    filenames = os.listdir(path)
    output = []
    for filename in filenames:
        dendrites = bdendrite_projections(path,filename)
        output.extend(dendrites)
    return output

def adendrite_compiler(path):
    filenames = os.listdir(path)
    output = []
    for filename in filenames:
        dendrites = adendrite_projections(path,filename)
        output.extend(dendrites)
    return output


##### Rad Estimation - ML Pipeline 

from sklearn.model_selection import train_test_split, cross_validate, GridSearchCV
from sklearn.neural_network import MLPRegressor
from sklearn.linear_model import ElasticNet
from sklearn import linear_model
from sklearn.metrics import classification_report, confusion_matrix
from sklearn import svm
from matplotlib import pyplot as plt
from sklearn import metrics
from sklearn import preprocessing
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import RandomizedSearchCV, cross_val_score, GridSearchCV

def accuracy(predictions):
    hits = 0
    divisor = len(predictions)
    for i in range(len(predictions)):
        upperbound = y_test.iloc[i] + y_std
        lowerbound = y_test.iloc[i] - y_std
        if predictions[i] > lowerbound:
            
            if upperbound > predictions[i]:
                hits += 1
    
    return float(hits) / divisor

def evaluate(model, test_features, test_labels):
    predictions = model.predict(test_features)
    errors = abs(predictions - test_labels)
    mape = 100 * np.mean(errors / test_labels)
    accuracy = 100 - mape
    print('Average Error:', format(np.mean(errors)))
    print('Accuracy', format(accuracy))
    # return accuracy


bdend = bdendrite_compiler(path)
Standard_X = preprocessing.StandardScaler()
df = pd.DataFrame(bdend, columns = ['BO', 'Rad','Coords', 'PL'])
x = df.drop(['Rad', 'Coords'], axis = 1)
y = df['Rad']
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size = 0.1)
x_train = Standard_X.fit_transform(x_train)
x_test = Standard_X.fit_transform(x_test)
y_std = np.std(y)


# lm = linear_model.LinearRegression()
# model1 = lm.fit(x_train, y_train) 
# predictions1 = model1.predict(x_test)


# ridge = linear_model.Ridge(alpha=.5)
# model2 = ridge.fit(x_train, y_train) 
# predictions2 = model2.predict(x_test)


lasso = linear_model.Lasso(alpha=0.1)
model3 = lasso.fit(x_train, y_train) 
predictions3 = model3.predict(x_test)


# bay_ridge = linear_model.BayesianRidge()
# model4 = bay_ridge.fit(x_train, y_train) 
# predictions4 = model4.predict(x_test)


# e_net = ElasticNet(random_state=0)
# model5 = e_net.fit(x_train, y_train) 
# predictions5 = model5.predict(x_test)


# print('Basal Model Performance')
# print("LM Score:", accuracy(predictions1))
# evaluate(model1, x_test, y_test)
# print("Ridge Score:", accuracy(predictions2))
# evaluate(model2, x_test, y_test)
# print("Lasso Score:", accuracy(predictions3))
# evaluate(model3, x_test, y_test)
# print("Bayesian Ridge Score:", accuracy(predictions4))
# evaluate(model4, x_test, y_test)
# print("Elastic Net Score:", accuracy(predictions5))
# evaluate(model5, x_test, y_test)


# # Number of trees in random forest
# n_estimators = [int(x) for x in np.linspace(start = 200, stop = 2000, num = 10)]
# # Number of features to consider at every split
# max_features = ['auto', 'sqrt']
# # Maximum number of levels in tree
# max_depth = [int(x) for x in np.linspace(10, 110, num = 11)]
# max_depth.append(None)
# # Method of selecting samples for training each tree
# bootstrap = [True, False]

# random_grid = {'n_estimators': n_estimators,
#                'max_features': max_features,
#                'max_depth': max_depth,
#                'bootstrap': bootstrap}

# RF = RandomForestRegressor()
# RF_random = RandomizedSearchCV(estimator = RF, param_distributions = random_grid, n_iter = 100, cv = 3, n_jobs = -1)
# model7 = RF_random.fit(x_train, y_train)
# best_params = model7.best_params_
# rfr = RandomForestRegressor(max_depth=best_params["max_depth"], n_estimators=best_params["n_estimators"], max_features = best_params['max_features'], bootstrap = best_params['bootsrap'], random_state=False, verbose=False)
# predictions7 = rfr.predict(x_test)
# RF_accuracy1 = evaluate(model7, x_test, y_test)
# print("RF Score:", accuracy(predictions7))


adend = adendrite_compiler(path)
Standard_X = preprocessing.StandardScaler()
df = pd.DataFrame(adend, columns = ['BO', 'Rad', 'Coords', 'PL'])
x = df.drop(['Rad', 'Coords'], axis = 1)
y = df['Rad']
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size = 0.1)
x_train = Standard_X.fit_transform(x_train)
x_test = Standard_X.fit_transform(x_test)
y_std = np.std(y)

# lm = linear_model.LinearRegression()
# model1 = lm.fit(x_train, y_train) 
# predictions1 = model1.predict(x_test)

# ridge = linear_model.Ridge(alpha=.5)
# model2 = ridge.fit(x_train, y_train) 
# predictions2 = model2.predict(x_test)

# lasso = linear_model.Lasso(alpha=0.1)
# model3 = lasso.fit(x_train, y_train) 
# predictions3 = model3.predict(x_test)

# bay_ridge = linear_model.BayesianRidge()
# model4 = bay_ridge.fit(x_train, y_train) 
# predictions4 = model4.predict(x_test)


e_net = ElasticNet(random_state=0)
model5 = e_net.fit(x_train, y_train) 
predictions5 = model5.predict(x_test)



# print('Apical Model Performance')
# print("LM Score:", accuracy(predictions1))
# evaluate(model1, x_test, y_test)
# print("Ridge Score:", accuracy(predictions2))
# evaluate(model2, x_test, y_test)
# print("Lasso Score:", accuracy(predictions3))
# evaluate(model3, x_test, y_test)
# print("Bayesian Ridge Score:", accuracy(predictions4))
# evaluate(model4, x_test, y_test)
# print("Elastic Net Score:", accuracy(predictions5))
# evaluate(model5, x_test, y_test)
