# SWC File -> Cell Cluster Distance Calculator 
# written by Jacob Brawer
# Allen Instite - 6/26/2019

import numpy as np
import pandas as pd
import os
from math import sin,cos, pi
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from neuron_morphology.validation.result import InvalidMorphology
from  Tkinter import *
import Tkinter, Tkconstants, tkFileDialog
from scipy import stats

from neuron_morphology.validation.result import NodeValidationError

def load_in_swc(path, filename):
    """loads in SWC file in the form list of lists"""
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
    SWC = output
    for i in range(len(SWC)):
        parent_node = SWC[i][-1]
        node_index = SWC[i][0]
        node_type = SWC[i][1]
        if parent_node == 1:
            if node_type == 1:
                SWC[i][1] = SWC[i+1][1]
    return SWC

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

def path_len(path, filename, proj_type):
    """calculates path length for a given projection"""
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

def axon_branch_order(path, filename):
    """Establishes branch order for all axons within an SWC"""
    SWC = load_in_swc(path, filename)
    axon_info = []
    soma = SWC[0]
    SWC = SWC[1:]
    counter = 0
    branch_index, branch_order, current_order = 1, 1, 1
    branch_change_dict = {}
    branch_dict = {}
    main_branches = []
    for node in SWC:        #finds the node locations for the start of each main branch
        if node[1] == 2:
            branch_index = node[-1]
            node_index = node[0]
            if branch_index == 1:
                main_branches.append(node_index)
                           
    branch_order = 1
    for node in SWC:
        if node[1] == 2:
               #searches all projections that branch from the main branch and establishes their node location
            branch_index = node[-1]
            node_index = node[0]
            if node_index in main_branches:
                axon_info.append([node_index, 1, branch_index])
                branch_change_dict[node_index] = branch_order 
            elif branch_index in main_branches:
                branch_order = 2 
                axon_info.append([node_index, branch_order, branch_index])
                branch_dict[node_index] = branch_order
            elif branch_index in branch_dict:
                branch_order = branch_dict[branch_index]
                axon_info.append([node_index, branch_order, branch_index])
                branch_dict[node_index] = branch_order

            elif node_index != branch_index - 1:
                branch_order = branch_dict[branch_index] + 1
                branch_change_dict[branch_index] = branch_order 
                axon_info.append([node_index, branch_order, branch_index])

    return axon_info 

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
    
def coords_calc(location,x,y,z):
    """returns the coordinates of a point as an array"""
    coords = np.array([x[location], y[location], z[location]]) #returns the coordinates based on a given index within the SWC dataset
    return coords

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
import Rad_ML
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import RandomizedSearchCV, cross_val_score, GridSearchCV


"""Pipeline for estimating radius using branch order + path length"""

bdend = Rad_ML.bdendrite_compiler(path)
Standard_X = preprocessing.StandardScaler()
df = pd.DataFrame(bdend, columns = ['BO', 'Rad', 'Coords', 'PL'])
x = df.drop(['Rad', 'Coords'], axis = 1)
y = df['Rad']
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size = 0.1)
x_train = Standard_X.fit_transform(x_train)
x_test = Standard_X.fit_transform(x_test)
y_std = np.std(y)

def accuracy(predictions):
    """calculates one measure of accuracy by testing whether a prediction is within one standard deviation of the mean radius"""
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
    """calculates another measure of accuracy by subtracting the MAPE from 100"""
    predictions = model.predict(test_features)
    errors = abs(predictions - test_labels)
    mape = 100 * np.mean(errors / test_labels)
    accuracy = 100 - mape
    print('Average Error:', format(np.mean(errors)))
    print('Accuracy', format(accuracy))

lasso = linear_model.Lasso(alpha=0.1)
b_model = lasso.fit(x_train, y_train) 
# predictions3 = model3.predict(x_test)

adend = Rad_ML.adendrite_compiler(path)
Standard_X = preprocessing.StandardScaler()
df = pd.DataFrame(adend, columns = ['BO', 'Rad', 'Coords', 'PL'])
x = df.drop(['Rad', 'Coords'], axis = 1)
y = df['Rad']
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size = 0.1)
x_train = Standard_X.fit_transform(x_train)
x_test = Standard_X.fit_transform(x_test)
y_std = np.std(y)

e_net = ElasticNet(random_state=0)
a_model = e_net.fit(x_train, y_train) 
predictions5 = model5.predict(x_test)

####### 

def axon_projections(parent_node, current_node, node_type, x, y, z, path, filename):
    """Returns a list of all axonal nodes and their coordinates based on cell number and which axonal projection it belongs to"""
    branch_list = axon_branch_order(path, filename)
    path_dict = path_len(path, filename, 2)
    axon_coords = []
    for i in range(len(node_type)):
        if node_type[i] == 2:
            for axon in branch_list:
                if axon[0] == current_node[i] and axon[-1] == parent_node[i]:
                    branch_order = axon[1]
                    proj_length = path_dict[current_node[i]][1]
                    rad = 0.1144
                    coords = coords_calc(i, x, y, z)
                    proj = [2, branch_order, rad, coords, proj_length]
                    axon_coords.append(proj)
    return axon_coords  # returns in the form (cell #, projection type, projection #, coordinates)

def bdendrite_projections(parent_node, current_node, node_type, x, y, z, path, filename):
    """Returns a list of all basal dendrite nodes and their coordinates based on cell number and which dendrite projection it belongs to"""
    branch_list = basal_branch_order(path, filename)
    path_dict = path_len(path, filename, 3)
    basal_coords = []
    for i in range(len(node_type)):
        if node_type[i] == 3:
            for basal in branch_list:
                if basal[0] == current_node[i] and basal[-1] == parent_node[i]:
                    branch_order = basal[1]
                    proj_length = path_dict[current_node[i]][1]
                    rad_info = pd.DataFrame([[branch_order, proj_length]], columns = ['BO', 'PL'])
                    rad = b_model.predict(rad_info)
                    coords = coords_calc(i, x, y, z)
                    proj = [3, branch_order, rad, coords, proj_length]
                    basal_coords.append(proj)
    return basal_coords  # returns in the form (cell #, projection type, projection #, coordinates)

def adendrite_projections(parent_node, current_node, node_type, x, y, z, path, filename):
    """Returns a list of all apical dendrite nodes and their coordinates based on cell number and which dendrite projection it belongs to"""
    branch_list = apical_branch_order(path, filename)
    path_dict = path_len(path, filename, 4)
    apical_coords = []
    for i in range(len(node_type)):
        if node_type[i] == 4:
            for apical in branch_list:
                if apical[0] == parent_node[i]:
                    branch_order = apical[1]
                    proj_length = path_dict[current_node[i]][1]
                    rad_info = pd.DataFrame([[branch_order, proj_length]], columns = ['BO', 'PL'])
                    rad = a_model.predict(rad_info)
                    coords = coords_calc(i,x,y,z)
                    proj = [4, branch_order, rad, coords, proj_length]
                    apical_coords.append(proj)
    return apical_coords  # returns in the form (cell #, projection type, projection #, coordinates)

def seg_dist(a0,a1,b0,b1,clampAll=False,clampA0=False,clampA1=False,clampB0=False,clampB1=False):
    ''' Given two lines defined by numpy.array pairs (a0,a1,b0,b1)
        Return the closest points on each segment and their distance
    '''
    # If clampAll=True, set all clamps to True
    if clampAll:
        clampA0=True
        clampA1=True
        clampB0=True
        clampB1=True
    
    # Calculate denomitator
    A = a1 - a0
    B = b1 - b0
    magA = np.linalg.norm(A)
    magB = np.linalg.norm(B)

    _A = A / magA
    _B = B / magB

    cross = np.cross(_A, _B);
    denom = np.linalg.norm(cross)**2

    # If lines are parallel (denom=0) test if lines overlap.
    # If they don't overlap then there is a closest point solution.
    # If they do overlap, there are infinite closest positions, but there is a closest distance
    if not denom:
        d0 = np.dot(_A,(b0-a0))

        # Overlap only possible with clamping
        if clampA0 or clampA1 or clampB0 or clampB1:
            d1 = np.dot(_A,(b1-a0))

            # Is segment B before A?
            if d0 <= 0 >= d1:
                if clampA0 and clampB1:
                    if np.absolute(d0) < np.absolute(d1):
                        return a0,b0,np.linalg.norm(a0-b0)
                    return a0,b1,np.linalg.norm(a0-b1)

            # Is segment B after A?
            elif d0 >= magA <= d1:
                if clampA1 and clampB0:
                    if np.absolute(d0) < np.absolute(d1):
                        return a1,b0,np.linalg.norm(a1-b0)
                    return a1,b1,np.linalg.norm(a1-b1)


        # Segments overlap, return distance between parallel segments
        return None,None,np.linalg.norm(((d0*_A)+a0)-b0)

    # Lines criss-cross: Calculate the projected closest points
    t = (b0 - a0);
    detA = np.linalg.det([t, _B, cross])
    detB = np.linalg.det([t, _A, cross])

    t0 = detA/denom;
    t1 = detB/denom;

    pA = a0 + (_A * t0) # Projected closest point on segment A
    pB = b0 + (_B * t1) # Projected closest point on segment B

    # Clamp projections
    if clampA0 or clampA1 or clampB0 or clampB1:
        if clampA0 and t0 < 0:
            pA = a0
        elif clampA1 and t0 > magA:
            pA = a1

        if clampB0 and t1 < 0:
            pB = b0
        elif clampB1 and t1 > magB:
            pB = b1

        # Clamp projection A
        if (clampA0 and t0 < 0) or (clampA1 and t0 > magA):
            dot = np.dot(_B,(pA-b0))
            if clampB0 and dot < 0:
                dot = 0
            elif clampB1 and dot > magB:
                dot = magB
            pB = b0 + (_B * dot)

        # Clamp projection B
        if (clampB0 and t1 < 0) or (clampB1 and t1 > magB):
            dot = np.dot(_A,(pB-a0))
            if clampA0 and dot < 0:
                dot = 0
            elif clampA1 and dot > magA:
                dot = magA
            pA = a0 + (_A * dot)

    return np.linalg.norm(pA-pB)

def min_dist(axon, dendrite):
    """calculates the minimum distance between the x,y,z coordinates of two nodes"""
    p1 = axon[4]
    p2 = dendrite[4]    
    dist = (np.linalg.norm(p1 - p2) - (dendrite[3] + axon[3])) / 0.1144 #pixel to um conversion (.1144 pixels / um)
    return dist

def min_proj_dist(axon, dendrite):
    """Iterates through every node for an axonal + dendritic pair and determines the node combination that results in the minimum distance"""
    upper_threshold = 50
    true_threshold = 20  #sets a threshold distance in microns
    a_to_d = []
    for i in range(len(axon)):
        for k in range(len(dendrite)):
            temp_min = min_dist(axon[i], dendrite[k])
            if temp_min < upper_threshold:
                if i < len(axon) - 1 and k < len(dendrite):
                    A0 = axon[i][4] 
                    A1 = axon[i+1][4]
                    B0 = dendrite[k][4]
                    B1 = dendrite[k+1][4]
                    new_min = seg_dist(A0,A1,B0,B1)
                    if new_min <= true_threshold:
                        a_to_d.append([axon[i], axon[i+1][4], dendrite[k], dendrite[k+1][4], new_min])
    return a_to_d

def cluster_dist(axons,bdend,adend):
    """Determines the apposition with the minimum distance between unique cell pairing and returns a list with the coordinates of each cell and the given distance"""
    #divides the list of each respective projection by cell
    final_output = []
    for i in range(len(axons)):
        for k in range(len(bdend)):
            if i != k:                  #ensures that projections from the same cell will not be compared
                
                a_to_d_syn = []
                a_to_bd = min_proj_dist(axons[i], bdend[k])
                if len(a_to_bd) != 0:
                    a_to_d_syn.extend([a_to_bd])          
                a_to_ad = min_proj_dist(axons[i], adend[k])
                if len(a_to_ad) != 0:
                    a_to_d_syn.extend([a_to_ad])
        
                if len(a_to_d_syn) != 0:  #ensures that the empty list a_to_ad_min defined at the beginning of the function is not added to the final output
                    final_output.extend(a_to_d_syn)
                    a_to_d_syn = []
    return final_output

def main():
    """Run this function to receive an output of the minimum distance axon-dendrite pairs between cells in a cluster"""
    axons, bdend, adend = [], [], []
    cell_counter = 0
    root = Tk()
    root.directory = tkFileDialog.askdirectory()
    path=root.directory
    #path = 'C:\Users\jacob.brawer\Downloads\Test_3cell'
    single_cells = os.listdir(path)
    # single_cells = single_cells[:-1]
    # rad_dict = {}
    output_dict = {}
    for cell in single_cells:
        # path = os.path.join(path, cell)
        cell_counter += 1
        SWC = re_org_swc(load_in_swc(path,cell))
        node_type = SWC[1]
        current_node = SWC[0]
        parent_node = SWC[-1]
        x = SWC[2]
        y = SWC[3]
        z = SWC[4]

        temp_axon = axon_projections(parent_node, current_node, node_type, x, y, z, path, cell)
        axon = []
        for ax in temp_axon:            #keeps track of cell number within each proction
            ax.insert(0,cell_counter)
            axon.append(ax)
        axons.append(axon)
        temp_bd = bdendrite_projections(parent_node, current_node, node_type, x, y, z, path, cell)
        bd = []
        for b in temp_bd:
            b.insert(0,cell_counter)
            bd.append(b)    
        bdend.append(bd)

        temp_ad = adendrite_projections(parent_node, current_node, node_type, x, y, z, path, cell)
        ad = []
        for a in temp_ad:
            a.insert(0,cell_counter)
            ad.append(a)
        adend.append(ad)
    return cluster_dist(axons,bdend,adend)

output = main()

def create_marker_file(markers):
	# get file name
	original_fn = str(path);
	new_fn = 'multipatch.marker'
	new_path = str(os.path.join(path, new_fn))
	new_path = remove_un(new_path)
	
	#create new file w/ header and add marker elements
	with open(new_path, 'w+') as f:
		# add header
		header = "##x,y,z,radius,shape,name,comment, color_r,color_g,color_b\n"
		f.writelines(header)
				
		# add markers 
		current_marker = 1
		while current_marker < len(markers):
			
			f.writelines(remove_un(str(markers[current_marker])))
			f.write('\n')
			current_marker = current_marker + 1
		
		f.close()
		print("A marker file has been generated ", new_fn)


def marker(connections):
"""Call this function to create a marker file for your output"""
#Put your output in 'connections'
    markers = []
    for c_to_c in connections:
        for synapse in c_to_c:
            syn = synapse[:-1]
            if len(syn) == 4:
                if syn[0][1] == 2:
                    markers.append([syn[0][4][0], syn[0][4][1], syn[0][4][2], 0, 0, 10, 0, 0, 0, 255])
                    markers.append([syn[1][0], syn[1][1], syn[1][2], 0, 0, 10, 0, 0, 0, 255])
                if syn[2][1] == 3 or syn[2][1] == 4:
                    markers.append([syn[2][4][0], syn[2][4][1], syn[2][4][2], 0, 0, 10, 0, 255, 0, 0])
                    markers.append([syn[3][0], syn[3][1], syn[3][2], 0, 0, 10, 0, 255, 0, 0])

            elif len(syn) == 3:
                if syn[1] == 2:
                    markers.append([syn[0][3][0], syn[0][3][1], syn[0][3][2], 0, 0, 10, 0, 255, 255, 255])
                if syn[1] == 3 or syn[1] == 4:
                    markers.append([syn[1][0], syn[1][1], syn[1][2], 0, 0, 10, 0, 255, 255, 255])
    create_marker_file(markers)
    return markers



def distribution(connections):
    """Returns descriptives about distances within the whole distribution (true and false positives)"""
    curve = []
    # for c_to_c in connections:
    for connection in connections:
        if len(connection) == 5:
            if int(connection[4]) > 0:
                curve.append(float(connection[4]))
        elif len(connection) == 3:
            if int(connection[2]) > 0:
                curve.append(float(connection[2]))
    print(stats.describe(curve))
    return curve
        
def true_distribution(connections, true1, true2):
    """Takes in the entire distribution (output), returns information about path length and synaptic distance for only true positive connections"""
    #note that cell numbers need to be changed (true1 and true2) in order to accomodate which cells are truly connected
    true_curve = []
    path_curve = []
    trues = []
    for c_to_c in connections:
        for connection in c_to_c:
            if len(connection) == 5:
                cell1 = connection[0]
                cell2 = connection[2]
                dist = float(connection[4])
                path1 = cell1[5]
                path2 = cell2[5]
            elif len(connection) == 3:
                cell1 = connection[0]
                cell2 = connection[1]
                dist = float(connection[2])
                path1 = cell1[5]
                path2 = cell2[5]
            if cell1[0] == true1 and cell2[0] == true2:
                if dist > 0:
                    true_curve.append(dist)
                    trues.append(connection)
                path_curve.append(path1)
                path_curve.append(path2)
    print(stats.describe(true_curve))
    print(stats.describe(path_curve))
    # print(np.mean(true_curve))
    # print(np.std(true_curve))
    # print(np.max(true_curve))
    # print(np.min(true_curve))
    # return [true_curve, path_curve]
    # return trues

def false_distribution(connections, true1, true2):
    """Takes in the entire distribution (output), returns information about path length and synaptic distance for only true positive connections"""
    #note that cell numbers need to be changed (true1 and true2) in order to accomodate which cells are truly connected
    false_curve = []
    rad_curve = []
    path_curve = []
    falses = []
    for c_to_c in connections:
        for connection in c_to_c:
            if len(connection) == 5:
                cell1 = connection[0]
                cell2 = connection[2]
                dist = float(connection[4])
                path1 = cell1[5]
                path2 = cell2[5]
            elif len(connection) == 3:
                cell1 = connection[0]
                cell2 = connection[1]
                dist = float(connection[2])
                path1 = cell1[5]
                path2 = cell2[5]
            if cell1[0] != true1 and cell2[0] != true2:
                if dist > 0:
                    false_curve.append(dist)
                    falses.append(connection)
                path_curve.append(path1)
                path_curve.append(path2)
    print(stats.describe(false_curve))
    print(stats.describe(path_curve))
    # print(np.mean(false_curve))
    # print(np.std(false_curve))
    # print(np.max(false_curve))
    # print(np.min(false_curve))
    # return [false_curve, path_curve]
    # return falses
