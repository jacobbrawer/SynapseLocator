
import os
import numpy
import sys

# Load swc file into multidimensional array
def load_in_swc(path, filename):
	r = []
	os.chdir(path)
	with open(os.path.join(path, filename), "r") as f:
		for line in f:
			if line.lstrip().startswith('#'):
				continue
			# read values. expected SWC format is:
			# ID, type, x, y, z, rad, parent
			# x, y, z and rad are floats. the others are ints
			toks = line.split()
			vals = [
						int(toks[0]),  			# node value
						int(toks[1]), 			# type
						float(toks[2]), 		# x coordinate
						float(toks[3]),			# y coordinate
						float(toks[4]),			# z coordinate
						float(toks[5]),			# radius
						int(toks[6].rstrip())	# parent node
					]
			# store this node
			r.append(vals)
	return r
	

# Generate swc file from neuron dictionary 
def to_swc(morpho, swc_file):
	with open(swc_file, "w") as f:
		f.write("#n,type,x,y,z,radius,parent\n")
		for n in morpho:
			vals = morpho[n]
			f.write("%d %d " % (vals['id'], vals['type']))
			f.write("%0.4f " % vals['x'])
			f.write("%0.4f " % vals['y'])
			f.write("%0.4f " % vals['z'])
			f.write("%0.4f " % vals['radius'])
			f.write("%d\n" % vals['parent'])

# Creates individual neuron dictionary for connected reconstruction from starting node n
def to_dict(n, morpho, fn):
	nodes = {}
	x = 0
	first = True

	while n < len(morpho):
		if morpho[n][6] == -1: 
			if first:
				pass
			else:
				break
		
		if first:
			parent = -1
			first = False
		else:
			parent = morpho[n][6] 
			
		node_dict = {
			'id' : morpho[n][0],
			'type' : morpho[n][1],
			'x' : morpho[n][2],
			'y' : morpho[n][3],
			'z' : morpho[n][4],
			'radius' : morpho[n][5],
			'parent' : parent
		}
		nodes[x] = node_dict
		x += 1
		n+= 1
	to_swc(nodes, fn)
	return n

#___________________________________________________________________________________________________________________________#

path = raw_input("Enter the path of the swc file(s):")

for file in os.listdir(path):
	if file.find(".swc") != -1:
		full_reconstruction = load_in_swc(path, file)
		
		original_fn = os.path.join(path, file)
		print "Separating:", file
		base = original_fn[:-4]
		n = 0
		cell_count = 0

		while n < len(full_reconstruction):
			if full_reconstruction[n][6] == -1: 
				cell_count += 1
				new_fn = base + "_separated" + str(cell_count) + ".swc"
				n = to_dict(n, full_reconstruction, new_fn)
			else:
				n += 1

		
		print "The swc file has been separated into %s files" % (cell_count)
		print

