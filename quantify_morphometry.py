# Python Script to quantify the pulmonary artery morphometry
#	Uses 1D input file for 3D segmented model.
#	Organizes morphometry based on a diameter-defined Strahler ordering model
#	Methods for morphometry based on Huang et al 1996
#		(W. Huang, R.T Yen, M McLaurine, and G. Bledsoe. Morphometry of the human pulmonary vasculature. Am Physiol Soc.1996)
# Melody Dong (06/2018)


# Import Header
import os
import sys
import csv
from collections import defaultdict
from os import listdir
from os.path import isfile, join
from tempfile import mkstemp
from shutil import move
from os import fdopen, remove
from glob import glob
import numpy as np
import matplotlib.pyplot

# Define Global Variables
nodeInfo = defaultdict(list) 		# { Node # : [X, Y, Z] }
jointNode = {} 						# {Seg # : Joint Node}
jointSeg = defaultdict(list) 		# {Inlet Seg # : [Outlet Seg #'s]}
segName = defaultdict(list) 		# {Vessel Name : [Seg #'s]'}
segNode = defaultdict(list) 		# {Seg # : [Node In, Node Out]}
segLength = {} 						# {Seg # : Seg Length}
segArea = defaultdict(list) 		# {Seg # : [Area In, Area Out]}
avgDiameters = defaultdict(float)
stdDev = defaultdict(list)


def main():	
	pname = r'/home/melody/PH/Scripts/Morphometry/quantify_morphometry/' # path name to 1D input file
	in_file = pname + 'SU0262_LPA_02.in' # input file

	# Get Relevant Segment, Node, and Joint Information from Input File
	nodeInfo = nodes(in_file)
	jointNode, jointSeg = joints(in_file)	
	segName, segNode, segLength, segArea = segments(in_file)

	# Get Morphometry-Based Segment Areas and Lengths
	bifLength, avgAreaSegment, huangSegments, avgLengthSegment, segmentsInHuangSegments = bifSegInfo(nodeInfo, jointSeg, segNode, segLength, segArea)

	# I don't think the bifDiameter is ever used?  -MD(8/1/18)
	# # Get dictionary of diameters
	# bifDiameter = defaultdict(list)
	# for name in avgAreaSegment:
	# 	diameter = np.mean(avgAreaSegment.get(name))/2*np.pi
	# 	bifDiameter[name] = diameter

	# print('jointSeg: ' + str(jointSeg))

	maxOrder = 15


	#Initial list
	#huangDiameters = (.020, .036, .056, .097, .15, .22, .34, .51, .77, 1.16, 1.75, 2.71, 4.16, 7.34, 14.80) #Average Diameters from Huang Paper
	huangDiameters = (.020, .036, .056, .097, .15, .22, .34, .51, .77, 1.16, 1.75, 2.71, 4.16, 7.34, 14.80) #Average Diameters from Huang Paper [mm]
	initDiameters = [i * 0.2 for i in huangDiameters] #Scale by scaling factor(0.518)
	huangStdDev = (.003, 0.005, 0.005, 0.012, 0.02, 0.02, 0.06, 0.04, 0.07, 0.10, 0.19, 0.35, 0.60, 1.14, 2.10) #Standard Deviations from Huang Paper
	initStdDev = [i * 0.2 for i in huangStdDev] #Scale by scaling factor(0.518)


	# Initialize average and stdev diameters for each order
	previousAverageDiameters = defaultdict(list)	# {Order #: avg Diameter}
	for i in range(1, 16):
		previousAverageDiameters[i] = 0

	previousStdDevs = defaultdict(list)
	for i in range(1, 16):
		previousStdDevs[i] = 0

	#print('Huang Diameters' + str(huangDiameters))	
	#print('Initial Diameters' + str(initDiameters))
	#print('Huang StdDev' + str(huangStdDev))
	#print('Initial StdDev' + str(initStdDev))


	#Convert Areas of each Huang Segment into Diameter of Each Huang Segment
	huangSegmentDiameters = defaultdict(list) # {Huang Segment # : Diameter}
	for vessel in avgAreaSegment:
		for index in range(0, len(avgAreaSegment.get(vessel))):
			huangSegmentDiameters[huangSegments.get(vessel)[index]] = 2*(np.sqrt(avgAreaSegment.get(vessel)[index]/np.pi))


	#Create dictionary of Lengths of each Huang Segment
	huangSegmentLengths = defaultdict(list) # {Huang Segment # : Length}
	for vessel in avgLengthSegment:
		for index in range(0, len(avgLengthSegment.get(vessel))):
			huangSegmentLengths[huangSegments.get(vessel)[index]] = avgLengthSegment.get(vessel)[index]

	# print('Huang Segment Diameters: ' + str(huangSegmentDiameters))
	

	segOrder = defaultdict(list)		# {Huang Seg # : Order}


	# Classify each segment into a diameter-based order using an initial order classification
	sortedDiameters = sorted(huangSegmentDiameters.items(), key=lambda x: x[1]) #sort diameters from smallest to largest
	currOrder = 1
	# for item in sortedDiameters:
	# 	for i in range(currOrder,16):
	# 		currOrder = i
	# 		if(item[1] < initDiameters[currOrder-1] + initStdDev[currOrder-1]):
	# 			segOrder[item[0]] = currOrder
	# 			break
	
	# for seg in sortedDiameters:
	# 	segOrder[seg[0]] = 7

	counter = 0
	currOrder = 1
	for seg in sortedDiameters:
		counter += 1
		if counter < len(sortedDiameters)/15:
			segOrder[seg[0]] = currOrder
		else:
			currOrder += 1
			counter = 0
			if currOrder > 15:
				currOrder = 15
			segOrder[seg[0]] = currOrder


	print('Seg Order: ' + str(segOrder))
	#print('Segment Diameters' + str(sortedDiameters))
	#print('Avg Area Segment' + str(avgAreaSegment))	
	#print('Previous Average Diameters ' + str(previousAverageDiameters))

	finished = False
	iterationNumber = 0
	iterations = []
	diametersPerIteration = defaultdict(list) #{Order # : [Avg Diameters]}

	

	while not finished:
		# Recalculate the average diameters and standard deviation of diameters in each order
		print('Seg Order: ' + str(segOrder))
		currOrder = 1
		currOrderSegments = []
		counter = 0
		iterationNumber += 1
		iterations.append(iterationNumber)
		print('Avg Diameters: ' + str(avgDiameters))
		print('Std Devs: ' + str(stdDev))
		for i in segOrder: #Iterate through all Huang Segments -> sorted from smallest diameter to largest
			print('i: ' + str(i))
			print('Seg Order: ' + str(segOrder.get(i)))
			counter+=1
			if(segOrder.get(i) == currOrder): #create np array of all diameters in current order
				currOrderSegments.append(huangSegmentDiameters.get(i))
			else:
				npCurrOrderSegments = np.asarray(currOrderSegments)
				print('npCurrOrderSegments: ' + str(npCurrOrderSegments))
				if not np.isnan(np.mean(npCurrOrderSegments)): #calculate average diameter for each order
					avgDiameters[currOrder] = np.mean(npCurrOrderSegments)
					diametersPerIteration[currOrder].append(np.mean(npCurrOrderSegments))
				if not np.isnan(np.std(npCurrOrderSegments)): #calculate average standard deviation for each order
					stdDev[currOrder] = np.std(npCurrOrderSegments)
				currOrderSegments = []
				currOrder = segOrder.get(i)
				currOrderSegments.append(huangSegmentDiameters.get(i))
			print('Current Order Segments: ' + str(currOrderSegments))

		npCurrOrderSegments = np.asarray(currOrderSegments)
		avgDiameters[currOrder] = np.mean(npCurrOrderSegments)
		stdDev[currOrder] = np.std(npCurrOrderSegments)

		print('Average Diameters ' + str(avgDiameters))
		print('Standard Deviations ' + str(stdDev))


		# Repeat classification based on new averages and standard deviation of order diameter
		#	Segment maintain order if:
		#		(D_(n-1) + SD_(n-1) + D_n - SD_n)/2 < Di < (D_n + SD_n + D_(n+1) - SD_(n+1))/2
		#			D_(n-1) = average diameter of order one less than current segment's order
		#			SD_(n-1) = standard deviation of order one less than current segment's order
		#			D_n = average diameter of current segment's order
		#			SD_n = standard deviation of current segment's order
		#			D_(n+1) = average diameter of order one greater than current segment's order
		#			SD_(n+1) = standard deviation of order one greater than current segment's order
		# Continue iterating until change in diameter and change in standard deviation are less than 1%
		currOrder = 1

		for segment in segOrder: #iterate through all Huang Segments
			currOrder = segOrder.get(segment)
			if currOrder > 1:
				#print(str(avgDiameters.get(currOrder-1, 0) + stdDev.get(currOrder-1, 0) + avgDiameters.get(currOrder) - stdDev.get(currOrder)/2) + " > " + str(huangSegmentDiameters.get(segment)))
			
				if ((avgDiameters.get(currOrder-1, 0) + stdDev.get(currOrder-1, 0) + avgDiameters.get(currOrder) - stdDev.get(currOrder))/2 > huangSegmentDiameters.get(segment)):
					print('True')
					segOrder[segment] = currOrder - 1
				
			if currOrder < maxOrder:
				#print(str(avgDiameters.get(currOrder) + stdDev.get(currOrder) + avgDiameters.get(currOrder+1, 0) - stdDev.get(currOrder+1, 0)) + " < " + str(huangSegmentDiameters.get(segment)))
			
				if ((avgDiameters.get(currOrder) + stdDev.get(currOrder) + avgDiameters.get(currOrder+1, 0) - stdDev.get(currOrder+1, 0)) < huangSegmentDiameters.get(segment)):
					print('True')
					segOrder[segment] = currOrder + 1
					if segOrder.get(segment) > maxOrder:
						maxOrder = segOrder.get(segment)
		
		#Check if change in diameter and change in standard deviation are less than 1%
		diameterWithin1Percent = True
		stdDevWithin1Percent = True
		for order in range(1, maxOrder+1): #Iterate through all orders
			divBy = 1
			if previousAverageDiameters.get(order, 0) != 0:
				divBy = previousAverageDiameters.get(order, 0)
			# print(str(previousAverageDiameters.get(order, 0) - avgDiameters.get(order, 0)))
			# print('Prev Avg Di: ' + str(previousAverageDiameters[order]))
			# print('Avg Di: ' + str(avgDiameters[order]))
			if (np.abs(previousAverageDiameters.get(order, 0) - avgDiameters.get(order, 0)))/divBy > .0001:
				diameterWithin1Percent = False
				break
			divBy = 1
			if(previousStdDevs.get(order, 0) != 0):
				divBy = previousStdDevs.get(order, 0)
			if (np.abs(previousStdDevs.get(order, 0) - stdDev.get(order, 0)))/divBy > .0001:
				stdDevWithin1Percent = False
				break

		if diameterWithin1Percent and stdDevWithin1Percent:
			finished = True

		previousAverageDiameters = avgDiameters.copy()
		previousStdDevs = stdDev.copy()

	#print('Std Devs: ' + str(stdDev))
	# print('Iterations: ' + str(iterations))
	# print('Diameters Per Iteration: ' + str(diametersPerIteration))
	currOrder = 1
	currOrderSegments = []
	counter = 0
	avgLengths = defaultdict(list)
	for iteration in diametersPerIteration:
		matplotlib.pyplot.plot(iterations, diametersPerIteration[iteration])
	matplotlib.pyplot.show()

	# print('Num Segments: ' + str(len(segName)))
	# print('Node Info: ' + str(len(nodeInfo)))

	# print('Segment Orders: ' + str(segOrder))
	# print('Average Diameter: ' + str(avgDiameters))
	# print('Std Devs: ' + str(stdDev))
	# print('Huang Segment Lengths' + str(huangSegmentLengths))

	#Calculate average length of vessels in each order
	isFirst = True
	length = 0
	currOrder = 1
	counter = 0
	for currSeg in segOrder:
		if isFirst:
			currOrder = segOrder.get(currSeg)
			isFirst = False
		if segOrder.get(currSeg) == currOrder:
			length += huangSegmentLengths.get(currSeg)
			counter += 1
		else:
			avgLengths[currOrder] = length/counter
			length = 0
			counter = 0
			currOrder += 1
			length += huangSegmentLengths.get(currSeg)
			counter +=1
	avgLengths[currOrder] = length/counter

	numSegmentsInOrder = defaultdict(list)
	currOrder = 1
	numTimes = 0
	for elem in segOrder:
		if segOrder.get(elem) == currOrder:
			numTimes+=1
		else:
			numSegmentsInOrder[currOrder] = numTimes
			numTimes = 1
			currOrder +=1
	#print(numSegmentsInOrder)
	maxOrder = len(numSegmentsInOrder)

	#print('Num Segments: ' + str(len(bifDiameter)))

	
	#print('Average Lengths: ' + str(avgLengths))
	#print('Num Segments in Order: ' + str(numSegmentsInOrder))
	#print('Avg Diameters: ' + str(avgDiameters))

	#Calculate number of Huang elements per Order
	isFirst = True
	counter = 0
	elementsPerOrder = defaultdict(list) #{Order # : # Elements in Order}
	for order in avgDiameters: #Iterate through all orders
		counter = 0
		for vessel in huangSegments: #Iterate through all vessels
			if not isFirst and currOrder == order:
				counter += 1
			isFirst = True
			for huangSeg in huangSegments.get(vessel): #Iterate through Huang Segments in each vessel
				if isFirst:
					currOrder = segOrder.get(huangSeg)
					isFirst = False
				if currOrder != segOrder.get(huangSeg) and currOrder == order:
					counter += 1
					currOrder = segOrder.get(huangSeg)
		elementsPerOrder[order] = counter

	#print('Elements Per Order: ' + str(elementsPerOrder))

	#Get Order of each SimVascular Segment
	svSegOrder = defaultdict(int) #{SV Seg # : Order #}
	for huangSegment in segOrder:
		for svSegment in segmentsInHuangSegments.get(huangSegment):
			svSegOrder[svSegment] = segOrder.get(huangSegment)



	#Create Connectivity Matrix
	connectivityMatrix = np.zeros((maxOrder+2, maxOrder+2))
	for i in range(0,maxOrder+2): #Set first row and first column to be labeled -> first column = parent artery orders; first row = child artery orders
		connectivityMatrix[0][i] = i
		connectivityMatrix[i][0] = i

	parentOrder = 1
	childOrder = 1

	# for vessel in segName: #iterate through each vessel
	# 	for segment in segName.get(vessel): #iterate through all SimVascular Segments (because nodes are based on SimVascular segments, not Huang segments)
	# 		numOut = len(jointSeg[segment])
	# 		if numOut > 1: #checks for bifurcation case
	# 			for childArtery in jointSeg[segment]:
	# 				if np.absolute(int(childArtery) - int(segment)) > 1: #checks if part of same vessel
	# 					connectivityMatrix[svSegOrder.get(int(segment))][svSegOrder.get(int(childArtery))] += 1
	# print(str(svSegOrder))
	#print(str(connectivityMatrix))
	# print(str(segOrder))
	# print(str(stdDev))


	#Divides each element of the connectivity matrix by the number of elements in the parent order
	for parentOrder in range(1,len(numSegmentsInOrder)+1):
		for childOrder in range(1,len(numSegmentsInOrder)+1):
			if(numSegmentsInOrder.get(parentOrder) > 0):
				connectivityMatrix[parentOrder][childOrder] /= elementsPerOrder.get(parentOrder, 1)
	
	# print(str(connectivityMatrix))
	# print('Average Diameters: ' + str(avgDiameters))
	# print('Std Devs: ' + str(stdDev))
	# print('Average Lengths: ' + str(avgLengths))

	#print('Sorted Diameters: ' + str(sortedDiameters))




# Because we only want to quantify segments between bifurcation points, 
# we need to find the length and areas of the segments in between bifurcations.
# This involves finding the segments in between bifurcations and adding them 
# together to create one long segment
def bifSegInfo(nodeInfo, jointSeg, segNode, segLength, segArea):
	# Define new, updated dictionaries + edge cases
	newNode = nodeInfo.copy() # make a copy of nodeInfo Dictionary
	nodeInd = 0
	newNode[int(nodeInd)] = newNode['0'] # edge case to get beginning node
	newNodeInfo = []
	prevOutNode = int('0')
	newSegInfo = []
	newSegInd = 0
	jointInd = 0
	avgAreaSegment = defaultdict(list) # {Vessel Name: [Avg Huang Segment Area ...]}
	avgLengthSegment = defaultdict(list) # {Vessel Name: [Avg Huang Segment Length ...]}
	avgLength = defaultdict(list)
	bifLength = defaultdict(list)
	huangSegments = defaultdict(list) #{Vessel Name: [Huang Segment #'s]}
	huangIndex = 0
	segmentsInHuangSegments = defaultdict(list) #{Huang Segment # : Segment #'s}


	# Find Bifurcation Cases
	for name in segName: #iterate through all vessels
		temp = segName[name]
		inSeg =  temp[0]
		outSeg = temp[-1]
		bifEndSeg = []
		bifBegSeg = [inSeg]
		bifFE = []
		bifArea = []
		avgBifArea = []
		avgBifLength = []
		huangSegmentsInVessel = []
		segments = []

		vesLength = 0
		vesFE = 0
		for seg in segName[name]: #iterate through all segments in vessel
			
			vesLength = vesLength + float(segLength[seg]) #sum segment length along vessel

			numOut = len(jointSeg[seg])
			if numOut > 1: #more than 1 outlet = bifurcation Case
				# Append Inlet and Outlet Segment of Vessel segments 
				bifEndSeg.append(seg) #append inlet segment number where bifurcation occurs
				bifBegSeg.append(str(int(seg)+1))

				# Sum the segment lengths until a bifurcation occurs, append non-bifurcating length
				bifLength[name] = vesLength
				avgBifLength.append(vesLength)
				vesLength = 0

		# Append edge segments 
		bifEndSeg.append(outSeg)
		bifLength[name] = vesLength
		bifFE.append(vesFE)
		avgBifLength.append(vesLength)
		avgLengthSegment[name] = avgBifLength


		# Construct the new segment card: new areas, nodes
		for ind in range(0,len(bifEndSeg)):

			# Fix Middle Nodes to join non-bifurcation vessel segments
			inNode = segNode[bifBegSeg[ind]]
			inNode = inNode[0] #define inlet node of abbrev segment
			outNode = segNode[bifEndSeg[ind]]
			outNode = outNode[-1] #define outlet node of abbrev segment

			newNode = fixMiddleNodes(inNode, outNode,newNode)

			node1 = nodeInd
			if not(int(inNode) == int(prevOutNode)): # Checks if previously defined segment connects with current segment
				nodeInd += 1
				newNode[int(nodeInd)] = newNode[inNode]
				node1 = nodeInd

			nodeInd += 1
			newNode[int(nodeInd)] = newNode[outNode] #defines keys as integer node #'s and abbreviates node list

			node2 = nodeInd
			
			for index in range(int(bifBegSeg[ind]), int(bifEndSeg[ind]) + 1):
				segments.append(index)
			segmentsInHuangSegments[huangIndex] = segments
			segments = []

			# New Areas
			Ain = segArea[bifBegSeg[ind]]
			Ain = Ain[0]
			Aout = segArea[bifEndSeg[ind]]
			Aout = Aout[1]

			bifArea.append([float(Ain), float(Aout)]) #append new areas to vessel segments' areas

			# Update Info
			prevOutNode = outNode
			newSegInd += 1
			huangSegmentsInVessel.append(huangIndex)
			huangIndex += 1
			avgBifArea.append(np.mean([float(Ain), float(Aout)]))


		avgAreaSegment[name] = avgBifArea
		huangSegments[name] = huangSegmentsInVessel
		# print(name)name
		# print('Segment Length' + str(bifLength))
		# print('Areas ' + str(bifArea))
		# print('Avg Area Segment ' + str(avgAreaSegment[name]))
	return(bifLength, avgAreaSegment, huangSegments, avgLengthSegment, segmentsInHuangSegments)





# Remove middle nodes of Abbreviated Segment - Maintain In and Out Node Coordinates and Joint Node
def fixMiddleNodes(inNode, outNode, newNode):
	midNodes = range(int(inNode)+1,int(outNode))

	for node in nodeInfo:
		if((int(node) in midNodes)):
			del newNode[node]

	return(newNode)


# Define Node INFO
def nodes(in_file):
	with open(in_file,'r') as old_in:
		for line in old_in:
			if line.find("NODE ")>=0 and not(line.find("# ")>=0):
				node = line.split(" ")
				nodeInfo[node[1]] = [float(node[2]), float(node[3]), float(node[4])] # { Node # : [X, Y, Z] }

	return nodeInfo


# Define Joint INFO
def joints(in_file):

	with open(in_file,'r') as inFile:
		for line in inFile:
			if line.find("JOINT ")>=0 and not(line.find("# ")>=0):
				temp_line1 = line.split() #Node info for joint
				line2 = next(inFile)
				temp_line2 = line2.split() #Inlet segment for joint
				line3 = next(inFile)
				temp_line3 = line3.split() #Outlet segments for joint

				jointNode[temp_line2[3]] = temp_line1[2] # {Seg # : Joint Node}
				jointSeg[temp_line2[3]] = temp_line3[3:len(temp_line3)+1] # {Inlet Seg # : Outlet Seg #}

	return jointNode, jointSeg


# Define Segment INFO
def segments(in_file):

	with open(in_file,'r') as inFile:
		for line in inFile:
			if line.find("SEGMENT ")>=0 and not(line.find("# ")>=0):
				temp_line = line.split() # segment info

				name = temp_line[1]
				name = name.split('_')
				ves_name = name[0]
				for n in name[1:-1]:
					ves_name = ves_name + "_" + n

				segName[ves_name].append(temp_line[2]) # {Seg Name : Seg #}
				segNode[temp_line[2]] = temp_line[5:7] # {Seg # : [Node In, Node Out]}
				segLength[temp_line[2]] = temp_line[3] # {Seg # : Seg Length}
				segArea[temp_line[2]] = temp_line[7:9] # {Seg # : [Area In, Area Out]}

	return segName, segNode, segLength, segArea




main()



