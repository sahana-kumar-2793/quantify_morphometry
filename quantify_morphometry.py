# Python Script to quantify the pulmonary artery morphometry
#	Uses 1D input file for 3D segmented model.
#	Organizes morphometry based on a diameter-defined Strahler ordering model
#	Methods for morphometry based on Huang et al 1996
#		(W. Huang, R.T Yen, M McLaurine, and G. Bledsoe. Morphometry of the human pulmonary vasculature. Am Physiol Soc.1996)
# Melody Dong (06/2018)
# Additions to Huang/Strahler Morphometry method by Sahana Kumar (08/2018)


# Import Header
import os
import System
import csv
from collections import defaultdict
from os import listdir
from os.path import isfile, join
from tempfile import mkstemp
from shutil import move
from os import fdopen, remove
from glob import glob
import numpy as np
# import matplotlib.pyplot
import csv

# Define Global Variables
nodeInfo = defaultdict(list) 		# { Node # : [X, Y, Z] }
jointNode = {} 						# {Seg # : Node #}
jointSeg = defaultdict(list) 		# {Inlet Seg # : [Outlet Seg #'s]}
segName = defaultdict(list) 		# {Vessel Name : [Seg #'s]'}
segNode = defaultdict(list) 		# {Seg # : [Node In, Node Out]}
segLength = {} 						# {Seg # : Seg Length}
segArea = defaultdict(list) 		# {Seg # : [Area In, Area Out]}
avgDiameters = defaultdict(float)	# {Order : Average Diameter}
stdDev = defaultdict(list)			# {Order : Std Dev for Diameter}


def main():	
	pname = r'C:/Users/melody/Documents/Marsden_Research/Scripts/Morphometry/quantify_morphometry/Input_Files/' # path name to 1D input file
	in_file = pname + 'N_15M.in' # input file


	# Get Relevant Segment, Node, and Joint Information from Input File
	nodeInfo = nodes(in_file)
	jointNode, jointSeg = joints(in_file)	
	segName, segNode, segLength, segArea = segments(in_file)

	# Get Morphometry-Based Segment Areas and Lengths
	bifLength, avgAreaSegment, huangSegments, avgLengthSegment, segmentsInHuangSegments = bifSegInfo(nodeInfo, jointSeg, segNode, segLength, segArea)

	maxOrder = 15


	# Modify 1D segments to adjust for actual bifurcation positions instead of grouping multiple into a simvascular segment
	# numSegments = 0
	# for inletSeg in segName:
	# 	if len(segName.get(inletSeg)) > 1:
	# 		numSegments += len(segName.get(inletSeg)) - 1

	# isFirst = True
	# currSegNum = 0
	# inletID = 0
	# nodeSeg = defaultdict(list) # {Vessel # : Node Seg #s}
	# nodeSegLength = defaultdict(list) # {Node Seg # : Node Seg Length}
	# nodeSegArea = defaultdict(float) # {Node Seg # : Node Seg Area}
	# nodeSegConnections = defaultdict(list) # {Inlet Node Seg # : Outlet Node Seg #'s'}
	# for vessel in segName: #Iterates through all vessels
	# 	isFirst = True
	# 	currInlet = []
	# 	currOutlet = []
	# 	vesNodeSeg = []
	# 	for inletSeg in segName.get(vessel): #Iterates through all SV Segments
	# 		if isFirst:
	# 			currInlet = nodeInfo.get(inletSeg)
	# 			currInletArea = segArea.get(inletSeg)[0]
	# 			isFirst = False
	# 			continue
	# 		if len(jointSeg.get(inletSeg)) > 1: 
	# 			for outletSeg in jointSeg.get(inletSeg):
	# 				vesNodeSeg.append(currSegNum)
	# 				currOutlet = nodeInfo.get(outletSeg)
	# 				currOutletArea = segArea.get(outletSeg)[1]
	# 				nodeSegArea[currSegNum] = (float(currInletArea) + float(currOutletArea))/2 #Calculates average area in nodeSeg
	# 				nodeSegLength[currSegNum] = np.sqrt(np.power(currOutlet[0] - currInlet[0], 2) + np.power(currOutlet[1] - currInlet[1], 2) + np.power(currOutlet[2] - currInlet[2], 2)) #Calculates distance between 2 nodes - new length
	# 				# nodeSegConnections[inletID] = currSegNum
	# 				# inletID = currSegNum
	# 				currSegNum += 1
	# 				currInlet = currOutlet
	# 		currOutlet = nodeInfo.get(inletSeg)
		

	# 	vesNodeSeg.append(currSegNum) # Append all nodeSegs in current vessel
	# 	currOutletArea = segArea.get(outletSeg)[1]
	# 	nodeSegArea[currSegNum] = (float(currInletArea) + float(currOutletArea))/2
	# 	nodeSegLength[currSegNum] = np.sqrt(np.power(currOutlet[0] - currInlet[0], 2) + np.power(currOutlet[1] - currInlet[1], 2) + np.power(currOutlet[2] - currInlet[2], 2))
	# 	# nodeSegConnections[inletID] = currSegNum
	# 	# inletID = currSegNum
	# 	currSegNum += 1 
	# 	nodeSeg[vessel] = vesNodeSeg

	# nodeSegDiameter = {k : 2*np.sqrt(v/(np.pi)) for k, v in nodeSegArea.items()}


	# Average diameters and stdev for Huang SEGMENTS [mm]
	huangDiameters = (.020, .036, .056, .097, .15, .22, .35, .51, .77, 1.17, 1.78, 2.81, 4.33, 7.31, 15.12) #Average Diameters of ELEMENTS from Huang Paper [mm]
	initDiameters = [i * 0.1 for i in huangDiameters] #convert to cm
	huangStdDev = (.003, 0.006, 0.005, 0.015, 0.02, 0.03, 0.09, 0.06, 0.10, 0.14, 0.25, 0.46, 0.68, 1.38, 1.81) #Standard Deviations from Huang Paper
	initStdDev = [i * 0.1 for i in huangStdDev] #convert to cm



	# Create dictionary of diameters of each Huang Segment
	huangSegmentDiameters = defaultdict(list) # {Huang Segment # : Diameter}
	for vessel in avgAreaSegment:
		for index in range(0, len(avgAreaSegment.get(vessel))):
			# Convert Areas of each Huang Segment into Diameter of Each Huang Segment
			huangSegmentDiameters[huangSegments.get(vessel)[index]] = 2*(np.sqrt(avgAreaSegment.get(vessel)[index]/np.pi))

	# Create dictionary of lengths of each Huang Segment
	huangSegmentLengths = defaultdict(list) # {Huang Segment # : Length}
	for vessel in avgLengthSegment:
		for index in range(0, len(avgLengthSegment.get(vessel))):
			huangSegmentLengths[huangSegments.get(vessel)[index]] = avgLengthSegment.get(vessel)[index]
	
	# Create dictionary of orders of each Huang Segment
	huangSegOrder = defaultdict(lambda: 0)		# {Huang Seg # : Order}



	#Initial Classification Schemes
	# (A) Classify all segments initially based on Huang's diameter order classification
	initHuang = False
	if initHuang:
		huangSegOrder = initHuangDiamOrder(huangSegmentDiameters, initDiameters, initStDev, huangSegOrder)

	# (B) Classify each segment into a diameter-based order from greatest to least
	# Keep order if within 15% of the greatest diameter in that order
	initLarge = False
	if initLarge:
		huangSegOrder = initLargestDiam(huangSegmentDiameters, huangSegOrder, segmentsInHuangSegments)
	
	# (C) Classify each segment into an initial order using the Strahler Ordering System
	initStrah = True
	if initStrah:
		huangSegOrder, svSegmentToHuangSegment = initStrahler(huangSegmentDiameters, segmentsInHuangSegments, huangSegOrder)

	# (D) Classify each segment into a random order from 1 to 15
	initRand = False
	if initRand:
		huangSegOrder = initRandom(huangSegmentDiameters, huangSegOrder)
	

	if initStrah:
		#Find first sv segment in each segmental
		svSegmentals = list()
		segmentals = list()

		for svSeg in segName.get('LPA'):
			if len(jointSeg.get(svSeg)) > 2:
				for outletSeg in jointSeg.get(svSeg):
					if np.abs(int(outletSeg)-int(svSeg)) > 1:
						svSegmentals.append(outletSeg)


		svSegmentalsToCut = len(svSegmentals) - 8

		for i in range(1, svSegmentalsToCut+1):
			vesselToRemove = svSegmentals[0]
			smallestHuangSeg = svSegmentToHuangSegment.get(svSegmentals[0])
			smallestDiameter = huangSegmentDiameters.get(smallestHuangSeg)
			for segmental in svSegmentals:
				currHuangSeg = svSegmentToHuangSegment.get(segmental)
				currDiameter = huangSegmentDiameters.get(currHuangSeg)
				if currDiameter < smallestDiameter:
					vesselToRemove = segmental
			svSegmentals.remove(vesselToRemove)


		svSegToVessel = defaultdict(str) #{SV Seg # : Vessel Name}
		for vessel in segName:
			for svSeg in segName.get(vessel):
				svSegToVessel[svSeg] = vessel



		#Find Segmentals
		for svSeg in svSegmentals:
			segmentals.append(svSegToVessel.get(svSeg))



		#Scale so segmentals have order 14
		greatestSegmentalOrder = 0
		for svSeg in svSegmentals:
			if huangSegOrder.get(int(svSegmentToHuangSegment.get(int(svSeg)))) > greatestSegmentalOrder:
				greatestSegmentalOrder = huangSegOrder.get(int(svSegmentToHuangSegment.get(int(svSeg))))

		add = 14 - greatestSegmentalOrder
		print('Huang Seg Order: ' + str(huangSegOrder))

		for huangSeg in huangSegOrder:
			huangSegOrder[huangSeg] += add




	finished = False
	iterationNumber = 0
	iterations = []

	# List of Huang Segments sorted by their Order from Smallest to Largest
	sortedSeg = sorted(huangSegOrder, key=huangSegOrder.get)
	sortedDiameters = sorted(huangSegmentDiameters.items(), key=lambda x: x[1])
	sortedSegOrder = sorted(huangSegOrder.items(), key=lambda x:x[1])

	
	# Calculate average order's diameters using initial ordering classification scheme
	currOrder = 1
	currOrderSegments = []
	for i in sortedSeg: #Iterate through all Huang Segments -> sorted from smallest to largest order
		if(huangSegOrder.get(i) == currOrder): #append diameters in current order to np array
			currOrderSegments.append(huangSegmentDiameters.get(i))
		else: #if segment not in current order, calculate the average diameters/stdev of all diameters in current order, move onto next order
			npCurrOrderSegments = np.asarray(currOrderSegments)
			if len(npCurrOrderSegments) >= 1: #ensures that array for current order is not empty
				avgDiameters[currOrder] = np.mean(npCurrOrderSegments) 
				stdDev[currOrder] = np.std(npCurrOrderSegments) #stdev diameter
			
			#Change current order to the order of the current segment that was different from previous
			currOrderSegments = []
			currOrder = huangSegOrder.get(i) 
			currOrderSegments.append(huangSegmentDiameters.get(i))

	# Edge case for last order in list
	npCurrOrderSegments = np.asarray(currOrderSegments)
	avgDiameters[currOrder] = np.mean(npCurrOrderSegments)
	stdDev[currOrder] = np.std(npCurrOrderSegments)

	previousAverageDiameters = avgDiameters.copy()
	previousStdDevs = stdDev.copy()	


	
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
	counter = 0
	while not finished:
		counter += 1
		currOrder = 6
		for segment in huangSegOrder: #iterate through all Huang Segments
			currOrder = huangSegOrder.get(segment)
			if currOrder > 1: #Check if Huang segment's diameter is smaller than the lower bound cutoff
				#if segment == 0:
					# print(str((avgDiameters.get(currOrder-1, 0) + stdDev.get(currOrder-1, 0) + avgDiameters.get(currOrder) - stdDev.get(currOrder))/2) + " > " + str(huangSegmentDiameters.get(segment)))
				if (avgDiameters.get(int(currOrder-1),0) + stdDev.get(int(currOrder-1),0) + avgDiameters.get(int(currOrder)) - stdDev.get(int(currOrder)))/2 > huangSegmentDiameters.get(segment):
					huangSegOrder[segment] = currOrder - 1
				
			if currOrder < maxOrder: #Check if Huang segment's diameter greater than the upper bound cutoff
				if ((avgDiameters.get(int(currOrder)) + stdDev.get(int(currOrder)) + avgDiameters.get(int(currOrder+1), 0) - stdDev.get(int(currOrder+1), 0)) < huangSegmentDiameters.get(segment)):
					huangSegOrder[segment] = currOrder + 1
					if huangSegOrder.get(segment) > maxOrder: #If order exceeds the max order, force the order to be the max Order
						huangSegOrder[segment] = maxOrder

		# Recalculate the average diameters and standard deviation of diameters in each order
		currOrder = 1
		currOrderSegments = []
		iterationNumber += 1
		iterations.append(iterationNumber)
		sortedSegOrder = sorted(huangSegOrder.items(), key=lambda x:x[1])
		orders = []
		avgDiameters.clear()
		for segment in sortedSegOrder: #Iterate through all Huang Segments -> sorted from smallest diameter to largest
			orders.append(currOrder)
			if(huangSegOrder.get(segment[0]) == currOrder): #append diameters in current to np array
				currOrderSegments.append(huangSegmentDiameters.get(segment[0]))
			else: #if segment not in current order, calculate the average diameters/stdev of all diameters in current order
				npCurrOrderSegments = np.asarray(currOrderSegments)
				if len(npCurrOrderSegments) >= 1: #ensures that array for current order is not empty
					avgDiameters[currOrder] = np.mean(npCurrOrderSegments) #avg diameter
					stdDev[currOrder] = np.std(npCurrOrderSegments) #stdev diameter
				currOrderSegments = []
				currOrder = huangSegOrder.get(segment[0]) #change current order to the order of the current segment that was different from previous
				currOrderSegments.append(huangSegmentDiameters.get(segment[0]))



		npCurrOrderSegments = np.asarray(currOrderSegments)
		avgDiameters[currOrder] = np.mean(npCurrOrderSegments)
		stdDev[currOrder] = np.std(npCurrOrderSegments)
		
		#Check if change in diameter and change in standard deviation are less than 1%
		diameterWithin1Percent = True
		stdDevWithin1Percent = True
		for order in range(1, maxOrder+1): #Iterate through all orders
			divBy = 1
			if previousAverageDiameters.get(order, 0) != 0:
				divBy = previousAverageDiameters.get(order, 0)
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



	# Calculate average length of vessels in each order
	avgLengths = defaultdict(float)
	isFirst = True
	length = 0
	currOrder = 1
	counter = 0
	sortedSegOrder = sorted(huangSegOrder.items(), key=lambda x:x[1])
	for currSeg in sortedSegOrder:
		if isFirst:
			currOrder = currSeg[1]
			isFirst = False
		if currSeg[1] == currOrder:
			length += huangSegmentLengths.get(currSeg[0])
			counter += 1
		else:
			avgLengths[currOrder] = length/counter
			length = 0
			counter = 0
			currOrder = currSeg[1]
			length += huangSegmentLengths.get(currSeg[0])
			counter +=1
	avgLengths[currOrder] = length/counter

	# maxOrder = len(numSegmentsInOrder)

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
					currOrder = huangSegOrder.get(huangSeg)
					isFirst = False
				if currOrder != huangSegOrder.get(huangSeg) and currOrder == order:
					counter += 1
					currOrder = huangSegOrder.get(huangSeg)
		elementsPerOrder[order] = counter


	#Get Order of each SimVascular Segment
	svSegOrder = defaultdict(int) #{SV Seg # : Order #}
	for huangSegment in huangSegOrder:
		for svSegment in segmentsInHuangSegments.get(huangSegment):
			svSegOrder[svSegment] = huangSegOrder.get(huangSegment)

	sortedDiameters = sorted(huangSegmentDiameters.items(), key=lambda x: x[1])
	sortedSegOrder = sorted(huangSegOrder.items(), key=lambda x:x[1])

	#Create Connectivity Matrix
	connectivityMatrix = np.zeros((maxOrder+2, maxOrder+2))
	for i in range(0,maxOrder+2): #Set first row and first column to be labeled -> first column = parent artery orders; first row = child artery orders
		connectivityMatrix[0][i] = i
		connectivityMatrix[i][0] = i

	parentOrder = 1
	childOrder = 1

	for vessel in segName: #iterate through each vessel
		for segment in segName.get(vessel): #iterate through all SimVascular Segments (because nodes are based on SimVascular segments, not Huang segments)
			numOut = len(jointSeg[segment])
			if numOut > 1: #checks for bifurcation case
				for childArtery in jointSeg[segment]:
					if np.absolute(int(childArtery) - int(segment)) > 1: #checks if part of same vessel
						connectivityMatrix[svSegOrder.get(int(segment))][svSegOrder.get(int(childArtery))] += 1



	#Divides each element of the connectivity matrix by the number of elements in the parent order
	for parentOrder in range(1,maxOrder+1):
		for childOrder in range(1,maxOrder+1):
			if elementsPerOrder.get(parentOrder, 1) != 0:
				connectivityMatrix[parentOrder][childOrder] /= elementsPerOrder.get(parentOrder, 1)
	
	print(str(huangSegOrder))


	# with open('connectivity_matrix.csv', 'w', newline='') as csvfile:
 # 		writer = csv.writer(csvfile, delimiter=' ', quotechar=',', quoting=csv.QUOTE_ALL)
 # 		for row in connectivityMatrix:
 # 			writer.writerow(row)




# (D) Classify each segment into a random order from 1 to 15
def initRandom(huangSegmentDiameters, huangSegOrder):
	sortedDiameters = sorted(huangSegmentDiameters.items(), key=lambda x: x[1]) #sort diameters from smallest to largest
	counter = 0
	currOrder = 1
	for seg in sortedDiameters:
		counter += 1
		if counter < len(sortedDiameters)/15:
			huangSegOrder[seg[0]] = currOrder
		else:
			currOrder += 1
			counter = 0
			if currOrder > 15:
				currOrder = 15
			huangSegOrder[seg[0]] = currOrder

	return(huangSegOrder)




# (C) Classify each segment into an initial order using the Strahler Ordering System
#	Terminal branches = order 1
#	When two branches of order n converge to 1 branch, converged branch = order n+1
#	If a higher order branch converges with a lower order branch, converged branch = higher order
def initStrahler(huangSegmentDiameters, segmentsInHuangSegments, huangSegOrder):
	# Correlate simvascular segment #'s to their corresponding huang segment #'s
	svSegmentToHuangSegment = defaultdict(int) #{SV Segment # : Huang Segment #}
	for huangSegment in segmentsInHuangSegments:
		for svSegment in segmentsInHuangSegments.get(huangSegment):
			svSegmentToHuangSegment[svSegment] = huangSegment

	
	# Define terminal branches as the end segment with no more joints (attached segment) as order 1
	for huangSegment in huangSegmentDiameters:
		svSegments = segmentsInHuangSegments.get(int(huangSegment))
		if len(svSegments) > 0:
			endSegment = svSegments[len(svSegments) - 1]
			if len(jointSeg.get(str(endSegment))) == 0:
				huangSegOrder[huangSegment] = 1


	finished = False
	counter = 0
	while not finished:
		finished = True
		counter += 1
		for parentSeg in huangSegmentDiameters: #iterate through all huang segments in a relatively random order
			if int(huangSegOrder[parentSeg]) == 0:
				finished = False
				svSegmentsInParentSeg = segmentsInHuangSegments.get(parentSeg)
				if len(svSegmentsInParentSeg) == 0:
					continue

				parentSVSeg = svSegmentsInParentSeg[len(svSegmentsInParentSeg) - 1]
				childSVSegments = jointSeg.get(str(parentSVSeg))
				childHuangOrders = []
				for childSVSeg in childSVSegments:
					childHuangSeg = svSegmentToHuangSegment.get(int(childSVSeg))
					if huangSegOrder[int(childHuangSeg)] == 0:
						childHuangOrders = []
						break
					else:
						childHuangOrders.append(huangSegOrder.get(childHuangSeg))

				if len(childHuangOrders) != 0:
					currOrder = childHuangOrders[0]
					same = True
					if len(childHuangOrders) == 1:
						huangSegOrder[parentSeg] = currOrder
					else:
						for i in range(1, len(childHuangOrders)):
							if childHuangOrders[i] > currOrder:
								currOrder = childHuangOrders[i]
								same = False
							if childHuangOrders[i] < currOrder:
								same = False
						if not same:
							huangSegOrder[parentSeg] = currOrder
						else:
							huangSegOrder[parentSeg] = currOrder + 1

		if counter == 1000:
			finished = True
			print('Timed Out')

	return(huangSegOrder, svSegmentToHuangSegment)




# (B) Classify each segment into a diameter-based order from greatest to least
#	keep order if within 15% of the greatest diameter in that order
def initLargestDiam(huangSegmentDiameters, huangSegOrder, segmentsInHuangSegments):
	sortedDiameters = sorted(huangSegmentDiameters.items(), key=lambda x: x[1], reverse=True) #sort diameters from largest to smallest
	currOrder = 15
	maxOrder = 15
	isFirst = True
	smallestOrder = 0
	greatestDiameter = 0 
	for segment in sortedDiameters:
		if isFirst:
			huangSegOrder[segment[0]] = currOrder
			greatestDiameter = segment[1]
			isFirst = False
			continue 
		if greatestDiameter*0.85 <= segment[1]:
			huangSegOrder[segment[0]] = currOrder
		else:
			currOrder -= 1
			huangSegOrder[segment[0]] = currOrder
			greatestDiameter = segment[1]
			smallestOrder = currOrder


	if smallestOrder < 1:
		add = 1 - smallestOrder
		for segment in huangSegOrder:
			huangSegOrder[segment] += add


	currOrder = 15
	finished = False
	counter = 0
	while not finished:
		counter += 1
		finished = True
		for huangSegment in range(0, len(huangSegmentDiameters)):
			if huangSegment == 0:
				huangSegOrder[huangSegment] = currOrder
				continue
			firstSeg = segmentsInHuangSegments.get(huangSegment)[0]
			parentSVSeg = 0
			for inletSeg in jointSeg:
				for outletSeg in jointSeg.get(inletSeg):
					if int(firstSeg) == int(outletSeg):
						parentSVSeg = inletSeg
			parentHuangSeg = 0
			for huangSeg in segmentsInHuangSegments:
				for seg in segmentsInHuangSegments.get(huangSeg):
					if int(seg) == int(parentSVSeg):
						parentHuangSeg = huangSeg
						break
			parentDiameter = huangSegmentDiameters.get(parentHuangSeg)
			if parentDiameter*.85 < huangSegmentDiameters.get(huangSegment) and parentDiameter*1.15 > huangSegmentDiameters.get(huangSegment):
				if huangSegOrder[parentHuangSeg] > 0:
					huangSegOrder[huangSegment] = huangSegOrder.get(parentHuangSeg)
				else:
					finished = False
			elif huangSegmentDiameters.get(huangSegment) > parentDiameter * 1.15:
				if huangSegOrder[parentHuangSeg] > 0:
					huangSegOrder[huangSegment] = huangSegOrder.get(parentHuangSeg) + 1
				else:
					finished = False
			else:
				if huangSegOrder[parentHuangSeg] > 0:
					huangSegOrder[huangSegment] = huangSegOrder.get(parentHuangSeg) - 1
				else:
					finished = False

	return(huangSegOrder)




# (A) Classify each segment into a diameter-based order using a Huang's order classification as the initial classification
# scheme starting at order 1
#	Sorts all segments into an initial classification based on the diameters of Huang's segment orders
#	If less than order n's mean diam + stdev diam, then defined as order n.  If not, keeps incrementing orders
#	until finds an order that is within range of the segment's diameter.
def initHuangDiamOrder(huangSegmentDiameters, initDiameters, initStDev, huangSegOrder):
	sortedDiameters = sorted(huangSegmentDiameters.items(), key=lambda x: x[1]) #sort diameters from smallest to largest
	currOrder = 1
	maxOrder = 15
	for item in sortedDiameters: # iterate through all Huang segments from smallest to largest
		for i in range(currOrder,16): # iterate through all orders starting with smallest order
			currOrder = i
			# Define segment of current order if diameter is within range
			if(item[1] < (initDiameters[currOrder-1] + initStdDev[currOrder-1])):
				huangSegOrder[item[0]] = currOrder
				break
			elif (item[1] > (initDiameters[maxOrder-1] + initStdDev[maxOrder-1])):
				huangSegOrder[item[0]] = maxOrder
				break

	return(huangSegOrder)




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
			
			tempArea = []
			# Track/organize simvascular segments into their respective huang segments
			for index in range(int(bifBegSeg[ind]), int(bifEndSeg[ind]) + 1):
				segments.append(index)
				area = segArea[str(index)]
				tempArea.append(float(area[0]))
			segmentsInHuangSegments[huangIndex] = segments
			tempArea.append(float(area[1]))
			segments = []

			# New Areas
			Ain = segArea[bifBegSeg[ind]]
			Ain = Ain[0]
			Aout = segArea[bifEndSeg[ind]]
			Aout = Aout[1]

			# MD 08/17/18
			# average through all segments in huang segments for huang area
			avgArea = np.mean(tempArea)

			bifArea.append([float(Ain), float(Aout)]) #append new areas to vessel segments' areas

			# Update Info
			prevOutNode = outNode
			newSegInd += 1
			huangSegmentsInVessel.append(huangIndex)
			huangIndex += 1
			avgBifArea.append(avgArea) #np.mean([float(Ain), float(Aout)]))


		# Track average area over all simvascular segments in Huang segment and the simvascular segment #'s in a huang segment
		avgAreaSegment[name] = avgBifArea
		huangSegments[name] = huangSegmentsInVessel

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