import time
import re
import sys
import os
import numpy as np

def findAnsysFileInFolder():
	textFileList = []
	currentDirectory = os.path.dirname(os.path.realpath(__file__))
	for file in os.listdir(currentDirectory):
		if file.endswith(".dat"):
			textFileList.append(str(file))

	modTextFileList = []
	ii = 1
	spaces = ' '*50

	numberOfFilesPrinted = 10

	for element in textFileList[0:numberOfFilesPrinted]:
		spacedLine = element + spaces[len(element):-len(str(ii))] + str(ii)
		modTextFileList.append(spacedLine)
		ii += 1

	textFiles = "\n".join(modTextFileList)

	if len(textFileList) > 0:
		print '\n--------- Textfiles found in current directory: ---------'
		headLine = 'Name:' + ' '*50
		print headLine[0:-(len('Name:')+len('FileNr:'))] + 'FileNr:'
		print textFiles
		#dispText = dispText + '\n' + textFiles + '\n\n'

		dispText = '\nFileNr or [Enter] :'
		inputNumber = raw_input(dispText)
		try:
			inputNumber = int(float(inputNumber))
		except:
			inputNumber = -1

		testCond = inputNumber < numberOfFilesPrinted and \
		 (inputNumber < len(textFileList)+1 and inputNumber > -1)

		if testCond:
			inputFile = textFileList[inputNumber-1]
		else:
			inputFile = raw_input('Inputfile from ANSYS: ')

	return (inputFile,currentDirectory)

def checkForOverwrite(currentDirectory,outputFile):
	foundFile = False
	for file in os.listdir(currentDirectory):
		if str(file) == outputFile:
			foundFile = True

	return foundFile


def splitList(listIn):
	outList = []
	subList = []
	for string in listIn:
		outList.append(string.split())
	return outList

def splitList2(listIn):
	outList = []
	dim = '2D'	#POSSIBLE TO EXTEND TO 3D
	for string in listIn:
		subList = []
		#print [float(i) for i in string.split(',')[1:]]
		coords = [float(i) for i in string.split(',')[1:]]
		if dim == '2D':
			coords = coords[0:3]			
		coords = [str(i) for i in coords]
		subList.extend(coords)
		#print coords
		#print subList
		outList.append(subList)
	return outList


def makeMatlabFile(fileType,prepList,fromPos,toPos):
	outputText = []
	if len(prepList) == 0:
		outputText.append(fileType+ ' = [ ' + '  '+ '];\n')
	elif len(prepList) == 1:
		outputText.append(fileType+ ' = [ ' + '  '.join(prepList[0][fromPos:toPos]) + '];\n')
	elif len(prepList)>1:
		outputText.append(fileType+ ' = [ ' + '  '.join(prepList[0][fromPos:toPos]))
		for element in prepList[1:-1]:
			tempLine = '\t  '  + '  '.join(element[fromPos:toPos])
			outputText.append(tempLine)
		outputText.append('\t  ' + '  '.join(prepList[-1][fromPos:toPos])+ '];\n')
	else:
		outputText = '[]\n'
	return outputText


(inputFile,currentDirectory) = findAnsysFileInFolder()

f = open(inputFile,'r')

recNodesStart = False
recNodesStop = False
recElementsStart = False
recElementsStop = False
recBoundaryStart = False
recBoundaryStop = False
nodeList = []
elementList = []
boundaryforceList = []
boundaryList = []
forceList = []
pressureList = []
MaterialCC = 0
MaterialID = 1
MaterialArray = np.zeros((50,4)) # Note - not using append - max 50 materials !

for line in f:
	#---------------FETCH NODES-------------
	if line[0:19] == '! Nodal coordinates':
		recNodesStart = True
	if  line[0:22] == '! Element connectivity':
		recNodesStop = True
	if recNodesStart and not recNodesStop:
		nodeList.append(line)

	#---------------FETCH ELEMENTS AND MATERIAL ID+DATA -------------
	if line[0:11] == '! Material ':
		#print int(line.split(' ')[2])
		MaterialID = int(line.split(' ')[2]);
		MaterialCC = MaterialCC + 1
	if line[0:3] == 'MP,':
		mline = line.strip().split(',')
		if mline[1].strip()=='EX':
			MaterialArray[int(MaterialID)-1,1] = mline[3]
		if mline[1].strip()=='PRXY':
			MaterialArray[int(MaterialID)-1,2] = mline[3]
		if mline[1].strip()=='DENS':
			MaterialArray[int(MaterialID)-1,3] = mline[3]
	if line[0:2] == 'R,':
		mline = line.strip().split(',')
		MaterialArray[int(MaterialID)-1,0] = mline[2]
	if line[0:14] == '! Connectivity':
		recElementsStart = True
	if  line[0:19] == '! Nodal diplacement':
		recElementsStop = True
	if recElementsStart and not recElementsStop:
		if line[0:3] == 'EN,':
			nline = line[4:-1] + '   ' +  str(MaterialID)
			elementList.append(nline)
		
	#---------------FETCH BOUNDARY-------------
	if line[0:19] == '! Nodal diplacement': #line[0:10] == 'ERESX,DEFA':
		recBoundaryStart = True
	if  line[0:6] == 'FINISH':
		recBoundaryStop = True
	if recBoundaryStart and not recBoundaryStop:
		if line[0:2]=='D,' or line[0:2]=='F,' or line[0:4]=='SFE,':
			#print line
			boundaryforceList.append(line)

f.close()

# Clean lists
nodeListNew = nodeList
nodeList = splitList2(nodeList[1:])

#print elementList

# clean material table
#print MaterialArray[0:10,:]
MaterialArray = MaterialArray[0:MaterialCC,:].tolist()
MaterialArray = [[str(y) for y in x] for x in MaterialArray]
#print MaterialArray

#---------------BLANK FIELDS NODELIST----------
tempList = []
for element in nodeList:
	if len(element) == 2:
		element.append('0')
	#elif len(element) == 3:
	#	del element[-1]
	tempList.append(element)
nodeList = tempList

elementList = splitList(elementList)
boundaryforceList = splitList(boundaryforceList)

#-------------- ADD NODE NUMBER ---------------
#for element1,element2 in zip(nodeList,nodeListNew[2:-1]):
#	nodeList.insert(0,element2[0])
#	print nodeList

def setBoundary(line):
	tempLine = line
	dof = line.pip(1).slplit(',')
	line.insert(1,dof)
	if line[0] == 'D':
		line[1]

#---------------BOUNDARY/FORCES----------------
def makeLine(line):
	outLine = []
	outLine.append(line[1])
	if line[2][1] == 'X':
		outLine.append('1')
	elif line[2][1] == 'Y':
		outLine.append('2')
	elif line[2][1] == 'Z':
		outLine.append('3')
	outLine.append(line[4])
	return outLine



lineTwoOfPressure = False
for line in boundaryforceList:
	if line[0] == 'SFE,':
		pressOpts = line[2].split(',')
		if pressOpts[1] == 'PRES':
			pressNode = line[1]
			lineTwoOfPressure = True
	elif lineTwoOfPressure:
		if pressOpts[2] == '1': 	#Uncertant what meaning this has. Only [1,PRES,1,R5.0] that has load value.
			pressLine = []
			pressLine.append(pressNode)
			pressLine.extend(line)
			pressureList.append(pressLine)
		lineTwoOfPressure = False
	else:
		if line[2]=='UX,': # SIMPLE
			dof = ['1']
		elif line[2]=='UY,': # CLAMPED
			dof = ['1','2','3']
		else:
			dof = ['1'] # SIMPLE IF NOTHING IS SPECIFIED !!!!!
		if line[0] == 'D,':
			#boundaryLine = makeLine(line)
			for i in dof:
				mline = [str((line[1].strip(',')) + "  " + i + "  " + str(float(line[3].strip(','))))]
				#print mline
				boundaryList.append(mline)
		elif line[0] == 'F,':
			mline = [str((line[1].strip(',')) + "  1  " + str(float(line[3].strip(','))))]
			forceList.append(mline)

#------------ CONSTRUCT OUTPUTFILE ---------------
outputText = []

if len(forceList) == 1:
	delLoadChoice = raw_input('One PointLoad found. Only use PointLoad as wdof? \"n\"/[Enter]:')
	outputText.append('% Inputfile generated from FlexTract-file. ')
	outputText.append('mesh.wdof = [' + forceList[0][0] + '];\n')
else:
	outputText.append('% Inputfile generated from FlexTract-file. \n')
        outputText.append('mesh.wdof = [' + forceList[0][0] + '];\n')

outputText.append('% Coordinates (x,y):')
outputText.extend(makeMatlabFile('mesh.X',nodeList,0,20))
#outputText.append('\n')
outputText.append('% Topology matrix IX(node1,node2,node3,node4,matID)')
outputText.extend(makeMatlabFile('mesh.IX',elementList,0,20))
#outputText.append('\n')
outputText.append('% Boundary conditions mat(node,ldof,disp)')
outputText.extend(makeMatlabFile('mesh.bound',boundaryList,0,20))
#outputText.append('\n')
outputText.append('% Prescribed point loads mat(node,ldof,force)')
outputText.extend(makeMatlabFile('mesh.PointLoads',forceList,0,20))
#outputText.append('\n')
outputText.append('% Prescribed pressure loads mat(element,ldof,force)')
outputText.extend(makeMatlabFile('mesh.PressureLoads',pressureList,0,20))
outputText.append('% Material table)')
outputText.extend(makeMatlabFile('mesh.Material',MaterialArray,0,20))

outputFile = inputFile.split('.dat')[0] + '.m'
foundFile = checkForOverwrite(currentDirectory,outputFile)
if foundFile:
	writeChoice = raw_input('File found with name: \"' + outputFile + '\"\nOverwrite? \"n\"/[Enter]: ')
else:
	writeChoice = 'w'

if not writeChoice or writeChoice == 'w':
	f = open(outputFile, 'w')
	for line in outputText:
		f.write(line)
		f.write('\n')
	if not writeChoice:
		print '\nFile overwritten: ' + outputFile
	else:
		print '\nFile written: ' + outputFile
else:
	print 'No file written: ' + outputFile

f.close()


