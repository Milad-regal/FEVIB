import time
import re
import sys
import os
import numpy as np

def findAnsysFileInFolder():
	textFileList = []
	currentDirectory = os.path.dirname(os.path.realpath(__file__))
	for file in os.listdir(currentDirectory):
		if file.endswith(".txt"):
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
		subList.append(string.split()[0])
		coords = re.findall('[\-\+]?\w+\.\w*[\-\+]?\d*',string)
		if dim == '2D':
			coords = coords[0:2]			
		subList.extend(coords)
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
MaterialID = 0
MaterialCC = -1 # count the four mp parameters since rho and thk have no names !
MaterialArray = np.zeros((50,4)) # Note - not using append - max 50 materials !

for line in f:
	#---------------FETCH NODES-------------
	if line[0:14] == 'NBLOCK,6,SOLID':
		recNodesStart = True
	if  line[0:10] == 'N,R5.3,LOC':
		recNodesStop = True
	if recNodesStart and not recNodesStop:
		nodeList.append(line)

	#---------------FETCH ELEMENTS-------------
	if line[0:15] == 'EBLOCK,19,SOLID':
		recElementsStart = True
	if  line[0:12] == 'EN,R5.5,ATTR' or line[0:12] == 'MPTEMP,R5.0,':
		recElementsStop = True
	if recElementsStart and not recElementsStop:
		elementList.append(line)

	# Get materials
	if line[0:2] == 'MP':
		# Keep count of the materials... hack !	
		MaterialCC=MaterialCC+1
		if MaterialCC>3:
			MaterialID=MaterialID+1
			MaterialCC=0
		
		# Check the material parametres - expects a certain order !!!!!	
		if MaterialCC==0:
			# thk
			MaterialArray[MaterialID,0] = float(line.split(',')[3])
		if MaterialCC==1:
			# E
			MaterialArray[MaterialID,1] = float(line.split(',')[6])
		if MaterialCC==2:
			# dens
			MaterialArray[MaterialID,3] = float(line.split(',')[3])
		if MaterialCC==3:
			# nu
			MaterialArray[MaterialID,2] = float(line.split(',')[6])

	#---------------FETCH BOUNDARY-------------
	if line[0:10] == 'ERESX,DEFA':
		recBoundaryStart = True
	if  line[0:6] == 'FINISH':
		recBoundaryStop = True
	if recBoundaryStart and not recBoundaryStop:
		boundaryforceList.append(line)

f.close()

nodeListNew = nodeList
nodeList = splitList2(nodeList[2:])

# reduce size of the material array
MaterialArray = MaterialArray[0:MaterialID+1,:].tolist()
MaterialArray = [[str(y) for y in x] for x in MaterialArray]

#---------------BLANK FIELDS NODELIST----------
tempList = []
for element in nodeList:
	if len(element) == 2:
		element.append('0')
	#elif len(element) == 3:
	#	del element[-1]
	tempList.append(element)
nodeList = tempList

elementList = splitList(elementList[2:-1])
boundaryforceList = splitList(boundaryforceList[2:-1])

#-------------- ADD NODE NUMBER ---------------
#for element1,element2 in zip(nodeList,nodeListNew[2:-1]):
#	nodeList.insert(0,element2[0])
#	print nodeList

#---------------APPEND MATERIAL---------------
tempElementList = []
for line in elementList:
	tempLine = line
	if not line:
		print 'no line 0'
		time.sleep(2)
	tempLine.append(line[0])
	tempElementList.append(tempLine)
elementList = tempElementList

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
		dof = line.pop(1).split(',')
		line.insert(1,dof[0])
		line.insert(2,dof[1])
		if line[0] == 'D,':
			boundaryLine = makeLine(line)
			boundaryList.append(boundaryLine)
		elif line[0] == 'F,':
			forceLine = makeLine(line)
			forceList.append(forceLine)

#------------ CONSTRUCT OUTPUTFILE ---------------
outputText = []

if len(forceList) == 1:
	delLoadChoice = raw_input('One PointLoad found. Only use PointLoad as wdof? \"n\"/[Enter]:')
	outputText.append('% Inputfile generated from ANSYS-file. ')
	outputText.append('mesh.wdof = [' + forceList[0][0] + ' ' + forceList[0][1]  +  '];\n')
else:
	outputText.append('% Inputfile generated from ANSYS-file. \n')

outputText.append('% Coordinates (x,y):')
outputText.extend(makeMatlabFile('mesh.X',nodeList,0,len(nodeList)))
#outputText.append('\n')
outputText.append('% Topology matrix IX(node1,node2,node3,node4)')
outputText.extend(makeMatlabFile('mesh.IX',elementList,10,16))
#outputText.append('\n')
outputText.append('% Boundary conditions mat(node,ldof,disp)')
outputText.extend(makeMatlabFile('mesh.bound',boundaryList,0,len(boundaryList)))
#outputText.append('\n')
outputText.append('% Prescribed point loads mat(node,ldof,force)')
outputText.extend(makeMatlabFile('mesh.PointLoads',forceList,0,len(boundaryList)))
#outputText.append('\n')
outputText.append('% Prescribed pressure loads mat(element,ldof,force)')
outputText.extend(makeMatlabFile('mesh.PressureLoads',pressureList,0,len(pressureList)))
outputText.append('% Material table)')
outputText.extend(makeMatlabFile('mesh.Material',MaterialArray,0,20))
#for element in outputText:
#	print element

outputFile = inputFile.split('.txt')[0] + '.m'
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
		print '\nFile overwritten'
	else:
		print '\nFile written'
else:
	print 'No file written.'

f.close()


