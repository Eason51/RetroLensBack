import copy
import json
import time
import random
from RXN import RXN
from rdkit import Chem


def GetRingSystems(mol, includeSpiro=False):
    ri = mol.GetRingInfo()
    systems = []
    for ring in ri.AtomRings():
        ringAts = set(ring)
        nSystems = []
        for system in systems:
            nInCommon = len(ringAts.intersection(system))
            if nInCommon and (includeSpiro or nInCommon>1):
                ringAts = ringAts.union(system)
            else:
                nSystems.append(system)
        nSystems.append(ringAts)
        systems = nSystems
    return systems


def markAISteps(graph):

	graph["isAI"] = True

	if("children" in graph):
		for child in graph["children"]:
			markAISteps(child)


def processAIRoute(AIRoute):

	failureCauseArr = ["Precursor is too expensive", 
		"Cannot find a route within the given steps",
		"Cannot find a route when excluding the molecule",
		"Cannot find a route when excluding the substructure"]

	if("isAvailable" in AIRoute and AIRoute["isAvailable"] == False):
		AIRoute["failureCause"] = random.choice(failureCauseArr)


	# if("children" not in AIRoute or len(AIRoute["children"]) == 0):
	# 	if(AIRoute["attributes"]["isExpandable"] == False 
	# 		and AIRoute["attributes"]["isCommercial"] == False):
	# 		AIRoute["isAvailable"] = False
	# 		AIRoute["failureCause"] = random.choice(failureCauseArr)
	# 		print(AIRoute["failureCause"])
	# 	else:
	# 		AIRoute["isAvailable"] = True
	
	if("children" in AIRoute):
		for child in AIRoute["children"]:
			processAIRoute(child)


def initialize_func(smiles, constraints):

	pathArr = []

	mol = Chem.MolFromSmiles(smiles)
	ringArr = GetRingSystems(mol)
	useAI = True
	for ring in ringArr:
		if(len(ring) > 12):
			useAI = False
			break

	# if(useAI):
	# 	try: 
	# 		pathArr = RXN.predictResult(smiles, constraints)
	# 	except KeyError as e:
	# 		print("key error")
	# 		print(e)


	returnedDict = {}
	returnedDict["smiles"] = smiles
	returnedDict["children"] = []
	
	if(len(pathArr) == 0):
		returnedDict["AIFailed"] = True
		# with open("failureSmiles.txt", "a") as failureFile:
		# 	failureFile.write(smiles)

	if(len(pathArr) > 0):
		if("children" in pathArr[0]):
			for child in pathArr[0]["children"]:
				returnedDict["children"].append(child)

			returnedDict["AIRoutes"] = pathArr
			returnedDict["handledByAI"] = True
			returnedDict["AIRoutesIndex"] = 0
			returnedDict["confidence"] = pathArr[0]["confidence"]

	for path in pathArr:
		markAISteps(path)

	print("initialize complete")

	return returnedDict

	# with open("output3.json", "r") as inputFile:
	# 	graph = json.load(inputFile)

	# processAIRoute(graph)

	# return graph

