import json
import time
import random
from RXN import RXN
import threading
from concurrent.futures import ThreadPoolExecutor
import globalVar
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


def executeRXN(graph, constraints, seconds):

	if(seconds != 0):
		time.sleep(seconds)
	pathArr = []

	globalVar.runningTime -= 20

	try:
		pathArr = RXN.predictResult(graph["smiles"], constraints)
	except KeyError as e:
		print("keyerror")
		print(e)
	except Exception as e:
		print("exception")
		print(e)

	# with open("output4.json", "w") as outputFile:
	# 	outputFile.write(json.dumps(pathArr))


	if(len(pathArr) == 0):
		graph["AIFailed"] = True

	if(len(pathArr) > 0):
		graph["children"] = []
		if("children" in pathArr[0]):
			for child in pathArr[0]["children"]:
				graph["children"].append(child)

			graph["AIRoutes"] = pathArr
			graph["AIRoutesIndex"] = 0
			graph["confidence"] = pathArr[0]["confidence"]


	for path in pathArr:
		markAISteps(path)



def checkRXN_func_helper(graph, constraints):

	nodeArr = []

	if("children" in graph):
		for child in graph["children"]:
			returnValue = checkRXN_func_helper(child, constraints)
			if(returnValue != None):
				nodeArr.append(returnValue)

	nodeTaskArr = []
	max_workers = 3
	nodeCount = 0
	with ThreadPoolExecutor(max_workers=max_workers) as executor:
		for node in nodeArr:
			task = None
			time.sleep(20)
			if(nodeCount <= max_workers):
				task = executor.submit(executeRXN, node, constraints, 0)
			else:
				task = executor.submit(executeRXN, node, constraints, 5)
			nodeTaskArr.append([node, task])


	for [node, task] in nodeTaskArr:
		if(task.exception() != None):
			print(node["smiles"])
			print(task.exception())



	if("handledByAI" in graph and graph["handledByAI"]
		and("children" not in graph or 
		("children" in graph and len(graph["children"]) == 0))):

		smiles = graph["smiles"]
		inputMol = Chem.MolFromSmiles(smiles)
		inputSmiles = Chem.MolToSmiles(inputMol)

		with open("failureSmiles.txt", "r") as failureFile:
			for failureSmiles in failureFile:
				if(failureSmiles != ""):
					failureMol = Chem.MolFromSmiles(failureSmiles)
					failureSmiles = Chem.MolToSmiles(failureMol)
					if(failureSmiles == inputSmiles):
						print("failure Smiles: ", smiles)
						graph["AIFailed"] = True
						return None

		ringArr = GetRingSystems(inputMol)
		for ring in ringArr:
			if(len(ring) > 12):
				print("failure Smiles: ", smiles)
				graph["AIFailed"] = True
				return None


		return graph

	return None




def checkRXN_func(graph, constraints):

	checkRXN_func_helper(graph, constraints)
	# time.sleep(1)

	print("checkRXN complete")
	return graph