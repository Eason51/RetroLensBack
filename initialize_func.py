import json
import time
import random

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


def initialize_func(smiles):

	with open("output3.json", "r") as inputFile:
		graph = json.load(inputFile)

	processAIRoute(graph)

	# time.sleep(1)

	return graph

