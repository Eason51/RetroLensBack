import json
import time
import random

def processAIRoute(AIRoute):

	failureCauseArr = ["Cannot find a route within the given price threshold for molecules", 
		"Cannot find a route within the given maximum retrosynthetic steps",
		"Cannot find a route when excluding the molecule",
		"Cannot find a route when excluding the substructure"]

	if(len(AIRoute["children"]) == 0):
		if(AIRoute["attributes"]["isExpandable"] == False 
			and AIRoute["attributes"]["isCommercial"] == False):
			AIRoute["isAvailable"] = False
			AIRoute["failureCause"] = random.choice(failureCauseArr)
		else:
			AIRoute["isAvailable"] = True
	
	for child in AIRoute["children"]:
		processAIRoute(child)


def checkRXN_func_helper(graph):

	if("children" in graph):
		for child in graph["children"]:
			checkRXN_func_helper(child)
	
	if("handledByAI" in graph and graph["handledByAI"]):
		with open("output.json", "r") as inputFile:
			AIRoutes = json.load(inputFile)

			for AIRoute in AIRoutes:
				processAIRoute(AIRoute)

			graph["AIRoutes"] = AIRoutes
			if("children" not in graph):
				graph["children"] = []

			graph["children"].append(AIRoutes[0])
			graph["AIRoutesIndex"] = 0



def checkRXN_func(graph):

	checkRXN_func_helper(graph)
	# time.sleep(1)


	return graph