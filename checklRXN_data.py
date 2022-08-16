import json
import time
import random

# def deleteID(AIRoute):

# 	if("children" in AIRoute):
# 		for child in AIRoute["children"]:
# 			checkRXN_func_helper(child)

# 	AIRoute.pop("id", None)


def markAISteps(graph):

	graph["isAI"] = True

	if("children" in graph):
		for child in graph["children"]:
			markAISteps(child)



def processAIRoute(AIRoute):

	AIRoute.pop("id", None)

	failureCauseArr = ["Cannot find a route within the given price threshold for molecules", 
		"Cannot find a route within the given maximum retrosynthetic steps",
		"Cannot find a route when excluding the molecule",
		"Cannot find a route when excluding the substructure"]

	if(len(AIRoute["children"]) == 0):
		if(AIRoute["isExpandable"] == False 
			and AIRoute["isCommercial"] == False):
			AIRoute["isAvailable"] = False
			AIRoute["failureCause"] = random.choice(failureCauseArr)
		else:
			AIRoute["isAvailable"] = True
	
	for child in AIRoute["children"]:
		processAIRoute(child)


def checkRXN_func_helper(graph):

	# print("haha")

	if("children" in graph):
		for child in graph["children"]:
			checkRXN_func_helper(child)
	
	if("handledByAI" in graph and graph["handledByAI"]
		and("children" not in graph or 
		("children" in graph and len(graph["children"]) == 0))):

		with open("output.json", "r") as inputFile:
			AIRoutes = json.load(inputFile)

			for AIRoute in AIRoutes:
				markAISteps(AIRoute)

			for AIRoute in AIRoutes:
				processAIRoute(AIRoute)

			graph["AIRoutes"] = AIRoutes
			if("children" not in graph):
				graph["children"] = []

			if("children" in AIRoutes[0]):
				for child in AIRoutes[0]["children"]:
					graph["children"].append(child)
			graph["AIRoutesIndex"] = 0
			graph["confidence"] = AIRoutes[0]["confidence"]





def checkRXN_data(graph, weights):

	checkRXN_func_helper(graph)
	# time.sleep(1)


	return graph