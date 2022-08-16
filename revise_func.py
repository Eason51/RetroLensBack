import itertools
import json
from criteria import confidence
from criteria import complexity_reduction
from criteria import convergence
from criteria import reactant_confidence
from scaffold import extractScaffold
from operator import itemgetter
from datetime import datetime


# rank = 1

# def revise_func_helper(graph):

# 	if("children" not in graph or len(graph["children"]) == 0):
# 		return

# 	global rank

# 	if(rank >= 10):
# 		return
		
# 	graph["rank"] = rank
# 	rank  = rank + 1
# 	graph["SAW"] = 10 - rank
# 	graph["influence"] = 0.2
# 	graph["reactionConfidence"] = 0.3
# 	graph["complexity"] = 0.4
# 	graph["convergence"] = 0.5
# 	graph["associatedSubtreeConfidence"] = 0.6


# 	if("children" in graph):
# 		for child in graph["children"]:
# 			revise_func_helper(child)


# 	# if("notAvailable" in graph):
# 	# 	graph["notAvailable"] = False


def removeRank(graph):

	if("rank" in graph):
		graph.pop("rank", None)
	
	if("children" in graph):
		for child in graph["children"]:
			removeRank(child)
	

def calculateWeights(weights):
	
	totalWeights = 0

	if(weights["influence"] == None):
		weights["influence"] = 0
		totalWeights += 0
	else:
		totalWeights += weights["influence"]
	
	if(weights["complexity"] == None):
		weights["complexity"] = 0
		totalWeights += 0
	else:
		totalWeights += weights["complexity"]
	
	if(weights["convergence"] == None):
		weights["convergence"] = 0
		totalWeights += 0
	else:
		totalWeights += weights["convergence"]
	
	if(weights["reactionConfidence"] == None):
		weights["reactionConfidence"] = 0
		totalWeights += 0
	else:
		totalWeights += weights["reactionConfidence"]
	
	if(weights["associatedSubtreeConfidence"] == None):
		weights["associatedSubtreeConfidence"] = 0
		totalWeights += 0
	else:
		totalWeights += weights["associatedSubtreeConfidence"]


	if(totalWeights == 0):
		weights["influence"] = 0.2
		weights["complexity"] = 0.2
		weights["convergence"] = 0.2
		weights["reactionConfidence"] = 0.2
		weights["associatedSubtreeConfidence"] = 0.2
	else:
		weights["influence"] /= totalWeights
		weights["complexity"] /= totalWeights
		weights["convergence"] /= totalWeights
		weights["reactionConfidence"] /= totalWeights
		weights["associatedSubtreeConfidence"] /= totalWeights


	

def countNode(graph):

	nodeCount = 1
	if("children" in graph):
		for child in graph["children"]:
			nodeCount += countNode(child)
	
	return nodeCount


def resetAISteps(graph):

	if("AIRegion" in graph):
		graph.pop("AIRegion", None)

	if("children" in graph):
		for child in graph["children"]:
			resetAISteps(child)


def resetProblemSteps(graph):

	if("problem" in graph):
		graph.pop("problem", None)
		
	if("candidate" in graph):
		graph.pop("candidate", None)

	if("SAW" in graph):
		graph.pop("SAW", None)

	if("aggregateConfidence" in graph):
		graph.pop("aggregateConfidence", None)

	if("childSubtreeConfidenceArr" in graph):
		graph.pop("childSubtreeConfidenceArr", None)

	if("complexity" in graph):
		graph.pop("complexity", None)

	if("convergence" in graph):
		graph.pop("convergence", None)

	if("influence" in graph):
		graph.pop("influence", None)

	if("normalizedInfluence" in graph):
		graph.pop("normalizedInfluence", None)

	if("rank" in graph):
		graph.pop("rank", None)

	if("reactionConfidence" in graph):
		graph.pop("reactionConfidence", None)

	if("children" in graph):
		for child in graph["children"]:
			resetProblemSteps(child)


def seperateHumanAISteps(graph):

	if(("handledByAI" in graph and graph["handledByAI"]
		and "AIRoutes" in graph and len(graph["AIRoutes"]) != 0)
		or ("AIRegion" in graph and graph["AIRegion"])):

		graph["AIRegion"] = True

		if("children" in graph):
			for child in graph["children"]:
				child["AIRegion"] = True

	if("children" in graph):
		for child in graph["children"]:
			seperateHumanAISteps(child)


def updateReactantConfidence_helper(graph):

	reactantStepArr = []

	if("children" in graph and len(graph["children"]) != 0):

		if("AIRegion" in graph and graph["AIRegion"]):
			if(("handledByAI" in graph and graph["handledByAI"]
				and "AIRoutes" in graph and len(graph["AIRoutes"]) != 0)):

				reactantStepArr.append(graph)
		else:
			reactantStepArr.append(graph)


		for child in graph["children"]:
			reactantStepArr += updateReactantConfidence_helper(child)
	
	return reactantStepArr


def updateReactantConfidence(graph):

	reactantStepArr = updateReactantConfidence_helper(graph)

	reactant_confidence(reactantStepArr)




def updateGraph(graph, weights, totalNodes):

	print("updateGraph")

	if("children" not in graph or len(graph["children"]) == 0):
		
		graph["influence"] = 1
		graph["normalizedInfluence"] = 1 / totalNodes

		if(("notAvailable" in graph and graph["notAvailable"])
			or ("isAvailable" in graph and not graph["isAvailable"])):
			

			graph["problem"] = True

		return

	print("handle children")

	if("children" in graph):
		for child in graph["children"]:
			updateGraph(child, weights, totalNodes)

		print("influnece")

		# influence
		influence = 1
		for child in graph["children"]:
			influence += child["influence"]
		
		graph["influence"] = influence
		graph["normalizedInfluence"] = influence / totalNodes


		reactantArr = []
		for child in graph["children"]:
			reactantArr.append(child["smiles"])

		print("reactionConfidence")

		if("AIRegion" in graph and graph["AIRegion"]):
			pass
		else:
			print("reaction confidence-------------------------------------------", end="\n\n\n")
			# reaction confidence
			now = datetime.now()
			current_time = now.strftime("%H:%M:%S")
			print("Start Time =", current_time)
			reactionConfidence = confidence(reactantArr, graph["smiles"])
			graph["reactionConfidence"] = reactionConfidence
			now = datetime.now()
			current_time = now.strftime("%H:%M:%S")
			print("Finish Time =", current_time)
			# graph["reactionConfidence"] = 0.5

		print("complexity")

		# complexity
		complexity = complexity_reduction(reactantArr, graph["smiles"])
		graph["complexity"] = complexity

		print("convergence")

		# convergence
		convergence_value = convergence(reactantArr, graph["smiles"])
		graph["convergence"] = convergence_value


	print("mark problem branches")

	# mark problem branches
	if("children" in graph):

		if(("notAvailable" in graph and graph["notAvailable"])
			or ("isAvailable" in graph and not graph["isAvailable"])):

			graph["problem"] = True

		for child in graph["children"]:

			# if ("notAvailable" in child and child["notAvailable"]):
			# 	print("hahaahah")

			if(("notAvailable" in child and child["notAvailable"])
				or ("isAvailable" in child and not child["isAvailable"])
				or ("problem" in child and child["problem"])):

				graph["problem"] = True


	print("aggregateConfidence")

	# aggregateConfidence
	if("AIRegion" in graph and graph["AIRegion"]):

		if("handledByAI" in graph and graph["handledByAI"]
			and "AIRoutes" in graph and len(graph["AIRoutes"]) != 0):

			# graph["aggregateConfidence"] = graph["reactionConfidence"]
			graph["aggregateConfidence"] = graph["confidence"]

	else:

		if("children" in graph and len(graph["children"]) != 0):
			graph["aggregateConfidence"] = graph["reactionConfidence"]
			for child in graph["children"]:
				if("aggregateConfidence" in child):
					graph["aggregateConfidence"] *= child["aggregateConfidence"]

		else:
			pass

	
	print("subtreeConfidence")

	# subtreeConfidence
	if("problem" in graph and graph["problem"] and
		"children" in graph and "AIRegion" not in graph):

		print("start subtreeConfidence")

		#[childIndex, subtreeConfidence]
		childSubtreeConfidenceArr = []
		childIndex = 0
		for child in graph["children"]:
			if("problem" in child and child["problem"]):
				
				totalInfluence = 0
				tempIndex = 0
				print(1)
				for child in graph["children"]:
					print(2)
					if(tempIndex == childIndex):
						tempIndex += 1
						continue
					
					print(3)
					totalInfluence += child["influence"]
					
					tempIndex += 1

				print("total influence", totalInfluence)
				subtreeConfidence = 0
				tempIndex = 0
				print(4)
				for child in graph["children"]:
					print(5)
					if(tempIndex == childIndex):
						tempIndex += 1
						continue

					if("problem" in child and child["problem"]):
						subtreeConfidence += 0
						tempIndex += 1
						continue

					print(6)
					if("aggregateConfidence" in child):
						print(child["aggregateConfidence"])
						subtreeConfidence += child["aggregateConfidence"] * (child["influence"] / totalInfluence)
					else:
						subtreeConfidence += reactant_confidence(child["smiles"]) * (child["influence"] / totalInfluence)

					tempIndex += 1

				childSubtreeConfidenceArr.append([childIndex, subtreeConfidence])

			childIndex += 1


		graph["childSubtreeConfidenceArr"] = childSubtreeConfidenceArr

		nodeSubtreeConfidence = 0
		totalInfluence = 0
		for [childIndex, subtreeConfidence] in graph["childSubtreeConfidenceArr"]:
			totalInfluence += graph["children"][childIndex]["influence"]
		
		for [childIndex, subtreeConfidence] in graph["childSubtreeConfidenceArr"]:
			nodeSubtreeConfidence += subtreeConfidence * (graph["children"][childIndex]["influence"] / totalInfluence)

		graph["associatedSubtreeConfidence"] = nodeSubtreeConfidence


	if("problem" in graph and graph["problem"] and
		"children" in graph and "AIRegion" not in graph):
		print()
		print(graph["smiles"])
		print(graph["childSubtreeConfidenceArr"])
		print()

	print("SAW Score")

	#SAW Score
	if("problem" in graph and graph["problem"] and
		"children" in graph and "AIRegion" not in graph
		and "childSubtreeConfidenceArr" in graph
		and len(graph["childSubtreeConfidenceArr"]) != 0):


		# [childIndex, sawScore]

		childSawArr = []
		for [childIndex, subtreeConfidence] in graph["childSubtreeConfidenceArr"]:

			sawScore = 0
			sawScore += (graph["influence"] / totalNodes) * weights["influence"]
			sawScore += graph["complexity"] * weights["complexity"]
			sawScore += graph["convergence"] * weights["convergence"]
			sawScore += graph["reactionConfidence"] * weights["reactionConfidence"]
			sawScore += subtreeConfidence * weights["associatedSubtreeConfidence"]

			childSawArr.append([childIndex, sawScore])					

		graph["childSawArr"] = childSawArr

		sawScore = 0
		totalInfluence = 0
		for [childIndex, childSawScore] in graph["childSawArr"]:
			totalInfluence += graph["children"][childIndex]["influence"]
		for [childIndex, childSawScore] in graph["childSawArr"]:
			sawScore += childSawScore * (graph["children"][childIndex]["influence"] / totalInfluence)

		graph["SAW"] = sawScore
		print("SAW score", graph["SAW"])

	
	print("node finished", end="\n")




def getCandidates(graph):

	candidateArr = []
	returnedSmiles = ""

	print("isavailable", "isAvailable" in graph and not graph["isAvailable"])
	print("notAvailable", ("notAvailable" in graph and graph["notAvailable"]))

	if(("isAvailable" in graph and not graph["isAvailable"]) or
		("notAvailable" in graph and graph["notAvailable"])):
		print("error", len(candidateArr))
		return [candidateArr, extractScaffold(graph["smiles"])]

	currentScaffold = extractScaffold(graph["smiles"])

	if("children" in graph and "problem" in graph and graph["problem"]):
		for child in graph["children"]:

			# print("before", len(candidateArr))
			# print("child")

			if("problem" in child and child["problem"]):

				print("child has problem")

				[childCandidateArr, problemSmiles] = getCandidates(child)
				if(len(problemSmiles) != 0):
					if(problemSmiles == currentScaffold):
						returnedSmiles = problemSmiles
					else:
						if("AIRegion" in graph and graph["AIRegion"]):
							returnedSmiles = ""
						else:
							if(graph not in candidateArr):
								candidateArr.append(graph)
								graph["candidate"] = True
							returnedSmiles = ""

				else:
					if("AIRegion" in graph and graph["AIRegion"]):
						candidateArr = []
					else:
						if(graph not in candidateArr):
							candidateArr.append(graph)
							graph["candidate"] = True
						candidateArr += childCandidateArr

			# else:
			# 	candidateArr.append(graph)
			# 	returnedSmiles = ""

			# print("after", len(candidateArr))
	
	return [candidateArr, returnedSmiles]



def getCandidateCombo(graph):

	candidateComboArr = []

	if("candidate" in graph and graph["candidate"]):
		
		candidateComboArr.append([{"id": graph["id"], "SAW": graph["SAW"]}])
		if("children" in graph):
			# [ [[combo1], [combo2]], [[combo1], [combo2]] ]
			tempCandidateComboArr = []
			for child in graph["children"]:
				if("candidate" in child and child["candidate"]):
					# [[combo1], [combo2], [combo3]]
					childCandidateComboArr = getCandidateCombo(child)
					if(len(childCandidateComboArr) != 0):
						tempCandidateComboArr.append(childCandidateComboArr)
				else:
					if("problem" in child and child["problem"]):
						return candidateComboArr


			if(len(tempCandidateComboArr) != 0):
				comboArr = list(itertools.product(*tempCandidateComboArr))
				processedComboArr = []
				for combo in comboArr:
					comboItem = []
					for item in combo:
						for step in item:
							comboItem.append(step)
					processedComboArr.append(comboItem)
				
				for combo in processedComboArr:
					candidateComboArr.append(combo)
	
	return candidateComboArr
			



def mergeSAW(graph):
	candidateComboArr = getCandidateCombo(graph)

	print("candidateComboArr", candidateComboArr)
	print("\n")
	
	comboSawArr = []
	for combo in candidateComboArr:

		combo = sorted(combo, key=itemgetter("SAW"))

		combination = {}
		combination["combo"] = []
		combination["SAW"] = 0
		for step in combo:
			combination["combo"].append(step["id"])
			combination["SAW"] += step["SAW"]
		
		comboSawArr.append(combination)
	
	sortedComboSawArr = sorted(comboSawArr, key=itemgetter("SAW"))
	print("sortedComboSawArr", sortedComboSawArr)
	return sortedComboSawArr
		



def rankSteps(graph):

	[candidateArr, problemSmiles] = getCandidates(graph)

	sortedCandidateArr = sorted(candidateArr, key=itemgetter("SAW"))

	rank = 1
	for candidate in sortedCandidateArr:
		print("udpate")
		candidate["rank"] = rank
		rank += 1

	for candidate in sortedCandidateArr:
		print("rank", candidate["rank"])
	

	


def revise_func(graph, weights):

	# global rank
	# rank = 1

	# revise_func_helper(graph)

	removeRank(graph)
	calculateWeights(weights)
	totalNodes = countNode(graph)
	resetAISteps(graph)
	resetProblemSteps(graph)
	seperateHumanAISteps(graph)
	# updateReactantConfidence(graph)
	updateGraph(graph, weights, totalNodes)
	with open("temp.json", "w") as outputFile:
		outputFile.write(json.dumps(graph))
	rankSteps(graph)
	comboSawArr = mergeSAW(graph)

	print("weights", weights)

	return {"graph": graph, "comboSawArr": comboSawArr}