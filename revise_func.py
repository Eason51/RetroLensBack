rank = 1

def revise_func_helper(graph):

	if("children" not in graph or len(graph["children"]) == 0):
		return

	global rank

	if(rank >= 10):
		return
		
	graph["rank"] = rank
	rank  = rank + 1
	graph["SAW"] = 10 - rank
	graph["influence"] = 0.2
	graph["reactionConfidence"] = 0.3
	graph["complexity"] = 0.4
	graph["convergence"] = 0.5
	graph["associatedSubtreeConfidence"] = 0.6


	if("children" in graph):
		for child in graph["children"]:
			revise_func_helper(child)


	# if("notAvailable" in graph):
	# 	graph["notAvailable"] = False
	

	



def revise_func(graph):

	global rank
	rank = 1

	revise_func_helper(graph)

	return graph