from flask import Flask, render_template
from flask import request
from flask import jsonify
from flask_cors import CORS
import time
import json
from checkRXN_func import checkRXN_func
from initialize_func import initialize_func
from revise_func import revise_func
from reconfigureConstraints_func import reconfigureConstraints_func
from checklRXN_data import checkRXN_data
import globalVar



def seperateHumanAISteps(graph):

	if(("handledByAI" in graph and graph["handledByAI"]
		and "AIRoutes" in graph and len(graph["AIRoutes"]) != 0)
		or ("isAI" in graph and graph["isAI"])):

		# graph["isAI"] = True

		if("children" in graph):
			for child in graph["children"]:
				child["isAI"] = True

	if("children" in graph):
		for child in graph["children"]:
			seperateHumanAISteps(child)


app = Flask(__name__)
CORS(app, supports_credentials=True)


@app.route('/', methods=["POST", "GET"])
def hello_world():
	# print(request, file=sys.stdout)
	
	print(request.data)

	return jsonify(username="haha", email="hehe")
	

@app.route("/initialize", methods=["POST"])
def initialize():

	inputMolecule = request.get_json(force=True, silent=False)

	globalVar.runningTime = 540
	graph = initialize_func(inputMolecule["smiles"], inputMolecule["constraints"])
	seperateHumanAISteps(graph)

	return jsonify(graph)


@app.route("/checkRXN", methods=["POST"])
def checkRXN():

	data = request.get_json(force=True, silent=False)

	# with open("output3.json", "w") as outputFile:
	# 	outputFile.write(json.dumps(graph))

	# with open("output.json", "r") as inputFile:
	# 	AIRoutes = json.load(inputFile)
	# 	graph["AIRoutes"] = AIRoutes
	# 	graph["children"] = AIRoutes[0]

	globalVar.runningTime = 540
	# graph = checkRXN_func(data["graph"], data["constraints"])
	graph = checkRXN_data(data["graph"], data["constraints"])

	seperateHumanAISteps(graph)

	# return jsonify(checkRXN_func(data["graph"], data["constraints"]))
	
	return jsonify(graph)


@app.route("/revise", methods=["POST"])
def revise():

	data = request.get_json(force=True, silent=False)
	# with open("reviseInput.json", "w") as outputFile:
	# 	outputFile.write(json.dumps(graph))

	# time.sleep(10)

	return jsonify(revise_func(data["graph"], data["weights"]))


@app.route("/reconfigureConstraints", methods=["POST"])
def reconfigureConstraints():

	data = request.get_json(force=True, silent=False)

	globalVar.runningTime = 540
	graph = reconfigureConstraints_func(data["graph"], data["constraints"])
	seperateHumanAISteps(graph)
	
	return jsonify(graph)


if __name__ == "__main__":
	# app.run(host="0.0.0.0", ssl_context="adhoc")
	app.run(host="0.0.0.0", port=5000)

