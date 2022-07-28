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



app = Flask(__name__)
CORS(app, supports_credentials=True)


@app.route('/', methods=["POST", "GET"])
def hello_world():
	# print(request, file=sys.stdout)
	
	print(request.data)

	return jsonify(username="haha", email="hehe")
	


@app.route("/checkRXN", methods=["POST"])
def checkRXN():

	graph = request.get_json(force=True, silent=False)
	# with open("output3.json", "w") as outputFile:
	# 	outputFile.write(json.dumps(graph))

	# with open("output.json", "r") as inputFile:
	# 	AIRoutes = json.load(inputFile)
	# 	graph["AIRoutes"] = AIRoutes
	# 	graph["children"] = AIRoutes[0]

	return jsonify(checkRXN_func(graph))


@app.route("/initialize", methods=["POST"])
def initialize():

	inputMolecule = request.get_json(force=True, silent=False)
	
	return initialize_func(inputMolecule["smiles"])


@app.route("/revise", methods=["POST"])
def revise():

	graph = request.get_json(force=True, silent=False)
	with open("reviseInput.json", "w") as outputFile:
		outputFile.write(json.dumps(graph))

	# time.sleep(10)

	return jsonify(revise_func(graph))


@app.route("/reconfigureConstraints", methods=["POST"])
def reconfigureConstraints():

	graph = request.get_json(force=True, silent=False)
	
	return jsonify(reconfigureConstraints_func(graph))


@app.route("/ketcher", methods=["GET"])
def getKetcher():
	
	return render_template("ketcher.html")


if __name__ == "__main__":
	# app.run(host="0.0.0.0", ssl_context="adhoc")
	app.run(host="0.0.0.0")

