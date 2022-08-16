from __future__ import annotations
from copyreg import pickle
import json
from operator import itemgetter
from typing import Dict, List
from numpy import product
from rxn4chemistry import RXN4ChemistryWrapper
from datetime import datetime
from requests.models import Response
import time
from datetime import datetime
import traceback
import globalVar
from rdkit import Chem



class RXN:

	rxn4chemistry_wrapper: RXN4ChemistryWrapper = None
	max_steps: int = 7
	ai_model: str = "12class-tokens-2021-05-14"
	fap: int = 0.65
	pruning_steps: int = 1
	nbeams: int = 30
	excludeSmiles = None
	excludeSubstructure = None
	price = 0

	@classmethod
	def initialize(cls)-> None:
		api_key = 'apk-331eb544c07538c4220d97524337f4466abf79a8410db0928c4476769753647f37b0009da2fa5a711e3960ce6bdf99e479124000f4dcf201870b687aba64788b8b82a59ce241270f8e5a19b39cf467e2'
		# api_key = "apk-0f058b2cd32807d420ae9bd9afd107274fb517646f6cbf9c532d91545ec7c5bba1194b20abfc8d7a654b382e5215ac1a706db21d02d6a8a8582578cad4b2c59546023c737f6e8f27ea25190eb435e08a"
		cls.rxn4chemistry_wrapper = RXN4ChemistryWrapper(api_key=api_key)
		current_time = datetime.now().strftime("%H:%M:%S")
		cls.rxn4chemistry_wrapper.create_project(current_time)


	@classmethod
	def predict(cls, smiles: str)-> Response:
		response: Response = cls.rxn4chemistry_wrapper.predict_automatic_retrosynthesis(
			product=smiles,
			max_steps=cls.max_steps,
			fap=cls.fap,
			nbeams=cls.nbeams,
			pruning_steps=cls.pruning_steps,
			ai_model=cls.ai_model,
			exclude_smiles=cls.excludeSmiles,
			exclude_substructures=cls.excludeSubstructure,
			availability_pricing_threshold = cls.price
		)

		while("response" in response and "error" in response["response"]):
			now = datetime.now()
			print(now.strftime("%H:%M:%S"))
			print(smiles)
			print("response error")
			print(response)
			time.sleep(61)
			response = cls.rxn4chemistry_wrapper.predict_automatic_retrosynthesis(
				product=smiles,
				max_steps=cls.max_steps,
				fap=cls.fap,
				nbeams=cls.nbeams,
				pruning_steps=cls.pruning_steps,
				ai_model=cls.ai_model,
				exclude_smiles=cls.excludeSmiles,
				exclude_substructures=cls.excludeSubstructure,
				availability_pricing_threshold = cls.price
			)   
			now = datetime.now()
			print(now.strftime("%H:%M:%S"))
			print("new response")
			print(smiles)
			print(response)
			print()

		results = None
		runningTime = 0
		while(results == None or ("status" in results and results["status"] != "SUCCESS") or
			"response" in results and "error" in results["response"]):


			if(runningTime >= globalVar.runningTime):
				raise KeyError("time")

			# print("response", response)

			time.sleep(60)

			try:
				results = cls.rxn4chemistry_wrapper.\
								get_predict_automatic_retrosynthesis_results(response["prediction_id"])
			except:
				print("rxn predict exception")
				traceback.print_exc()


			# with open("output3.json", "w") as outputFile:
			# 	outputFile.write(json.dumps(results))

			# PROCESSING\
			now = datetime.now()
			print(now.strftime("%H:%M:%S"))
			print(smiles)
			if("status" in results):
				print(results["status"])
				print()
			else:
				print("result error")
				print(results)
				print()

			errorCount = 0
			if("status" in results):
				if(results["status"] == "PROCESSING"):
					runningTime += 60
			else:
				errorCount += 1
				time.sleep(2)

			if(errorCount == 3):
				raise KeyError("status not found")

		with open("initialize.json", "w") as outputFile:
			outputFile.write(json.dumps(results))

		return results


	@classmethod
	def processAIRoute(cls, AIRoute):
		failureCauseArr = ["Cannot find a route within the given price threshold for molecules", 
		"Cannot find a route within the given maximum retrosynthetic steps",
		"Cannot find a route when excluding the molecule",
		"Cannot find a route when excluding the substructure"]

		if(len(AIRoute["children"]) == 0):
			if(AIRoute["isExpandable"] == False 
				and AIRoute["isCommercial"] == False):
				AIRoute["isAvailable"] = False
				AIRoute["failureCause"] = ""
			else:
				AIRoute["isAvailable"] = True

			if(AIRoute["isExpandable"] ==  True
				and AIRoute["isCommerical"] == False ):
				AIRoute["failureCause"] = failureCauseArr[1]
		
		for child in AIRoute["children"]:
			RXN.processAIRoute(child)


	@classmethod
	def processResult(cls, resultJson):

		pathArr = []

		if("retrosynthetic_paths" in resultJson):
			for path in resultJson["retrosynthetic_paths"]:
				if(path["confidence"] != 0):
					pathArr.append(path)

		for path in pathArr:
			RXN.processAIRoute(path)

		pathArr = sorted(pathArr, key=itemgetter("confidence"), reverse=True)

		return pathArr


	# O1[C@@H]2OC3C=C(C4C5CCC(C=5C(=O)OC=4C=3[C@@H]2C=C1)=O)OC
	@classmethod
	def predictResult(cls, smiles, constraints):

		inputMol = Chem.MolFromSmiles(smiles)
		inputSmiles = Chem.MolToSmiles(inputMol)

		with open("failureSmiles.txt", "r") as failureFile:
			for failureSmiles in failureFile:
				if(failureSmiles != ""):
					failureMol = Chem.MolFromSmiles(failureSmiles)
					failureSmiles = Chem.MolToSmiles(failureMol)
					if(failureSmiles == inputSmiles):
						print("failure Smiles: ", smiles)
						return []

		
		if(constraints["mssr"] == None):
			cls.max_steps = 7
		else:
			cls.max_steps = constraints["mssr"]
		cls.price = constraints["price"]
		# cls.price = 0
		cls.excludeSmiles = constraints["excludeSmiles"]
		cls.excludeSubstructure = constraints["excludeSubstructure"]


		# print(cls.max_steps)
		# print(cls.price)

		pathArr = []
		now = datetime.now()
		print(now.strftime("%H:%M:%S"))
		print("initiailize", smiles)
		print()

		try: 
			resultDict = RXN.predict(smiles)
			pathArr = RXN.processResult(resultDict)
		except Exception as e:
			now = datetime.now()
			print(now.strftime("%H:%M:%S"))
			print("RXN exception")
			print(smiles)
			print(e)


		# with open("temp.json", "w") as outputFile:
		# 	json.dump(pathArr, outputFile)

		return pathArr


	
RXN.initialize()











	


