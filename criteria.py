from rdkit import DataStructs
from rdkit import Chem
import pandas as pd
from standalone_model_numpy import SCScorer
import os
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import re
import math
from datetime import datetime
import time
import numpy as np
from concurrent.futures import ThreadPoolExecutor,  ProcessPoolExecutor



# 影响面积
def impact_area(subtree):
    if 'children' not in subtree:
        return 1
    else:
        sub_impact = 0
        for i in subtree['children']:
            sub_impact = sub_impact + impact_area(i)
        return sub_impact + 1

# 输入smiles
def similarity(mole_A, mole_B):
    mole_A = Chem.MolFromSmiles(mole_A)
    mole_B = Chem.MolFromSmiles(mole_B)
    fps1 = Chem.RDKFingerprint(mole_A)
    fps2 = Chem.RDKFingerprint(mole_B)
    return DataStructs.FingerprintSimilarity(fps1, fps2)

# https://deepchemdata.s3.us-west-1.amazonaws.com/datasets/USPTO_FULL.csv
# 传进smiles

# def reactant_confidence2(reactantStepArr):
# 	print("len of reactantConfidenceArr: ", len(reactantStepArr))
# 	similar_molecule = []
# 	article_count = []
# 	count = []

# 	for reactant in reactantStepArr:
# 		similar_molecule.append([])
# 		article_count.append([])
# 		count.append(0)

# 	data = pd.read_csv('USPTO_50K.csv')
# 	reactions = data['reactions']

# 	# reactionCount = 0
# 	# start = time.time()
# 	for i in reactions:

# 		# if(reactionCount % 1000 == 0):
# 		# 	print(reactionCount)
# 		# 	print("time: ", time.time() - start)
# 		# reactionCount += 1
# 		first = i.find('>')+1
# 		second = i[first:].find('>')+1
# 		product = i[first:][second:]

# 		reactantCount = 0
# 		for reactant in reactantStepArr:
# 			if(similarity(reactant["smiles"], product) >= 0.8):
				
# 				count[reactantCount] += 1
# 				if(product in similar_molecule[reactantCount]):
# 					article_count[reactantCount][similar_molecule[reactantCount].index(product)] \
# 						= article_count[reactantCount][similar_molecule[reactantCount].index(product)] + 1
# 				else:
# 					similar_molecule[reactantCount].append(product)
# 					article_count[reactantCount].append(1)
			
# 			reactantCount += 1

	
# 	confidence = []
# 	for reactant in reactantStepArr:
# 		confidence.append(0)

# 	reactantCount = 0
# 	for reactant in reactantStepArr:
# 		for i in range(len(similar_molecule[reactantCount])):
# 			confidence[reactantCount] += \
# 				similarity(reactant["smiles"], similar_molecule[reactantCount][i]) * article_count[reactantCount][i]
			
# 		reactantCount += 1

# 	reactantCount = 0
# 	for reactant in reactantStepArr:
# 		if(count[reactantCount] != 0):
# 			reactant["reactantConfidence"] = confidence[reactantCount] / count[reactantCount]
# 		else:
# 			reactant["reactantConfidence"] = 0

# 		reactantCount += 1


def reactant_confidence_helper(reactant, reactions):
    similar_molecule = []
    article_count = []
    count = 0
    print(1)
    # data = pd.read_csv('USPTO_FULL.csv')
    # data = pd.read_csv('USPTO_50K.csv')
    # reactions = data['reactions']

    for i in reactions:
        
        first = i.find('>')+1
        second = i[first:].find('>')+1
        product = i[first:][second:]
        if similarity(reactant,product)>=0.8:
            print(3)
            count = count + 1
            if product in similar_molecule:
                article_count[similar_molecule.index(product)] = article_count[similar_molecule.index(product)]+1
            else:
                similar_molecule.append(product)
                article_count.append(1)
            
	
    confidence = 0
    for i in range(len(similar_molecule)):
        # print(4)
        confidence = confidence + similarity(reactant, similar_molecule[i])*article_count[i]

    # if(count != 0):
    #     return confidence/count
    # else:
    #     return 0
    return {"similar_molecule": similar_molecule, "article_count": article_count, "count": count}



def reactant_confidence_multithreading(reactant):
	data = pd.read_csv("USPTO_50K.csv")
	reactionsData = data["reactions"]
	reactionsArr = np.array_split(reactionsData, 4)

	taskArr = []
	with ProcessPoolExecutor(max_workers=4) as executor:

		for reactions in reactionsArr:
			task = executor.submit(reactant_confidence_helper, reactant, reactions)
			taskArr.append(task)

	count = 0
	similar_molecule = []
	article_count = []
	for task in taskArr:
		result = task.result()
		count += result["count"]
		moleculeIndex = 0
		for molecule in result["similar_molecule"]:
			if(molecule in similar_molecule):
				article_count[similar_molecule.index(molecule)] += result["article_count"][moleculeIndex]
			else:
				similar_molecule.append(molecule)
				article_count.append(result["article_count"][moleculeIndex])			
			
			moleculeIndex += 1

	
	confidence = 0
	for i in range(len(similar_molecule)):
		confidence = confidence + similarity(reactant, similar_molecule[i])*article_count[i]

	if(count != 0):
		return confidence/count
	else:
		return 0




def reactant_confidence(reactant):
    similar_molecule = []
    article_count = []
    count = 0

    # data = pd.read_csv('USPTO_FULL.csv')
    data = pd.read_csv('USPTO_50K.csv')
    reactions = data['reactions']

    for i in reactions:
        
        first = i.find('>')+1
        second = i[first:].find('>')+1
        product = i[first:][second:]
        if similarity(reactant,product)>=0.8:
            count = count + 1
            if product in similar_molecule:
                article_count[similar_molecule.index(product)] = article_count[similar_molecule.index(product)]+1
            else:
                similar_molecule.append(product)
                article_count.append(1)
            
	
    confidence = 0
    for i in range(len(similar_molecule)):
        # print(4)
        confidence = confidence + similarity(reactant, similar_molecule[i])*article_count[i]

    if(count != 0):
        return confidence/count
    else:
        return 0

def atom_number(molecule):
    mol = Chem.MolFromSmiles(molecule)
    formula = CalcMolFormula(mol)
    pattern = re.compile(r'\d+')  # 查找数字
    result = pattern.findall(formula)
    num = 0
    for i in result:
        num = num + int(i)

    return num

# 反应的confidence (传入所有reactant的list和产物)
def confidence(reactant, product):
    all_atom_num = 0
    # 这里能不能只考虑骨架？
    for i in reactant:
        all_atom_num = all_atom_num + atom_number(i)

    confidence = 0
    for i in reactant:
        confidence = confidence + atom_number(i)/all_atom_num*reactant_confidence_multithreading(i)

    return confidence

# 反应的confidence (传入所有的reactant的list和产物)
def complexity_reduction(reactant, product):
    project_root = ""
    model = SCScorer()
    model.restore(
        os.path.join(project_root, 'models', 'full_reaxys_model_1024bool', 'model.ckpt-10654.as_numpy.json.gz'))
    score_reactant = 1
    for i in reactant:
        SC = model.get_score_from_smi(i)
        SC = SC[1]
        simi = 1 - (SC-1)/4
        score_reactant = score_reactant * simi

    SC_product = model.get_score_from_smi(product)
    SC_product = SC_product[1]
    score_product = 1 - (SC_product-1)/4

    return score_reactant/score_product

def convergence(reactant, product):
    sum = 0
    for i in reactant:
        sum = sum + abs(atom_number(product)/len(reactant) - atom_number(i))
    mae = sum/len(reactant)

    return 1/(1+mae)
