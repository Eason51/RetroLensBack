from rdkit import Chem, rdBase
from rdkit.Chem import Draw,rdDepictor
from rdkit.Chem.Draw import IPythonConsole, rdMolDraw2D
from rdkit.Chem.Scaffolds import rdScaffoldNetwork
import math
import os
import pyvis
from pyvis.network import Network
import inspect

def moltosvg(mol,molSize=(450,250),kekulize=True):
    mc=rdMolDraw2D.PrepareMolForDrawing(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    opts = drawer.drawOptions()
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg.replace('svg:','')



def extractScaffold(smiles):

	params = rdScaffoldNetwork.ScaffoldNetworkParams()
	params.includeScaffoldsWithoutAttachments=False

	molecule = Chem.MolFromSmiles(smiles)

	net1 = rdScaffoldNetwork.CreateScaffoldNetwork([molecule],params)

	nodemols1 = [Chem.MolFromSmiles(x) for x in net1.nodes]
	
	# print("scaffold")
	# for mol in nodemols1:
	# 	print(Chem.MolToSmiles(mol))
	# print()

	if(len(nodemols1) >= 2):
		return Chem.MolToSmiles(nodemols1[1])
	elif(len(nodemols1) == 1):
		return Chem.MolToSmiles(nodemols1[0])
	else:
		return ""

