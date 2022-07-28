from __future__ import annotations
from copyreg import pickle
import json
from typing import Dict, List
from numpy import product
from rxn4chemistry import RXN4ChemistryWrapper
from datetime import datetime
from requests.models import Response

import time


class Molecule:

    def __init__(self, smiles: str)-> None:
        self.smiles: str = smiles


class UserMolecule(Molecule):
    
    def __init__(self, smiles: str)-> None:
        super().__init__(smiles)

        self.initScaffold()

    def initScaffold(self) -> None:
        pass


class AIMolecule(Molecule):

    def __init__(self, smiles: str) -> None:
        super().__init__(smiles)


class RetroTreeNode:

    def __init__(self, molecule: Molecule, parent: RetroTreeNode)-> None:
        self.molecule: Molecule = molecule
        self.parent: RetroTreeNode = parent
        self.childrenArr: List[RetroTreeNode] = []
        self.errorCase = False

    def checkError(self):
        if(len(self.childrenArr) > 0):
            for child in self.childrenArr:
                if(child.checkError):
                    self.errorCase = True

        return self.errorCase



class UserRetroTreeNode(RetroTreeNode):

    def __init__(self, molecule: Molecule, parent: RetroTreeNode) -> None:
        super().__init__(molecule, parent)


class AIRetroTreeNode(RetroTreeNode):

    def __init__(self, 
        molecule: Molecule, 
        parent: RetroTreeNode,
        confidence: float,
        isExpandable: bool,
        isCommercial: bool
    ) -> None:
        super().__init__(molecule, parent)
        self.confidence = confidence
        self.isExpandable = isExpandable
        self.isCommercial = isCommercial 

    
    def checkError(self):
        if(len(self.childrenArr) == 0):
            if(self.isCommercial == False):
                self.errorCase = True
                return self.errorCase
        else:
            return super().checkError()





class RetroTreeRoute:

    def __init__(self, root: RetroTreeNode)-> None:
        self.root: RetroTreeNode = root
        self.routeIndexArr: List[int] = [0]

        self.initRetroTreeRoute()
    
    def initRetroTreeRoute(self)-> None:
        pass


class RXN:

    rxn4chemistry_wrapper: RXN4ChemistryWrapper = None
    max_steps: int = 100
    ai_model: str = "12class-tokens-2021-05-14"
    fap: int = 0.6
    pruning_steps: int = 0
    nbeams: int = 0

    @classmethod
    def initialize(cls)-> None:
        api_key = 'apk-331eb544c07538c4220d97524337f4466abf79a8410db0928c4476769753647f37b0009da2fa5a711e3960ce6bdf99e479124000f4dcf201870b687aba64788b8b82a59ce241270f8e5a19b39cf467e2'
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
            ai_model=cls.ai_model
        )

        while("response" in response and "error" in response["response"]):
            print("response error")
            time.sleep(61)
            response = cls.rxn4chemistry_wrapper.predict_automatic_retrosynthesis(
                product=smiles,
                max_steps=cls.max_steps,
                fap=cls.fap,
                nbeams=cls.nbeams,
                pruning_steps=cls.pruning_steps,
                ai_model=cls.ai_model
            )   

        results = None
        while(results == None or results["status"] != "SUCCESS"):

            time.sleep(15)
            results = cls.rxn4chemistry_wrapper.\
                            get_predict_automatic_retrosynthesis_results(response["prediction_id"])

        return results


    @classmethod
    def generateRetroTree(cls, parent: RetroTreeNode, routeDict: Dict)-> RetroTreeNode:
        currentMolecule: AIMolecule = AIMolecule(routeDict["smiles"])
        currentTreeNode: AIRetroTreeNode = AIRetroTreeNode(
            currentMolecule, 
            parent, 
            routeDict["confidence"],
            routeDict["isExpandable"],
            routeDict["isCommercial"]
            )

        if(parent != None):
            parent.childrenArr.append(currentTreeNode)

        if(len(routeDict["children"]) != 0):
            for childrenRouteDict in routeDict["children"]:
                cls.generateRetroTree(currentTreeNode, childrenRouteDict)

        return currentTreeNode




def main():

    """Input A molecule to RXN"""
    RXN.initialize()
    smiles: str = "C1C2(O)C3C4OC(O)C5C(O)CCC5C4CC4(CCCC(=C2C=C1)C34)OC"
    results: Response = RXN.predict(smiles)
    
    with open("result4.json", "w") as outputFile:
        json.dump(results, outputFile)


if __name__ == "__main__":
    main()
