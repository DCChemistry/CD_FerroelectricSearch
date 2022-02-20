import matplotlib.pyplot as plt
from Util import *
from pymatgen.core.composition import Composition
#Smidt imports
from pymatgen.symmetry.groups import SYMM_DATA, sg_symbol_from_int_number

import numpy as np #added just in case - may add numpy objects later
import matplotlib.pyplot as plt
import os
import re
from pymatgen.ext.matproj import MPRester #used for NoPolarVar()
import numpy as np
import pymatgen.analysis.local_env as localEnv
import pandas as pd

class Analysis:

    def __init__(self, searchFileName, orderOfFilters, elementList):
        #orderOfFilters is the order of the keys from 'filters' dictionary
        self.searchFileName = searchFileName
        self.elementList = elementList

        #keys are the 'codes' appended to a certain file
        filters = {
                    "NP": Analysis.NonPolar,
                    "oneCDSite": Analysis.OneCDSite,
                    "BAndT": Analysis.BinAndTern,
                    "lteq30sites": Analysis.LTorEQ30Sites,
                    "noTox": Analysis.NoToxicElements,
                    "onlyOxy": Analysis.KeepOnlyOxyAnion,
                    "chosenCDElem": self.OnlyChosenCDElements,
                    "specOS": self.RedSearchSpecificOS,
                    "ME": self.MutuallyExclusiveElements,
                    "NoPolarVar": Analysis.NoPolarVar,
                    "SiteEquiv": Analysis.SiteEquivalence
        }

        for counter, filter in enumerate(orderOfFilters):
            if(counter==0):
                self.ReadAnalyseWrite(filters[filter], "CD", filter, counter)
            else:
                self.ReadAnalyseWrite(filters[filter], orderOfFilters[counter-1], filter, counter)


    def ReadAnalyseWrite(self, analysisType, prevAnalysisTag, newAnalysisTag, numberInQueue): #numberInQueue is to show the order each filter was applied in
        """AnalysisType is the name of the method used to analyse the data, e.g. NonPolar.
           prevAnalysisTag is the text appeneded to the end of the analysis file you want to load.
           newAnalysisTag that will be appended to the end of the analysis file you want to create."""
        
        if(not os.path.isfile(f"{self.searchFileName}{newAnalysisTag}_{numberInQueue+1}.json")):
            print(f"\nStarting {newAnalysisTag} analysis:")
            results = ReadJSONFile(f"{self.searchFileName}{prevAnalysisTag}_{numberInQueue}")
            analysisResults = analysisType(results)
            SaveDictAsJSON(f"{self.searchFileName}{newAnalysisTag}_{numberInQueue+1}", analysisResults)
            print(f"{newAnalysisTag} analysis complete.")
            # ^ numberInQueue+1 starts from 1, hence numberInQueue without the +1 is the previous numberInQueue
            numOfMatInPrevAnal = len(list(results.keys()))
            numOfMatInCurrentAnal = len(list(analysisResults.keys()))
            print(f"{numOfMatInCurrentAnal} materials identified.")
            print(f"{numOfMatInPrevAnal-numOfMatInCurrentAnal} materials removed from previous analysis ({prevAnalysisTag}).")
        else:
            print(f"{newAnalysisTag} analysis has already been done for search {self.searchFileName}.")


    @staticmethod
    def NonPolar(results): #add in stuff to check if NP file already exists
        print(f"Starting non polar analysis.")
        #taken from blondegeek (Tess Smidt) ferroelectric_search_site repo:
        # https://github.com/blondegeek/ferroelectric_search_site/blob/82a2b939922bfd7d5f7e40ed3616dab9f0065a83/example_notebooks/nonpolar-polar_structure_pairs_for_SIMPLE_group_subgroup_relations_not_isotropy.ipynb
        
        # This is a list of the point groups as noted in pymatgen
        point_groups = []
        for i in range(1,231):
            symbol = sg_symbol_from_int_number(i)
            point_groups.append(SYMM_DATA['space_group_encoding'][symbol]['point_group'])

        # Note that there are 40 of them, rather than 32.

        # This is because multiple conventions are used for the same point group.
        # This dictionary can be used to convert between them.
        point_group_conv = {'321' :'32', '312': '32', '3m1' :'3m', '31m': '3m',
                            '-3m1' : '-3m', '-31m': '-3m', '-4m2': '-42m', '-62m': '-6m2' }

        # Using this dictionary we can correct to the standard point group notation.
        corrected_point_groups = [point_group_conv.get(pg, pg) for pg in point_groups]
        # Which produces the correct number of point groups. 32.


        # polar_point_groups = ['1', '2', 'm', 'mm2', '4', '4mm', '3', '3m', '3m1', '31m','6', '6mm']
        # There are 10 polar point groups
        polar_point_groups = ['1', '2', 'm', 'mm2', '4', '4mm', '3', '3m', '6', '6mm']

        # Polar spacegroups have polar point groups.
        polar_spacegroups = []
        # There are 230 spacegroups
        for i in range(1,231):
            symbol = sg_symbol_from_int_number(i)
            pg = SYMM_DATA['space_group_encoding'][symbol]['point_group']
            if point_group_conv.get(pg, pg) in polar_point_groups:
                polar_spacegroups.append(i)
        # 68 of the 230 spacegroups are polar.


        nonPolarResults = {}
        for material in results:
            if(results[material]["spacegroup.number"] not in polar_spacegroups):
                nonPolarResults[material] = results[material]
        
        print(f"NP analysis complete. {len(nonPolarResults.keys())} materials found.")
        return nonPolarResults #consider incorporating file reading and writing directly into this function, rather than doing it at initialisation
    
    @staticmethod
    def DisplayElemAndNo(formula):
        """Takes a formula, and returns a dictionary of the elements present as keys, with respective values being the no. of atoms of that element."""
        comp = Composition(formula)
        elementsFromFormulaWithNo = comp.formula #breaks formula into a form that is digestible by the program

        #elements in material
        alphabetRegex = re.compile('[a-zA-Z]+') #need this for removing numbers from string
        elementsFromFormula = alphabetRegex.findall(elementsFromFormulaWithNo) #returns list of elements present

        #number of each element in material
        numberRegex = re.compile('[0-9]+')
        noOfElemFromFormula = numberRegex.findall(elementsFromFormulaWithNo)
        noOfElemFromFormula = [int(num) for num in noOfElemFromFormula] #entries are strings otherwise - typecasting to int necessary.

        materialElemAndNo = dict(zip(elementsFromFormula, noOfElemFromFormula))

        return materialElemAndNo

    @staticmethod
    def OneCDSite(results):
        print("Starting one CD site filter.")
        primCellResults = {}
        for material in results:
            formula = results[material]["pretty_formula"]
            comp = Composition(formula)
            reducedFormulaAtomNo = comp.num_atoms
            noOfAtomsInStruct = results[material]["nsites"]
            noOfPrimCells = noOfAtomsInStruct/reducedFormulaAtomNo
            CDelement = results[material]["CDelement"]
            materialElemAndNo = Analysis.DisplayElemAndNo(formula)
            noOfCDSites = materialElemAndNo[CDelement]
            print(noOfCDSites)
            if(noOfPrimCells == 1 and noOfCDSites == 1):
                primCellResults[material] = results[material]
            elif(noOfAtomsInStruct%reducedFormulaAtomNo!=0):
                #^ Identifies materials that have a fractional no. of primitive cell in the calculated cell.
                # There are very few of these, and the ones found in my run were H2 and N2 only. All other materials have an int no. of prim. cells.
                print(f"Fractional no. of primitive cells found for material: {results[material]['pretty_formula']}.\nPrim. cell ratio = {noOfPrimCells}\n")
        
        print(f"Primitive cell analysis complete. {len(primCellResults.keys())} materials found.")
        return primCellResults
    
    @staticmethod
    def BinAndTern(results):
        print("Starting binary and ternary compound filtering.")
        binAndTernResults = {}
        for material in results:
            noOfElements = results[material]["nelements"]
            if(noOfElements in [2,3]): #this may look scuffed, but it works
                binAndTernResults[material] = results[material]
        
        print(f"Binary and ternary compounds found. {len(binAndTernResults.keys())} materials identified.")
        return binAndTernResults

    @staticmethod
    def LTorEQ30Sites(results):
        print("Starting less than or equal to 30 sites (atoms in the cell) filtering.")
        ltOrEq30Results = {}
        for material in results:
            noOfSites = results[material]["nsites"]
            if(noOfSites <= 30):
                ltOrEq30Results[material] = results[material]
        
        print(f"Site filtering (<= 30) filtering complete. {len(ltOrEq30Results.keys())} materials found.")
        return ltOrEq30Results

    @staticmethod
    def DisplayElements(formula):
        """Returns a list of elements present from a given formula."""
        comp = Composition(formula)
        elementsFromFormulaWithNo = comp.formula
        alphabetRegex = re.compile('[a-zA-Z]+') #need this for removing numbers from string
        elementsFromFormula = alphabetRegex.findall(elementsFromFormulaWithNo) #returns list of elements present
        return elementsFromFormula

    @staticmethod
    def NoToxicElements(results): #removing Cd and Pb - they're toxic, might as well.
        toxicElements = ["Pb", "Cd", "As"]
        print(f"Starting removal of compounds that contain: {toxicElements}.")
        noToxicElemResults = {}
        for material in results:
            formula = results[material]["pretty_formula"]
            listOfElemPresent = Analysis.DisplayElements(formula)
            if not any(toxicElem in listOfElemPresent for toxicElem in toxicElements):
                #^ this just works - if I had to explain it, if there is no overlap between the two lists, continue, otherwise do nothing
                # ^ using this solution https://stackoverflow.com/questions/62115746/can-i-check-if-a-list-contains-any-item-from-another-list
                noToxicElemResults[material] = results[material]

        print(f"Toxic elements {toxicElements} removed. {len(noToxicElemResults.keys())} materials remain.")
        return noToxicElemResults


    def MutuallyExclusiveElements(self, results):
        mutuallyExclusiveResults = {}
        for material in results:
            formula = results[material]["pretty_formula"]
            listOfElemPresent = Analysis.DisplayElements(formula)
            counter = 0 #the counter should never exceed 1 (not mutually exclusive if >1)
            for elem in self.elementList:
                if(elem in listOfElemPresent):
                    counter += 1 #one of the mutually exlusive elements is present
                    if(counter > 1): #a second element from the mutually exclusive list is in this material - don't want it, move on
                        break
            if(counter == 1): #only if one of the materials from self.elementList is present can the material move on to the next stage
                mutuallyExclusiveResults[material] = results[material]
        
        return mutuallyExclusiveResults

    def OnlyChosenCDElements(self, results):
        chosenCDElemResults = {}
        for material in results:
            CDElem = results[material]["CDelement"]
            if(CDElem in self.elementList):
                chosenCDElemResults[material] = results[material]
        return chosenCDElemResults

    #specific filters, no longer general - used for ReducedSearch1 and ReducedSearch2
    @staticmethod
    def KeepOnlyOxyAnion(results):
        oxyResults = {}
        for material in results:
            oxStates = results[material]["OxStates"]
            #need to then identify the oxidation states that contain a zero, then remove the material
            nonOxyAnions = 0
            for oxState in oxStates:
                oxyRegex = re.compile(f'[O-]')
                checkingForOxy = oxyRegex.findall(oxState)
                if("-" in checkingForOxy and "O" not in checkingForOxy):
                    nonOxyAnions += 1
            if(nonOxyAnions == 0): #if only oxy (O^z- anions present), save material
                oxyResults[material] = results[material]
        return oxyResults

    def RedSearchSpecificOS(self, results): #requires OnlyChosenCDElements to be applied first
        specificOSResults = {}
        specOSRejects = {}
        specificOxStates = {"Sn": ["Sn2+", "Sn4+"],
                            "Pb": ["Pb2+", "Pb4+"],
                            "Sb": ["Sb3+", "Sb5+"],
                            "Bi": ["Bi3+", "Bi5+"]}
        for material in results:
            CDElem = results[material]["CDelement"]
            oxStates = results[material]["OxStates"]
            alphabetRegex = re.compile('[a-zA-Z]+')
            oxStatesAlphaRegex = [alphabetRegex.findall(element)[0] for element in oxStates]
            CDinstances = oxStatesAlphaRegex.count(CDElem)

            expectedOxStates = specificOxStates[CDElem]
            if(all(elem in oxStates for elem in expectedOxStates) and len(specificOxStates[CDElem])==CDinstances):
                specificOSResults[material] = results[material]
            else:
                specOSRejects[material] = results[material]
        SaveDictAsJSON(f"{self.searchFileName}specOSRejects", specOSRejects) #saving a file of the rejects
        return specificOSResults

    @staticmethod
    def NoPolarVar(results): #ONLY USE WHEN THERE AREN'T MANY MATERIALS LEFT - USE AS LAST FILTER
        #Code to get list of polar space groups taken from BlondeGeek (Smidt) on GitHub; also used in NonPolar() above:

        # This is a list of the point groups as noted in pymatgen
        point_groups = []
        for i in range(1,231):
            symbol = sg_symbol_from_int_number(i)
            point_groups.append(SYMM_DATA['space_group_encoding'][symbol]['point_group'])

        # Note that there are 40 of them, rather than 32.

        # This is because multiple conventions are used for the same point group.
        # This dictionary can be used to convert between them.
        point_group_conv = {'321' :'32', '312': '32', '3m1' :'3m', '31m': '3m',
                            '-3m1' : '-3m', '-31m': '-3m', '-4m2': '-42m', '-62m': '-6m2' }

        # Using this dictionary we can correct to the standard point group notation.
        corrected_point_groups = [point_group_conv.get(pg, pg) for pg in point_groups]
        # Which produces the correct number of point groups. 32.


        # polar_point_groups = ['1', '2', 'm', 'mm2', '4', '4mm', '3', '3m', '3m1', '31m','6', '6mm']
        # There are 10 polar point groups
        polar_point_groups = ['1', '2', 'm', 'mm2', '4', '4mm', '3', '3m', '6', '6mm']

        # Polar spacegroups have polar point groups.
        polar_spacegroups = []
        # There are 230 spacegroups
        for i in range(1,231):
            symbol = sg_symbol_from_int_number(i)
            pg = SYMM_DATA['space_group_encoding'][symbol]['point_group']
            if point_group_conv.get(pg, pg) in polar_point_groups:
                polar_spacegroups.append(i)
        # 68 of the 230 spacegroups are polar.



        noPolarVarResults = {}
        for material in results:
            formula = results[material]["pretty_formula"]

            #query structure copied from DatabaseSearch.py
            APIkey = None #done so that APIkey is not lost in the scope of the with block
            with open("APIkey.txt", "r") as f:
                APIkey= f.read()

            variants = None #done so that variants exists outside the scope of the with block
            with MPRester(APIkey) as mpr:
    
                criteria = {"pretty_formula": {"$eq": formula}}

                properties = ['material_id', 'spacegroup.number']
                variants = mpr.query(criteria, properties, chunk_size=10000) #looking for all entries of a material
                
            polarCounter = 0
            for i in range(len(variants)):
                spacegroupNumber = variants[i]["spacegroup.number"]
                if(spacegroupNumber in polar_spacegroups):
                    polarCounter += 1
            if(polarCounter == 0): #if a material doesn't have any polar variants, save it
                noPolarVarResults[material] = results[material]            
        
        return noPolarVarResults

    @staticmethod
    def SiteEquivalence(results, tolerance=0.01): #tolerance: checking that bond lengths are within this margin of error
        siteEquivalenceResults = {}
        for material in results:
            
            print()
            print(material)
            
            CDElem = results[material]["CDelement"]

            APIkey = None #done so that APIkey is not lost in the scope of the with block
            with open("APIkey.txt", "r") as f:
                APIkey= f.read()

            structure = None #done so that structure exists outside the scope of the with block
            with MPRester(APIkey) as mpr:
                structure = mpr.get_structure_by_material_id(material)
            neighbourInfoList = []
            noOfAtoms = len(structure.as_dict()["sites"])
            for i in range(noOfAtoms): #for each atom in the material
                elem = structure.as_dict()["sites"][i]["species"][0]["element"]
                if(elem==CDElem): #only looking for the neighbours of CD atoms
                    neighbourInfo = {"atomType": [], "nnDistance": []}
                    noOfNeighbours = len(localEnv.BrunnerNN_real().get_nn_info(structure, i))
                    for j in range(noOfNeighbours): #for each of atom i's neighbours
                        atomType = localEnv.BrunnerNN_real().get_nn_info(structure, i)[j]["site"]._species
                        nnDistance = localEnv.BrunnerNN_real().get_nn_info(structure, i)[j]["site"].nn_distance
                        neighbourInfo["atomType"].append(str(atomType)) #atomType is normally some non-standard object
                        neighbourInfo["nnDistance"].append(nnDistance)
                    neighbourInfoList.append(neighbourInfo)

            neighbourDF = pd.DataFrame(data=neighbourInfoList)

            sameNoOfNeighboursList = []
            for i in range(neighbourDF.shape[0]):
                for j in range(neighbourDF.shape[0]):
                    if(i!=j and i<j): #if lists 1 and 2 have been compared, don't want to compare 2 and 1 (so i<j) and don't want to
                                    #compare a list to itself (i!=j)
                        sameNoOfNeighbours = len(list(neighbourDF["atomType"])[i]) == len(list(neighbourDF["atomType"])[j])
                        sameNoOfNeighboursList.append(sameNoOfNeighbours)
            print(f"Same no. of neighbours?: {sameNoOfNeighboursList}")

            if(False not in sameNoOfNeighboursList): #if the no. of neighbours is the same in each list, False shouldn't be in the list
                neighbourIdentityLists = list(neighbourDF["atomType"])

                neighbourComparisonList = []
                for i in range(neighbourDF.shape[0]):
                    for j in range(neighbourDF.shape[0]):
                        if(i!=j and i<j):
                            neighbourComparison = set(neighbourIdentityLists[i]) == set(neighbourIdentityLists[j])
                            neighbourComparisonList.append(neighbourComparison)
                print(f"Same neighbouring atoms?: {neighbourComparisonList}")


                if(False not in neighbourComparisonList):
                    nnDistanceLists = list(neighbourDF["nnDistance"])
                    nnDistanceLists = [list(np.sort(x)) for x in nnDistanceLists]
                    # ^ lists contain numpy floats - need to use np.sort then convert back to list (don't want numpy array)
                    # ^ want to compare distances one by one - if sites identical, the values should be sorted into the same order
                    #(this avoids problems with distances to neighbours from a different site being counted in a different order)

                    siteComparisonList = []
                    for i in range(neighbourDF.shape[0]):
                        for j in range(neighbourDF.shape[0]):
                            if(i!=j and i<j):
                                distanceComparisonList = []
                                for k in range(len(nnDistanceLists[0])):
                                    # ^ we've already established no. of neighbours equal for each site - can use first list to get no. of
                                    #neighbours
                                    distanceDiff = abs(nnDistanceLists[i][k]-nnDistanceLists[j][k])
                                    distanceComparison = distanceDiff < tolerance
                                    distanceComparisonList.append(distanceComparison)
                                if(False not in distanceComparisonList):
                                    siteComparisonList.append(True)
                                elif(False in distanceComparisonList):
                                    siteComparisonList.append(False)

                    print(f"Same (within tolerance) distance to neighbours?: {siteComparisonList}")

                    if(False not in siteComparisonList):
                        print(f"Candidate found: {material}")
                        siteEquivalenceResults[material] = results[material]
                    else:
                        print(f"{material} failed on nnDistance criterion.")
                else:
                    print(f"{material} failed on the atomType criterion.")
            else:
                print(f"{material} failed on the same no. of neighbours criterion.")

        return siteEquivalenceResults