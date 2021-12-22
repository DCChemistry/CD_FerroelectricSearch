import matplotlib.pyplot as plt
from Util import *
from pymatgen.core.composition import Composition
#Smidt imports
from pymatgen.symmetry.groups import SYMM_DATA, sg_symbol_from_int_number

import numpy as np #added just in case - may add numpy objects later
import matplotlib.pyplot as plt
import os

class Analysis:

    def __init__(self, searchFileName):
        self.searchFileName = searchFileName
        #from CD results to NP results
        self.ReadAnalyseWrite("NonPolar", "CDCandidates", "NP")

        #from NP results to one primitive cell results
        self.ReadAnalyseWrite("OnePrimCell", "NP", "onePC")


    def ReadAnalyseWrite(self, AnalysisType, prevAnalysisTag, newAnalysisTag):
        """AnalysisType is the name of the method used to analyse the data, e.g. NonPolar.
           prevAnalysisTag is the text appeneded to the end of the analysis file you want to load.
           newAnalysisTag that will be appended to the end of the analysis file you want to create."""
        
        if(not os.path.isfile(f"{self.searchFileName}{newAnalysisTag}.json")):
            print(f"Starting {newAnalysisTag} analysis.")
            results = ReadJSONFile(f"{self.searchFileName}{prevAnalysisTag}")
            analysisResults = eval(f"Analysis.{AnalysisType}(results)")
            SaveDictAsJSON(f"{self.searchFileName}{newAnalysisTag}", analysisResults)
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
    def OnePrimCell(results):
        print("Starting primitive cell analysis.")
        primCellResults = {}
        for material in results:
            comp = Composition(results[material]["pretty_formula"])
            reducedFormulaAtomNo = comp.num_atoms
            noOfAtomsInStruct = int(results[material]["nsites"])
            noOfPrimCells = noOfAtomsInStruct/reducedFormulaAtomNo
            if(noOfPrimCells == 1):
                primCellResults[material] = results[material]
            elif(noOfAtomsInStruct%reducedFormulaAtomNo!=0):
                #^ Identifies materials that have a fractional no. of primitive cell in the calculated cell.
                # There are very few of these, and the ones found in my run were H2 and N2 only. All other materials have an int no. of prim. cells.
                print(f"Fractional no. of primitive cells found for material: {results[material]['pretty_formula']}.\nPrim. cell ratio = {noOfPrimCells}\n")
        
        print(f"Primitive cell analysis complete. {len(primCellResults.keys())} materials found.")
        return primCellResults