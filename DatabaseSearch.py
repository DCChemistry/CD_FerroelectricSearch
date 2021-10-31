from pymatgen.ext.matproj import MPRester
from pymatgen.core.periodic_table import Element
import numpy as np
import matplotlib.pyplot as plt
from json_tricks import dumps, loads #the json module doesn't support non-standard types (such as the output from MAPI),
                                     #but json_tricks does
from ChargeDisproportation import * #this imports all functions, variables and classes within the
#ChargeDisproportionation.py file
import os
    


def ListOfTheElements(elementsExcluded=None):
    noOfElements = 118 #as of 2021
    atomicNos = np.arange(1, noOfElements+1) #the function stops just before the given value (default step = 1)
    
    if(elementsExcluded != None):
        atomicNos = [z for z in atomicNos if z not in elementsExcluded]
    
    symbolsTypeElement = [Element.from_Z(z) for z in atomicNos]
    symbols = [str(symbol) for symbol in symbolsTypeElement]
    
    return symbols

def SaveDictAsJSON(fileName, dictionary, indent=None):
    with open(fileName, "w") as f:
        f.write(dumps(dictionary, indent=indent)) #don't need to read this since it's just a 'checkpoint'

def ReadJSONFile(fileName):
    with open(fileName, "r") as f:
        return loads(f.read()) #loads() returns the string from f.read() as dict


def main():

    if(not os.path.isfile("NonRadSearch.json")): #if given file doesn't exist, then run the search
        

        #Now it's time for the search
        radElements = [43, 61]+list(np.arange(84, 118+1)) #atomic nos.: 43, 61, 84-118. https://www.epa.gov/radiation/radioactive-decay
        radElementSymbols = [str(Element.from_Z(radE)) for radE in radElements] #converting the atomic nos. into Element-type
                                                                                #variables, and then converting them into strings
        nonRadElements = ListOfTheElements(radElements) #list of elements to search with

        APIkey = None #done so that APIkey is not lost in the scope of the with block
        with open("APIkey.txt", "r") as f:
            APIkey= f.read()

        results = None #done so that results exists outside the scope of the with block
        with MPRester(APIkey) as mpr:
            criteria = {"elements": {"$in": nonRadElements}} #want to find materials that contain any of the listed elements, hence $in
            properties = ["material_id", "pretty_formula", "spacegroup.symbol", "spacegroup.crystal_system"]
            results = mpr.query(criteria, properties, chunk_size=10000) #it's ok to change chunk_size since I'm asking for small things -
                                                                        #it speeds things up tremendously (asking for larger amounts of
                                                                        #data less often - latency - sending info back and forth takes
                                                                        # time)

        SaveDictAsJSON("NonRadSearch.json", results)

    else:
        results = ReadJSONFile("NonRadSearch.json")
    
    resultsCD = CheckForCD(results)

    SaveDictAsJSON("NonRadSearchCDCandidates.json", resultsCD, indent=4)

if __name__ == "__main__": #if this file is run, call main function
    main()