from pymatgen.ext.matproj import MPRester
from pymatgen.core.periodic_table import Element
import numpy as np
import matplotlib.pyplot as plt
from Util import *
from ChargeDisproportation import * #this imports all functions, variables and classes within the
#ChargeDisproportionation.py file
import os
import time
from Analysis import*
    


def ListOfTheElements(elementsExcluded=None):
    noOfElements = 118 #as of 2021
    atomicNos = np.arange(1, noOfElements+1) #the function stops just before the given value (default step = 1)
    
    if(elementsExcluded != None):
        atomicNos = [z for z in atomicNos if z not in elementsExcluded]
    
    symbolsTypeElement = [Element.from_Z(z) for z in atomicNos]
    symbols = [str(symbol) for symbol in symbolsTypeElement]
    
    return symbols


def NonRadElements():
    """Returns all of the non-radioactive elements, and the radioactive ones as a separate result"""
    radElements = [43, 61]+list(np.arange(84, 118+1)) #atomic nos.: 43, 61, 84-118. https://www.epa.gov/radiation/radioactive-decay
    radElementSymbols = [str(Element.from_Z(radE)) for radE in radElements] #converting the atomic nos. into Element-type
                                                                            #variables, and then converting them into strings
    nonRadElements = ListOfTheElements(radElements) #list of elements to search with
    return nonRadElements, radElementSymbols


def DatabaseSearch(searchFileName, elementList, excludeList, orderOfFilters, noOfTasks=1024):
    print("Hello, and welcome to your database search!")
    
    if(not os.path.isfile(f"{searchFileName}.json")): #if given file doesn't exist, then run the search
        print(f"{searchFileName}.json not found. Creating file and querying the Materials Project Database.")



        APIkey = None #done so that APIkey is not lost in the scope of the with block
        with open("APIkey.txt", "r") as f:
            APIkey= f.read()

        results = None #done so that results exists outside the scope of the with block
        with MPRester(APIkey) as mpr:
 
            criteria = {"elements": {"$in": elementList, "$nin": excludeList}}
            # ^ want to find materials that contain any of the listed elements, hence $in, $nin excludes elements in given list,
            # and $gt is simply 'greater than' - ferroelectrics are insulators, and DFT underestimates band gaps greatly,
            # so if the band gap is > 0, that means the band gap is sizeable (likely insulator). NEW, RUN THIS SOON

            properties = ['material_id', 'pretty_formula', 'spacegroup.number', 'band_gap','nsites', "e_above_hull", "nelements"]
            results = mpr.query(criteria, properties, chunk_size=10000) #it's ok to change chunk_size since I'm asking for small things -
                                                                        #it speeds things up tremendously (asking for larger amounts of
                                                                        #data less often - latency - sending info back and forth takes
                                                                        # time)

        SaveDictAsJSON(searchFileName, results)

    else:
        print(f"File {searchFileName}.json found. Loading file.")
        results = ReadJSONFile(searchFileName)

    t1 = time.time()
    CheckForCD(results, searchFileName, noOfTasks)
    #MultiProcessing(searchFileName, "NP", Analysis.NonPolar)
    Analysis(searchFileName, orderOfFilters)
    t2 = time.time()
    print(f"Task took {t2-t1:.2f} s")




def main():
    nonRadElements, radElements = NonRadElements()
    #DatabaseSearch("NonRadSearch2", nonRadElements, radElements)
    DatabaseSearch("NonRadSearch2", nonRadElements, radElements, ["NP", "oneCDSite"])

if __name__ == "__main__": #if this file is run, call the chosen function below
    #import cProfile
    #cProfile.run("main()") #does the same thing as a profiler - helpful for CheckForCD which takes a long time to run
    main()
