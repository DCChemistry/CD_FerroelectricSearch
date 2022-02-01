from pymatgen.ext.matproj import MPRester
from pymatgen.core.periodic_table import Element
import numpy as np
import matplotlib.pyplot as plt
from Util import *
from ChargeDisproportation import * #this imports all functions, variables and classes within the
                                    #ChargeDisproportionation.py file
import os
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

def AddElementLists(*elementLists):
    """Combines any number of lists of element symbols, removing duplicates"""
    combinedElemList = []
    for elementList in elementLists:
        combinedElemList.extend(elementList)
    
    #now removing duplicates
    combinedElemList = list(dict.fromkeys(combinedElemList)) #this removes duplicates - make list into dict (no duplicate keys), then back into list
    return combinedElemList

def AtomicSymbols(listOfAtomicNumbers, removeFromPeriodicTable = False): #this is a more generic version of the NonRadElements function
    """Takes a list of atomic numbers and returns a list of the respective atomic symbols.
    If removeFromPeriodicTable = True, then the elements given will be removed from the periodic table using ListOfTheElements,
    and return the remaining elements."""

    if(removeFromPeriodicTable == True):
        tableWithElemExcluded = ListOfTheElements(listOfAtomicNumbers)
        return tableWithElemExcluded
    else:
        elementSymbols = [str(Element.from_Z(elem)) for elem in listOfAtomicNumbers]
        return elementSymbols


def DatabaseSearch(searchFileName, elementList, excludeList, orderOfFilters=None, noOfTasks=1024):
    print("\n\nHello, and welcome to your database search!")
    
    if(not os.path.isfile(f"{searchFileName}.json")): #if given file doesn't exist, then run the search
        print(f"{searchFileName}.json not found. Creating file and querying the Materials Project Database.")



        APIkey = None #done so that APIkey is not lost in the scope of the with block
        with open("APIkey.txt", "r") as f:
            APIkey= f.read()

        results = None #done so that results exists outside the scope of the with block
        with MPRester(APIkey) as mpr:
 
            criteria = {"elements": {"$in": elementList, "$nin": excludeList}, "nsites": {"$lt": 50}, "nelements": {"$eq": 3}} #added "nsites" onwards for ReducedSearch1
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


    CheckForCD(results, searchFileName, noOfTasks)

    Analysis(searchFileName, orderOfFilters, elementList)




def main():
    nonRadElements, radElements = NonRadElements()
    transitionMetalNos = list(np.arange(21, 31))+list(np.arange(39, 49))+list(np.arange(72, 81))+list(np.arange(104, 113))
    transitionMetalSymbols = AtomicSymbols(transitionMetalNos)
    excludedElementsRedSearch1 = AddElementLists(transitionMetalSymbols, radElements)
    notableAnions = ["F", "Cl", "Br", "I", "N"]
    excludedElementsRedSearch2 = AddElementLists(transitionMetalSymbols, radElements, notableAnions)
    #DatabaseSearch("NonRadSearch2", nonRadElements, radElements)
    #DatabaseSearch("NonRadSearch2", nonRadElements, radElements, orderOfFilters=["NP", "oneCDSite"])
    DatabaseSearch("ReducedSearch1", ["Sn", "Sb", "Pb", "Bi"], excludedElementsRedSearch1, orderOfFilters=["NP", "chosenCDElem", "ME", "specOS"], noOfTasks=300)
    DatabaseSearch("ReducedSearch2", ["Bi"], excludedElementsRedSearch2, orderOfFilters=["NP", "onlyOxy", "chosenCDElem", "specOS"], noOfTasks=300)

if __name__ == "__main__": #if this file is run, call the chosen function below
    #import cProfile
    #cProfile.run("main()") #does the same thing as a profiler - helpful for CheckForCD which takes a long time to run
    main()
