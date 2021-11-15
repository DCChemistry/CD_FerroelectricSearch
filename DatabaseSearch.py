from pymatgen.ext.matproj import MPRester
from pymatgen.core.periodic_table import Element
import numpy as np
import matplotlib.pyplot as plt
from json_tricks import dumps, loads #the json module doesn't support non-standard types (such as the output from MAPI),
                                     #but json_tricks does
from ChargeDisproportation import * #this imports all functions, variables and classes within the
#ChargeDisproportionation.py file
import os
import time
    


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


def SaveDictAsJSON(fileName, dictionary, indent=None):
    with open(fileName, "w") as f:
        f.write(dumps(dictionary, indent=indent)) #don't need to read this since it's just a 'checkpoint'

def ReadJSONFile(fileName):
    with open(fileName, "r") as f:
        return loads(f.read()) #loads() returns the string from f.read() as dict


def DatabaseSearch(searchFileName, elementList, excludeList):

    if(not os.path.isfile(f"{searchFileName}.json")): #if given file doesn't exist, then run the search
        



        APIkey = None #done so that APIkey is not lost in the scope of the with block
        with open("APIkey.txt", "r") as f:
            APIkey= f.read()

        results = None #done so that results exists outside the scope of the with block
        with MPRester(APIkey) as mpr:
 
            criteria = {"elements": {"$in": elementList, "$nin": excludeList}, "band_gap": {"$gt": 0.0}}
            # ^ want to find materials that contain any of the listed elements, hence $in, $nin excludes elements in given list,
            # and $gt is simply 'greater than' - ferroelectrics are insulators, and DFT underestimates band gaps greatly,
            # so if the band gap is > 0, that means the band gap is sizeable (likely insulator). NEW, RUN THIS SOON

            properties = ["material_id", "pretty_formula", "spacegroup.symbol", "spacegroup.crystal_system", "e_above_hull"]
            results = mpr.query(criteria, properties, chunk_size=10000) #it's ok to change chunk_size since I'm asking for small things -
                                                                        #it speeds things up tremendously (asking for larger amounts of
                                                                        #data less often - latency - sending info back and forth takes
                                                                        # time)

        SaveDictAsJSON(f"{searchFileName}.json", results)

    else:
        results = ReadJSONFile(f"{searchFileName}.json") #[:1000] #just want to run program with the first 1000 results - remove slice later
    
    t1 = time.time()
    CheckForCD(results, searchFileName) #this is a little scuffed, but it works. it used to be resultsCD = ...
    t2 = time.time()
    print(f"Task took {t2-t1:.2f} s")


def HistogramMaker(CDResults):
    """Takes CheckForCD results in JSON format (dict), and returns histograms for each property that was requested in the original
    MAPI query."""
    
    properties = list(CDResults[list(CDResults)[0]].keys())
    #^ ok, so, there are a few things going on here that will likely make no sense in future, but I'll try to explain them now while
    # it's fresh:
    # 1) list(CDResults) returns a list of the keys of CDResults (I guess like using .keys(), but apparently you don't need
    # that part).
    # 2) I want to just look at the first element inside the CDResults dict, so I use list(CDResults)[0], which gives me the first
    # key (material id string).
    # 3) Now, then the key I got can be used to get the associated value (dictionary of the properties I ordered in the MAPI query),
    # using CDResults[list(CDResults)[0]].
    # 4) I then want the keys of the dictionary with all the properties inside ('property tags') and turn them into a list, from
    # which I can use to create individual histograms for each of the ordered properties; the resultant line of code is what is
    # above: list(CDResults[list(CDResults)[0]].keys())

    #ok, gotta do some more stuff here - for loop using the properties list to create histograms. Can probably make a separate
    # program to use the currently existing CD results file.

def main():
    nonRadElements, radElements = NonRadElements()
    #DatabaseSearch("NonRadSearch", nonRadElements, radElements)
    popularElements = ["Mn", "V", "Fe", "Ni", "Co", "Cu", "Bi", "Ti", "Eu", "Sm", "Yb"]
    DatabaseSearch("popularCOElemSearch", popularElements, radElements)
    #transitionMetals = list(np.arange(21, 30+1))+list(np.arange(39, 48+1))+list(np.arange(72, 80+1))+list(np.arange(104, 112+1))
    #fBlock = list(np.arange(57, 71+1))+list(np.arange(89, 103+1))
    #semiMetals = 
    #DatabaseSearch("FullSearch", , radElements)

if __name__ == "__main__": #if this file is run, call the chosen function below
    import cProfile
    cProfile.run("main()")
    #main()
