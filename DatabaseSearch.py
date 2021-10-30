from pymatgen import Composition
from pymatgen import MPRester
from pymatgen.core.periodic_table import Element
import numpy as np
import matplotlib.pyplot as plt
from json_tricks import dumps, loads #the json module doesn't support non-standard types (such as the output from MAPI),
                                     #but json_tricks does

ironOxides = ["FeO", "Fe3O4", "LuFe2O4", "Fe2O3", "BiAlO3", "BiInO3"]

def OxidationStateCalc(formula):
    """Returns a pymatgen list-like object of elements with their respective oxidation states (OS) from a given chemical
    formula."""
    chemical = Composition(formula) #converting formula to pymatgen Composition object
    oxStates = list(chemical.add_charges_from_oxi_state_guesses())
    oxStates = [str(element) for element in oxStates]
    return oxStates #returns the number of elements, each with a charge assigned
                                                         #(multiple instances of an element with different charges will be
                                                         #shown if the OS isn't 'normal' - this is the basis of this program)

def SiteCentredCO(material):
    """Returns a material that is likely to undergo charge disproportionation (CD) (as well as the element that is likely
    to be undergoing CD) from a given formula."""
    oxStates = OxidationStateCalc(material)
    print(f"Elements and their OS: {oxStates}")
    oxStates = [element[:-2] for element in oxStates] #removes the charge from each element, e.g. Fe2+ and Fe3+ -> Fe and Fe
    
    print(f"\noxStates: {oxStates}")###############
    elements = list(set(oxStates)) #done to remove duplicates - this is done to avoid analysing an element for CD more than once
    for element in elements:
        print(f"Elements: {elements}")##############
        print(f"Testing element: {element}")#################
        instances = oxStates.count(element) #originally, this was "elements.count(element)", which was using the non-duplicate
                                            #list. Things why the try below always failed, since no CO materials could ever be
                                            #found
        print(f"Instances of {element}: {instances}")
        if(instances>1):
            print(f"CD element: {element}")############
            return material, element

def CheckForCD(listOfMaterials):
    """Takes in a list and returns a dictionary of materials with they key as the material that seems to have undergone charge
    disproportionation (CD), and the value as the element undergoing CD."""
    
    siteCOmaterials = {}
    for material in listOfMaterials:
        print(f"\n\nMaterial: {material}")
        try:
            COmaterial, CDelement = SiteCentredCO(material)
            siteCOmaterials[COmaterial] = CDelement

        except TypeError:
            pass
            #TypeError has occurred - cannot unpack non-iterable NoneType object. Material was {material}"

    return siteCOmaterials       


def ListOfTheElements(elementsExcluded=None):
    noOfElements = 118 #as of 2021
    atomicNos = np.arange(1, noOfElements+1) #the function stops just before the given value (default step = 1)
    
    if(elementsExcluded != None):
        atomicNos = [z for z in atomicNos if z not in elementsExcluded]
    
    symbolsTypeElement = [Element.from_Z(z) for z in atomicNos]
    symbols = [str(symbol) for symbol in symbolsTypeElement]
    
    return symbols

def CleanUpResults(results):
    """Will take a MAPI output, take only the values of each dictionary element, and make the elements into lists of
    those values - it makes the results easier to deal with."""
    resultsClean = [list(results[results.index(i)].values()) for i in results]
    #Regarding the above line - i refers to an element in the results list, so results.index(i) means 'return the index
    #from the results list using the element i', which is then being used to acquire the values of each element (where
    #each element in the list is a dictionary), since the values give me the actual results - in this case, the formula
    #and space group symbol. That value then needs to be converted to a list, otherwise the type would be "dict_values".
    
    resultsCleanTranspose = list(zip(*resultsClean)) #separates properties into separate lists, e.g.
    #one element will be a list of "pretty_formula"s and the other would be "spacegroup.symbol"s.
    
    return resultsCleanTranspose



#Now it's time for the search
radElements = [43, 61]+list(np.arange(84, 118+1)) #atomic nos.: 43, 61, 84-118. https://www.epa.gov/radiation/radioactive-decay
radElementSymbols = [str(Element.from_Z(radE)) for radE in radElements] #converting the atomic nos. into Element-type
                                                                        #variables, and then converting them into strings
print(f"Radioactive elements: {radElementSymbols}")

nonRadElements = ListOfTheElements(radElements) #list of elements to search with
print(nonRadElements)
APIkey = input("Please input your API key: ")





results = None #done so that results exists outside the scope of the with block
with MPRester(APIkey) as mpr:
    criteria = {"elements": {"$in": nonRadElements}} #want to find materials that contain any of the listed elements, hence $in
    properties = ["pretty_formula", "spacegroup.symbol"]
    results = mpr.query(criteria, properties)
    
    #^ currently outputting 0, since my criteria states that all of the elements I listed should be in a compound
    #- impossible - need to look at MongoDB operator syntax again (fixed - wrote '$all' instead of '$in' - see above)
    #resultsCD = CheckForCD(results) - currently won't work since 'results' isn't just a list.
    #First, I need to make a list of 'pretty_formula's, and then I can use the CD function
    #- maybe create a function that does this? Issue is that I'll have different queries every time, so there's not much point
    #in making a general function.

    with open("NonRadSearch.txt", "w") as f: #this is the thing I'm working on now
        f.write(dumps())
    
#new resultsCD using 'CleanUpResults':
resultsClean = CleanUpResults(results)
formulas = resultsClean.pop(0) #removing and saving the first element from resultsClean
spaceGroupSymbols = resultsClean.pop(0) #will use later for finding most common space groups for CD FE candidates
# ^ regarding the last line, the pop index has to be 0 since the previous pop removed the item which originally had index 0
resultsCD = CheckForCD(formulas)

#will want to include space group plotting thingy since this process takes so long. Actually, may want to write results to a
#file to do further (and faster) analysis later (I don't want to lok through the entire database everytime I want to analyse
#my results).