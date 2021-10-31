from pymatgen.core.composition import Composition
import numpy as np #added just in case - may add numpy objects later
import re

alphabetRegex = re.compile('[a-zA-Z]+') #need this for removing numbers from oxidation states (see below)
#- a-zA-Z gets all letters (upper and lowercase), and the subsequent '+' reads those letters as a group (otherwise
#you get hellish results), e.g. ["Z", "n"]

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

    #v removes the charge from each element, e.g. Fe2+ and Fe3+ -> Fe and Fe
    oxStates = [alphabetRegex.findall(element)[0] for element in oxStates]
    #^ alphabetRegex.findall() looks for the occurence of letters in each OS and returns an array with the string.
    #Need [0] to pop the result (the string) out of the array (of length 1).
    
    elements = list(set(oxStates)) #done to remove duplicates - this is done to avoid analysing an element for CD more than once
    for element in elements:
        instances = oxStates.count(element) #originally, this was "elements.count(element)", which was using the non-duplicate
                                            #list. Hence why the try below always failed, since no CO materials could ever be
                                            #found
        if(instances>1):
            return material, element

def CheckForCD(results):
    """Takes in a list and returns a dictionary of materials with they key as the material that seems to have undergone charge
    disproportionation (CD), and the value as the element undergoing CD."""
    
    siteCOmaterials = {}
    for i, material in enumerate(results):
        try:
            formula = material["pretty_formula"]
            COmaterial, CDelement = SiteCentredCO(formula)
            siteCOmaterials[COmaterial] = CDelement

        except TypeError:
            pass
            #TypeError has occurred - cannot unpack non-iterable NoneType object. Material was {material}"

        #v this is done to show that the program is currently on this function, and that it works
        if(i%100 == 0):   
            print(i)

        if(i==500):
            break

    return siteCOmaterials   