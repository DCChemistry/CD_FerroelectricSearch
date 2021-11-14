from os import cpu_count
from pymatgen.core.composition import Composition
import numpy as np #added just in case - may add numpy objects later
import re
import concurrent.futures
import multiprocessing
import math



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

    alphabetRegex = re.compile('[a-zA-Z]+') #need this for removing numbers from oxidation states (see below)
    #- a-zA-Z gets all letters (upper and lowercase), and the subsequent '+' reads those letters as a group (otherwise
    #you get hellish results), e.g. ["Z", "n"]

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

    siteCOmaterials = {}

    #i=0

    for material in results:
        try:
            formula = material["pretty_formula"]
            COmaterial, CDelement = SiteCentredCO(formula)
            print(f"{formula} is FE.") #printing FE material to check for signs of life, & to confirm all is working as it should
            siteCOmaterials[material["material_id"]] = {
                    "pretty_formula": COmaterial,
                    "spacegroup.symbol": material["spacegroup.symbol"],
                    "spacegroup.crystal_system": material["spacegroup.crystal_system"],
                    "CDelement": CDelement
            }

            #i+=1
            #print(i)

        except TypeError:
            pass
            #TypeError has occurred - cannot unpack non-iterable NoneType object. Material was {material}"

        #v this is done to show that the program is currently on this function, and that it works


    return siteCOmaterials 


def MultiThreadedCheckForCD(results):
    
    processor_count = multiprocessing.cpu_count()

    noOfTasks = 16*processor_count

    materialsPerTask = math.floor(len(results)/noOfTasks)

    tasks = []
    for i in range(noOfTasks):

        tasks.append(results[i*materialsPerTask: min((i+1)*materialsPerTask, len(results))])
        #regarding min(), it returns the lower of the two: either (i+1)*tasksPerProcessor, or the length of the results array (using this
        #also takes the 'remaining tasks' into account).

    with concurrent.futures.ProcessPoolExecutor() as executor:

        futures = []
        for i in range(noOfTasks): #and thus, the length of futures = noOfTasks
            futures.append(executor.submit(CheckForCD, tasks[i]))
            #^ calls the CheckForCD function with each task created above from the executor thread that the process pool is on.

        siteCOmaterials = {}
        for future in futures:
            siteCOmaterials.update(future.result())
            print(f"Task {future} complete!")
            #dict.update is analogous to append for a list, and future.result() returns the results from each task (essentially squashing
            #together all the individual task results into one place).
        
        return siteCOmaterials

