import os
from pymatgen.core.composition import Composition
import numpy as np #added just in case - may add numpy objects later
import re
import concurrent.futures
import multiprocessing
import math
from Util import SaveDictAsJSON, ReadJSONFile


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

class CheckForCD:

    def __init__(self, results, fileName, noOfTasks):
        self.fileName = fileName
        processor_count = multiprocessing.cpu_count()
        self.noOfTasks = noOfTasks
        self.MultiThreadedCheckForCD(results)

    def FileMerger(self):
        """Merges the task files together into one, and deletes the task files."""
        mergedFileDict = {}
        for i in range(self.noOfTasks):
            taskResults = ReadJSONFile(f"{self.fileName}_task_{i}")
            mergedFileDict.update(dict(taskResults))

        SaveDictAsJSON(f"{self.fileName}CDCandidates.json", mergedFileDict, indent=4)
        for i in range(self.noOfTasks):
            os.remove(f"{self.fileName}_task_{i}.json")

    def CheckForCDTaskMaster(self, results, taskNo):
        SaveDictAsJSON(f"{self.fileName}_task_{taskNo}.json", self.CheckForCD(results))
        print(f"Task {taskNo}/{self.noOfTasks} complete!")

    def CheckForCD(self, results):

        siteCOmaterials = {}

        for material in results:
            try:
                formula = material["pretty_formula"]
                COmaterial, CDelement = SiteCentredCO(formula)
                siteCOmaterials[material["material_id"]] = {
                        "pretty_formula": COmaterial,
                        "spacegroup.symbol": material["spacegroup.symbol"],
                        "spacegroup.crystal_system": material["spacegroup.crystal_system"],
                        "e_above_hull": material["e_above_hull"],
                        "CDelement": CDelement
                }

            except TypeError:
                pass
                #TypeError has occurred - cannot unpack non-iterable NoneType object. Material was {material}"
        
            #v this is done to show that the program is currently on this function, and that it works

            #if((i+1)%100 == 0):   WOULD LIKE TO GET THIS TO WORK - NEED TASKNO TO MATCH THE TASKNO IN THE MULTITHREADED VERSION (SEE BELOW)
            #    print(f"{i+1}/{len(results)} of task {taskNo} done.")  

        return siteCOmaterials 


    def MultiThreadedCheckForCD(self, results):
        


        materialsPerTask = math.floor(len(results)/self.noOfTasks)

        tasks = []
        for i in range(self.noOfTasks):

            tasks.append(results[i*materialsPerTask: min((i+1)*materialsPerTask, len(results))])
            #regarding min(), it returns the lower of the two: either (i+1)*tasksPerProcessor, or the length of the results array (using this
            #also takes the 'remaining tasks' into account).

        with concurrent.futures.ProcessPoolExecutor() as executor:

            for i in range(self.noOfTasks): #and thus, the length of futures = noOfTasks
                executor.submit(self.CheckForCDTaskMaster, tasks[i], i)
                #^ calls the CheckForCD function with each task created above from the executor thread that the process pool is on.
        self.FileMerger() #merges the task files into a new, single file, and deletes the task files.