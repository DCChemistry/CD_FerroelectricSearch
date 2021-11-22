# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 18:29:23 2021

@author: Dan
"""

import re
import glob
from json_tricks import loads, dumps
import numpy as np
import os
import math
import time
from pymatgen.core.composition import Composition

def ReadJSONFile(fileName):
    with open(fileName, "r") as f:
        return loads(f.read()) #loads() returns the string from f.read() as dict

def SaveDictAsJSON(fileName, dictionary, indent=None):
    with open(fileName, "w") as f:
        f.write(dumps(dictionary, indent=indent)) #don't need to read this since it's just a 'checkpoint'

def RemainingTaskIndices(fileName):
    cwd = os.getcwd()
    

    path = cwd+f"/{fileName}_task_*"
    tasks = glob.glob(path)
    taskIndices = []
    for i in range(len(tasks)):
        noRegex = re.compile('[0-9]+')
        # ^ should get the first instance of a 'number section'
        taskIndex = int(re.findall(noRegex, tasks[i])[0])
        taskIndices.append(taskIndex)
    assignedTaskIndices = list(np.arange(1024))
    remainingTaskIndices = list(set(assignedTaskIndices)-set(taskIndices))
    # ^ sets allow this comparison to take place using the minus operator
    return remainingTaskIndices

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

    for material in results:
        print(f"\nNow working on {material['pretty_formula']}...")
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
            print(f"{material['pretty_formula']} was ferroelectric.")
        except TypeError:
            print(f"{material['pretty_formula']} was not ferroelectric.")
            #TypeError has occurred - cannot unpack non-iterable NoneType object. Material was {material}"
    
        #v this is done to show that the program is currently on this function, and that it works

        #if((i+1)%100 == 0):   WOULD LIKE TO GET THIS TO WORK - NEED TASKNO TO MATCH THE TASKNO IN THE MULTITHREADED VERSION (SEE BELOW)
        #    print(f"{i+1}/{len(results)} of task {taskNo} done.")  

    return siteCOmaterials 

searchFileName = "NonRadSearch"
remainingTaskIndices = RemainingTaskIndices(searchFileName)
print(remainingTaskIndices)

results = ReadJSONFile(f"{searchFileName}.json")

noOfTasks = 1024
materialsPerTask = math.floor(len(results)/noOfTasks)

tasks = []
for i in range(noOfTasks):

    tasks.append(results[i*materialsPerTask: min((i+1)*materialsPerTask, len(results))])
    #regarding min(), it returns the lower of the two: either (i+1)*tasksPerProcessor, or the length of the results array (using this
    #also takes the 'remaining tasks' into account).

longTaskIndex = remainingTaskIndices[0]
longTask = tasks[longTaskIndex]
longTaskFormulas = [formula["pretty_formula"] for formula in longTask]
#print(longTaskFormulas)
print(f"No. of elements in task: {len(longTaskFormulas)}")
theWhopper = longTask.pop(1) #the big boy
SaveDictAsJSON("Problem Child.json", theWhopper, indent=4)
t1 = time.time()
longTaskCDCandidates = CheckForCD(longTask)
t2 = time.time()
print(f"Task took {t2-t1} s!")
SaveDictAsJSON(f"LongTask{longTaskIndex}.json", longTaskCDCandidates)
#%%
def FileMerger(noOfTasks, searchFileName):
    """Merges the task files together into one, and deletes the task files."""
    mergedFileDict = {}
    for i in range(noOfTasks):
        taskResults = ReadJSONFile(f"{searchFileName}_task_{i}.json")
        mergedFileDict.update(dict(taskResults))

    SaveDictAsJSON(f"{searchFileName}CDCandidates.json", mergedFileDict, indent=4)
    for i in range(noOfTasks):
        os.remove(f"{searchFileName}_task_{i}.json")

FileMerger(noOfTasks, searchFileName)