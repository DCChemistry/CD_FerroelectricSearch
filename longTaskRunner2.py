# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 18:29:23 2021

@author: Dan
"""

import re
import glob
import numpy as np
import os
import math
import time
from pymatgen.core.composition import Composition
from Util import *
from datetime import datetime

def ResumeFrom2(fileName, noOfTasks):

    generalTaskFileName = f"{fileName}_task_*"
    tasks = glob.glob(generalTaskFileName)
    if(len(tasks) != 0): #if task files for this run exist, do the following
        taskIndices = []
        for i in range(len(tasks)):
            taskNoPart = tasks[i].replace(fileName, "")
            noRegex = re.compile('[0-9]+')
            # ^ should get the first instance of a 'number section'
            taskIndex = int(re.findall(noRegex, taskNoPart)[0])
            taskIndices.append(taskIndex)
        assignedTaskIndices = list(np.arange(noOfTasks))
        remainingTaskIndices = list(set(assignedTaskIndices)-set(taskIndices))
        print(remainingTaskIndices)
        # ^ sets allow this comparison to take place using the minus operator to get the remaining indices
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

def CheckForCD(results): #simplified for this purpose
    for material in results:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print(f"\n[{current_time}]: Now working on {material['pretty_formula']}...")
        try:
            formula = material["pretty_formula"]
            COmaterial, CDelement = SiteCentredCO(formula)
            # ^ need this since this is the part that actually does the analysis, even if nothing is done with the returned values
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            print(f"[{current_time}]: {material['pretty_formula']} was ferroelectric.")
        except TypeError:
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            print(f"[{current_time}]: {material['pretty_formula']} was not ferroelectric.")
            #TypeError has occurred - cannot unpack non-iterable NoneType object. Material was {material}"

def LongTaskRunner(searchFileName):
    noOfTasks = 1024
    remainingTaskIndices = ResumeFrom2(searchFileName, noOfTasks)

    results = ReadJSONFile(searchFileName)


    materialsPerTask = math.floor(len(results)/noOfTasks)

    tasks = []
    for i in range(noOfTasks):

        tasks.append(results[i*materialsPerTask: min((i+1)*materialsPerTask, len(results))])
        #regarding min(), it returns the lower of the two: either (i+1)*tasksPerProcessor, or the length of the results array (using this
        #also takes the 'remaining tasks' into account).

    for i in remainingTaskIndices:
        print(f"Working on task {i}")
        longTask = tasks[i]
        longTaskFormulas = [formula["pretty_formula"] for formula in longTask]
        #print(longTaskFormulas)
        print(f"No. of elements in task: {len(longTaskFormulas)}")
        t1 = time.time()
        CheckForCD(longTask)
        t2 = time.time()
        print(f"Task took {t2-t1} s!")
        #SaveDictAsJSON(f"LongTask{longTaskIndex}", longTaskCDCandidates)

def FileMerger(noOfTasks, searchFileName):
    """Merges the task files together into one, and deletes the task files."""
    mergedFileDict = {}
    for i in range(noOfTasks):
        taskResults = ReadJSONFile(f"{searchFileName}_task_{i}.json")
        mergedFileDict.update(dict(taskResults))

    SaveDictAsJSON(f"{searchFileName}CDCandidates.json", mergedFileDict, indent=4)
    for i in range(noOfTasks):
        os.remove(f"{searchFileName}_task_{i}.json")

#FileMerger(noOfTasks, searchFileName)


def Main():
    searchFileName = "NonRadSearch2_oneAtom"
    LongTaskRunner(searchFileName)

Main()