# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 16:30:18 2022

@author: Dan
"""

from pymatgen.core.composition import Composition
import re

def DisplayElements(formula):
    """Returns a list of elements present from a given formula."""
    comp = Composition(formula)
    elementsFromFormulaWithNo = comp.formula
    alphabetRegex = re.compile('[a-zA-Z]+') #need this for removing numbers from string
    elementsFromFormula = alphabetRegex.findall(elementsFromFormulaWithNo) #returns list of elements present
    return elementsFromFormula

def MutuallyExclusiveElements(results, mutuallyExclusiveElementList):
    mutuallyExclusiveResults = []
    for material in results:
        listOfElemPresent = DisplayElements(material)
        counter = 0 #the counter should never exceed 1 (not mutually exclusive if >1)
        for elem in mutuallyExclusiveElementList:
            if(elem in listOfElemPresent):
                counter += 1 #one of the mutually exlusive elements is present
                if(counter > 1): #a second element from the mutually exclusive list is in this material - don't want it, move on
                    break
        if(counter == 1):
            mutuallyExclusiveResults.append(material)
    
    return mutuallyExclusiveResults

formulas = ["Sn2BiO7", "BiO2", "SbPbF6", "Sn", "SbO5Bi2", "Bi"]
MEList = ["Sn", "Pb", "Bi", "Sb"]
mutuallyExclusiveResults = MutuallyExclusiveElements(formulas, MEList)
print(mutuallyExclusiveResults)