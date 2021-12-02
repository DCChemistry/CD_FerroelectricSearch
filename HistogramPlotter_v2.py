# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 16:31:20 2021

@author: Dan
"""
import matplotlib.pyplot as plt
from Util import SaveDictAsJSON, ReadJSONFile
import numpy as np
import os
    

def HistogramMaker(CDResults):
    """Takes CheckForCD results in JSON format (dict), and returns histograms for each property that was requested in the original
    MAPI query."""
   
    folderName = "CDResults_plots"
    if(not os.path.exists(folderName)):
        os.mkdir(folderName)
    #fig.suptitle("Important properties", fontsize=18) #currently this overlaps with the top subplot title
    elements = [r["CDelement"] for r in CDResults.values()] #getting a list of the values for a certain property, e.g. space group, crystal system, etc.
    counter = {} #saves the type of item, e.g. space group symbols as keys and their respective values are the frequency at which
                    #they appeared in the items list
    noOfItems = len(elements)
    for i, element in enumerate(elements):
        print(f"{i+1}/{noOfItems}")
        if(element in counter): #if item is already in the counter dictionary
            counter[element] += 1 #then increase the count by 1 for that item
        else: #if item (key) isn't already 
            counter[element] = 1
    counter = {k: v for k, v in sorted(counter.items(), key=lambda item: item[1], reverse=True)} #reverse=True inverts the sorted order.
    # ^ https://stackoverflow.com/questions/613183/how-do-i-sort-a-dictionary-by-value
    
    #Charge disproportionation element frequency plot
    plt.figure(figsize=(12, 5))
    plt.bar(range(len(counter)), counter.values(), align="center")
    plt.xticks(range(len(counter)))
    plt.title("CDelement Frequency")
    plt.ylabel("No. of materials")
    plt.xlabel("Charge disproportionation elements")

    plt.xticks(range(len(counter)), [key[0:5] for key in list(counter.keys())])
    plt.tight_layout()
    plt.savefig(f"{folderName}/NonRadSearchFreq.png", dpi=300, bbox_inches="tight")

    #Energy above hull plot
    blankLists = []
    for i in range(len(elements)):
        blankLists.append([])

    eAboveHullForCDElem = dict(zip(counter.keys(), blankLists))
    # ^ for each element (using counter.keys() because I want to use the same element order), there is initially a blank list attached as a value.
    # I will be appending e_above_hull values to these lists, then using np.mean to get average e_above_hull values.
    # By checking the CDElement parameter for each result, can append to the relevant list. This needs to be done here instead of above with counter, 
    # since I want to have the same order of elements as the sorted counter dictionary.#
    nullCases= {} #contains entries that have an e_above_hull = null
    for id in CDResults: #iterating over a dict gives its keys, so the keys in this case are the material IDs
        CDElem = CDResults[id]["CDelement"] #gets me the correct key for eAboveHullForCDElem
        eAboveHull = CDResults[id]["e_above_hull"] #gives me the e_above_hull I want to append to the respective list
        if(eAboveHull==None):
            nullCases.update({id: CDResults[id]})
            # ^ for some reason, some entries have 'null' as their e_above_hull value, despite the web site entries having values
        else:
            eAboveHullForCDElem[CDElem].append(eAboveHull) #appending the e_above_hull to the appropriate list
    SaveDictAsJSON("Energy_above_hull_is_null", nullCases, indent=4)
    
    eAboveHullForCDElem = {k:float(f"{np.mean(v):.2e}") for (k,v) in eAboveHullForCDElem.items()}
    # ^ converting the e_above_hull lists into averages in scientific (exponent) notation
    plt.figure(figsize=(12, 5))
    plt.bar(range(len(eAboveHullForCDElem)), eAboveHullForCDElem.values(), align="center")
    plt.xticks(range(len(counter)))
    plt.title("Energy Above Convex Hull Per CD Element")
    plt.ylabel("Energy above hull / eV")
    plt.xlabel("Charge dissproportionation elements")
    
    plt.xticks(range(len(counter)), [key[0:5] for key in list(counter.keys())])
    plt.tight_layout()
    plt.savefig(f"{folderName}/NonRadSearchHull.png", dpi=300, bbox_inches="tight")

        
        
                


"""with MPRester() as mpr:
    # Collect all the entries for SiO2
    results = mpr.query(criteria='SiO2', 
                 properties=["material_id", "spacegroup.crystal_system"])

    # Collect the crystal systems into a list and count them
    systems = [r["spacegroup.crystal_system"] for r in results]
    count_dict = {}
    for system in systems:
        if system in count_dict:
            count_dict[system] += 1
        else:
            count_dict[system] = 1
            
    # Plot the distribution
    plt.bar(range(len(count_dict)), count_dict.values(), align='center')
    plt.xticks(range(len(count_dict)), [key[0:5] for key in list(count_dict.keys())])
    plt.title('Distribution of crystal systems within SiO2 compounds')
    plt.show()"""

CDResults = ReadJSONFile("NonRadSearchCDCandidates")
HistogramMaker(CDResults)