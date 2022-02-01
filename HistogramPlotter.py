# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 16:31:20 2021

@author: Dan
"""
import matplotlib.pyplot as plt
from json_tricks import loads
import os

def ReadJSONFile(fileName):
    with open(fileName, "r") as f:
        return loads(f.read()) #loads() returns the string from f.read() as dict
    

def HistogramMaker(CDResults):
    """Takes CheckForCD results in JSON format (dict), and returns histograms for each property that was requested in the original
    MAPI query."""
    
    #properties = list(CDResults[list(CDResults)[0]].keys())
    #properties.remove("pretty_formula")
    
    
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
   
    fig = plt.figure(figsize=(8, 6))
    
    #fig.suptitle("Important properties", fontsize=18) #currently this overlaps with the top subplot title
    propertiesToPlot = ["spacegroup.symbol", "spacegroup.crystal_system", "CDelement"]
    for i, property in enumerate(propertiesToPlot):
        items = [r[property] for r in CDResults.values()] #getting a list of the values for a certain property, e.g. space group, crystal system, etc.
        counter = {} #saves the type of item, e.g. space group symbols as keys and their respective values are the frequency at which
                     #they appeared in the items list
        for item in items:
            if(item in counter): #if item is already in the counter dictionary
                counter[item] += 1 #then increase the count by 1 for that item
            else: #if item (key) isn't already 
                counter[item] = 1
    
        ax = fig.add_subplot(len(propertiesToPlot), 1, i+1)
        ax.bar(range(len(counter)), counter.values(), align="center")
        ax.set_xticks(range(len(counter)))
        
        if(property == "spacegroup.symbol"):
            ax.set_xticklabels([key[0:5] for key in list(counter.keys())], fontsize=3, rotation=90)
        if(property == "spacegroup.crystal_system"):
            ax.set_xticklabels([key[0:5] for key in list(counter.keys())])
        if(property == "CDelement"):
            ax.set_xticklabels([key[0:5] for key in list(counter.keys())], fontsize=7, rotation=90)
        ax.set_title(f"{property}")

    plt.tight_layout()
    plt.savefig("PopularElemSearch.png", dpi=300, bbox_inches="tight")

        
        
                


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

CDResults = ReadJSONFile("popularCOElemSearchCDCandidates.json")
HistogramMaker(CDResults)