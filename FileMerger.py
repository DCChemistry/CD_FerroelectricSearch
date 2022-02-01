# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 11:49:21 2021

@author: Dan
"""

from json_tricks import dumps, loads

def SaveDictAsJSON(fileName, dictionary, indent=None):
    with open(fileName, "w") as f:
        f.write(dumps(dictionary, indent=indent)) #don't need to read this since it's just a 'checkpoint'

def ReadJSONFile(fileName):
    with open(fileName, "r") as f:
        return loads(f.read()) #loads() returns the string from f.read() as dict
    
siteCOMaterials = {}
for i in range(64):
    taskResults = ReadJSONFile(f"popularCOElemSearch_task_{i+1}.json")
    siteCOMaterials.update(dict(taskResults))

SaveDictAsJSON("popularCOElemSearchCDCandidates.json", siteCOMaterials, indent=4)