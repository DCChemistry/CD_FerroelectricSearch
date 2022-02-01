# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 11:12:11 2021

@author: Dan
"""

from json_tricks import dumps, loads
import matplotlib.pyplot as plt
import seaborn as sn



def SaveDictAsJSON(fileName, dictionary, indent=None):
    with open(fileName, "w") as f:
        f.write(dumps(dictionary, indent=indent)) #don't need to read this since it's just a 'checkpoint'

def ReadJSONFile(fileName):
    with open(fileName, "r") as f:
        return loads(f.read()) #loads() returns the string from f.read() as dict

noOfValues = []
for i in range(64):
    results = ReadJSONFile(f"popularCOElemSearch_task_{i+1}.json")
    noOfValues.append(len(results.keys()))

sn.set_context("paper")
plt.bar(range(len(noOfValues)), noOfValues)
#plt.xticks(range(len(counter)), [key[0:5] for key in list(counter.keys())], rotation=90)
plt.xlabel("Task No.")
plt.ylabel("No. of values in task")
plt.title(f"Task Size Distribution")
plt.show()