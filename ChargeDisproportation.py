import os
from pymatgen.core.composition import Composition
import numpy as np #added just in case - may add numpy objects later
import re
import concurrent.futures
import math
from Util import SaveDictAsJSON, ReadJSONFile
import matplotlib.pyplot as plt
import glob #for unix style searching


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
        self.noOfTasks = noOfTasks
        if(not os.path.isfile(f"{self.fileName}CD_0.json")):
            print("Starting CD analysis.")
            self.MultiThreadedCheckForCD(results)
            self.CDResultsPlotter() #after results are obtained, frequency and e_above_hull plots should be made
        else:
            print(f"CD analysis has already been done for search {self.fileName}")

    def FileMerger(self):
        """Merges the task files together into one, and deletes the task files."""
        mergedFileDict = {}
        for i in range(self.noOfTasks):
            taskResults = ReadJSONFile(f"{self.fileName}_task_{i}")
            mergedFileDict.update(dict(taskResults))

        SaveDictAsJSON(f"{self.fileName}CD_0", mergedFileDict, indent=4)
        for i in range(self.noOfTasks):
            os.remove(f"{self.fileName}_task_{i}.json")

    def CheckForCDTaskMaster(self, results, taskNo):
        SaveDictAsJSON(f"{self.fileName}_task_{taskNo}", self.CheckForCD(results))
        print(f"Task {taskNo}/{self.noOfTasks} complete!")

    def CheckForCD(self, results):

        siteCOmaterials = {}

        for material in results:
            try:
                formula = material["pretty_formula"]
                COmaterial, CDelement = SiteCentredCO(formula)
                siteCOmaterials[material["material_id"]] = {
                        "pretty_formula": COmaterial,
                        "spacegroup.number": material["spacegroup.number"],
                        "band_gap": material["band_gap"],
                        "nsites": material["nsites"],
                        "e_above_hull": material["e_above_hull"],
                        "nelements": material["nelements"],
                        "CDelement": CDelement
                }

            except TypeError:
                pass
                #TypeError has occurred - cannot unpack non-iterable NoneType object. Material was {material}"
        
            #v this is done to show that the program is currently on this function, and that it works

            #if((i+1)%100 == 0):   WOULD LIKE TO GET THIS TO WORK - NEED TASKNO TO MATCH THE TASKNO IN THE MULTITHREADED VERSION (SEE BELOW)
            #    print(f"{i+1}/{len(results)} of task {taskNo} done.")  

        return siteCOmaterials 

    def ResumeFrom(self):

        generalTaskFileName = f"{self.fileName}_task_*" #general form of task file names
        tasks = glob.glob(generalTaskFileName)
        if(len(tasks) != 0): #if task files for this run exist, do the following
            taskIndices = []
            for i in range(len(tasks)):
                noRegex = re.compile('[0-9]+')
                # ^ should get the first instance of a 'number section'
                taskNoPart = tasks[i].replace(self.fileName, "") #each task is a task file name, and thus follows the expected general format
                # ^ remove self.fileName from the general task file name, leaving _task_x.json, where x is the task index
                taskIndex = int(re.findall(noRegex, taskNoPart)[0]) #gotta add in [0] to get a good result from regex
                # ^ applying the number regex to taskNoPart, and converting the resultant string into an int
                taskIndices.append(taskIndex)
            assignedTaskIndices = list(np.arange(self.noOfTasks))
            remainingTaskIndices = list(set(assignedTaskIndices)-set(taskIndices))
            # ^ sets allow this comparison to take place using the minus operator to get the remaining indices
            return remainingTaskIndices

    def MultiThreadedCheckForCD(self, results):
        
        remainingTaskIndices = self.ResumeFrom()
        if(remainingTaskIndices == None): #if task files for this run don't exist
            indices = range(self.noOfTasks)
            print(f"Initiating fresh run. No task files for {self.fileName} found.")
        else: #task files already exist
            indices = remainingTaskIndices
            print("Task files found. Continuing from where we left off.")
            print(f"Number of remaining tasks: {len(indices)}")

        materialsPerTask = math.floor(len(results)/self.noOfTasks)
        tasks = []
        for i in range(self.noOfTasks):

            tasks.append(results[i*materialsPerTask: min((i+1)*materialsPerTask, len(results))])
            #regarding min(), it returns the lower of the two: either (i+1)*tasksPerProcessor, or the length of the results array (using this
            #also takes the 'remaining tasks' into account).

        with concurrent.futures.ProcessPoolExecutor() as executor:

            for i in indices: #indices will differ depending on whether task files already exist (see start of func def)
                executor.submit(self.CheckForCDTaskMaster, tasks[i], i)
                #^ calls the CheckForCD function with each task created above from the executor thread that the process pool is on.
        self.FileMerger() #merges the task files into a new, single file, and deletes the task files.


    def CDResultsPlotter(self): #formerly HistogramMaker()
        """Takes CheckForCD results in JSON format (dict), and returns histograms for each property that was requested in the original
        MAPI query."""

        CDResults = ReadJSONFile(f"{self.fileName}CD_0")
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
        plt.savefig(f"{folderName}/{self.fileName}Freq.png", dpi=300, bbox_inches="tight")

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
                print(f"{id} has e_above_hull=null")
                nullCases.update({id: CDResults[id]})
                # ^ for some reason, some entries have 'null' as their e_above_hull value, despite the web site entries having values
            else:
                eAboveHullForCDElem[CDElem].append(eAboveHull) #appending the e_above_hull to the appropriate list
        if(len(nullCases.keys()) > 0): #it isn't always the case that there are null values for e_above_hull in a search
            SaveDictAsJSON(f"Energy_above_hull_is_null_{self.fileName}", nullCases, indent=4)
        
        eAboveHullForCDElem = {k:np.mean(v) for (k,v) in eAboveHullForCDElem.items()}

        plt.figure(figsize=(12, 5))
        plt.bar(range(len(eAboveHullForCDElem)), eAboveHullForCDElem.values(), align="center")
        plt.xticks(range(len(counter)))
        plt.title("Average Energy Above Convex Hull Per CD Element")
        plt.ylabel("Average Energy above hull / eV")
        plt.xlabel("Charge dissproportionation elements")
        
        plt.xticks(range(len(counter)), [key[0:5] for key in list(counter.keys())])
        plt.tight_layout()
        plt.savefig(f"{folderName}/{self.fileName}Hull.png", dpi=300, bbox_inches="tight")