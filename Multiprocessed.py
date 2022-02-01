import os
import numpy as np #added just in case - may add numpy objects later
import re
import concurrent.futures
import math
from Util import *
import glob #for unix style searching
import multiprocessing


class MultiProcessing: #used for any analysis after CheckForCD - may introduce this for CheckForCD if it's been shown to work

    def __init__(self, fileName, analysisTag, func, noOfTasks=1024, cores=4):

        self.fileName = fileName
        self.analysisTag = analysisTag
        # ^ analysisTag will be appended to the resultant file, 
        # e.g. the resultant file from CheckForCD has 'CDCandidates' at the end of the file name.

        if(not os.path.isfile(f"{self.fileName}{self.analysisTag}.json")):
            
            print(f"Starting {self.analysisTag} analysis.")
            self.func = func
            # func is the function that is being 'multiprocessed' - appears in TaskMaster
            self.noOfTasks = noOfTasks
            self.cores = cores

            results = ReadJSONFile(fileName)
            self.MultiThreading(results)

        else:
            print(f"{self.analysisTag} analysis has already been done for search {self.fileName}")
        

    def FileMerger(self):
        """Merges the task files together into one, and deletes the task files."""
        mergedFileDict = {}
        for i in range(self.noOfTasks):
            taskResults = ReadJSONFile(f"{self.fileName}_task_{i}")
            mergedFileDict.update(dict(taskResults))

        SaveDictAsJSON(f"{self.fileName}{self.analysisTag}_({len(mergedFileDict.keys())})", mergedFileDict)
        for i in range(self.noOfTasks):
            os.remove(f"{self.fileName}_task_{i}.json")

    def TaskMaster(self, results, taskNo):
        SaveDictAsJSON(f"{self.fileName}_task_{taskNo}", self.func(results))
        print(f"Task {taskNo}/{self.noOfTasks} complete!")

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


    def MultiThreading(self, results):
        
        remainingTaskIndices = self.ResumeFrom()
        if(remainingTaskIndices == None): #if task files for this run don't exist
            indices = range(self.noOfTasks)
            print(f"Initiating fresh run. No task files for the given file name, {self.fileName}{self.analysisTag}, found.")
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

        with concurrent.futures.ProcessPoolExecutor(max_workers = min(multiprocessing.cpu_count(), self.cores)) as executor:

            for i in indices: #indices will differ depending on whether task files already exist (see start of func def)
                executor.submit(self.TaskMaster, tasks[i], i)
                #^ calls the CheckForCD function with each task created above from the executor thread that the process pool is on.
        self.FileMerger() #merges the task files into a new, single file, and deletes the task files.