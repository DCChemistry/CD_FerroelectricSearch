import os
import glob
import re
import numpy as np

def ResumeFrom(fileName, noOfTasks):
        cwd = os.getcwd()

        path = f"{fileName}_task_*"
        tasks = glob.glob(path)
        print(tasks)
        if(len(tasks) != 0): #if task files for this run exist, do the following
            taskIndices = []
            for i in range(len(tasks)):
                noRegex = re.compile('[0-9]+')
                # ^ should get the first instance of a 'number section'
                taskIndex = int(re.findall(noRegex, tasks[i])[0])
                taskIndices.append(taskIndex)
            assignedTaskIndices = list(np.arange(noOfTasks))
            remainingTaskIndices = list(set(assignedTaskIndices)-set(taskIndices))
            # ^ sets allow this comparison to take place using the minus operator to get the remaining indices
            return remainingTaskIndices

remainingTaskIndices = ResumeFrom("NonRadSearch2_band_gap0", 1024) #0 is a number, the regex stops on the first number. Heck.
#print(len(remainingTaskIndices))

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

remainingTaskIndices = ResumeFrom2("NonRadSearch2_band_gap0", 1024) #0 is a number, the regex stops on the first number. Heck.
print(len(remainingTaskIndices))