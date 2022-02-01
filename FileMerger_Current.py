from Util import *
import os

class CheckForCD:

    def __init__(self, fileName, noOfTasks):
        self.fileName = fileName
        self.noOfTasks = noOfTasks

    def FileMerger(self):
        """Merges the task files together into one, and deletes the task files."""
        mergedFileDict = {}
        for i in range(self.noOfTasks):
            taskResults = ReadJSONFile(f"{self.fileName}_task_{i}")
            mergedFileDict.update(dict(taskResults))

        SaveDictAsJSON(f"{self.fileName}CDCandidates", mergedFileDict, indent=4)
        for i in range(self.noOfTasks):
            os.remove(f"{self.fileName}_task_{i}.json")

dataInstance = CheckForCD("NonRadSearch2", 1024)
dataInstance.FileMerger()