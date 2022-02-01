from Util import *

def MergeDictFiles(mergedFileName, *fileNames): #can take any number of file names
    mergedFileDict = {}
    for fileName in fileNames:
        results = ReadJSONFile(fileName)
        mergedFileDict.update(dict(results))
    
    SaveDictAsJSON(mergedFileName, mergedFileDict)

def MergeListFiles(mergedFileName, *fileNames): #can take any number of file names
    mergedFileList = []
    for fileName in fileNames:
        results = ReadJSONFile(fileName)
        print(type(results))
        mergedFileList.extend(results)
    
    SaveDictAsJSON(mergedFileName, mergedFileList)


MergeDictFiles("NonRadSearch2", "NonRadSearch2_band_gap0", "NonRadSearch2_band_gap_gt0")
# ^ only works if both json files are in dictionary format. If they're in list format, an error is raised (and not an intuitive one).