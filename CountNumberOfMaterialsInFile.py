from Util import *
# from pymatgen.ext.matproj import MPRester
# from pymatgen.core.periodic_table import Element
# import numpy as np

# def ListOfTheElements(elementsExcluded=None):
#     noOfElements = 118 #as of 2021
#     atomicNos = np.arange(1, noOfElements+1) #the function stops just before the given value (default step = 1)
    
#     if(elementsExcluded != None):
#         atomicNos = [z for z in atomicNos if z not in elementsExcluded]
    
#     symbolsTypeElement = [Element.from_Z(z) for z in atomicNos]
#     symbols = [str(symbol) for symbol in symbolsTypeElement]
    
#     return symbols

# elementList = ListOfTheElements()

# APIkey = None #done so that APIkey is not lost in the scope of the with block
# with open("APIkey.txt", "r") as f:
#     APIkey= f.read()

# results = None #done so that results exists outside the scope of the with block
# with MPRester(APIkey) as mpr:

#     criteria = {"elements": {"$in": elementList}}
#     # ^ want to find materials that contain any of the listed elements, hence $in, $nin excludes elements in given list,
#     # and $gt is simply 'greater than' - ferroelectrics are insulators, and DFT underestimates band gaps greatly,
#     # so if the band gap is > 0, that means the band gap is sizeable (likely insulator). NEW, RUN THIS SOON

#     properties = ['material_id', 'pretty_formula']
#     results = mpr.query(criteria, properties, chunk_size=10000) #it's ok to change chunk_size since I'm asking for small things -
#                                                                 #it speeds things up tremendously (asking for larger amounts of
#                                                                 #data less often - latency - sending info back and forth takes
#                                                                 # time)
# print(len(results))

origResults = ReadJSONFile("NonRadSearch2")
print(len(origResults))

CDResults = ReadJSONFile("NonRadSearch2CD_0")
print(len(list(CDResults.keys())))

NPResults = ReadJSONFile("NonRadSearch2NP_1")
print(len(list(NPResults.keys())))