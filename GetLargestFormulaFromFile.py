from Util import *
from pymatgen.core.composition import Composition

print("Starting program.")
CDResults = ReadJSONFile("NonRadSearch2CD_0")
print("File loaded.")
noOfMaterials = len(list(CDResults.keys())) #maybe loading keys is faster than loading values? Probably won't matter for a one-time call.

noAtomsList = []

for i, material in enumerate(CDResults.values()):
    comp = Composition(material["pretty_formula"])
    noAtomsList.append(comp.num_atoms)

noAtomsList.sort() #just sorts the list, I guess - can't even save to a variable, it just modifies the original list
largestMaterialSize = noAtomsList.pop(-1)
print(largestMaterialSize)