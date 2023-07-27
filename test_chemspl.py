import chemspl
from lxml import etree as ET
from rdkit.Chem import rdChemReactions

et, element_chem_map = chemspl.read_spl("reaction-2875_spl.xml")


spl = chemspl.createNewSPLFromSingleReaction([val for val in element_chem_map.values() if isinstance(val, rdChemReactions.ChemicalReaction)][0])
print(ET.tostring(spl,pretty_print=True))
spl_et = ET.ElementTree(spl)
spl_et.write("reaction-2875_spl1.xml", pretty_print=True)

et, element_chem_map = chemspl.read_spl("reaction-2875_spl1.xml")


spl = chemspl.createNewSPLFromSingleReaction([val for val in element_chem_map.values() if isinstance(val, rdChemReactions.ChemicalReaction)][0])
print(ET.tostring(spl,pretty_print=True))
spl_et = ET.ElementTree(spl)
spl_et.write("reaction-2875_spl2.xml", pretty_print=True)