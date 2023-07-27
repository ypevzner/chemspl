import chemspl
from rdkit.Chem import rdChemReactions
from lxml import etree as ET

et, element_chem_map = chemspl.read_spl("reaction-isoxanzolines-2_yp.xml")

reaction_objects = [val for val in element_chem_map.values() if isinstance(val, rdChemReactions.ChemicalReaction)]
spl = chemspl.createNewSPLFromMultipleReactions(reaction_objects)
#spl = chemspl.createNewSPLFromMultipleReactions(list(element_chem_map.values())[-2:])
print(ET.tostring(spl,pretty_print=True))
spl_et = ET.ElementTree(spl)
spl_et.write("reaction-isoxanzolines-2_yp_1.xml", pretty_print=True)

et, element_chem_map = chemspl.read_spl("reaction-isoxanzolines-2_yp_1.xml")


reaction_objects = [val for val in element_chem_map.values() if isinstance(val, rdChemReactions.ChemicalReaction)]
spl = chemspl.createNewSPLFromMultipleReactions(reaction_objects)
print(ET.tostring(spl,pretty_print=True))
spl_et = ET.ElementTree(spl)
spl_et.write("reaction-isoxanzolines-2_yp_2.xml", pretty_print=True)