import chemspl
from lxml import etree as ET
from rdkit.Chem import rdchem
from rdkit.Chem import AllChem as Chem

# m = Chem.MolFromMolFile('alanine.mol')
# spl = chemspl.createNewSPLFromSingleSubstance(m)
# spl_et = ET.ElementTree(spl)
# spl_et.write("alanine.xml", pretty_print=True)


et, element_chem_map = chemspl.read_spl("alanine.xml")


spl = chemspl.createNewSPLFromSingleSubstance([val for val in element_chem_map.values() if isinstance(val, rdchem.Mol)][0])
print(ET.tostring(spl,pretty_print=True))
spl_et = ET.ElementTree(spl)
spl_et.write("alanine_spl1.xml", pretty_print=True)

et, element_chem_map = chemspl.read_spl("alanine_spl1.xml")


spl = chemspl.createNewSPLFromSingleSubstance([val for val in element_chem_map.values() if isinstance(val, rdchem.Mol)][0])
print(ET.tostring(spl,pretty_print=True))
spl_et = ET.ElementTree(spl)
spl_et.write("alanine_spl2.xml", pretty_print=True)