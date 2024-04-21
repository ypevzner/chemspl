import chemspl
from lxml import etree as ET
from rdkit.Chem import rdchem
from rdkit.Chem import AllChem as Chem

suppl=Chem.SDMolSupplier("amino_acids.mol")
#for mol in suppl:
spl = chemspl.createNewSPLFromMultipleSubstances(suppl)

print(ET.tostring(spl,pretty_print=True))
spl_et = ET.ElementTree(spl)
spl_et.write("amino_acids.xml", pretty_print=True)