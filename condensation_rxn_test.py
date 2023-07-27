from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdChemReactions
import chemspl
from lxml import etree as ET

#create reaction objecdt from file
rxn_from_rxn = Chem.ReactionFromRxnFile("condencation_reaction.rxn")
rxn_from_rxn.Initialize()

#create spl document and save it as file
spl = chemspl.createNewSPLFromSingleReaction(rxn_from_rxn)
print(ET.tostring(spl,pretty_print=True))
spl_et = ET.ElementTree(spl)
spl_et.write("condencation_reaction_spl.xml", pretty_print=True)

#read the spl that was jus written
et, element_chem_map = chemspl.read_spl("condencation_reaction_spl.xml")
rxn = [val for val in element_chem_map.values() if isinstance(val, rdChemReactions.ChemicalReaction)][0]
rxn.Initialize()
aa1 = Chem.MolFromMolFile('lysine.mol')

#print(Chem.MolToMolBlock(Chem.rdmolops.AddHs(aa1)),file=open('lysine_hs.mol','w+'))
aa2 = Chem.MolFromMolFile('alanine.mol')
aa3 = Chem.MolFromMolFile('glycine.mol')
#Chem.rdmolops.RemoveHs(
print(Chem.MolToSmiles(aa1))
products=rxn.RunReactants((aa1,aa2))
print(len(products))
print(len(products[0]))
print(Chem.MolToSmiles(products[0][0]))
print(Chem.MolToSmiles(products[0][1]))
products_step2=rxn.RunReactants((products[0][0],aa3))

Chem.rdCoordGen.AddCoords(products_step2[0][0])
print(Chem.MolToMolBlock(products_step2[0][0]),file=open('3AA.mol','w+'))