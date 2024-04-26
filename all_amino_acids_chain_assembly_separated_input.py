import chemspl as chemspl
from lxml import etree as ET
from rdkit.Chem import rdchem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdChemReactions
import os

splns = "urn:hl7-org:v3"
nsmap = {"h": splns, "xsi": "http://www.w3.org/2001/XMLSchema-instance"}


def get_polypeptide_sequence(element_chem_map):
    seq=""
    try:
        for element in element_chem_map:
            seq=element.xpath("./h:moiety/h:subjectOf/h:characteristic/h:value[@mediaType='application/x-aa-seq']", namespaces=nsmap)[0].text
    except:
        pass
    return seq

def get_reactant_by_code(element_chem_map,code_system, code):
    for obj in element_chem_map.values():
        if isinstance(obj, rdchem.Mol) and obj.GetProp('code:' + code_system)==code:
            return obj
    return None

def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

output_objects=[]
polypeptide_et, polypeptide_element_chem_map = chemspl.read_spl("all_amino_acids_chain.xml")
amino_acids_et, amino_acids_element_chem_map = chemspl.read_spl("amino_acids.xml")
reaction_et, reaction_element_chem_map = chemspl.read_spl("condensation_reaction.xml")
polypeptide_sequence=get_polypeptide_sequence(polypeptide_element_chem_map)
rxn = [val for val in reaction_element_chem_map.values() if isinstance(val, rdChemReactions.ChemicalReaction)][0]
rxn.Initialize()
if len(polypeptide_sequence)>1:
    for i in range(len(polypeptide_sequence)):
        if i==0:
            products=rxn.RunReactants((get_reactant_by_code(amino_acids_element_chem_map,"1.3.6.1.4.1.32366.1.1.17",polypeptide_sequence[0]),get_reactant_by_code(amino_acids_element_chem_map,"1.3.6.1.4.1.32366.1.1.17",polypeptide_sequence[1])))
        elif i>1:
            products=rxn.RunReactants((products[0][0],get_reactant_by_code(amino_acids_element_chem_map,"1.3.6.1.4.1.32366.1.1.17",polypeptide_sequence[i])))

Chem.rdCoordGen.AddCoords(products[0][0])

input_tree = ET.parse("all_amino_acids_chain.xml")
moiety = input_tree.xpath("//h:identifiedSubstance/h:identifiedSubstance/h:moiety", namespaces=nsmap)[0]
moiety = chemspl.appendSubstanceCharacteristics(moiety,products[0][0])

# this is to ensure indentation after appending new elements
output_root=ET.fromstring(ET.tostring(input_tree))
indent(output_root)

output_tree = ET.ElementTree(output_root)
output_tree.write("all_amino_acids_chain_with_molfile.xml", pretty_print=True)



