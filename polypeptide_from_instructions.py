import chemspl
from lxml import etree as ET
from rdkit.Chem import rdchem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdChemReactions

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

def get_reactant_by_code(element_chem_map,prop_name, code):
    for obj in element_chem_map.values():
        if isinstance(obj, rdchem.Mol) and obj.GetProp(prop_name)==code:
            return obj
    return None
    #return [val for val in element_chem_map.values() if (isinstance(val, rdchem.Mol) and val.GetProp(prop_name)==code)][0]


et, element_chem_map = chemspl.read_spl("polypeptide2.xml")
polypeptide_sequence=get_polypeptide_sequence(element_chem_map)
rxn = [val for val in element_chem_map.values() if isinstance(val, rdChemReactions.ChemicalReaction)][0]
rxn.Initialize()
if len(polypeptide_sequence)>1:
    for i in range(len(polypeptide_sequence)):
        if i==0:
            products=rxn.RunReactants((get_reactant_by_code(element_chem_map,"One Letter Code",polypeptide_sequence[0]),get_reactant_by_code(element_chem_map,"One Letter Code",polypeptide_sequence[1])))
        elif i>1:
            products=rxn.RunReactants((products[0][0],get_reactant_by_code(element_chem_map,"One Letter Code",polypeptide_sequence[i])))

Chem.rdCoordGen.AddCoords(products[0][0])
print(Chem.MolToMolBlock(products[0][0]),file=open('from_polypeptide2.mol','w+'))

