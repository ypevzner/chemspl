from rdkit.Chem import AllChem as Chem
import rdkit
from lxml import etree as ET
from lxml import objectify
from lxml.builder import E 
from lxml.builder import ElementMaker
import uuid
import numpy as np
from datetime import datetime
import re

#builds RDKit reaction object from reaction SPL file
def read_reaction_spl(spl_file):

    et=ET.parse(spl_file,parser=ET.XMLParser(encoding="UTF-8"))
    root=et.getroot()

    #remove namespace so I don't have to bother with it for now
    for elem in root.getiterator():
        if not (
            isinstance(elem, ET._Comment)
            or isinstance(elem, ET._ProcessingInstruction)
        ):
            # Remove a namespace URI in the element's name
            elem.tag = ET.QName(elem).localname
            
    root_xml=ET.tostring(et)
    for reaction in root.xpath("./component/structuredBody/component/section/subject2/specification/component/processStep/component/processStep"):
        interactors=reaction.xpath("./interactor")
        
    reactants_molfiles=[]
    products_molfiles=[]
    for interactor in interactors:
        interactor_molfile = interactor.xpath("./identifiedSubstance/identifiedSubstance/moiety/subjectOf/characteristic/value[@mediaType='application/x-mdl-molfile']")[0].text
        if interactor.get("typeCode") == "CSM":
            reactants_molfiles.append(interactor_molfile)
        elif interactor.get("typeCode") == "PRD":
            products_molfiles.append(interactor_molfile)

    rxn = Chem.ChemicalReaction()
    for reactant_molfile in reactants_molfiles:
        rxn.AddReactantTemplate(Chem.MolFromMolBlock(reactant_molfile))

    for product_molfile in products_molfiles:
        rxn.AddProductTemplate(Chem.MolFromMolBlock(product_molfile))

    return rxn


urn="urn:hl7-org:v3"
xsi = "http://www.w3.org/2001/XMLSchema-instance" 
schemaLocation="urn:hl7-org:v3 https://www.accessdata.fda.gov/spl/schema/spl.xsd"

def numbers_only(val):
    result = re.sub('[^0-9]','', val)
    return result

def make_zero_based_index(val):
    return int(numbers_only(val))-1

#returns a tuple representing atoms in an inchi canonical order
def get_inchi_atom_order(mol):
    inchi_aux_info=Chem.inchi.MolToInchiAndAuxInfo(mol)[1]
    atom_order_string = re.search("/N:([0-9,;]+)",inchi_aux_info).group()
    return tuple([make_zero_based_index(i) for i in atom_order_string.split(",")])

def create_interactor_spl(reaction, type_code, index_in_reaction,code):
    if code == "reactant":
        template = reaction.GetReactantTemplate(index_in_reaction)
    elif code == "product":
        template = reaction.GetProductTemplate(index_in_reaction-reaction.GetNumReactantTemplates())

    canonically_ordered_template = Chem.RenumberAtoms(template, get_inchi_atom_order(template))
    nsmap={None: urn, 'xsi': xsi}
    qname=ET.QName(xsi,"type")
    
    molfile_value_element = ET.Element('value',{qname: "ED"})
    molfile_value_element.set("mediaType","application/x-mdl-molfile")
    molfile_value_element.text = ET.CDATA(Chem.MolToMolBlock(canonically_ordered_template))

    inchi_value_element = ET.Element('value',{qname: "ED"})
    inchi_value_element.set("mediaType","application/x-inchi")
    inchi_value_element.text = Chem.MolToInchi(canonically_ordered_template)

    inchikey_value_element = ET.Element('value',{qname: "ED"})
    inchikey_value_element.set("mediaType","application/x-inchi-key")
    inchikey_value_element.text = Chem.MolToInchiKey(canonically_ordered_template)

    return E("interactor",E("priorityNumber",value=str(index_in_reaction)),
                                  E("functionCode",code=code,codeSystem="1.3.6.1.4.1.32366.1.1"),
                                  E("identifiedSubstance",
                                    E("identifiedSubstance",
                                      E("moiety",
                                        E("partMoiety"),
                                        E("subjectOf",
                                          E("characteristic",
                                            E("code",displayName="Chemical Structure",codeSystem="2.16.840.1.113883.3.26.1.1",code="C103240"),
                                            molfile_value_element
                                          )
                                        ),
                                        E("subjectOf",
                                          E("characteristic",
                                            E("code",displayName="Chemical Structure",codeSystem="2.16.840.1.113883.3.26.1.1",code="C103240"),
                                            inchi_value_element
                                          )
                                        ),
                                        E("subjectOf",
                                          E("characteristic",
                                            E("code",displayName="Chemical Structure",codeSystem="2.16.840.1.113883.3.26.1.1",code="C103240"),
                                            inchikey_value_element
                                          )
                                        )))),
                                  typeCode=type_code)

#creates a string representing an entire SPL document for a given reaction
def make_reaction_spl(rxn):
    interactors=[]
    for interactor_index in range(rxn.GetNumReactantTemplates()):
        interactors.append(create_interactor_spl(rxn, "CSM", interactor_index,"reactant"))

    for interactor_index in range(rxn.GetNumReactantTemplates(), rxn.GetNumReactantTemplates()+rxn.GetNumProductTemplates()):
        interactors.append(create_interactor_spl(rxn, "PRD", interactor_index,"product"))

    reaction_spl = E("component",E("processStep",*interactors))

    document_begin = """<?xml version="1.0" encoding="UTF-8"?>
    <?xml-stylesheet href="https://www.accessdata.fda.gov/spl/stylesheet/spl.xsl" type="text/xsl"?>
    <document xmlns="urn:hl7-org:v3" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="urn:hl7-org:v3 https://www.accessdata.fda.gov/spl/schema/spl.xsd">
    <id root=\"""" + str(uuid.uuid4()) + """\"/>
    <effectiveTime value=\"""" + datetime.now().strftime("%Y%m%d%H%M%S") + """\"/>
    <setIf root=\"""" + str(uuid.uuid4()) + """\"/>
    <versionNumber value="1"/>
    <component>
        <structuredBody>
        <component>
            <section>
            <id root=\"""" + str(uuid.uuid4()) + """\"/>
            <code/>
            <title/>
            <text/>
                <effectiveTime value="20210928133013"/>
                <subject2>
                <specification>
                    <code nullFlavor="NI"/>
                    <component>
                    <processStep>
                        """

    document_end = """                </processStep>
                </component>
                </specification>
            </subject2>
            </section>
        </component>
        </structuredBody>
    </component>
    </document>"""
    ET.cleanup_namespaces(reaction_spl)
    ET.indent(reaction_spl, space="  ", level=10)
    complete_spl_document = document_begin + ET.tostring(reaction_spl,pretty_print=True).decode().replace('value xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"',"value") + document_end

    return complete_spl_document