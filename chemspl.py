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

<<<<<<< Updated upstream
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

=======
splns = "urn:hl7-org:v3"
nsmap = {"h": splns, "xsi": "http://www.w3.org/2001/XMLSchema-instance"}
xsiType = ET.QName("xsi", "type")
xsiSchemaLocation = ET.QName("xsi", "schemaLocation")
E = ElementMaker(namespace=splns, nsmap=nsmap)
elementChemMap = {} # a map of Element objects to Chem objects

# builds RDKit reaction and substance objects found in a reaction SPL file
# There may be multiple reaction objects, so we
# - find all chemical things in the file 
# - generate substance or reaction objects from them
# - associate them with the document element
def read_spl(spl_file): 
    elementTree = ET.parse(spl_file, parser=ET.XMLParser(encoding="UTF-8"))
    rootElement = elementTree.getroot()
    # go over all the reactions and interactors, but you can instead go over all substances first
    # and create the required object and associate them on the elementChem map,
    for substanceElement in rootElement.xpath(".//h:identifiedSubstance/h:identifiedSubstance",namespaces=nsmap):
        molfileText = substanceElement.xpath("./h:moiety/h:subjectOf/h:characteristic/h:value[@mediaType='application/x-mdl-molfile']", namespaces=nsmap)[0].text
        substanceObject = Chem.MolFromMolBlock(molfileText)
        elementChemMap[substanceElement] = substanceObject
    # now go over all reactions and create the reaction objects.
    for reactionElement in rootElement.xpath(".//h:processStep/h:component/h:processStep", namespaces=nsmap):
        reactionObject = Chem.ChemicalReaction()
        for reactantElement in reactionElement.xpath("./h:interactor[@typeCode = 'CSM']/h:identifiedSubstance/h:identifiedSubstance", namespaces=nsmap):
            reactionObject.AddReactantTemplate(elementChemMap[reactantElement])
        for productElement in reactionElement.xpath("./h:interactor[@typeCode = 'PRD']/h:identifiedSubstance/h:identifiedSubstance", namespaces=nsmap):
            reactionObject.AddProductTemplate(elementChemMap[productElement])
        elementChemMap[reactionElement] = reactionObject
    return elementTree, elementChemMap

# this read_spl function returns two things:
#  1. an elementTree
#  2. the elementChemMap populated with rdkit Chem objects 
# you can iterate over the map to find those Chem objects if you don't care about the document structure
# or you can always go back to the document structure and find any Chem object for any element that represents one

>>>>>>> Stashed changes
def numbers_only(val):
    result = re.sub('[^0-9]','', val)
    return result

def make_zero_based_index(val):
    return int(numbers_only(val))-1

# returns a tuple representing atoms in an inchi canonical order
def get_inchi_atom_order(mol):
    inchi_aux_info=Chem.inchi.MolToInchiAndAuxInfo(mol)[1]
    atom_order_string = re.search("/N:([0-9,;]+)",inchi_aux_info).group()
    return tuple([make_zero_based_index(i) for i in atom_order_string.split(",")])

<<<<<<< Updated upstream
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
=======
# returns a molecule with canonically renumbered atom numbers
def canonicalize(template): #XXX: I am not sure what a "template" is different from the actual substance object, if it's a substance object, we should call it that
    return Chem.RenumberAtoms(template, get_inchi_atom_order(template))

def createSubstanceElement(substanceObject):
    substanceObject = canonicalize(substanceObject) # we just do the canonical ordering and now this is our substanceObject
    return E("identifiedSubstance", # there is going to be a lot more stuff in the future, like ids, etc.
              E("identifiedSubstance",
               E("moiety",
                 E("partMoiety"),
                 E("subjectOf",
                   E("characteristic",
                     E("code", displayName="Chemical Structure", codeSystem="2.16.840.1.113883.3.26.1.1", code="C103240"),
                     E("value", {xsiType: "ED", "mediaType": "application/x-mdl-molfile"},
                     ET.CDATA(Chem.MolToMolBlock(substanceObject))))),
                 E("subjectOf",
                  E("characteristic",
                     E("code", displayName="Chemical Structure", codeSystem="2.16.840.1.113883.3.26.1.1", code="C103240"),
                     E("value", {xsiType: "ED", "mediaType": "application/x-inchi"},
                     Chem.MolToInchi(substanceObject)))),
                 E("subjectOf",
                   E("characteristic",
                     E("code",displayName="Chemical Structure",codeSystem="2.16.840.1.113883.3.26.1.1",code="C103240"),
                     E("value", {xsiType: "ED", "mediaType": "application/x-inchi-key"}, 
                     Chem.MolToInchiKey(substanceObject)
                      )
                    )
                  )
                )
              )
            )

def createReactionElement(reactionObject):                             
    reactionElement = E("processStep")
    priorityNumber = 1
    for i in range(reactionObject.GetNumReactantTemplates()):
        substanceObject = reactionObject.GetReactantTemplate(i)
        reactionElement.append(
            E("interactor", 
                E("priorityNumber", value=str(priorityNumber)),
                E("functionCode",code="reactant",codeSystem="1.3.6.1.4.1.32366.1.1"),
                createSubstanceElement(substanceObject)))
        priorityNumber+=1

    for i in range(reactionObject.GetNumProductTemplates()):
        substanceObject = reactionObject.GetProductTemplate(i) #XXX: I don't think we need to subtract the number of reactants here, else we go in with the priorityNumber directly
        reactionElement.append(
            E("interactor", 
                E("priorityNumber", value=str(priorityNumber)),
                E("functionCode",code="product",codeSystem="1.3.6.1.4.1.32366.1.1"),
                createSubstanceElement(substanceObject)))
        priorityNumber+=1
    return reactionElement                        
        
# this creates a complete SPL with only a single reaction object, this isn't very useful but it is as far as it has been done
def createNewSPLFromSingleReaction(reactionObject):
    # IMPORTANT: don't manually hack XML with string concatenation, because this doesn't take care of proper escaping etc.
    root = E("document", {xsiSchemaLocation: "urn:hl7-org:v3 https://www.accessdata.fda.gov/spl/schema/spl.xsd"})
    #root.addprevious(ET.ProcessingInstruction('xml-stylesheet href="https://www.accessdata.fda.gov/spl/stylesheet/spl.xsl" type="text/xsl"'))
    root.addprevious(ET.ProcessingInstruction('xml-stylesheet', text='href="https://www.accessdata.fda.gov/spl/stylesheet/spl.xsl" type="text/xsl"'))    
    documentId = uuid.uuid4() 
    documentSetId = documentId
    versionNumber = 1        
    root.append(E("id", root=str(documentId)))
    root.append(E("effectiveTime", value=datetime.now().strftime("%Y%m%d%H%M%S")))
    root.append(E("setId", root=str(documentSetId)))
    root.append(E("versionNumber", value=str(versionNumber)))
    root.append(E("component",
          E("structuredBody",
            E("component",
              E("section",
                E("id", str(uuid.uuid4())),
                E("code"),
                E("title"),
                E("text"),
                E("effectiveTime", value=datetime.now().strftime("%Y%m%d%H%M%S")),
                  E("subject2",
                    E("specification",
                      E("code", nullFlavor="NI"),
                      E("component",
                        createReactionElement(reactionObject)))))))))
    return root
>>>>>>> Stashed changes
