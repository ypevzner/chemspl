from rdkit.Chem import AllChem as Chem
from lxml import etree as ET
from lxml.builder import ElementMaker # not E, you need the ElementMaker where you set a namespace and nsmap!
import uuid
from datetime import datetime
import re

splns = "urn:hl7-org:v3"
nsmap = {"h": splns, "xsi": "http://www.w3.org/2001/XMLSchema-instance"}
nsmap_none = {None: splns, "xsi": "http://www.w3.org/2001/XMLSchema-instance"}
xsiType = ET.QName("xsi", "type")
xsiSchemaLocation = ET.QName("xsi", "schemaLocation")
E = ElementMaker(namespace=splns, nsmap=nsmap_none)
elementChemMap = {} # a map of Element objects to Chem objects

# builds RDKit reaction and substance objects found in a reaction SPL file
# There may be multiple reaction objects, so we
# - find all chemical things in the file 
# - generate substance or reaction objects from them
# - associate them with the document element
def read_spl(spl_file): 
    elementChemMap = {}
    elementTree = ET.parse(spl_file, parser=ET.XMLParser(encoding="UTF-8"))
    rootElement = elementTree.getroot()
    # go over all the reactions and interactors, but you can instead go over all substances first
    # and create the required object and associate them on the elementChem map,
    for substanceElement in rootElement.xpath(".//h:identifiedSubstance/h:identifiedSubstance",namespaces=nsmap):
        molfileText = substanceElement.xpath("./h:moiety/h:subjectOf/h:characteristic/h:value[@mediaType='application/x-mdl-molfile']", namespaces=nsmap)[0].text
        substanceObject = Chem.MolFromMolBlock(molfileText)
        try:
            for code in substanceElement.xpath("./h:code",namespaces=nsmap):
                print(code)
                substanceObject.SetProp(code.attrib["displayName"],code.attrib["code"])
        except:
            pass
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

# returns a molecule with canonically renumbered atom numbers
def canonicalize(substance):
    frags = Chem.GetMolFrags(substance,asMols=True)
    canonicalized_substance = Chem.Mol()
    #in case of multi-fragment molecules, canonicalize each fragment individually then combine back into a single mol
    for frag in frags:
        #No need to reorder fragment atoms if fragment only consists of one atom
        if frag.GetNumAtoms() ==1:
            reordered_frag = frag
        else:
            if Chem.MolToSmiles(frag)!='[H]O[H]':
                frag_atom_order = get_inchi_atom_order(frag)
                reordered_frag = Chem.RenumberAtoms(frag,frag_atom_order)
            else:
                reordered_frag=frag
        canonicalized_substance = Chem.CombineMols(canonicalized_substance,reordered_frag)  
    return canonicalized_substance

def createSubstanceElement_old(substanceObject):
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

def createSubstanceElement(substanceObject):
    identified_substance_element = E("identifiedSubstance")
    for prop_name in substanceObject.GetPropNames():    
        identified_substance_element.append(E("code",displayName=prop_name, codeSystem="1.3.6.1.4.1.32366.1.1.17",code=substanceObject.GetProp(prop_name)))
    
    substanceObject = canonicalize(substanceObject) # we just do the canonical ordering and now this is our substanceObject
    identified_substance_element.append(E("moiety",
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

    return E("identifiedSubstance", # there is going to be a lot more stuff in the future, like ids, etc.
              identified_substance_element
            )

def createReactionElement(reactionObject):                             
    reactionElement = E("processStep")
    priorityNumber = 0
    for i in range(reactionObject.GetNumReactantTemplates()):
        substanceObject = reactionObject.GetReactantTemplate(i)
        reactionElement.append(
            E("interactor", 
                E("priorityNumber", value=str(priorityNumber)),
                E("functionCode",code="reactant",codeSystem="1.3.6.1.4.1.32366.1.1"),
                createSubstanceElement(substanceObject),
              typeCode="CSM"))
        priorityNumber+=1

    for i in range(reactionObject.GetNumProductTemplates()):
        substanceObject = reactionObject.GetProductTemplate(i) #XXX: I don't think we need to subtract the number of reactants here, else we go in with the priorityNumber directly
        reactionElement.append(
            E("interactor", 
                E("priorityNumber", value=str(priorityNumber)),
                E("functionCode",code="product",codeSystem="1.3.6.1.4.1.32366.1.1"),
                createSubstanceElement(substanceObject),
              typeCode="PRD"))
        priorityNumber+=1
    return reactionElement


# this creates a complete SPL with only a single reaction object
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
                E("id", root=str(uuid.uuid4())),
                E("code"),
                E("title"),
                E("text"),
                E("effectiveTime", value=datetime.now().strftime("%Y%m%d%H%M%S")),
                  E("subject2",
                    E("specification",
                      E("code", nullFlavor="NI"),
                      E("component",
                        E("processStep",
                          E("component",
                          createReactionElement(reactionObject)))))))))))
    return root

def createNewSPLFromMultipleReactions(reactionObjects):
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
    reactions_component=E("component")
    for reactionObject in reactionObjects:
      reactions_component.append(createReactionElement(reactionObject))

    root.append(E("component",
      E("structuredBody",
        E("component",
          E("section",
            E("id", root=str(uuid.uuid4())),
            E("code"),
            E("title"),
            E("text"),
            E("effectiveTime", value=datetime.now().strftime("%Y%m%d%H%M%S")),
              E("subject2",
                E("specification",
                  E("code", nullFlavor="NI"),
                    E("component",
                      E("processStep",
                        reactions_component)))))))))
    return root

# this creates a complete SPL with only a single chemical substance object, this isn't very useful but it is as far as it has been done
def createNewSPLFromSingleSubstance(molObject):
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
                E("id", root=str(uuid.uuid4())),
                E("code"),
                E("title"),
                E("text"),
                E("effectiveTime", value=datetime.now().strftime("%Y%m%d%H%M%S")),
                  E("subject",
                    createSubstanceElement(molObject)))))))
    return root