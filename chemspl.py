from rdkit.Chem import AllChem as Chem
from lxml import etree as ET
from lxml.builder import ElementMaker # not E, you need the ElementMaker where you set a namespace and nsmap!
import uuid
from datetime import datetime
import re

splns = "urn:hl7-org:v3"
#nsmap = {None: splns, "xsi": "http://www.w3.org/2001/XMLSchema-instance"}
nsmap = {"h": splns, "xsi": "http://www.w3.org/2001/XMLSchema-instance"}
#nsmap = {"": splns, "xsi": "http://www.w3.org/2001/XMLSchema-instance"}
xsiType = ET.QName("xsi", "type") #XXX: this was called "qname" but you are misleading the reader, it's not any qName, it's the name xsi:type
xsiSchemaLocation = ET.QName("xsi", "schemaLocation")
E = ElementMaker(namespace=splns, nsmap=nsmap) #Note: see now your E is operating in the correct namespace

elementChemMap = {} # a map of Element objects to Chem objects

#builds RDKit reaction object from reaction SPL file
#FIXME: what if there are multiple reaction objects?
#FIXME: what if we have the SPL file with only substances, not reactions?
#We should just go through the document as is, and associate the rdkit.Chem 
#This means this function is not "read_reaction" SPL, but rather 
# - find all chemical things in the file 
# - and generate substance or reaction objects from them
# - and associate them with the document element
#So, instead of calling this method:
#def read_reaction_spl(spl_file): we should say:
def read_spl(spl_file): 
    # and we would simply return the ET's result while also updating the elementChemMap with the Chem objects.

    elementTree = ET.parse(spl_file, parser=ET.XMLParser(encoding="UTF-8"))
    rootElement = elementTree.getroot()
    #elementTree.write("test1.xml")
    #remove namespace so I don't have to bother with it for now
    #FIXME: this should not be necessary, this is a hack which ultimately will fail.
    #just use an nsmap {None: "urn:hl7-org:v3"} as the second argument of the XPath function
    #for elem in root.getiterator():
    #    if not (
    #        isinstance(elem, ET._Comment)
    #        or isinstance(elem, ET._ProcessingInstruction)
    #    ):
    #        # Remove a namespace URI in the element's name
    #        elem.tag = ET.QName(elem).localname
    #just delete the above, and use xpath("...", nsmap)
    
    #root_xml=ET.tostring(et) #FIXME: unused variable, delete it
   
    # now you go over all the reactions and interactors, but you can instead go over all substances first
    # and create the required object and associate them on the elementChem map,
    # and then you go over all reactions and create the reaction objects.
    #print(nsmap)
    for substanceElement in rootElement.xpath(".//h:identifiedSubstance/h:identifiedSubstance",namespaces=nsmap):
        molfileText = substanceElement.xpath("./h:moiety/h:subjectOf/h:characteristic/h:value[@mediaType='application/x-mdl-molfile']", namespaces=nsmap)[0].text
        substanceObject = Chem.MolFromMolBlock(molfileText)
        elementChemMap[substanceElement] = substanceObject
    
    # now we can do the reactions (if any)
    #YP debug
    #xpath_find = rootElement.xpath(".//processStep/component/processStep", namespaces=nsmap)
    #print(ET.tostring(xpath_find, encoding='utf8').decode('utf8'))
    for reactionElement in rootElement.xpath(".//h:processStep/h:component/h:processStep", namespaces=nsmap):
        #XXX: when you called this "rxn" does that mean this is an RXN file structure or just a cute abbreviation?
        #XXX: we should always use full words that mean exactly what they say, so replace "rxn" with "reaction" eveywhere
        reactionObject = Chem.ChemicalReaction()
        for reactantElement in reactionElement.xpath("./h:interactor[@typeCode = 'CSM']/h:identifiedSubstance/h:identifiedSubstance", namespaces=nsmap):
            reactionObject.AddReactantTemplate(elementChemMap[reactantElement])
        for productElement in reactionElement.xpath("./h:interactor[@typeCode = 'PRD']/h:identifiedSubstance/h:identifiedSubstance", namespaces=nsmap):
            reactionObject.AddProductTemplate(elementChemMap[productElement])
        elementChemMap[reactionElement] = reactionObject

    return elementTree, elementChemMap

# now after this read_spl function executes, you have two things
# 1. an elementTree
# 2. the elementChemMap populated with rdkit Chem objects 
# you can iterate over the map to find those Chem objects if you don't care about the document structure
# or you can always go back to the document structure and find any Chem object for any element that represents one

# To generate SPL documents, there is certainly more to it than just making such simple reaction files.

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

#returns get a molecule with canonically renumbered atom numbers
def canonicalize(template): #XXX: I am not sure what a "template" is different from the actual substance object, if it's a substance object, we should call it that
    return Chem.RenumberAtoms(template, get_inchi_atom_order(template))

#XXX: this is not a great name, you want to have this function do one thing: create substance representations
#XXX: whether they are interactors in reactions doesn't matter, nor how you access those by "index_in_reaction"
#def create_interactor_spl(reaction, type_code, index_in_reaction,code):
#    if code == "reactant":
#       template = reaction.GetReactantTemplate(index_in_reaction)
#   elif code == "product":
#       template = reaction.GetProductTemplate(index_in_reaction-reaction.GetNumReactantTemplates())
#XXX: I am not sure what a "template" is different from the actual substance object, if it's a substanceObject, we should call it that

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
# see how clean this now looks? 

#creates a string representing an entire SPL document for a given reaction
#def make_reaction_spl(rxn):
#XXX: no, that would be wrong, a document may have multiple reactions, do each reaction at a time
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
    # root.append(
    #     E("id", root=str(documentId)),
    #     E("effectiveTime", value=datetime.now().strftime("%Y%m%d%H%M%S")),
    #     E("setId", root=str(documentSetId)),
    #     E("versionNumber", value=str(versionNumber)),
    #     E("component",
    #       E("structuredBody",
    #         E("component",
    #           E("section",
    #             E("id", str(uuid.uuid4())),
    #             E("code"),
    #             E("title"),
    #             E("text"),
    #             E("effectiveTime", value=datetime.now().strftime("%Y%m%d%H%M%S")),
    #               E("subject2",
    #                 E("specification",
    #                   E("code", nullFlavor="NI"),
    #                   E("component",
    #                     createReactionElement(reactionObject)))))))))
    return root