## ChemSPL

  

ChemSPL is a python module for interconverting between SPL files and RDKit chemistry objects such as chemical substances and reactions

SPL (Structured Product Labeling) file represents a document markup standard approved by Health Level Seven (HL7) and adopted by FDA as a mechanism for exchanging product and facility information.
SPL file is in the XML format.

#### Main Methods:
#### _read_spl(spl_filename)_ 
Reads a specified SPL file.
Returns an ElementTree object (from Python's lxml library) that contains the entire SPL document, as well as a dictionary representing a map between main moieties in the SPL file (substances, reactions, etc.) and their corresponding RDKit objects in the order they appear in the input SPL file.

#### _createNewSPLFromSingleReaction(reaction_object)_ 
Given an RDKit chemical reaction object, returns an ElementTree object containing the entire SPL document that describes the input reaction.

#### _createNewSPLFromMultipleReactions(list of reaction_objects)_ 
Given a list of RDKit chemical reaction objects, returns an ElementTree object containing the entire SPL document that describes the input multi-reaction process