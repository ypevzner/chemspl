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

#### _createNewSPLFromSingleReaction(rdkit mol_object)_ 
Given an RDKit Mol object, returns an ElementTree object containing the entire SPL document that describes the input molecule.

#### _createNewSPLFromMultipleReactions(list of rdkit mol_objects)_ 
Given a list of RDKit Mol objects, returns an ElementTree object containing the entire SPL document that describes the input multi-Mol set

#### _createNewSPLFromMixedObjects(rdkit objects)_ 
Given a mixed list of RDKit Mol and reaction objects, returns an ElementTree object containing the entire SPL document that describes the input Mol and reaction objects

#### _appendSubstanceCharacteristics(ElementTree subjectOf_element, rdkit object)_ 
Given an input XML element, returns an ElementTree object containing the input ElementTree but with addition of molfile, InChI and InChiKey characteristics

#### Example: create_substances_spl.py
This script takes as input a molfile containing all 20 standard amino acids with their one letter codes (amino_acids.mol) and writes out an SPL file containing the amino acids as substances (amino_acids.xml).

#### Example: all_amino_acids_chain_assembly.py
The script reads in three SPL files. First - substance describing a polypeptide consisting of all 20 standard amino acids as defined by their one letter codes (all_amino_acid_chain.xml). Second - definition of the 20 standard amino acids (amino_acids.xml). Third - the condensation reaction definition that describes how to attach amino acids (condensation_reaction.xml). Then the script performs the condensation reaction according to the SPL instructions using the RDKit reaction methods according to the input polypeptide definition to stich together the 20 peptide chain. Then a molfile, InChi and InChIKey are generated for the resulting molecule and, together with the input polypeptide definition, written out to output file all_amino_acids_chain_with_molfile.xml
