import chemspl

rxn = chemspl.read_reaction_spl("reaction-2875_spl.xml")
spl = chemspl.make_reaction_spl(rxn)
with open("reaction-2875_spl1.xml",'w') as f:
    f.write(spl)

rxn = chemspl.read_reaction_spl("reaction-2875_spl1.xml")
spl = chemspl.make_reaction_spl(rxn)
with open("reaction-2875_spl2.xml",'w') as f:
    f.write(spl)