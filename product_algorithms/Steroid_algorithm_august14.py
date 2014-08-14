# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Process book for the steroid algorithm code.

# <markdowncell>

# Let's get all the relevant smiles for the proteins we're going to use. 

# <codecell>

import sys                           
sys.path.append('/Home')

from SteroidFunctions import *

# <headingcell level=1>

# New Anabolic Steroids

# <codecell>

TestRxns = pd.read_table("Keggsteroids.csv")
TestRxns.head(11)

# <codecell>

test = pcp.get_properties('CanonicalSMILES', 'Cholesterol', 'name')
print str(test[0]['CanonicalSMILES'])
Chem.MolFromSmiles( str(test[0]['CanonicalSMILES']) )

# <codecell>

start = timeit.default_timer()

anabolic = 'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C'
substrates = []
products = []
targets = []
results = []
enzymes = []
for i in range(len(TestRxns)):
    try:
        results.append( combinatorial_rxn(anabolic, TestRxns['Substrates'][i], TestRxns['Products'][i]) )
        substrates.append( TestRxns['Substrates'][i] )
        products.append( TestRxns['Products'][i] )
        targets.append( anabolic )
        enzymes.append( TestRxns['Enzymes'][i] )
    except:
        pass
    
stop = timeit.default_timer()   
unique_prods, unique_substrates, unique_products, unique_enzymes = filter_results(results, substrates, products, TestRxns, enzymes)

#search pubchem for unique hits
novel_compounds = []
for i in unique_prods:
    searches = pcp.get_compounds('CanonicalSMILES', str(Chem.MolToSmiles(i)), 'smiles')
    if str(searches) == '[Compound()]':
        novel_compounds.append( i )
        
#put things in a df
unique_prods_smiles = []
novel_smiles = []
for i in unique_prods:
    unique_prods_smiles.append( Chem.MolToSmiles(i) )
for i in novel_compounds:
    novel_smiles.append( Chem.MolToSmiles(i) )

novels =[ ]includeAtomNumbers=False
for i in unique_prods_smiles:
    if i in novel_smiles:
        novels.append( 'Novel' )
    else:
        novels.append( 'Found in Pubchem' )
Results = pd.DataFrame({'Native Substrate': unique_substrates, 'Native Product': unique_products, 'Products': unique_prods_smiles, 'Novel Compound?': novels, 'Enzymes': unique_enzymes})
        
print 'Novel Compounds found...'   
print len(novel_compounds)
print "Runtime..."
print stop - start    
Draw.MolsToGridImage(unique_prods,molsPerRow=4, includeAtomNumbers=False)

# <codecell>

len(list(set(Results['Products'].tolist())))

# <codecell>

#Results.to_csv("Desmosterol_one_iteration.csv")
#Results.to_csv("epoxycholesterol_one_iteration.csv")
#Results.to_csv('Lophenol_one_iteration.csv')
Results.to_csv("26-hydroxycholesterol_one_iteration.csv")

# <markdowncell>

# Code for two iterations.

# <codecell>

substrates = []
products = []
targets = []
results = []
enzyme1 = []
enzyme2 = []
#start = timeit.default_timer()
for i in range(len(Results)):
    for j in range(len(Results)):
        try:
            results.append( combinatorial_rxn( Results['Products'][j], Results['Native Substrate'][i], Results['Native Product'][i]) )
            substrates.append( Results['Native Substrate'][i] )
            products.append( Results['Native Product'][i] )
            targets.append( Results['Products'][j] )
            enzyme1.append( Results['Enzymes'][i] )
            enzyme2.append( Results['Enzymes'][j] )
        except:
            pass

# <codecell>

resultsmiles = []
for i in results:
    resultsmiles.append( Chem.MolToSmiles(i) )
It2 = pd.DataFrame({'Native Substrate': substrates, 'Native Product': products, 'Results': resultsmiles, 'Enzyme 1': enzyme1, 'Enzyme 2': enzyme2})
It2.head(10)

filtered_results = filter_results2(It2['Results'])


sf = []
pf = []
rsmf = []
ef1 = []
ef2 = []

for i in filtered_results:
    for j in range(len(It2)):
        if i == It2['Results'].irow(j):
            sf.append( It2['Native Substrate'].irow(j) )
            pf.append( It2['Native Product'].irow(j) )
            rsmf.append( It2['Results'].irow(j) )
            ef1.append( It2['Enzyme 1'].irow(j) )
            ef2.append( It2['Enzyme 2'].irow(j) )
It2F = pd.DataFrame({'Native Substrate': sf, 'Native Product': pf, 'Results': rsmf, 'Enzyme 1': ef1, 'Enzyme 2': ef2})

It2filtered = It2F.drop_duplicates(cols=['Results'])

#search pubchem for unique hits
novel_compounds = []
for i in It2filtered['Results']:
    try:
        searches = pcp.get_compounds('CanonicalSMILES', str(i), 'smiles')
    except:
        pass
    if str(searches) == '[Compound()]':
        novel_compounds.append( i ) 
        
novel_id = []
for i in It2filtered['Results']:
    if i in novel_compounds:
        novel_id.append( "Novel" )
    else:
        novel_id.append( "In Pubchem" )
It2filtered['Novelty'] = novel_id
       
print 'Novel Compounds found...'   
print len(novel_compounds)
print "Runtime..."
print stop - start

It2filtered.head(10)

# <codecell>

test = Results
vis = []
for i in test:
    if Chem.MolFromSmiles(i) != None:
        vis.append( Chem.MolFromSmiles( i ) )
Draw.MolsToGridImage(vis, molsPerRow=8, includeAtomNumbers=False)

# <codecell>

vis = []
for i in It2filtered['Results'].tolist():
    vis.append( Chem.MolFromSmiles(i) )

# <headingcell level=1>

# Troubleshooting

# <codecell>

#sz = Chem.MolFromSmiles( TestRxns['Substrates'].irow(0) )
#pz = Chem.MolFromSmiles( TestRxns['Products'].irow(0) )
#Draw.MolsToGridImage([t, sz, pz, metastructure],molsPerRow=4, includeAtomNumbers=False)
#product_modifications

# <codecell>

TestRxns = TestRxns[TestRxns['Enzymes'] == 'ARK1C3']
TestRxns

# <codecell>

def modify_substrate(substrate_modifications, steroid1, steroid2, steroids):
    m1 = Chem.MolFromSmiles(steroid1)
    m2 = Chem.MolFromSmiles(steroid2)
    ms = Chem.MolFromSmiles(steroids)
    ############################################################
    #removes atoms that are removed via lyase activity, first must find atoms that are removed from the native substrate to the product
    #then we have to compare those atoms to our non-native substrate then systematically remove them the tricky thing here will be indexing (as always)
    native_lexicon = create_lexicon(ms, m2)
    substrates_lexicon = create_lexicon(m1, ms)
    
    ms_atoms = []
    ms_m2_atoms = []
    m2matchingatoms = []
    for atom in ms.GetAtoms():
        ms_atoms.append( atom.GetIdx() )
    for i in native_lexicon:
        ms_m2_atoms.append( i[0] )
        m2matchingatoms.append( i[1] )
    #see which atoms don't have overlap i.e. things that need to be deleted
    ms_cleaved = []
    for i in ms_atoms:
        if i not in ms_m2_atoms:
            ms_cleaved.append( i )
            
    mp_atoms = []
    for atom in m2.GetAtoms():
        mp_atoms.append( atom )
     
    mp_unique_atoms = []
    for i in mp_atoms:
        if i.GetIdx() not in m2matchingatoms:
            mp_unique_atoms.append( i.GetIdx() )
    
    #aromatic atoms will screw this code up, we need to make sure the atoms we're going to delete are due to aromaticitiy         
    aromatic_check = []
    aromatic_idxi = []
    aromatic_idxj = []
    for i in mp_unique_atoms:
        for j in mp_unique_atoms:
            try:
                aromatic_check.append( str( m2.GetBondBetweenAtoms(i, j).GetIsAromatic() ) )
                aromatic_idxi.append( i )
                aromatic_idxj.append( j )
            except:
                pass

    #translate ms_cleaved to our target
    target_specifics = []
    for i in ms_cleaved:
        for j in substrates_lexicon:
            if i == j[1]:
                target_specifics.append( j[0] )
           
    if 'True' not in aromatic_check:   
        AROMATIC_FLAG = None
        temp_lex = create_lexicon( ms, m1 )
        em = Chem.EditableMol(m1)
        for i in range(len(ms_cleaved)):
            temp_lex = create_lexicon( ms, m1 )
            for j in temp_lex:
                if j[0] == ms_cleaved[i]:
                    print j[1]
                    deletion_atom = j[1]
                    em.RemoveAtom(deletion_atom)
                    m1 = em.GetMol()
        m1smiles = Chem.MolToSmiles( m1 )
        m1smiles = clean_smiles( m1smiles )
        try:
            m1 = Chem.MolFromSmiles( m1smiles )
        except:
            pass
    else:
        AROMATIC_FLAG = 'GO'
        
    ##########################################    
    native_lexicon = create_lexicon(ms, m2)
    substrates_lexicon = create_lexicon(m1, ms)
    
    ms_atoms = []
    ms_m2_atoms = []
    m2matchingatoms = []
    for atom in ms.GetAtoms():
        ms_atoms.append( atom.GetIdx() )
    for i in native_lexicon:
        ms_m2_atoms.append( i[0] )
        m2matchingatoms.append( i[1] )
    #see which atoms don't have overlap i.e. things that need to be deleted
    ms_cleaved = []
    for i in ms_atoms:
        if i not in ms_m2_atoms:
            ms_cleaved.append( i )
            
    mp_atoms = []
    for atom in m2.GetAtoms():
        mp_atoms.append( atom )
     
    mp_unique_atoms = []
    for i in mp_atoms:
        if i.GetIdx() not in m2matchingatoms:
            mp_unique_atoms.append( i.GetIdx() )
    
    #aromatic atoms will screw this code up, we need to make sure the atoms we're going to delete are due to aromaticitiy         
    aromatic_check = []
    aromatic_idxi = []
    aromatic_idxj = []
    for i in mp_unique_atoms:
        for j in mp_unique_atoms:
            try:
                aromatic_check.append( str( m2.GetBondBetweenAtoms(i, j).GetIsAromatic() ) )
                aromatic_idxi.append( i )
                aromatic_idxj.append( j )
            except:
                pass

    #translate ms_cleaved to our target
    target_specifics = []
    for i in ms_cleaved:
        for j in substrates_lexicon:
            if i == j[1]:
                target_specifics.append( j[0] )
           
    if 'True' not in aromatic_check:   
        AROMATIC_FLAG = None
        temp_lex = create_lexicon( ms, m1 )
        em = Chem.EditableMol(m1)
        for i in range(len(ms_cleaved)):
            temp_lex = create_lexicon( ms, m1 )
            for j in temp_lex:
                if j[0] == ms_cleaved[i]:
                    print j[1]
                    deletion_atom = j[1]
                    em.RemoveAtom(deletion_atom)
                    m1 = em.GetMol()
        m1smiles = Chem.MolToSmiles( m1 )
        m1smiles = clean_smiles( m1smiles )
        try:
            m1 = Chem.MolFromSmiles( m1smiles )
        except:
            pass
    else:
        AROMATIC_FLAG = 'GO'   
        
    
    
    #create lexicon to compare atom indices
    lexicon = create_lexicon(m1, m2)
    ############################################################
    #This will be the double bonds specific to the substrate 
    #I have to do this because double bonds seem to be more specific than single bonds in RDkit, 
    #By knowing the exact position of the double bonds I need to remove and add, I can more accurately transform the molecule
    
    m2bondtypessub = []
    m2bondidxsub = []
    m2bondstartsub = []
    m2bondendsub = []
    
    m1bondtypessub = []
    m1bondidxsub = []
    m1bondstartsub = []
    m1bondendsub = []
      
    for i in lexicon:
        idx1 = i[0]
        m1bondidxsub.append( idx1 )
        m1bondtypessub.append( m1.GetBondWithIdx(idx1).GetBondType() )
    for i in m1bondidxsub:
        m1bondstartsub.append( m1.GetBondWithIdx(int(i)).GetBeginAtomIdx() )
        m1bondendsub.append( m1.GetBondWithIdx(int(i)).GetEndAtomIdx() )
    
    for i in m1bondstartsub:
        for j in lexicon:
            if i == j[0]:
                m2bondstartsub.append( j[1] )
    for i in m1bondendsub:
        for j in lexicon:
            if i == j[0]:
                m2bondendsub.append( j[1] )
    for i in m1bondidxsub:
        for j in lexicon:
            if i == j[0]:
                m2bondidxsub.append( j[1] ) 
    for i in m2bondidxsub:
        m2bondtypessub.append( m2.GetBondWithIdx(int(i)).GetBondType() )
        
    bondindices = ''
    bondindicessub = ''
    
    em = Chem.EditableMol(m1)    
    if len(BondIndex( ms, m2 )) != 0 and AROMATIC_FLAG != 'GO': #compares bonds between both the substrate, product, and target to determine if anything needs to be added
        bondindicessub = BondIndex( m1, m2 )
        bondindicessub = bondindicessub[bondindicessub['M1 Bond Type'] != bondindicessub['M2 Bond Type']]
        #bondindices = bondindices[bondindices['M1 Bond Index'] != bondindices['M1 Bond Start']]
        em = Chem.EditableMol(m1)
        if len(BondIndex(m1, m2)) != 0: ###This logic gate passes the bond formation if the substrate > product chemistry has no bond changes. I may actually destroy this step entirely...
            for i in range(len(bondindicessub)):
                if str(bondindicessub['M1 Bond Type'].irow(i)) == str('SINGLE') and str(bondindicessub['M2 Bond Type'].irow(i)) == str('DOUBLE'):
                    em.RemoveBond(int(bondindicessub['M1 Bond Start'].irow(i)), int(bondindicessub['M1 Bond End'].irow(i)))
                    em.AddBond(int(bondindicessub['M1 Bond Start'].irow(i)), int(bondindicessub['M1 Bond End'].irow(i)), Chem.BondType.DOUBLE)
                else:
                    str(bondindicessub['M1 Bond Type'].irow(i)) == str('DOUBLE')
                    em.RemoveBond(int(bondindicessub['M1 Bond Start'].irow(i)), int(bondindicessub['M1 Bond End'].irow(i)))
                    em.AddBond(int(bondindicessub['M1 Bond Start'].irow(i)), int(bondindicessub['M1 Bond End'].irow(i)), Chem.BondType.SINGLE)
    else:
        pass
        
    m1 = em.GetMol()
        
    ############################################################iterate through bonds in both molecules to see if we need to delete any
    #This will be the double bonds specific to the product
    m2bondtypes = []
    m2bondidx = []
    m2bondstart = []
    m2bondend = []
    
    m1bondtypes = []
    m1bondidx = []
    m1bondstart = []
    m1bondend = []
    
    for i in lexicon:
        idx2 = i[1]
        m2bondidx.append( idx2 )
        m2bondtypes.append( m2.GetBondWithIdx(idx2).GetBondType() )
    for i in m2bondidx:
        m2bondstart.append( m2.GetBondWithIdx(int(i)).GetBeginAtomIdx() )
        m2bondend.append( m2.GetBondWithIdx(int(i)).GetEndAtomIdx() )
        
    m1bondtypes = []
    m1bondidx = []
    m1bondstart = []
    m1bondend = []
    
    for i in m2bondstart:
        for j in lexicon:
            if i == j[1]:
                m1bondstart.append( j[0] )
    for i in m2bondend:
        for j in lexicon:
            if i == j[1]:
                m1bondend.append( j[0] )
    for i in m2bondidx:
        for j in lexicon:
            if i == j[1]:
                m1bondidx.append( j[0] ) 
    for i in m1bondidx:
        m1bondtypes.append( m1.GetBondWithIdx(int(i)).GetBondType() )
                
    #exceptions are atoms that extend from the start bond that actually are irrelevant because where they should go don't exist in the starting molecule
    #this is problematic because they could add bonds to places they shouldn't be
    exceptions = []    
    for i in range(len(m2bondend)):
        if m2bondend[i] not in [int(j[1]) for j in lexicon]:
            exceptions.append( m2bondend[i] )
    
    try: #this try skips the deconstruction phase if it doesn't need to happen
        bondindices = BondIndex( m2, m1 )
        patt1 = Chem.MolFromSmarts(MCS.FindMCS([Chem.MolFromSmiles(steroid1), Chem.MolFromSmiles(steroid2)]).smarts)
        #remove specififed atoms from the substrate
        to_delete_atoms = substrate_modifications['Substrate Unique Atoms'].tolist()
    
        #deletes atoms that are in the substrate but not the product, could be from lyases or whatnot
        for i in range(len(substrate_modifications['Substrate Unique Atoms'])):
            try:
                patt1 = Chem.MolFromSmarts(MCS.FindMCS([m1, Chem.MolFromSmiles(steroid2)]).smarts)
                matching1 = m1.GetSubstructMatch(patt1)
                matching1 = list(matching1)
                #below is indices in m, ordered as patt‘s atoms
                index1 = range(m1.GetNumAtoms())
                #these are the atoms in the substrate that are NOT in the product
                substrate_specific_atoms = list(set(index1) - set(matching1))
                
                #delete just the first substrate_specific_atom because the whole molecule will reindex 
                em = Chem.EditableMol(m1)
                em.RemoveAtom(to_delete_atoms[0])
            
                #need to fix the valences of the atoms we deleted earlier, for whatever reason a radical or a hydrogen is thrown on and the valence is incorrect.
                #Using bondindicessub because it's coordinates are 100% reliable
                for j in range(len(bondindicessub)):
                    if int(bondindicessub['M1 Bond End'].irow(j)) == int(substrate_specific_atoms[0]):
                        deletion_target_num = int(bondindicessub['M1 Bond Start'].irow(j))
                        deletion_target = m1.GetAtomWithIdx( int(bondindicessub['M1 Bond Start'].irow(j)) )
                        deletion_neighbors = [x.GetIdx() for x in deletion_target.GetNeighbors()] 
                        #get the bonds that are connected to indexed atom but not the ones that are in the 'meta-structure'
                        em = Chem.EditableMol(m1)
                        newidx = em.AddAtom(Chem.Atom(6))
                        if len(deletion_neighbors) == 1:
                            em.AddBond(newidx, deletion_neighbors[0], Chem.BondType.SINGLE)
                        elif len(deletion_neighbors) == 2:
                            em.AddBond(newidx, deletion_neighbors[0], Chem.BondType.SINGLE)
                            em.AddBond(newidx, deletion_neighbors[1], Chem.BondType.SINGLE)
                        elif len(deletion_neighbors) == 3:
                            em.AddBond(newidx, deletion_neighbors[0], Chem.BondType.SINGLE)
                            em.AddBond(newidx, deletion_neighbors[1], Chem.BondType.SINGLE)
                            em.AddBond(newidx, deletion_neighbors[2], Chem.BondType.SINGLE)
                        em.RemoveAtom(deletion_target_num)
            except:
                pass
    except:
        bondindices = ''
        bondindicessub = ''
                
    m1 = em.GetMol()
    #Chem.SanitizeMol(m1)
    return m1, bondindices, bondindicessub, AROMATIC_FLAG

# <codecell>

TestRxns = TestUni
comblist = []
i = 56
    
target = TestRxns['Substrates'].irow(i)
native_substrate = TestRxns['Substrates'].irow(i)
native_product = TestRxns['Products'].irow(i)

t = Chem.MolFromSmiles( target )
s = Chem.MolFromSmiles( native_substrate )
p = Chem.MolFromSmiles( native_product )

unique = unique_atoms(t, s)[0]
exceptions = [x for x in unique]
unique2 = unique_atoms(p, t)[0] #relative to the product and substrate, this is so you don't add too much stuff to your target
uniquekeep = unique_atoms(p, s)[0]
exceptions2 = list(set(unique2) - set(uniquekeep))

#deconstruction phase
substrate_modifications = pattern_findersub(target, native_substrate, exceptions)
metastructure, chart, sub, flag = modify_substrate(substrate_modifications, target, native_product, native_substrate)

if flag == 'GO':
    exceptions2 = [] #list(set(unique2) - set(uniquekeep))
    product_modifications = pattern_finderprod(p, target, exceptions2, native_substrate)
    new_product = modify_metastructure(product_modifications, p, target)
else:
    product_modifications = pattern_finderprod(metastructure, native_product, exceptions2, native_substrate)
    new_product = modify_metastructure(product_modifications, metastructure, native_product)
 
metastructure

# <codecell>

create_lexicon(metastructure, p)

# <codecell>

Draw.MolsToGridImage([t, s, p, metastructure, new_product],molsPerRow=8, includeAtomNumbers=True)

# <codecell>

Draw.MolsToGridImage([t, s, p, metastructure, new_product],molsPerRow=8, includeAtomNumbers=False)

# <headingcell level=1>

# Piecing together the algorithm

# <codecell>

from SteroidFunctions import *

test = pcp.get_properties('CanonicalSMILES', 'Androstenedione', 'name')
print str(test[0]['CanonicalSMILES'])
Chem.MolFromSmiles( str(test[0]['CanonicalSMILES']) )

# <codecell>

a = 56

Draw.MolsToGridImage([Chem.MolFromSmiles( TestUni['Substrates'].irow(a)), Chem.MolFromSmiles( TestUni['Products'].irow(a))], includeAtomNumbers=False)

# <codecell>

TestRxns = pd.read_table("Keggsteroids.csv")
TestUni = TestRxns
TestUni = TestUni.drop_duplicates(cols=['Enzymes'])
print len(TestUni)
TestUni.head(4)

# <codecell>

chol = 'CC12CC(=O)C3C(C1CCC2(C(=O)CO)O)CCC4=CC(=O)C=CC34C'

tally = []

for z in range(len(TestUni)):
    try:
        t = Chem.MolFromSmiles( chol ) 
        s = Chem.MolFromSmiles( TestUni['Substrates'].irow(z) )
        p = Chem.MolFromSmiles( TestUni['Products'].irow(z) )
        tally.append( scan_atoms(t, s, p) )
        if z % 10 == 0:
                print z  
    except:
        tally.append( 'do not proceed' )

# <codecell>

#Draw.MolsToGridImage([t, Chem.MolFromSmiles( TestUni['Substrates'].irow(0)), Chem.MolFromSmiles( TestUni['Products'].irow(0))], includeAtomNumbers=True)

# <codecell>

TestUni['Anabolic Compatible'] = tally 
print len(TestUni)
TestUni2 = TestUni[TestUni['Anabolic Compatible'] == 'proceed']
print len(TestUni2)

# <codecell>

TestUni2

# <codecell>

start = timeit.default_timer()

anabolic = chol
substrates = []
products = []
targets = []
results = []
enzymes = []
for x in range(len(TestUni2)):
    try:
        results.append( explore_substrate(anabolic, x, TestUni2) )
        substrates.append( TestUni2['Substrates'].irow(x) )
        products.append( TestUni2['Products'].irow(x) )
        targets.append( anabolic )
        enzymes.append( TestUni2['Enzymes'].irow(x) )
    except:
        pass
    
stop = timeit.default_timer()
#unique_prods, unique_substrates, unique_products, unique_enzymes = filter_results(results, substrates, products, TestUni, enzymes)
unique_prods, unique_substrates, unique_products, unique_enzymes = results, substrates, products, enzymes

print "Finished reactions... searching hits on PubChem... "

#search pubchem for unique hits
novel_compounds = []
for i in unique_prods:
    searches = pcp.get_compounds('CanonicalSMILES', str(Chem.MolToSmiles(i)), 'smiles')
    if str(searches) == '[Compound()]':
        novel_compounds.append( i )
        
#put things in a df
unique_prods_smiles = []
novel_smiles = []
for i in unique_prods:
    unique_prods_smiles.append( Chem.MolToSmiles(i) )
for i in novel_compounds:
    novel_smiles.append( Chem.MolToSmiles(i) )

novels =[ ]
for i in unique_prods_smiles:
    if i in novel_smiles:
        novels.append( 'Novel' )
    else:
        novels.append( 'Found in Pubchem' )
Results = pd.DataFrame({'Native Substrate': unique_substrates, 'Native Product': unique_products, 'Products': unique_prods_smiles, 'Novel Compound?': novels, 'Enzymes': unique_enzymes})
Results = Results.drop_duplicates(cols=['Products'])
        
print 'Novel Compounds found...'   
print len(novel_compounds)
print "Runtime..."
print stop - start    
Draw.MolsToGridImage(unique_prods,molsPerRow=8, includeAtomNumbers=False)

# <codecell>

Chem.MolFromSmiles( Results['Native Substrate'][0] )

# <codecell>

Chem.MolFromSmiles( 

# <codecell>

Results2, product_pictures2 = prediction_algorithm(str(Results['Products'][4]), TestUni)

# <codecell>

#Results.to_csv("Lathosterol_hits.csv")
product_pictures2

# <codecell>

Results2

# <codecell>

a = 56

asd = Draw.MolsToGridImage([Chem.MolFromSmiles( TestUni['Substrates'].irow(a)), Chem.MolFromSmiles( TestUni['Products'].irow(a))], includeAtomNumbers=False)
print TestUni['Enzymes'].irow(56)
asd

# <markdowncell>

# Now that we have a bunch of TestRxns, we need to run the target as the substrate and see if the product is made. 

# <codecell>

#matches = []
#new_products = []

#for i in range(len(TestRxns)):
#    try:
#        product = linear_exploration(TestRxns['Substrates'].irow(i), i, TestRxns)
#        product = Chem.MolToSmiles( product, canonical=True )
#        new_products.append( product )
#        
#        p1 = Chem.MolFromSmiles( TestRxns['Products'].irow(i) ) 
#        p = Chem.MolToSmiles( p1, canonical=True  )
#        
#        if str(product) == str(p):
#            matches.append( 'matches' )
#        else:
#            matches.append( 'no match' )
#    except:
#        matches.append( 'glitch' )
#        new_products.append( 'glitch' )

# <codecell>

TestRxns['matches'] = matches
TestRxns['New Products'] = new_products
count = 0
for i in matches:
    if i == 'matches':
        count = count + 1
print count    

# <codecell>

TestRxns = TestRxns[TestRxns['matches'] != 'matches']
print len(TestRxns)

# <codecell>

def prediction_algorithm(target_smiles, df):

    start = timeit.default_timer()
    TestUni = df
    chol = str(target_smiles)
    tally = []
    for z in range(len(TestUni)):
        try:
            t = Chem.MolFromSmiles( chol ) 
            s = Chem.MolFromSmiles( TestUni['Substrates'].irow(z) )
            p = Chem.MolFromSmiles( TestUni['Products'].irow(z) )
            tally.append( scan_atoms(t, s, p) )
            if z % 10 == 0:
                    print z  
        except:
            tally.append( 'do not proceed' )
            
    TestUni['Anabolic Compatible'] = tally 
    TestUni2 = TestUni[TestUni['Anabolic Compatible'] == 'proceed']
    print "Potential hits..."
    print len(TestUni2)

    anabolic = chol
    substrates = []
    products = []
    targets = []
    results = []
    enzymes = []
    for x in range(len(TestUni2)):
        try:
            results.append( explore_substrate(anabolic, x, TestUni2) )
            substrates.append( TestUni2['Substrates'].irow(x) )
            products.append( TestUni2['Products'].irow(x) )
            targets.append( anabolic )
            enzymes.append( TestUni2['Enzymes'].irow(x) )
        except:
            pass
        
    stop = timeit.default_timer()
    #unique_prods, unique_substrates, unique_products, unique_enzymes = filter_results(results, substrates, products, TestUni, enzymes)
    unique_prods, unique_substrates, unique_products, unique_enzymes = results, substrates, products, enzymes
    
    print "Finished reactions... searching hits on PubChem... "
    
    #search pubchem for unique hits
    novel_compounds = []
    for i in unique_prods:
        searches = pcp.get_compounds('CanonicalSMILES', str(Chem.MolToSmiles(i)), 'smiles')
        if str(searches) == '[Compound()]':
            novel_compounds.append( i )
            
    #put things in a df
    unique_prods_smiles = []
    novel_smiles = []
    for i in unique_prods:
        unique_prods_smiles.append( Chem.MolToSmiles(i) )
    for i in novel_compounds:
        novel_smiles.append( Chem.MolToSmiles(i) )
    
    novels =[ ]
    for i in unique_prods_smiles:
        if i in novel_smiles:
            novels.append( 'Novel' )
        else:
            novels.append( 'Found in Pubchem' )
    Results = pd.DataFrame({'Native Substrate': unique_substrates, 'Native Product': unique_products, 'Products': unique_prods_smiles, 'Novel Compound?': novels, 'Enzymes': unique_enzymes})
    Results = Results.drop_duplicates(cols=['Products'])
            
    print 'Novel Compounds found...'   
    print len(novel_compounds)
    print "Runtime..."
    print stop - start    
    product_pictures = Draw.MolsToGridImage(unique_prods,molsPerRow=8, includeAtomNumbers=False)
    
    return Results, product_pictures

# <codecell>


