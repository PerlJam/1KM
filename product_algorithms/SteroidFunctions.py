# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import FunctionalGroups
import os
import gzip,random
from rdkit.Chem import MCS
import pandas as pd
import timeit
import pubchempy as pcp
from bioservices import Kegg
from bioservices import *
from compiler.ast import flatten
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import Chem
import numpy as np
from rdkit.Chem import AllChem

# <codecell>

def create_lexicon(molecule1, molecule2):
    #Chem.Kekulize(molecule1)
    #Chem.Kekulize(molecule2)
    patt1 = Chem.MolFromSmarts(MCS.FindMCS([molecule2, molecule1], matchValences=True).smarts)
    matching1 = molecule2.GetSubstructMatch(patt1)
    matching2 = molecule1.GetSubstructMatch(patt1)
    
    #below is indices in m, ordered as patt‘s atoms
    index1 = range(molecule2.GetNumAtoms())
    #these are the atoms in the product that are NOT in the metastructure
    product_specific_atoms = list(set(index1) - set(matching1))

    matching1 = zip(matching1, range(molecule2.GetNumAtoms()) )
    matching2 = zip(matching2, range(molecule1.GetNumAtoms()) )
    
    #lexicon for what values equal what. this is a bit confusing but it's the product's substructure that's similar with the meta-metastructure's
    #then the corresponding atom on the meta-metastructure to the metastructure
    
    lexicon = sorted(zip([int(i[0]) for i in matching2], [int(i[0]) for i in matching1]))
    return lexicon

def unique_atoms(s1, s2):
    lexicon = create_lexicon(s1, s2)
    
    s1Idx = []
    s2Idx = []
    for atom in s1.GetAtoms():
        s1Idx.append( atom.GetIdx() )
    for atom in s2.GetAtoms():
        s2Idx.append( atom.GetIdx() )
        
    substrate_unique_atoms = list(set(s1Idx) - set([x[0] for x in lexicon]))
    product_unique_atoms = list(set(s2Idx) - set([x[1] for x in lexicon]))
    
    return substrate_unique_atoms, product_unique_atoms

def protect(s1):
    '''returns the atoms that need to be protected, this is done by subtracting all atoms from the unique atoms'''
    s1list = []
    for atom in s1.GetAtoms():
        s1list.append( atom.GetIdx() )
    
        
    substrate_unique_atoms, product_unique_atoms = unique_atoms(s1, s2)
    print substrate_unique_atoms
    for i in substrate_unique_atoms:
        atom = s1.GetAtomWithIdx(i)
        neighbors = [x.GetIdx() for x in atom.GetNeighbors()] 
    for i in neighbors:
        substrate_unique_atoms.append( i )
        atom = s1.GetAtomWithIdx(i)
        alphaneighbors = [x.GetIdx() for x in atom.GetNeighbors()] 
    for i in alphaneighbors:
        substrate_unique_atoms.append( i )
    
    protect = list(set(s1list) - set(substrate_unique_atoms))
    print substrate_unique_atoms
    for atom in s1.GetAtoms():
        for i in protect:
            if atom == i:
                atom.SetProp('_protected', '1')
    return s1

# <codecell>

def pattern_findersub(steroid1, steroid2, exceptions):
    m1 = Chem.MolFromSmiles(steroid1)
    m2 = Chem.MolFromSmiles(steroid2)
    patt1 = Chem.MolFromSmarts(MCS.FindMCS([Chem.MolFromSmiles(steroid1), Chem.MolFromSmiles(steroid2)]).smarts)
    
    matching1 = m1.GetSubstructMatch(patt1)
    matching1 = list(matching1)
    
    ####################important exception line
    for i in exceptions:
        matching1.append( i )
    
    #below is indices in m, ordered as patt‘s atoms
    index1 = range(Chem.MolFromSmiles(steroid1).GetNumAtoms())
    #these are the atoms in the substrate that are NOT in the product
    substrate_specific_atoms = list(set(index1) - set(matching1))
    
    del_bonds = []
    add_connections = []
    add_bonds = []
    del_connections = []
    
    for i in substrate_specific_atoms:
        atom = m1.GetAtomWithIdx(i)
        
        #get the bonds that are connected to indexed atom but not the ones that are in the 'meta-structure'
        neighbors = [x.GetIdx() for x in atom.GetNeighbors()] 
        extra_bonds = list(set(neighbors) & set(substrate_specific_atoms))
        
        #get the bonds of these atoms that need to be deleted 
        bondtype = []
        for bond in extra_bonds:
            bond = str(m1.GetBondBetweenAtoms(i, bond).GetBondType())
            bond = bond.replace('rdkit.Chem.rdchem.BondType.', '')
            bondtype.append( bond )

        del_connections.append( extra_bonds )
        del_bonds.append( bondtype )
   
    substrate_modifications = pd.DataFrame({'Substrate Unique Atoms': substrate_specific_atoms, 'Connections to be deleted': del_connections, 'Bonds to be deleted': del_bonds})
    return substrate_modifications

# <codecell>

def pattern_finderprod(metastructure, steroid2, exceptions, steroid1):
    ms = Chem.MolFromSmiles(steroid1)
    mp = Chem.MolFromSmiles(steroid2)
    m2 = Chem.MolFromSmiles(steroid2)
    
    native_lexicon = create_lexicon(ms, mp)
    
    mp_atoms = []
    ms_mp_atoms = []
    for atom in mp.GetAtoms():
        mp_atoms.append( atom.GetIdx() )
    for i in native_lexicon:
        ms_mp_atoms.append( i[1] )
    
    product_specific_atoms = []
    for i in mp_atoms:
        if i not in ms_mp_atoms:
            product_specific_atoms.append( i )
    
    neighbors = []
    del_connections = []
    del_bonds = []
    add_connections = []
    add_bonds = []
    atom_type = []
    neighboring_atoms = []
    
    for i in product_specific_atoms:
        atom = m2.GetAtomWithIdx(i)
        neighbors = [x.GetIdx() for x in atom.GetNeighbors()] 
        neighboring_atoms.append ( neighbors )
        extra_bonds = list(set(neighbors) & set(product_specific_atoms))
        #atoms to be added
        adjacent_atoms = []
        for adjacent in neighbors:
            adjacent_atoms.append( m2.GetBondBetweenAtoms(i, adjacent).GetBondType() )
        
        #get the bonds of these atoms that need to be deleted 
        bondtype = []
        for bond in extra_bonds:
            bond = str(m2.GetBondBetweenAtoms(i, bond).GetBondType())
            bond = bond.replace('rdkit.Chem.rdchem.BondType.', '')
            bondtype.append( bond )
            
        #get the type of atom to be added
        atom_type.append( atom.GetAtomicNum() )
            
        del_connections.append( extra_bonds )
        del_bonds.append( bondtype )
        add_connections.append( neighbors )
        add_bonds.append( adjacent_atoms )
        
        
        
    product_modifications = pd.DataFrame({'Product Unique Atoms': product_specific_atoms, 'Connections to be deleted': del_connections, 'Bonds to be deleted': del_bonds, 
'Connections to be added': add_connections, 'Bonds to be added': add_bonds, 'Atomic Number': atom_type, 'Neighbors': neighboring_atoms})
    
    return product_modifications

# <codecell>

def BondIndex(mol1, mol2):
    
    lexicon = create_lexicon(mol1, mol2)
    
    Chem.Kekulize(mol1)
    Chem.Kekulize(mol2)
    
    bonds = []
    bondstart = []
    bondend = []
    bonds2 = []
    bondstart2 = []
    bondend2 = []
    m1lex = []
    m2lex = []
    
    for i in lexicon:
        m1lex.append( i[0] )
        
    for i in m1lex:
        for j in lexicon:
            if i == j[0]:
                m2lex.append( j[1] ) 
                
    for i in m1lex:
        for j in m1lex:
            if j < i:
                try:
                    bonds.append( str(mol1.GetBondBetweenAtoms(i,j).GetBondType()) )
                except:
                    bonds.append( 'none' )
                bondstart.append( i )
                bondend.append( j )

    bondsdf = pd.DataFrame({'M1 Bond Type': bonds, 'M1 Bond Start': bondstart, 'M1 Bond End': bondend})#, 'Bond Type2': bonds2, 'Bond Start2': bondstart2, 'Bond End2': bondend2})
    #bondsdf = bondsdf[bondsdf['Bond Type'] !=  bondsdf['Bond Type2']]
    bondsdf = bondsdf[bondsdf['M1 Bond Type'] != 'none']
    bondtrans2b = bondsdf['M1 Bond Start'].tolist()
    bondtrans2e = bondsdf['M1 Bond End'].tolist()
    
    ###translate bond start and endings to molecule two, find equivalent atoms
    bond2start = []
    bond2end = []
    for i in bondtrans2b:
        for j in lexicon:
            if i == j[0]:
                bond2start.append( j[1] )
    for i in bondtrans2e:
        for j in lexicon:
            if i == j[0]:
                bond2end.append( j[1] )
                 
                
    divis = len(bondtrans2b)
    test = []
    for i in bondtrans2b:
        for j in lexicon:
            if i == j[0]:
                test.append( j[1] )
            else:
                test.append( 'failure' )
    test2 = chunkIt(test, divis)
    for i in test2:
        i = list(set(i))
        for j in i:
            if len(i) == 2:
                if j != 'failure':
                    bondstart2.append( j )
            elif len(i) ==1:
                bondstart2.append( j )
                
    divis = len(bondtrans2b)
    test = []
    for i in bondtrans2e:
        for j in lexicon:
            if i == j[0]:
                test.append( j[1] )
            else:
                test.append( 'failure' )
    test2 = chunkIt(test, divis)
    for i in test2:
        i = list(set(i))
        for j in i:
            if len(i) == 2:
                if j != 'failure':
                    bondend2.append( j )
            elif len(i) ==1:
                bondend2.append( j )
                
    for i in range(len(bondend2)):
        try:
            bonds2.append( str(mol2.GetBondBetweenAtoms(bond2start[i], bond2end[i]).GetBondType()) )
        except:
            bonds2.append('failure')
                
    bondsdf['M2 Bond Type'] = bonds2
    bondsdf['M2 Bond Start'] = bond2start
    bondsdf['M2 Bond End'] = bond2end
    
    #bondsdf = pd.DataFrame({'Bond Type': bonds, 'Bond Start': bondstart, 'Bond End': bondend, 'Bond Type2': bonds2, 'Bond Start2': bondstart2, 'Bond End2': bondend2})
    bondsdf = bondsdf[bondsdf['M1 Bond Type'] !=  bondsdf['M2 Bond Type']]
    bondsdf = bondsdf[bondsdf['M1 Bond Type'] != 'failure']
    bondsdf = bondsdf[bondsdf['M2 Bond Type'] != 'failure']
    
    return bondsdf

def chunkIt(seq, num):
  avg = len(seq) / float(num)
  out = []
  last = 0.0

  while last < len(seq):
    out.append(seq[int(last):int(last + avg)])
    last += avg

  return out

# <codecell>

def AddBonds(new_product, product):
    '''Adds bonds if there is an analogous carbon in both the m1, and m2 molecules, does not yet have any exception filtering.'''
    bondindices = BondIndex2(product, new_product)
    bondindices = bondindices[bondindices['M1 Bond Type'] != 'SINGLE']
    lexicon = create_lexicon(new_product, product)
    
    for i in lexicon:
        if i[1] == bondindices['M1 Bond End'].irow(0):
            bondend = i[0]
        if i[1] == bondindices['M1 Bond Start'].irow(0):
            bondstart = i[0]
    
    ########################################################
    #edits bonds to look like product
    
    em = Chem.EditableMol(new_product)
    for i in range(len(bondindices)):
        if str(bondindices['M1 Bond Type'].irow(i)) == str('DOUBLE'):
            em.RemoveBond(int(bondend), int(bondstart))
            em.AddBond(int(bondend), int(bondstart), Chem.BondType.DOUBLE)
        else:
            str(bondindices['M1 Bond Type'].irow(i)) == str('SINGLE')
            em.RemoveBond(int(bondend), int(bondstart))
            em.AddBond(int(bondend), int(bondstart), Chem.BondType.SINGLE)
    m1 = em.GetMol()        
    #########################################################
    return m1

# <codecell>

def BondIndex2(mol1, mol2):
    lexicon = create_lexicon(s, p)
    bonds = []
    bondstart = []
    bondend = []
    bonds2 = []
    bondstart2 = []
    bondend2 = []
    m1lex = []
    m2lex = []
    
    for i in lexicon:
        m1lex.append( i[0] )
        
    for i in m1lex:
        for j in lexicon:
            if i == j[0]:
                m2lex.append( j[1] ) 
                
    for i in m1lex:
        for j in m1lex:
            if j < i:
                try:
                    bonds.append( str(mol1.GetBondBetweenAtoms(i,j).GetBondType()) )
                except:
                    bonds.append( 'none' )
                bondstart.append( i )
                bondend.append( j )

    bondsdf = pd.DataFrame({'M1 Bond Type': bonds, 'M1 Bond Start': bondstart, 'M1 Bond End': bondend})#, 'Bond Type2': bonds2, 'Bond Start2': bondstart2, 'Bond End2': bondend2})
    #bondsdf = bondsdf[bondsdf['Bond Type'] !=  bondsdf['Bond Type2']]
    bondsdf = bondsdf[bondsdf['M1 Bond Type'] != 'none']
    return bondsdf

# <codecell>

def combinatorial_rxn(target, native_substrate, native_product):
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
    metastructure, chart, sub = modify_substrate(substrate_modifications, target, native_product, native_substrate)
    product_modifications = pattern_finderprod(metastructure, native_product, exceptions2, native_substrate)
    new_product = modify_metastructure(product_modifications, metastructure, native_product)
    #new_product = cleanup(new_product, native_substrate, native_product)
    #if len(BondIndex(s, p)) != 0:
    #    try:
    #        new_product = AddBonds(new_product, p)
    #    except:
    #        new_product = new_product

    return new_product

# <codecell>

def filter_results(results, sucsubstrates, sucproducts, df, enzymes):
    rsmiles = []
    smiles = []
    for i in results:
        rsmiles.append( Chem.MolToSmiles( i ) )
        
    rsmiles2 = []
    for i in rsmiles:
        if str(i[:8]) == 'C.C.O.O.':
            rsmiles2.append( i.replace('C.C.O.O.', ''))
        else:
            rsmiles2.append( i )
    rsmiles3 = []
    for i in rsmiles2:
        if str(i[:2]) == 'O.':
            rsmiles3.append( i.replace('O.', ''))
        else:
            rsmiles3.append( i )
    
    products = []
    filtered_products = []
    filtered_substrates = []
    for i in rsmiles3:
        products.append( Chem.MolFromSmiles(i) )
    products2 = []
    for i in range(len(products)):
        products2.append( products[i] )
        filtered_products.append( sucproducts[i] )
        filtered_substrates.append( sucsubstrates[i] )
    #for i in products2:
    #    for atom in i.GetAtoms():
    #        atom.GetNumRadicalElectrons(0)
    #    i = Chem.AddHs(i)
    #    Chem.SanitizeMol(i)
    
    rsmilesfinal = []
    for i in products2:
        try:
            rsmilesfinal.append( Chem.MolToSmiles(i) )
        except:
            pass
        
        
    background_pool = [] 
    for i in range(len(df)):
        tp = Chem.MolFromSmiles( df['Products'].irow(i) )
        ts = Chem.MolFromSmiles( df['Substrates'].irow(i) )
        background_pool.append( Chem.MolToSmiles( tp ) )
        background_pool.append( Chem.MolToSmiles( ts ) )

        
    rsmiles4 = []
    filtered_products2 = []
    filtered_substrates2 = []
    filtered_enzymes = []
    for i in range(len(rsmilesfinal)):
        if rsmilesfinal[i] not in background_pool:
            rsmiles4.append( Chem.MolFromSmiles( rsmilesfinal[i] ) )
            filtered_products2.append( filtered_products[i] )
            filtered_substrates2.append( filtered_substrates[i] )
            filtered_enzymes.append( enzymes[i] )
        
    return rsmiles4, filtered_substrates2, filtered_products2, filtered_enzymes

def filter_results2(results):
    rsmiles = results
    rsmiles2 = []
    rsmiles3 = []
    for i in rsmiles:
        if str(i[:8]) == 'C.C.O.O.':
            rsmiles2.append( i.replace('C.C.O.O', ''))
        else:
            rsmiles2.append( i )
    for i in rsmiles2:
        if str(i[:2]) == 'O.':
            rsmiles3.append( i.replace('O.', ''))
        else:
            rsmiles3.append( i )
    return rsmiles3

# <codecell>

def linear_exploration(inputsmiles, ith, TestRxns):
    comblist = []
    i = ith
    
    try: #try starting from substrate
        target = TestRxns['Substrates'][i]
        native_substrate = TestRxns['Substrates'][i]
        native_product = TestRxns['Products'][i]
        
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
            #unique2 = [] #unique_atoms(p, p)[0] #relative to the product and substrate, this is so you don't add too much stuff to your target
            #uniquekeep = unique_atoms(p, p)[0]
            exceptions2 = [] #list(set(unique2) - set(uniquekeep))
            product_modifications = pattern_finderprod(p, target, exceptions2, native_substrate)
            new_product = modify_metastructure(product_modifications, p, target)
        else:
            product_modifications = pattern_finderprod(metastructure, native_product, exceptions2, native_substrate)
            new_product = modify_metastructure(product_modifications, metastructure, native_product)
    
    except: #try starting from product
        target = TestRxns['Substrates'][i]
        native_substrate = TestRxns['Substrates'][i]
        native_product = TestRxns['Products'][i]
        
        t = Chem.MolFromSmiles( target )
        s = Chem.MolFromSmiles( native_substrate )
        p = Chem.MolFromSmiles( native_product )
        
        Chem.Kekulize(p, clearAromaticFlags=True)
        
        unique = unique_atoms(p, p)[0]
        exceptions = [x for x in unique]
        unique2 = unique_atoms(p, p)[0] #relative to the product and substrate, this is so you don't add too much stuff to your target
        uniquekeep = unique_atoms(p, p)[0]
        exceptions2 = list(set(unique2) - set(uniquekeep))
        
        #deconstruction phase
        #substrate_modifications = pattern_findersub(native_product, native_product, exceptions)
        #metastructure, chart, sub, flag = modify_substrate(substrate_modifications, native_product, native_product, native_substrate)
        
        product_modifications = pattern_finderprod(p, target, exceptions2, native_substrate)
        new_product = modify_metastructure(product_modifications, p, target)

    print Chem.MolToSmiles(new_product)
    return new_product

# <codecell>

def explore_substrate(inputsmiles, ith, TestRxns):
    comblist = []
    i = ith
    
    #try: #try starting from substrate
    target = inputsmiles
    native_substrate = TestRxns['Substrates'].irow(i)
    native_product = TestRxns['Products'].irow(i)
    
    t = Chem.MolFromSmiles( inputsmiles )
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
        #unique2 = [] #unique_atoms(p, p)[0] #relative to the product and substrate, this is so you don't add too much stuff to your target
        #uniquekeep = unique_atoms(p, p)[0]
        exceptions2 = [] #list(set(unique2) - set(uniquekeep))
        product_modifications = pattern_finderprod(p, target, exceptions2, native_substrate)
        new_product = modify_metastructure(product_modifications, p, target)
    else:
        product_modifications = pattern_finderprod(metastructure, native_product, exceptions2, native_substrate)
        new_product = modify_metastructure(product_modifications, metastructure, native_product)

    '''except: #try starting from product
        target = TestRxns['Substrates'][i]
        native_substrate = TestRxns['Substrates'][i]
        native_product = TestRxns['Products'][i]
        
        t = Chem.MolFromSmiles( target )
        s = Chem.MolFromSmiles( native_substrate )
        p = Chem.MolFromSmiles( native_product )
        
        Chem.Kekulize(p, clearAromaticFlags=True)
        
        unique = unique_atoms(p, p)[0]
        exceptions = [x for x in unique]
        unique2 = unique_atoms(p, p)[0] #relative to the product and substrate, this is so you don't add too much stuff to your target
        uniquekeep = unique_atoms(p, p)[0]
        exceptions2 = list(set(unique2) - set(uniquekeep))
        
        #deconstruction phase
        #substrate_modifications = pattern_findersub(native_product, native_product, exceptions)
        #metastructure, chart, sub, flag = modify_substrate(substrate_modifications, native_product, native_product, native_substrate)
        
        product_modifications = pattern_finderprod(p, target, exceptions2, native_substrate)
        new_product = modify_metastructure(product_modifications, p, target)'''

    return new_product

# <codecell>

def get_nearby_atoms(mol, targetidx, radius):
    idxs = []
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    for i in mol.GetAtoms():
        at1Coords = np.array(conf.GetAtomPosition(i.GetIdx()))
        at2Coords = np.array(conf.GetAtomPosition(targetidx))
        if np.linalg.norm(at1Coords-at2Coords) < float(radius):
            idxs.append( i.GetIdx() )
    return idxs
            
def clean_smiles(smiles):
    smiles = str(smiles) 
    for i in range(len(smiles)):
        if smiles[-i] == '.':
            i = i - 1
            return smiles[-i:]

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
 
    #create lexicon to compare atom indices'''
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

def modify_metastructure(product_modifications, metastructure, steroid2):
    m2 = Chem.MolFromSmiles(steroid2)
    #msubstrate = Chem.MolFromSmiles(steroid1)
    patt1 = Chem.MolFromSmarts(MCS.FindMCS([metastructure, m2], matchValences=True).smarts)
    
    #convert number of product specifc atom to our metastructure
    anchors = [] #anchors are in the MCS, they will ultimately be deleted but are important for figuring out where to add bonds
    anchortype = []
    lexicon = create_lexicon(metastructure, m2)
    for i in product_modifications['Connections to be added'].tolist():
        for k in i: #connections in list if there are multiple
            for j in lexicon:
                if k == j[1]:
                    anchors.append( j[0] ) 
                    atom = metastructure.GetAtomWithIdx(int(j[0]))
                    anchortype.append( atom.GetAtomicNum()) 

    neighbors = [] #for every index gives the neighbors of the same index in the pandas DF earlier             
    for i in anchors:
        adjacent_atoms = []
        atom = metastructure.GetAtomWithIdx(i)
        adjacent_atoms = [x.GetIdx() for x in atom.GetNeighbors()] 
        neighbors.append( adjacent_atoms )
    
    #add the product-specific atoms
    em = Chem.EditableMol(metastructure)
    newindexes = []
    newanchors = []
    
    ###KEEP THESE DELETED####
    '''for i in product_modifications['Atomic Number'].tolist():
        newidx = em.AddAtom(Chem.Atom( int(i) ))
        newindexes.append( newidx )
    for i in range(len(anchors)):
        newanchor = em.AddAtom(Chem.Atom( anchortype[i] ))
        newanchors.append( newanchor )'''
        
    
    #####logic gate for if a carboxyl like addition is going on
    similar_indices = []
    for i in range(len(product_modifications)):
        for j in range(len(product_modifications)):
            if product_modifications['Connections to be added'].irow(i) == product_modifications['Connections to be added'].irow(j):
                similar_indices.append( i ) 
            
    if len(similar_indices) > 2:
        #translate neighbor number to what it corresponds to in m1
        m1 = em.GetMol()
        for i in lexicon:
            if product_modifications['Neighbors'].irow( similar_indices[0] )[0] == i[1]:
                neighbor = i[0]
                atom = m1.GetAtomWithIdx(neighbor)
                neighbortype = atom.GetAtomicNum()
                neighbor_of_neighbor = [x.GetIdx() for x in atom.GetNeighbors()]
                neighboranchor = em.AddAtom(Chem.Atom( int(neighbortype) ))
        new_atoms = []
        for i in range(len(product_modifications)):
            new_atom = em.AddAtom( Chem.Atom( int(product_modifications['Atomic Number'].irow(i)) ))
            new_atoms.append( new_atom )
            
        for i in range(len(product_modifications)):
                if str(product_modifications['Bonds to be added'].irow(i)[0]) == 'DOUBLE':
                    em.AddBond( int(neighboranchor),int(new_atoms[i]), Chem.BondType.DOUBLE)
                elif str(product_modifications['Bonds to be added'].irow(i)[0]) == 'SINGLE':
                    em.AddBond(int(neighboranchor),int(new_atoms[i]), Chem.BondType.SINGLE)
        for i in neighbor_of_neighbor:
            em.AddBond(int(i), int(neighboranchor), Chem.BondType.SINGLE)
            
        #get rid of old anchor
        for i in list(set(anchors)):
            em.RemoveAtom(i)
      
    else:
        #add the product-specific atoms
        em = Chem.EditableMol(metastructure)
        newindexes = []
        newanchors = []
        
        for i in product_modifications['Atomic Number'].tolist():
            newidx = em.AddAtom(Chem.Atom( int(i) ))
            newindexes.append( newidx )
        for i in range(len(anchors)):
            newanchor = em.AddAtom(Chem.Atom( anchortype[i] ))
            newanchors.append( newanchor )
            
        mref = em.GetMol()
        
        #combine the new atom with it's new anchor
        for i in range(len(newindexes)):
            if str(product_modifications['Bonds to be added'][i][0]) == 'DOUBLE':
                em.AddBond(int(newindexes[i]), int(newanchors[i]), Chem.BondType.DOUBLE)
            elif str(product_modifications['Bonds to be added'][i][0]) == 'SINGLE':
                em.AddBond(int(newindexes[i]), int(newanchors[i]), Chem.BondType.SINGLE)
                
        #combine new structure (newanchor + new atom) to the neighbors of the old anchor
        for i in range(len(anchors)):
            for j in range(len(neighbors[i])):
                atom = mref.GetAtomWithIdx( int(neighbors[i][j]) )
                em.AddBond(int(newanchors[i]), int(neighbors[i][j]), Chem.BondType.SINGLE)
                    
        #get rid of old anchor
        for i in anchors:
            em.RemoveAtom(i)
            
    m1 = em.GetMol()
    for atom in m1.GetAtoms():
        atom.SetNumRadicalElectrons(0)
    Chem.SanitizeMol(m1)
    return m1

# <codecell>

def scan_atoms(t, s, p):
    #scheme: take things that differ from the product and substrate and the target and substrate and compare those differences
    #if they're close in space, don't proceed with reaction 
    ##############
    '''native reaction'''
    native_rxn = create_lexicon(s,p)
    p_s_atoms = []
    patoms = []
    for atom in p.GetAtoms():
        patoms.append( atom.GetIdx() )
    for i in native_rxn:
        p_s_atoms.append( i[1] )
    product_specific_atoms = []
    for i in patoms:
        if i not in p_s_atoms:
            product_specific_atoms.append( i )
            
    '''native reaction bonds'''
    native_bonds = BondIndex(p, s)
    if len(native_bonds) > 0:
        for i in range(len(native_bonds)):
            product_specific_atoms.append( int(native_bonds['M1 Bond Start'].irow(i)) )
            product_specific_atoms.append( int(native_bonds['M1 Bond End'].irow(i)) )
    
    ##############
    '''synthetic reaction'''
    target_specific_atoms = []
    synthetic_rxn = create_lexicon(s,t)
    s_t_atoms = []
    t_atoms = []
    for atom in t.GetAtoms():
        t_atoms.append( atom.GetIdx() )
    for i in synthetic_rxn:
        s_t_atoms.append( i[1] )
    target_specific_atoms = []
    for i in t_atoms:
        if i not in s_t_atoms:
            target_specific_atoms.append( i )
            
    '''synthetic reaction bonds'''
    synthetic_bonds = BondIndex(t, s)
    if len(synthetic_bonds) > 0:
        for i in range(len(synthetic_bonds)):
            target_specific_atoms.append( int(synthetic_bonds['M1 Bond Start'].irow(i)) )
            target_specific_atoms.append( int(synthetic_bonds['M1 Bond End'].irow(i)) )
     
    if len(product_specific_atoms) > 0:
        for j in product_specific_atoms:
            tight_enzyme_specificity = get_nearby_atoms(p, j, 2.5) 
    else:
        tight_enzyme_specificity = []
    substrate_to_product = create_lexicon( s, p )
    equivalent_sub_atoms = []
    for i in tight_enzyme_specificity: #transferring the like atoms from the product to substrate...
        for j in substrate_to_product:
            if i == j[1]:
                equivalent_sub_atoms.append( j[0] )

    target_to_substrate = create_lexicon( t, s )
    atoms_encoding_specificity = []
    for j in equivalent_sub_atoms: #get the atoms in the substrate that should encode specificity (i.e. those around the ones in the product)
        atoms_encoding_specificity = get_nearby_atoms(s, j, 1.5)
    #atoms_encoding_specificity = [item for sublist in atoms_encoding_specificity for item in sublist] #flatten
        
    equivalent_target_atoms = []
    for i in atoms_encoding_specificity:
        for j in target_to_substrate:
            if i == j[1]:
                equivalent_target_atoms.append( j[0] )
    
    if len(equivalent_target_atoms) == 0:
        return "do not proceed"
    #if len(equivalent_sub_atoms) != len(equivalent_target_atoms):
    #    return 'do not proceed'
    for i in equivalent_target_atoms:
        if i in target_specific_atoms:
            return 'do not proceed'
        else:
            return 'proceed'

def naive_find_alternative_substrates(input_smiles_list, df, df_row, threshold):
    '''uses the naive approach to finding an alternative substrate by comparing the area around where an enzyme performs transformation and areas that are unique to your target'''
    input_smiles_list = tanimoto_candidates(df['Substrates'].irow(df_row), list(set(input_smiles_list)), threshold) ###this function quickly narrows down potential candidates for doing the naive alternative substrate search
    tally = []
    print "Progress..."
    for z in range(len(input_smiles_list)):
        try:
            t = input_smiles_list[z][1]
            s = Chem.MolFromSmiles( df['Substrates'].irow(df_row) )
            p = Chem.MolFromSmiles( df['Products'].irow(df_row) )
            tally.append( scan_atoms(t, s, p) )
            if z % 10 == 0 and z != 0:
                    print z  
        except:
            tally.append( 'do not proceed' )
    #mol_smiles = [Chem.MolToSmiles(x, canonical=True) for x in input_smiles_list[1]]
    mol_smiles = []
    for i in range(len(input_smiles_list)):
        mol_smiles.append( Chem.MolToSmiles(input_smiles_list[i][1]) )
    mol_scores = ['%.4f'%x[0] for x in input_smiles_list]
    compatible_mols = pd.DataFrame({'Alternative Substrate SMILES': mol_smiles, 'Tanimoto Similarity': mol_scores, 'Compatible': tally})
    compatible_mols = compatible_mols[compatible_mols['Compatible'] == 'proceed']
    #print compatible_mols['Alternative Substrate SMILES'].tolist()
    mol_pictures = []
    for i in compatible_mols['Alternative Substrate SMILES'].tolist():
        mol_pictures.append( Chem.MolFromSmiles(i) )
    mol_pics = Draw.MolsToGridImage(mol_pictures, molsPerRow=8) 
    print "Finished"
    return compatible_mols, mol_pics

def tanimoto_candidates(target, steroidlist, threshold):
    '''given a list of compounds, will compare to your target---can be used for quick candidate truncation'''
    steroidmols = [Chem.MolFromSmiles(i) for i in steroidlist]

    for m in steroidmols: AllChem.Compute2DCoords(m)
    steroidlist_fps=[AllChem.GetMorganFingerprintAsBitVect(x,2) for x in steroidmols]
    
    #recombine endogenous structures with their scores
    sims = DataStructs.BulkTanimotoSimilarity(steroidlist_fps[0],steroidlist_fps)
    nbrs = sorted(zip(sims,steroidmols),reverse=True)
    
    #grab bottom 10% of matches
    negative_structures = [x[1] for x in nbrs[:20]]
    negative_smiles = []
    for i in negative_structures:
        negative_smiles.append( Chem.MolToSmiles(i) )
        
    nbrs_filtered = []
    for i in nbrs:
        if i[0] > threshold:
            nbrs_filtered.append( i )
    #Draw.MolsToGridImage([x[1] for x in nbrs_filtered[:]],legends=['%.4f'%x[0] for x in nbrs_filtered])        
    return nbrs_filtered


def product_prediction_algorithm(target_smiles_list, df, row):

    start = timeit.default_timer()
    targets = []
    substrates = []
    products = []
    results = []
    enzymes = []
    for x in xrange(len(df)):
        try:
            results.append( explore_substrate(target_smiles_list[x], row, df) )
            substrates.append( df['Substrates'].irow(row) )
            products.append( df['Products'].irow(row) )
            targets.append( target_smiles_list[x] )
            enzymes.append( df['Enzymes'].irow(row) )
        except:
            pass
        
    stop = timeit.default_timer()
    #unique_prods, unique_substrates, unique_products, unique_enzymes = filter_results(results, substrates, products, TestUni, enzymes)
    unique_prods, unique_substrates, unique_products, unique_enzymes = results, substrates, products, enzymes
    
    print "Finished reactions... searching hits on PubChem... "
    
    #search pubchem for unique hits
    novel_compounds = []
    for i in unique_prods:
	try:
	    searches = pcp.get_compounds('CanonicalSMILES', str(Chem.MolToSmiles(i)), 'smiles')
	    if str(searches) == '[Compound()]':
	        novel_compounds.append( i )
	except:
	    pass
            
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



def prediction_algorithm2(target_smiles, df):
    '''iterates through an entire df instead of a target smiles list'''
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
            if z % 10 == 0 and z != 0:
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

