# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Process book for the steroid algorithm code.

# <markdowncell>

# Let's get all the relevant smiles for the proteins we're going to use. 

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
import pickle
from bioservices import Kegg
from bioservices import *
from compiler.ast import flatten
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import matplotlib.pyplot as plt

# <headingcell level=1>

# Bioservices Magic

# <markdowncell>

# Upload steroids.

# <codecell>

import sys                           
sys.path.append('/Home')
from RF_functions import *

steroidlist = pickle.load(open('steroidlist', 'rb'))
Keggsteroids_df = pd.read_table("Keggsteroids.csv")
endogenous_steroids = Keggsteroids_df['Substrates'].tolist()
endogenous_steroids2 = Keggsteroids_df['Products'].tolist()
for i in endogenous_steroids2:
    endogenous_steroids.append( i )
endogenous_steroids = list(set(endogenous_steroids))
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import gzip
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw

# <markdowncell>

# Search for compound on KEGG.

# <codecell>

retrieve_enzyme_substrates('ec:1.14.15.6')

# <codecell>

s = Kegg()
compound_search = (s.get(str('ec:1.14.15.6')))
p = KeggParser()
d = p.parse(compound_search)

# <codecell>

def retrieve_enzyme_substrates(enzyme_name):
    s = Kegg()
    compound_search = (s.get(str(enzyme_name)))
    p = KeggParser()
    d = p.parse(compound_search)
    
    ####################### grabs all the reactions that are dealt with by this enzyme
    all_rxns = []
    
    if d['all_reac'][0] == 'R': ####filter to make sure that we do have a list, if the first thing is R then there's only one reaction and no need to iterate
        all_rxns.append( d['all_reac'] )
    else:
        for i in d['all_reac']:
            for j in i.split():
                if j[0] == 'R':
                    if j[-1] == ';':
                        j = j[:-1]
                    all_rxns.append( j )
          
    ####################### grabs all the smiles of the substrates that are involved with these reactions            
    smiles = []
    real_name = []
    
    for i in all_rxns:
        reactants = [] ##build reactants
        reaction = s.get(i)
        for i in reaction.split('\n'):
            if 'EQUATION' in i:
                for j in range(len( i.split() )):
                    if i.split()[j][0] == '<':
                        eq = j
                for j in range(len( i.split() )):
                    if i.split()[j][0] == "C" and j < eq:
                        reactants.append( i.split()[j] ) 
        reactant = max(reactants)

        reactantinfo = s.get(reactant)
        for i in reactantinfo.split('\n'):
            if 'NAME' in i:
                reactantname = i.split()[1]
        if reactantname[-1] == ';':
            reactantname = reactantname[:-1]
        if reactantname != 'Reduced':
            real_name.append( reactantname )
        try:
            test = pcp.get_properties('CanonicalSMILES', reactantname, 'name')
            smiles.append( str(test[0]['CanonicalSMILES']) )  
        except:
            pass
    
    return real_name, smiles

Keggsteroids_df = pd.read_table("Keggsteroids.csv")

# <headingcell level=1>

# Workflow

# <headingcell level=4>

# Retrieve substrates given an enzyme

# <codecell>

s = Kegg()
s.find('enzyme', 'cyp7a1')

# <codecell>

compound_names, positive_smiles = retrieve_enzyme_substrates('ec:1.14.13.17')

print compound_names
print positive_smiles

# <headingcell level=4>

# Compute the coordinates of both the positive hits and all the endogenous steroids

# <codecell>

positive_structures = []
for i in positive_smiles:
    positive_structures.append( Chem.MolFromSmiles(i) )
endogenous_structures = []
for i in endogenous_steroids:
    endogenous_structures.append( Chem.MolFromSmiles(i) )

for m in positive_structures: AllChem.Compute2DCoords(m)
for m in endogenous_structures: AllChem.Compute2DCoords(m)
    
positive_fps=[AllChem.GetMorganFingerprintAsBitVect(x,2) for x in positive_structures]
endogenous_fps=[AllChem.GetMorganFingerprintAsBitVect(x,2) for x in endogenous_structures]

# <codecell>

#recombine endogenous structures with their scores
sims = DataStructs.BulkTanimotoSimilarity(positive_fps[0],endogenous_fps)
nbrs = sorted(zip(sims,endogenous_structures),reverse=False)

# <codecell>

#grab bottom 10% of matches
negative_structures = [x[1] for x in nbrs[:12]]
negative_smiles = []
for i in negative_structures:
    negative_smiles.append( Chem.MolToSmiles(i) )

# <codecell>

#Draw.MolsToGridImage([x[1] for x in nbrs[:12]],legends=['%.4f'%x[0] for x in nbrs])

# <headingcell level=4>

# Let's code the cross validation section

# <codecell>

activity = []
smiles = []


for i in positive_smiles:
    smiles.append( i )
for i in negative_smiles:
    smiles.append( i )
    
for i in range(len(positive_smiles)):
    activity.append( 1 )
for i in range(len(negative_smiles)):
    activity.append( 0 )
    
for i in negative_smiles:
    try:
        results = pcp.get_compounds(str(i), 'smiles')
        cpn = pcp.Compound.from_cid( int(str(results[0])[9:-1]) )
        compound_names.append( cpn.iupac_name )   
    except:
        compound_names.append( 'error parsing' )

# <codecell>

training_set = pd.DataFrame({'Smiles': smiles, 'Activity': activity, 'Compound Name': compound_names})

# <codecell>

import numpy as np
import pandas as pd
from random import sample

# <headingcell level=4>

# Scrambles up the training_set dataframe so I can test the rest of it

# <codecell>

percentages = list(np.arange(0, 1.05, 0.05))

validation_per_percent = []
for iteration in percentages:
    cross_validation_scores = []
    for x in xrange(100):
        # given data frame df
        # create random index
        rindex =  np.array(sample(xrange(len(training_set)), int( len(training_set)*.25 )))
        
        #get some random rows from df
        pulled_data = training_set.ix[rindex]
        pulled_smiles = pulled_data['Smiles'].tolist()
        training_set_culled = training_set.drop(rindex)
        
        predicted_values = [] 
        actual_values = pulled_data['Activity'].tolist()
        mols = training_set_culled['Smiles'].tolist()
        activities = training_set_culled['Activity'].tolist()
        
        try:
            for j in range(len(pulled_smiles)):
                if rf_validate(mols, activities, pulled_smiles[j]) > iteration:
                    if pulled_smiles[j] not in mols:
                        predicted_values.append( 1 )
                else:
                    predicted_values.append( 0 ) 
            #cross_validation_scores.append( predicted_values == actual_values )
            cross_validation_scores.append( predicted_values == actual_values )
        except:
            pass   
    #validation_per_percent.append( float(cross_validation_scores.count(True))/len(cross_validation_scores) )
    validation_per_percent.append( cross_validation_scores )

# <codecell>

validation_per_percent

# <codecell>

import numpy.numarray as na

from pylab import *

labels = ["Baseline", "System"]
data = validation_per_percent

xlocations = na.array(range(len(data)))+0.5
width = 0.5
bar(xlocations, data, width=width)
yticks(range(0, 1))
xticks(xlocations+ width/2, labels)
xlim(0, xlocations[-1]+width*2)
title("Average Ratings on the Training Set")
gca().get_xaxis().tick_bottom()
gca().get_yaxis().tick_left()

show()

# <codecell>

rf_classifier(mols, activities, pulled_smiles[2])

# <codecell>

vis = []
for i in values:
    vis.append( Chem.MolFromSmiles(i) )
Draw.MolsToGridImage(vis,molsPerRow=8, includeAtomNumbers=False)

# <codecell>

def find_alternative_substrates(enzyme):
    
    ###grabs the substrates of the given enzyme and their smiles
    compound_names, positive_smiles = retrieve_enzyme_substrates(enzyme) 
    
    ###generate 3D coordinates for each molecule
    positive_structures = []
    for i in positive_smiles:
        positive_structures.append( Chem.MolFromSmiles(i) )
    endogenous_structures = []
    for i in endogenous_steroids:
        endogenous_structures.append( Chem.MolFromSmiles(i) )
    
    for m in positive_structures: AllChem.Compute2DCoords(m)
    for m in endogenous_structures: AllChem.Compute2DCoords(m)
        
    positive_fps=[AllChem.GetMorganFingerprintAsBitVect(x,2) for x in positive_structures]
    endogenous_fps=[AllChem.GetMorganFingerprintAsBitVect(x,2) for x in endogenous_structures]
    
    #recombine endogenous structures with their scores
    sims = DataStructs.BulkTanimotoSimilarity(positive_fps[0],endogenous_fps)
    nbrs = sorted(zip(sims,endogenous_structures),reverse=False) 
    
    #grab bottom 10% of matches
    negative_structures = [x[1] for x in nbrs[:12]]
    negative_smiles = []
    for i in negative_structures:
        negative_smiles.append( Chem.MolToSmiles(i) )
        
    #Draw.MolsToGridImage([x[1] for x in nbrs[:12]],legends=['%.4f'%x[0] for x in nbrs])
    
    ###gets the information for each molecule that we've grabbed thus far
    activity = []
    smiles = []
    
    for i in positive_smiles:
        smiles.append( i )
    for i in negative_smiles:
        smiles.append( i )
        
    for i in range(len(positive_smiles)):
        activity.append( 1 )
    for i in range(len(negative_smiles)):
        activity.append( 0 )
        
    for i in negative_smiles:
        try:
            results = pcp.get_compounds(str(i), 'smiles')
            cpn = pcp.Compound.from_cid( int(str(results[0])[9:-1]) )
            compound_names.append( cpn.iupac_name )   
        except:
            compound_names.append( 'error parsing' )
    
    ###our model
    training_set = pd.DataFrame({'Smiles': smiles, 'Activity': activity, 'Compound Name': compound_names})
    
    ###identity is hard coded in the other file as 70%
    values = [] 
    mols = training_set['Smiles'].tolist()
    activities = training_set['Activity'].tolist()
    
    for i in range(len(steroidlist)):
        if str( rf_classifier(mols, activities, steroidlist[i]) ) == '[1]':
            if steroidlist[i] not in mols:
                values.append( steroidlist[i] )
                
    return values, training_set

# <codecell>

s = Kegg()
s.find('enzyme', 'STE1')

# <codecell>

values, training_set = find_alternative_substrates('ec:1.1.1.181')

# <codecell>

valuesi = []
for i in values:
    i = Chem.MolFromSmiles(i)
    i = Chem.MolToSmiles( i, canonical=True ) 
    if '.' not in i:
        valuesi.append( i )

vis = []
for i in list(set(valuesi)):
    vis.append( Chem.MolFromSmiles(i) )
Draw.MolsToGridImage(vis,molsPerRow=8, includeAtomNumbers=False)

# <headingcell level=1>

# Quick candidate truncation.

# <codecell>

def tanimoto_candidates(target, steroidlist):
    '''given a list of compounds, will compare to your target'''
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
        if i[0] > .40:
            nbrs_filtered.append( i )
    #Draw.MolsToGridImage([x[1] for x in nbrs_filtered[:]],legends=['%.4f'%x[0] for x in nbrs_filtered])        
    return nbrs_filtered

# <codecell>

Draw.MolsToGridImage([x[1] for x in nbrs_filtered[:]],legends=['%.4f'%x[0] for x in nbrs_filtered])

# <headingcell level=1>

# Grabbing enzymes that are associated with a molecule.

# <codecell>

s.find('compound', 'cholesterol')
print s.get(str('cpd:C00001'))

# <codecell>


