# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# $Id$
#
#  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met: 
#
#     * Redistributions of source code must retain the above copyright 
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following 
#       disclaimer in the documentation and/or other materials provided 
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
#       nor the names of its contributors may be used to endorse or promote 
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Created by Sereina Riniker, Aug 2013


from rdkit import Chem
from rdkit import RDConfig
from rdkit import DataStructs
from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit.Chem import rdmolops
from rdkit.Chem import Draw
import numpy
import math
import copy
from matplotlib import cm


def GetAtomicWeightsForFingerprint(refMol, probeMol, fpFunction, metric=DataStructs.DiceSimilarity):
  """
  Calculates the atomic weights for the probe molecule 
  based on a fingerprint function and a metric.

  Parameters:
    refMol -- the reference molecule
    probeMol -- the probe molecule
    fpFunction -- the fingerprint function
    metric -- the similarity metric

  Note:
    If fpFunction needs additional parameters, use a lambda construct
  """
  if hasattr(probeMol, '_fpInfo'): delattr(probeMol, '_fpInfo')
  if hasattr(refMol, '_fpInfo'): delattr(refMol, '_fpInfo')
  refFP = fpFunction(refMol, -1)
  probeFP = fpFunction(probeMol, -1)
  baseSimilarity = metric(refFP, probeFP)
  # loop over atoms
  weights = []
  for atomId in range(probeMol.GetNumAtoms()):
    newFP = fpFunction(probeMol, atomId)
    newSimilarity = metric(refFP, newFP)
    weights.append(baseSimilarity - newSimilarity)
  if hasattr(probeMol, '_fpInfo'): delattr(probeMol, '_fpInfo')
  if hasattr(refMol, '_fpInfo'): delattr(refMol, '_fpInfo')
  return weights


def GetAtomicWeightsForModel(probeMol, fpFunction, predictionFunction):
  """
  Calculates the atomic weights for the probe molecule based on 
  a fingerprint function and the prediction function of a ML model.

  Parameters:
    probeMol -- the probe molecule
    fpFunction -- the fingerprint function
    predictionFunction -- the prediction function of the ML model
  """
  if hasattr(probeMol, '_fpInfo'): delattr(probeMol, '_fpInfo')
  probeFP = fpFunction(probeMol, -1)
  baseProba = predictionFunction(probeFP)
  # loop over atoms
  weights = []
  for atomId in range(probeMol.GetNumAtoms()):
    newFP = fpFunction(probeMol, atomId)
    newProba = predictionFunction(newFP)
    weights.append(baseProba - newProba)
  if hasattr(probeMol, '_fpInfo'): delattr(probeMol, '_fpInfo')
  return weights


def GetStandardizedWeights(weights):
  """
  Normalizes the weights,
  such that the absolute maximum weight equals 1.0.

  Parameters:
    weights -- the list with the atomic weights
  """
  tmp = [math.fabs(w) for w in weights]
  currentMax = max(tmp)
  if currentMax > 0:
    return [w/currentMax for w in weights], currentMax
  else:
    return weights, currentMax


def GetSimilarityMapFromWeights(mol, weights, colorMap=cm.PiYG, scale=-1, size=(250, 250), sigma=None, 
                                coordScale=1.5, step=0.01, colors='k', contourLines=10, alpha=0.5, **kwargs):
  """
  Generates the similarity map for a molecule given the atomic weights.

  Parameters:
    mol -- the molecule of interest
    colorMap -- the matplotlib color map scheme
    scale -- the scaling: scale < 0 -> the absolute maximum weight is used as maximum scale
                          scale = double -> this is the maximum scale
    size -- the size of the figure
    sigma -- the sigma for the Gaussians
    coordScale -- scaling factor for the coordinates
    step -- the step for calcAtomGaussian
    colors -- color of the contour lines
    contourLines -- if integer number N: N contour lines are drawn
                    if list(numbers): contour lines at these numbers are drawn
    alpha -- the alpha blending value for the contour lines
    kwargs -- additional arguments for drawing
  """
  if mol.GetNumAtoms() < 2: raise ValueError("too few atoms")
  fig = Draw.MolToMPL(mol, coordScale=coordScale, size=size, **kwargs)
  if sigma is None:
    if mol.GetNumBonds() > 0:
      bond = mol.GetBondWithIdx(0)
      idx1 = bond.GetBeginAtomIdx()
      idx2 = bond.GetEndAtomIdx()
      sigma = 0.3 * math.sqrt(sum([(mol._atomPs[idx1][i]-mol._atomPs[idx2][i])**2 for i in range(2)]))
    else:
      sigma = 0.3 * math.sqrt(sum([(mol._atomPs[0][i]-mol._atomPs[1][i])**2 for i in range(2)]))
    sigma = round(sigma, 2)
  x, y, z = Draw.calcAtomGaussians(mol, sigma, weights=weights, step=step)
  # scaling
  if scale <= 0.0: maxScale = max(math.fabs(numpy.min(z)), math.fabs(numpy.max(z)))
  else: maxScale = scale
  # coloring
  fig.axes[0].imshow(z, cmap=colorMap, interpolation='bilinear', origin='lower', extent=(0,1,0,1), vmin=-maxScale, vmax=maxScale)
  # contour lines
  # only draw them when at least one weight is not zero
  if len([w for w in weights if w != 0.0]):
      fig.axes[0].contour(x, y, z, contourLines, colors=colors, alpha=alpha, **kwargs)
  return fig


def GetSimilarityMapForFingerprint(refMol, probeMol, fpFunction, metric=DataStructs.DiceSimilarity, **kwargs):
  """
  Generates the similarity map for a given reference and probe molecule, 
  fingerprint function and similarity metric.

  Parameters:
    refMol -- the reference molecule
    probeMol -- the probe molecule
    fpFunction -- the fingerprint function
    metric -- the similarity metric.
    kwargs -- additional arguments for drawing
  """
  weights = GetAtomicWeightsForFingerprint(refMol, probeMol, fpFunction, metric)
  weights, maxWeight = GetStandardizedWeights(weights)
  fig = GetSimilarityMapFromWeights(probeMol, weights, **kwargs)
  return fig, maxWeight


def GetSimilarityMapForModel(probeMol, fpFunction, predictionFunction, **kwargs):
  """
  Generates the similarity map for a given ML model and probe molecule, 
  and fingerprint function.

  Parameters:
    probeMol -- the probe molecule
    fpFunction -- the fingerprint function
    predictionFunction -- the prediction function of the ML model
    kwargs -- additional arguments for drawing
  """
  weights = GetAtomicWeightsForModel(probeMol, fpFunction, predictionFunction)
  weights, maxWeight = GetStandardizedWeights(weights)
  fig = GetSimilarityMapFromWeights(probeMol, weights, **kwargs)
  return fig, maxWeight
  

apDict = {}
apDict['normal'] = lambda m, bits, minl, maxl, bpe, ia: rdMD.GetAtomPairFingerprint(m, minLength=minl, maxLength=maxl, ignoreAtoms=ia)
apDict['hashed'] = lambda m, bits, minl, maxl, bpe, ia: rdMD.GetHashedAtomPairFingerprint(m, nBits=bits, minLength=minl, maxLength=maxl, ignoreAtoms=ia)
apDict['bv'] = lambda m, bits, minl, maxl, bpe, ia: rdMD.GetHashedAtomPairFingerprintAsBitVect(m, nBits=bits, minLength=minl, maxLength=maxl, nBitsPerEntry=bpe, ignoreAtoms=ia)

# usage:   lambda m,i: GetAPFingerprint(m, i, fpType, nBits, minLength, maxLength, nBitsPerEntry)
def GetAPFingerprint(mol, atomId=-1, fpType='normal', nBits=2048, minLength=1, maxLength=30, nBitsPerEntry=4):
  """
  Calculates the atom pairs fingerprint with the torsions of atomId removed.

  Parameters:
    mol -- the molecule of interest
    atomId -- the atom to remove the pairs for (if -1, no pair is removed)
    fpType -- the type of AP fingerprint ('normal', 'hashed', 'bv')
    nBits -- the size of the bit vector (only for fpType='bv')
    minLength -- the minimum path length for an atom pair
    maxLength -- the maxmimum path length for an atom pair
    nBitsPerEntry -- the number of bits available for each pair
  """
  if fpType not in ['normal', 'hashed', 'bv']: raise ValueError("Unknown Atom pairs fingerprint type")
  if atomId < 0:
    return apDict[fpType](mol, nBits, minLength, maxLength, nBitsPerEntry, 0)
  if atomId >= mol.GetNumAtoms(): raise ValueError("atom index greater than number of atoms")
  return apDict[fpType](mol, nBits, minLength, maxLength, nBitsPerEntry, [atomId])


ttDict = {}
ttDict['normal'] = lambda m, bits, ts, bpe, ia: rdMD.GetTopologicalTorsionFingerprint(m, targetSize=ts, ignoreAtoms=ia)
ttDict['hashed'] = lambda m, bits, ts, bpe, ia: rdMD.GetHashedTopologicalTorsionFingerprint(m, nBits=bits, targetSize=ts, ignoreAtoms=ia)
ttDict['bv'] = lambda m, bits, ts, bpe, ia: rdMD.GetHashedTopologicalTorsionFingerprintAsBitVect(m, nBits=bits, targetSize=ts, nBitsPerEntry=bpe, ignoreAtoms=ia)

# usage:   lambda m,i: GetTTFingerprint(m, i, fpType, nBits, targetSize)
def GetTTFingerprint(mol, atomId=-1, fpType='normal', nBits=2048, targetSize=4, nBitsPerEntry=4):
  """
  Calculates the topological torsion fingerprint with the pairs of atomId removed.

  Parameters:
    mol -- the molecule of interest
    atomId -- the atom to remove the torsions for (if -1, no torsion is removed)
    fpType -- the type of TT fingerprint ('normal', 'hashed', 'bv')
    nBits -- the size of the bit vector (only for fpType='bv')
    minLength -- the minimum path length for an atom pair
    maxLength -- the maxmimum path length for an atom pair
    nBitsPerEntry -- the number of bits available for each torsion
  """
  if fpType not in ['normal', 'hashed', 'bv']: raise ValueError("Unknown Topological torsion fingerprint type")
  if atomId < 0:
    return ttDict[fpType](mol, nBits, targetSize, nBitsPerEntry, 0)
  if atomId >= mol.GetNumAtoms(): raise ValueError("atom index greater than number of atoms")
  return ttDict[fpType](mol, nBits, targetSize, nBitsPerEntry, [atomId])


# usage:   lambda m,i: GetMorganFingerprint(m, i, radius, fpType, nBits, useFeatures)
def GetMorganFingerprint(mol, atomId=-1, radius=2, fpType='bv', nBits=2048, useFeatures=False):
  """
  Calculates the Morgan fingerprint with the counts of atomId removed.

  Parameters:
    mol -- the molecule of interest
    radius -- the maximum radius
    fpType -- the type of Morgan fingerprint: 'count' or 'bv'
    atomId -- the atom to remove the counts for (if -1, no count is removed)
    nBits -- the size of the bit vector (only for fpType = 'bv')
    useFeatures -- if false: ConnectivityMorgan, if true: FeatureMorgan
  """
  if fpType not in ['bv', 'count']: raise ValueError("Unknown Morgan fingerprint type")
  if not hasattr(mol, '_fpInfo'):
    info = {}
    # get the fingerprint
    if fpType == 'bv': molFp = rdMD.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits, useFeatures=useFeatures, bitInfo=info)
    else: molFp = rdMD.GetMorganFingerprint(mol, radius, useFeatures=useFeatures, bitInfo=info)
    # construct the bit map
    if fpType == 'bv': bitmap = [DataStructs.ExplicitBitVect(nBits) for x in range(mol.GetNumAtoms())]
    else: bitmap = [[] for x in range(mol.GetNumAtoms())]
    for bit, es in info.iteritems():
      for at1, rad in es:
        if rad == 0: # for radius 0
          if fpType == 'bv': bitmap[at1][bit] = 1
          else: bitmap[at1].append(bit)
        else: # for radii > 0
          env = Chem.FindAtomEnvironmentOfRadiusN(mol, rad, at1)
          amap = {}
          submol = Chem.PathToSubmol(mol, env, atomMap=amap)
          for at2 in amap.keys():
            if fpType == 'bv': bitmap[at2][bit] = 1
            else: bitmap[at2].append(bit)
    mol._fpInfo = (molFp, bitmap)

  if atomId < 0:
    return mol._fpInfo[0]
  else: # remove the bits of atomId
    if atomId >= mol.GetNumAtoms(): raise ValueError("atom index greater than number of atoms")
    if len(mol._fpInfo) != 2: raise ValueError("_fpInfo not set")
    if fpType == 'bv':
      molFp = mol._fpInfo[0] ^ mol._fpInfo[1][atomId] # xor
    else: # count
      molFp = copy.deepcopy(mol._fpInfo[0])
      # delete the bits with atomId
      for bit in mol._fpInfo[1][atomId]:
        molFp[bit] -= 1
    return molFp

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

# <codecell>

def rf_classifier(enzymes, success, candidate):
    # generate fingeprints: Morgan fingerprint with radius 2
    mols = []
    for i in enzymes:
        mols.append( Chem.MolFromSmiles(i) )
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in mols]  
    
    # convert the RDKit explicit vectors into numpy arrays
    np_fps = []
    for fp in fps:
      arr = numpy.zeros((1,))
      DataStructs.ConvertToNumpyArray(fp, arr)
      np_fps.append(arr)
    
    # get a random forest classifiert with 100 trees
    rf = RandomForestClassifier(n_estimators=100, random_state=1123)
    
    # train the random forest
    # with the first two molecules being actives (class 1) and
    # the last two being inactives (class 0)
    ys_fit = success
    rf.fit(np_fps, ys_fit)
    
    # use the random forest to predict a new molecule
    m6 = Chem.MolFromSmiles(candidate)
    fp = numpy.zeros((1,))
    DataStructs.ConvertToNumpyArray(AllChem.GetMorganFingerprintAsBitVect(m6, 2), fp)
    
    if rf.predict_proba(fp)[0][1] > 0.95:
        return '[1]'
    else:
        return '[0]'

def rf_validate(enzymes, success, candidate):
    # generate fingeprints: Morgan fingerprint with radius 2
    mols = []
    for i in enzymes:
        mols.append( Chem.MolFromSmiles(i) )
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in mols]  
    
    # convert the RDKit explicit vectors into numpy arrays
    np_fps = []
    for fp in fps:
      arr = numpy.zeros((1,))
      DataStructs.ConvertToNumpyArray(fp, arr)
      np_fps.append(arr)
    
    # get a random forest classifiert with 100 trees
    rf = RandomForestClassifier(n_estimators=100, random_state=1123)
    
    # train the random forest
    # with the first two molecules being actives (class 1) and
    # the last two being inactives (class 0)
    ys_fit = success
    rf.fit(np_fps, ys_fit)
    
    # use the random forest to predict a new molecule
    m6 = Chem.MolFromSmiles(candidate)
    fp = numpy.zeros((1,))
    DataStructs.ConvertToNumpyArray(AllChem.GetMorganFingerprintAsBitVect(m6, 2), fp)
    
    return rf.predict_proba(fp)[0][1] 


def rf_classifier2(enzymes, success, candidate):
    # generate fingeprints: Morgan fingerprint with radius 2
    mols = []
    for i in enzymes:
        mols.append( Chem.MolFromSmiles(i) )
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in mols]  
    
    # convert the RDKit explicit vectors into numpy arrays
    np_fps = []
    for fp in fps:
      arr = numpy.zeros((1,))
      DataStructs.ConvertToNumpyArray(fp, arr)
      np_fps.append(arr)
    
    # get a random forest classifiert with 100 trees
    rf = RandomForestClassifier(n_estimators=100, random_state=1123)
    
    # train the random forest
    # with the first two molecules being actives (class 1) and
    # the last two being inactives (class 0)
    ys_fit = success
    rf.fit(np_fps, ys_fit)
    
    # use the random forest to predict a new molecule
    m6 = Chem.MolFromSmiles(candidate)
    fp = numpy.zeros((1,))
    DataStructs.ConvertToNumpyArray(AllChem.GetMorganFingerprintAsBitVect(m6, 2), fp)
    
    #return rf.predict_proba(fp)
    #return rf.predict(fp)

    # helper function
    def getProba(fp, predictionFunction):
      return predictionFunction(fp)[0][1]
    
    fig, maxweight = GetSimilarityMapForModel(m6, GetMorganFingerprint, lambda x: getProba(x, rf.predict_proba))
    return fig

# <codecell>

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestClassifier
import numpy
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
import pandas as pd

# <codecell>

def rf_classifier2(enzymes, success, candidate):
    # generate fingeprints: Morgan fingerprint with radius 2
    mols = []
    for i in enzymes:
        mols.append( Chem.MolFromSmiles(i) )
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in mols]  
    
    # convert the RDKit explicit vectors into numpy arrays
    np_fps = []
    for fp in fps:
      arr = numpy.zeros((1,))
      DataStructs.ConvertToNumpyArray(fp, arr)
      np_fps.append(arr)
    
    # get a random forest classifiert with 100 trees
    rf = RandomForestClassifier(n_estimators=100, random_state=1123)
    
    # train the random forest
    # with the first two molecules being actives (class 1) and
    # the last two being inactives (class 0)
    ys_fit = success
    rf.fit(np_fps, ys_fit)
    
    # use the random forest to predict a new molecule
    m6 = Chem.MolFromSmiles(candidate)
    fp = numpy.zeros((1,))
    DataStructs.ConvertToNumpyArray(AllChem.GetMorganFingerprintAsBitVect(m6, 2), fp)
    
    #return rf.predict_proba(fp)
    #return rf.predict(fp)

    # helper function
    def getProba(fp, predictionFunction):
      return predictionFunction(fp)[0][1]
    
    fig, maxweight = GetSimilarityMapForModel(m6, GetMorganFingerprint, lambda x: getProba(x, rf.predict_proba))
    return fig

