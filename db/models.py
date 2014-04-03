from __future__ import unicode_literals

from django.db import models
from django.db.models.query import QuerySet
from django.db.models.sql.query import Query
from django.db.models.sql.compiler import SQLCompiler

import logging

logger = logging.getLogger(__name__)

class SmallMolecule(models.Model):
    sm_id = models.CharField(max_length=8, unique=True)
    name = models.TextField()
    smiles = models.TextField()
    inchi = models.TextField()
    inchi_key = models.TextField()
    IUPAC_Name = models.TextField()
    alternative_names = models.TextField()
    Metabolic_Network_ids = models.TextField()
    pubchem_cid = models.CharField(max_length=15)
    chebi_id = models.CharField(max_length=15)
    kegg_id = models.CharField(max_length=15)
    hmdb_id = models.CharField(max_length=15)
    drugbank_id = models.CharField(max_length=15)
    
class Gene(models.Model):
    gene_id = models.CharField(max_length=8, unique=True)
    name = models.TextField()
    alternate_name = models.TextField()
    entrez_id = models.TextField()
    entrez_symbol = models.TextField()
    

class Reaction(models.Model):
    name = models.TextField()
    r_id = models.TextField(max_length=8, unique=True)
    chemical_equation = models.TextField()
    #substrate(s)
    #product(s)
    directionality = models.TextField()
    #enzyme(s)
    #subsystem(s)
    confidence_score = models.TextField()
    kegg_id = models.TextField()
    ec_number = models.TextField()
    alternate_name = models.TextField()
    reference = models.TextField()
    comment = models.TextField()

    