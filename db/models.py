from __future__ import unicode_literals
import logging

from django.db import models
# from django.db.models.query import QuerySet
# from django.db.models.sql.query import Query
# from django.db.models.sql.compiler import SQLCompiler
from django.db import connection

logger = logging.getLogger(__name__)

from reports.utils.gray_codes import create_substance_id

def create_id():
    # TODO: There is no way to pass an argument to
    # to the default on a model, so this is why we have to grope around for the
    # sequence value.  The sequence here is used with the model idfield  
    try:
        cursor = connection.cursor()
        cursor.execute("SELECT nextval('db_substance_id_seq')")
        row = cursor.fetchone();
        val = row[0]
    
        new_id = create_substance_id(val)
        if(logger.isEnabledFor(logging.DEBUG)):
            logger.debug(str(('created_id', val, new_id)))
    #     if val > 1:
    #         cursor.execute("SELECT setval('db_substance_id_seq', %s)", [val-1])
        return new_id
    except Exception, e:
        logger.warn(str(('create_substance_id fails', e)))
        return None
        
# create a substance table, as an easy way of creating the db_substance_id_seq
class Substance(models.Model):
    comment = models.TextField()
    
    
class SmallMolecule(models.Model):
    name = models.TextField()
    
    sm_id = models.CharField(
        max_length=8, unique=True, 
        default=create_id)
    
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
    
    gene_id = models.CharField(
        max_length=8, unique=True, 
        default=create_id)
    
    name = models.TextField()
    alternate_name = models.TextField()
    entrez_id = models.TextField()
    entrez_symbol = models.TextField()
    sequence = models.TextField()
        
class Protein(models.Model):
    
    protein_id = models.CharField(
        max_length=8, unique=True, 
        default=create_id)
    
    uniprot_id = models.CharField(max_length=15)
    amino_acid_sequence = models.TextField()
    gene = models.ForeignKey('Gene')

class Reaction(models.Model):
        
    reaction_id = models.CharField(
        max_length=8, unique=True, 
        default=create_id)
    
    name = models.TextField()
    chemical_equation = models.TextField()
    #substrate(s)
    #product(s)
    directionality = models.TextField()
    #enzyme(s)
    #subsystem(s)
    confidence_score = models.TextField(blank=True)
    kegg_id = models.TextField(blank=True)
    ec_number = models.TextField(blank=True)
    alternate_name = models.TextField(blank=True)
    reference = models.TextField(blank=True)
    comment = models.TextField(blank=True)

class ExpressionHost(models.Model):
    host_id = models.CharField(
        max_length=8, unique=True, 
        default=create_id)

    name = models.TextField()
    taxonomic_id = models.TextField()
    taxonomic_branch = models.TextField()
    genotype = models.CharField(max_length=32)
    genome = models.TextField()
    expression_data = models.TextField(blank=True)
    metabolic_model = models.TextField(blank=True)
    reference = models.IntegerField(null=True, blank=True)
    
