from __future__ import unicode_literals

from django.db import models
from django.db.models.query import QuerySet
from django.db.models.sql.query import Query
from django.db.models.sql.compiler import SQLCompiler

import logging

logger = logging.getLogger(__name__)


from lims.hms.gray_codes import create_substance_id
from django.db import connection
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
        logger.info(str(('created_id', val, new_id)))
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

    