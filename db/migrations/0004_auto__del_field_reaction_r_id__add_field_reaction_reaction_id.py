# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Deleting field 'Reaction.r_id'
        db.delete_column(u'db_reaction', 'r_id')

        # Adding field 'Reaction.reaction_id'
        db.add_column(u'db_reaction', 'reaction_id',
                      self.gf('django.db.models.fields.CharField')(default='UWPRGWTM', unique=True, max_length=8),
                      keep_default=False)


    def backwards(self, orm):

        # User chose to not deal with backwards NULL issues for 'Reaction.r_id'
        raise RuntimeError("Cannot reverse this migration. 'Reaction.r_id' and its values cannot be restored.")
        
        # The following code is provided here to aid in writing a correct migration        # Adding field 'Reaction.r_id'
        db.add_column(u'db_reaction', 'r_id',
                      self.gf('django.db.models.fields.TextField')(max_length=8, unique=True),
                      keep_default=False)

        # Deleting field 'Reaction.reaction_id'
        db.delete_column(u'db_reaction', 'reaction_id')


    models = {
        u'db.gene': {
            'Meta': {'object_name': 'Gene'},
            'alternate_name': ('django.db.models.fields.TextField', [], {}),
            'entrez_id': ('django.db.models.fields.TextField', [], {}),
            'entrez_symbol': ('django.db.models.fields.TextField', [], {}),
            'gene_id': ('django.db.models.fields.CharField', [], {'default': "'UWPRGWTQ'", 'unique': 'True', 'max_length': '8'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.TextField', [], {})
        },
        u'db.reaction': {
            'Meta': {'object_name': 'Reaction'},
            'alternate_name': ('django.db.models.fields.TextField', [], {}),
            'chemical_equation': ('django.db.models.fields.TextField', [], {}),
            'comment': ('django.db.models.fields.TextField', [], {}),
            'confidence_score': ('django.db.models.fields.TextField', [], {}),
            'directionality': ('django.db.models.fields.TextField', [], {}),
            'ec_number': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'kegg_id': ('django.db.models.fields.TextField', [], {}),
            'name': ('django.db.models.fields.TextField', [], {}),
            'reaction_id': ('django.db.models.fields.CharField', [], {'default': "'UWPRGWTP'", 'unique': 'True', 'max_length': '8'}),
            'reference': ('django.db.models.fields.TextField', [], {})
        },
        u'db.smallmolecule': {
            'IUPAC_Name': ('django.db.models.fields.TextField', [], {}),
            'Meta': {'object_name': 'SmallMolecule'},
            'Metabolic_Network_ids': ('django.db.models.fields.TextField', [], {}),
            'alternative_names': ('django.db.models.fields.TextField', [], {}),
            'chebi_id': ('django.db.models.fields.CharField', [], {'max_length': '15'}),
            'drugbank_id': ('django.db.models.fields.CharField', [], {'max_length': '15'}),
            'hmdb_id': ('django.db.models.fields.CharField', [], {'max_length': '15'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'inchi': ('django.db.models.fields.TextField', [], {}),
            'inchi_key': ('django.db.models.fields.TextField', [], {}),
            'kegg_id': ('django.db.models.fields.CharField', [], {'max_length': '15'}),
            'name': ('django.db.models.fields.TextField', [], {}),
            'pubchem_cid': ('django.db.models.fields.CharField', [], {'max_length': '15'}),
            'sm_id': ('django.db.models.fields.CharField', [], {'default': "'UWPRGWTK'", 'unique': 'True', 'max_length': '8'}),
            'smiles': ('django.db.models.fields.TextField', [], {})
        },
        u'db.substance': {
            'Meta': {'object_name': 'Substance'},
            'comment': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        }
    }

    complete_apps = ['db']